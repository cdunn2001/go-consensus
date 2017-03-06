package main

/*
#include "common.h"
#include "foo.h"
#include <stdio.h>
#include <stdlib.h>
void foo() {
	fputs("there\n", stderr);
}
*/
import "C"

import "unsafe"

import "regexp"
import "sync"

//import "flag"
import "bufio"

//import "bytes"
import "io"
import "fmt"
import "github.com/jessevdk/go-flags"
import "os"
import "sort"

//import "strconv"
import "strings"

func Use(vals ...interface{}) {
	for _, val := range vals {
		_ = val
	}
}

//type NtSlice []byte
type NtSlice string // since we have to convert to null-terminated c-array

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
func get_longest_reads(seqs []NtSlice, max_n_read, max_cov_aln int) []NtSlice {
	longest_n_reads := max_n_read
	if max_cov_aln > 0 {
		longest_n_reads = 1
		seed_len := len(seqs[0])
		read_cov := 0
		for _, seq := range seqs[1:] {
			if read_cov/seed_len > max_cov_aln {
				break
			}
			longest_n_reads += 1
			read_cov += len(seq)
		}
		longest_n_reads = min(longest_n_reads, max_n_read)
	}
	// In Python, the slice is stopped at end. In Go, we must test.
	longest_n_reads = min(longest_n_reads, len(seqs))
	return seqs[:longest_n_reads]
}

type ByLongestNtSlice []NtSlice

func (slice ByLongestNtSlice) Len() int {
	return len(slice)
}
func (slice ByLongestNtSlice) Less(i, j int) bool {
	return len(slice[i]) > len(slice[j])
}
func (slice ByLongestNtSlice) Swap(i, j int) {
	slice[i], slice[j] = slice[j], slice[i]
}
func get_longest_sorted_reads(seqs []NtSlice, max_n_read, max_cov_aln int) []NtSlice {
	n := len(seqs)
	my_seqs := make([]NtSlice, n)
	copy(my_seqs, seqs)
	sort.Sort(ByLongestNtSlice(my_seqs[1:n]))
	return get_longest_reads(my_seqs, max_n_read, max_cov_aln)
}
func get_consensus_with_trim(datum SeqDatum) (string, string) {
	seqs := datum.seqs
	seed_id := datum.seed_id
	config := datum.config
	Use(seqs, seed_id, config)
	trim_seqs := make([]NtSlice, 0)
	seed := seqs[0]
	trim_seqs = append(trim_seqs, seed)
	for _, seq := range seqs[1:] {
		s1, e1, s2, e2, aln_size, aln_score, c_status := get_alignment(seq, seed, config.edge_tolerance)
		Use(s2, e2, aln_size)
		if c_status == "none" {
			continue
		}
		if (aln_score > 1000) && (e1-s1 > 500) {
			e1 -= config.trim_size
			s1 += config.trim_size
			//trim_seqs.append( (e1-s1, seq[s1:e1]) )
			// We do not need (e1-s1) because that is len(seq[s1:e1])
			trim_seqs = append(trim_seqs, seq[s1:e1])
		}
	}
	sort.Sort(ByLongestNtSlice(trim_seqs[1:len(trim_seqs)]))
	if len(trim_seqs[1:]) > config.max_n_read { // Is this correct? -cdunn
		// seqs already sorted, dont' sort again -- jason
		trim_seqs = get_longest_reads(trim_seqs, config.max_n_read, config.max_cov_aln)
	}
	cseqs := copy_seq_ptrs(trim_seqs)
	consensus_data_ptr := C.generate_consensus(&cseqs[0],
		C.uint(len(trim_seqs)),
		C.uint(config.min_cov),
		C.uint(config.K),
		C.double(config.min_idt))
	free_seq_ptrs(cseqs)
	consensus := C.GoString(consensus_data_ptr.sequence)
	C.free_consensus_data(consensus_data_ptr)
	return consensus, seed_id
}
func get_alignment(seq1 NtSlice, seq0 NtSlice, edge_tolerance int) (int, int, int, int, int, int, string) {
	var K uint32 = 8
	lk_ptr := C.allocate_kmer_lookup(C.seq_coor_t(1 << (K * 2)))
	sa_ptr := C.allocate_seq(C.seq_coor_t(len(seq0)))
	sda_ptr := C.allocate_seq_addr(C.seq_coor_t(len(seq0)))
	cseq0 := C.CString(string(seq0))
	cseq1 := C.CString(string(seq1))
	defer C.free(unsafe.Pointer(cseq0))
	defer C.free(unsafe.Pointer(cseq1))
	C.add_sequence(0, C.uint(K), cseq0, C.seq_coor_t(len(seq0)), sda_ptr, sa_ptr, lk_ptr)
	C.mask_k_mer(1<<(K*2), lk_ptr, 16)
	kmer_match_ptr := C.find_kmer_pos_for_seq(cseq1, C.seq_coor_t(len(seq1)), C.uint(K), sda_ptr, lk_ptr)
	defer C.free_kmer_match(kmer_match_ptr)
	aln_range_ptr := C.find_best_aln_range2(kmer_match_ptr, C.seq_coor_t(K), C.seq_coor_t(K*50), 25)
	defer C.free_aln_range(aln_range_ptr)
	//x,y = zip( * [ (kmer_match.query_pos[i], kmer_match.target_pos[i]) for i in range(kmer_match.count )] )
	cs1, ce1, cs0, ce0, km_score := aln_range_ptr.s1, aln_range_ptr.e1, aln_range_ptr.s2, aln_range_ptr.e2, aln_range_ptr.score
	s0 := int(cs0)
	s1 := int(cs1)
	e0 := int(ce0)
	e1 := int(ce1)
	e1 += int(K + K/2)
	e0 += int(K + K/2)
	len_1 := len(seq1)
	len_0 := len(seq0)
	if e1 > len_1 {
		e1 = len_1
	}
	if e0 > len_0 {
		e0 = len_0
	}
	aln_size := 1
	aln_score := -1
	aln_q_s := -1
	aln_q_e := -1
	aln_t_s := -1
	aln_t_e := -1
	if e1-s1 > 500 {
		aln_size = max(e1-s1, e0-s0)
		aln_score = int(km_score * 48)
		aln_q_s = s1
		aln_q_e = e1
		aln_t_s = s0
		aln_t_e = e0
	}

	C.free_seq_addr_array(sda_ptr)
	C.free_seq_array(sa_ptr)
	C.free_kmer_lookup(lk_ptr)
	if s1 > edge_tolerance && s0 > edge_tolerance {
		return 0, 0, 0, 0, 0, 0, "none"
	}
	if len_1-e1 > edge_tolerance && len_0-e0 > edge_tolerance {
		return 0, 0, 0, 0, 0, 0, "none"
	}

	if e1-s1 > 500 && aln_size > 500 {
		return s1, s1 + aln_q_e - aln_q_s, s0, s0 + aln_t_e - aln_t_s, aln_size, aln_score, "aln"
	}
	return 0, 0, 0, 0, 0, 0, "none"
}
func copy_seq_ptrs(seqs []NtSlice) []*C.char {
	result := make([]*C.char, len(seqs))
	for i, seq := range seqs {
		result[i] = (*C.char)(C.CString(string(seq)))
	}
	return result
}
func free_seq_ptrs(cseqs []*C.char) {
	for _, cseq := range cseqs {
		C.free(unsafe.Pointer(cseq))
	}
}
func get_consensus_without_trim(datum SeqDatum) (string, string) {
	seqs := datum.seqs
	seed_id := datum.seed_id
	config := datum.config
	Use(seqs, seed_id, config)
	if len(seqs) > config.max_n_read {
		seqs = get_longest_sorted_reads(seqs, config.max_n_read, config.max_cov_aln)
	}
	cseqs := copy_seq_ptrs(seqs)
	consensus_data_ptr := C.generate_consensus(&cseqs[0],
		C.uint(len(seqs)),
		C.uint(config.min_cov),
		C.uint(config.K),
		C.double(config.min_idt))
	free_seq_ptrs(cseqs)
	consensus := C.GoString(consensus_data_ptr.sequence)
	C.free_consensus_data(consensus_data_ptr)
	return consensus, seed_id
}

type SeqConfig struct {
	K                                              int
	min_cov, max_n_read, edge_tolerance, trim_size int
	min_cov_aln, max_cov_aln                       int
	min_idt                                        float32
}

type SeqDatum struct {
	seqs    []NtSlice
	seed_id string
	config  *SeqConfig
}

func Get_seq_data(config *SeqConfig, min_n_read int, min_len_aln int) <-chan SeqDatum {
	data_ch := make(chan SeqDatum) // xfer whole seq
	go func() {
		get_seq_data(config, min_n_read, min_len_aln, data_ch)
		close(data_ch)
	}()
	return data_ch
}
func get_seq_data(config *SeqConfig, min_n_read int, min_len_aln int, data_ch chan SeqDatum) {
	const max_len = 100000
	seqs := make([]NtSlice, 0)
	seed_id := ""
	seed_len := 0
	read_cov := 0
	read_ids := make(map[string]bool)
	reader := bufio.NewReader(os.Stdin)
	for {
		text, err := reader.ReadString('\n') // err if no trailing newline
		if err == io.EOF {
			break
		} else if err != nil {
			panic(err)
		}
		text = strings.TrimSpace(text)
		if len(text) == 0 {
			continue
		}
		start := text[0]
		if start == '+' {
			//fmt.Fprintf(os.Stderr, "+++ %d\n", seed_len)
			if len(seqs) >= min_n_read && read_cov/seed_len >= config.min_cov_aln {
				seqs = get_longest_sorted_reads(seqs, config.max_n_read, config.max_cov_aln)
				//yield (seqs, seed_id, config)
				datum := SeqDatum{
					seqs:    seqs,
					seed_id: seed_id,
					config:  config,
				}
				data_ch <- datum
			}
			seqs = make([]NtSlice, 0)
			read_ids = make(map[string]bool)
			seed_id = ""
			read_cov = 0
		} else if start == '*' {
			seqs = make([]NtSlice, 0)
			read_ids = make(map[string]bool)
			seed_id = ""
			read_cov = 0
		} else if start == '-' {
			break
		} else {
			if len(text) > max_len {
				text = text[:max_len-1]
			}
			parts := strings.Fields(text)
			if len(parts) != 2 {
				panic(text) // unexpected
			}
			read_id := parts[0]
			seq := NtSlice(parts[1])
			//fmt.Fprintf(os.Stderr, "%s -> %d\n", read_id, len(seq))
			if len(seq) >= min_len_aln {
				if len(seqs) == 0 {
					seqs = append(seqs, seq) //the "seed"
					seed_len = len(seq)
					seed_id = read_id
				}
				if _, ok := read_ids[read_id]; !ok {
					//avoidng using the same read twice. seed is used again here by design
					seqs = append(seqs, seq)
					read_ids[read_id] = true
					read_cov += len(seq)
				}
			}
		}
	}
	Use(seed_id, seed_len, read_cov, read_ids)
	fmt.Fprintln(os.Stderr, "CHRIS")
	fmt.Fprintln(os.Stderr, "len(seqs)", len(seqs))
}

func format_seq(seq string, col int) string {
	lines := make([]string, 0)
	for i := 0; i < len(seq); i += col {
		bound := i + col
		if bound > len(seq) {
			bound = len(seq)
		}
		lines = append(lines, seq[i:bound])
	}
	return strings.Join(lines, "\n")
}
func findall_good_regions(cns string) []string {
	return regexp.MustCompile("[ACGT]+").FindAllString(cns, -1)
}

type ByShortestString []string

func (slice ByShortestString) Len() int {
	return len(slice)
}
func (slice ByShortestString) Less(i, j int) bool {
	return len(slice[i]) < len(slice[j])
}
func (slice ByShortestString) Swap(i, j int) {
	slice[i], slice[j] = slice[j], slice[i]
}
func main() {
	C.fputs(C.CString("In main()\n"), C.stderr)
	C.foo()
	C.poo()
	println("hello!")
	var opts struct {
		N_core      int `long:"n_core" default:"24" description:"number of processes used for generating consensus; 0 for main process only"`
		Min_cov     int `long:"min_cov" default:"6" description:"minimum coverage to break the consensus"`
		Min_cov_aln int `long:"min_cov_aln" default:"10" description:"minimum coverage of alignment data; a seed read with less than MIN_COV_ALN average depth of coverage will be completely ignored"`
		Max_cov_aln int `long:"max_cov_aln" default:"0" description:"maximum coverage of alignment data; a seed read with more than MAX_COV_ALN average depth of coverage of the longest alignments will be capped, excess shorter alignments will be ignored"`
		// 0 to emulate previous behavior

		Min_len_aln int `long:"min_len_aln" default:"0" description:"minimum length of a sequence in an alignment to be used in consensus; any shorter sequence will be completely ignored"`
		// 0 to emulate previous behavior

		Min_n_read     int     `long:"min_n_read" default:"10" description:"1 + minimum number of reads used in generating the consensus; a seed read with fewer alignments will be completely ignored"`
		Max_n_read     int     `long:"max_n_read" default:"500" description:"1 + maximum number of reads used in generating the consensus"`
		Trim           bool    `long:"trim" description:"trim the input sequence with k-mer spare dynamic programming to find the mapped range"`
		Output_full    bool    `long:"output_full" description:"output uncorrected regions too"`
		Output_multi   bool    `long:"output_multi" description:"output multi correct regions"`
		Min_idt        float32 `long:"min_idt" default:"0.70" description:"minimum identity of the alignments used for correction"`
		Edge_tolerance int     `long:"edge_tolerance" default:"1000" description:"for trimming, the there is unaligned edge leng > edge_tolerance, ignore the read"`
		Trim_size      int     `long:"trim_size" default:"50" description:"the size for triming both ends from initial sparse aligned region"`
	}
	//args, err := flags.ParseArgs(&opts, os.Args)
	//args, err := flags.Parse(&opts)
	parser := flags.NewParser(&opts, flags.Default)
	parser.Usage = "[OPTIONS]\n\nA simple multi-processor consensus sequence generator."
	args, err := parser.Parse()
	if err != nil {
		//panic(err)
		os.Exit(2)
	}
	fmt.Fprintf(os.Stderr, "%v %v %v\n", opts, args, err)
	//fmt.Println(opts, args, err)
	type GetConsensusFunc func(SeqDatum) (string, string)
	var get_consensus GetConsensusFunc
	if opts.Trim {
		get_consensus = get_consensus_with_trim
	} else {
		get_consensus = get_consensus_without_trim
	}
	K := 8
	Use(get_consensus, K)
	//config = args.min_cov, K, args.max_n_read, args.min_idt, args.edge_tolerance, args.trim_size, args.min_cov_aln, args.max_cov_aln
	config := SeqConfig{
		K:              K,
		min_cov:        opts.Min_cov,
		max_n_read:     opts.Max_n_read,
		min_idt:        opts.Min_idt,
		edge_tolerance: opts.Edge_tolerance,
		trim_size:      opts.Trim_size,
		min_cov_aln:    opts.Min_cov_aln,
		max_cov_aln:    opts.Max_cov_aln,
	}
	fmt.Fprintf(os.Stderr, "%v\n", config)
	//fmt.Fprintf(os.Stderr, "len(data) %d\n", (len(data)))
	//for _, datum := range data {
	done_ch := make(chan bool)
	write_ch := make(chan string)
	go func() {
		for blob := range write_ch {
			fmt.Print(blob)
		}
		done_ch <- true
	}()
	var wg sync.WaitGroup
	for datum := range Get_seq_data(&config, 1, 2) {
		go func() {
			//defer wg.Done()
			defer func() {
				println("Done")
				wg.Done()
			}()
			cns, seed_id := get_consensus(datum)
			//println(len(cns), seed_id)
			if len(cns) < 500 {
				return
			}
			if opts.Output_full {
				result := fmt.Sprintf("> %s_f\n%s\n", seed_id, cns)
				write_ch <- result
				return
			}
			cns_goods := findall_good_regions(cns)
			if len(cns_goods) == 0 {
				return
			}
			if opts.Output_multi {
				seq_i := 0
				for _, cns_seq := range cns_goods {
					if len(cns_seq) < 500 {
						return
					}
					if seq_i >= 10 {
						break
					}
					result := fmt.Sprintf(">prolog/%s%01d/%d_%d\n%s\n",
						seed_id, seq_i, 0, len(cns_seq),
						format_seq(cns_seq, 80),
					)
					write_ch <- result
					seq_i += 1
				}
			} else {
				sort.Sort(ByShortestString(cns_goods))
				result := fmt.Sprintf(">%s\n%s\n", seed_id, cns[len(cns)-1])
				write_ch <- result
			}
		}()
		wg.Add(1)
		println("Added")
	}
	println("Waiting...")
	wg.Wait()
	close(write_ch)
	<-done_ch // Wait for writer to finish.
}
