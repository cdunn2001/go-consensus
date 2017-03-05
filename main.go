package main

//import "flag"
import "fmt"
import "github.com/jessevdk/go-flags"
import "os"
import "sort"

func Use(vals ...interface{}) {
	for _, val := range vals {
		_ = val
	}
}

type NtSlice []byte

func min(a, b int) int {
	if a < b {
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
func get_consensus_with_trim(c_input int) (string, int) {
	return "hi", 0
}
func get_consensus_without_trim(c_input int) (string, int) {
	return "bye", 1
}
func main() {
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
	fmt.Println(opts, args, err)
	type GetConsensusFunc func(int) (string, int)
	var get_consensus GetConsensusFunc
	if opts.Trim {
		get_consensus = get_consensus_with_trim
	} else {
		get_consensus = get_consensus_without_trim
	}
	K := 8
	Use(get_consensus, K)
	//config = args.min_cov, K, args.max_n_read, args.min_idt, args.edge_tolerance, args.trim_size, args.min_cov_aln, args.max_cov_aln
	type SeqData struct {
		K                                              int
		min_cov, max_n_read, edge_tolerance, trim_size int
		min_cov_aln, max_cov_aln                       int
		min_idt                                        float32
	}
	config := SeqData{
		K:              K,
		min_cov:        opts.Min_cov,
		max_n_read:     opts.Max_n_read,
		min_idt:        opts.Min_idt,
		edge_tolerance: opts.Edge_tolerance,
		trim_size:      opts.Trim_size,
		min_cov_aln:    opts.Min_cov_aln,
		max_cov_aln:    opts.Max_cov_aln,
	}
	fmt.Println(config)
	Use(config)
}
