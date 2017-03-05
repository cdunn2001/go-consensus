package main

//import "flag"
import "fmt"
import "github.com/jessevdk/go-flags"
import "os"

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
	fmt.Println(args, err)
}
