package main

import (
	"fmt"
	"os"
	"strings"
	"sync"
	"time"

	"github.com/koeng101/pichia_toolkit/fragment"
	"github.com/koeng101/poly"
)

func main() {
	start := time.Now()
	finalTable := poly.ReadGbk("data/pichia_chr1.gb").GetOptimizationTable(poly.GetCodonTable(11))
	finalTable = poly.AddCodonTable(finalTable, poly.ReadGbk("data/pichia_chr2.gb").GetOptimizationTable(poly.GetCodonTable(11)))
	finalTable = poly.AddCodonTable(finalTable, poly.ReadGbk("data/pichia_chr3.gb").GetOptimizationTable(poly.GetCodonTable(11)))
	finalTable = poly.AddCodonTable(finalTable, poly.ReadGbk("data/pichia_chr4.gb").GetOptimizationTable(poly.GetCodonTable(11)))

	// Write codon table to data/pichiaTable.json
	poly.WriteCodonJSON(finalTable, "data/pichiaTable.json")
	tableFinish := time.Since(start)
	fmt.Println(fmt.Sprintf("Finished codon table generation in %fs", tableFinish.Seconds()))
	start = time.Now()

	enzymesChan := make(chan poly.Fasta, 100)
	go poly.ReadFASTAConcurrent("data/enzymes.fasta", enzymesChan)
	enzymeFragments := make(map[string][]string)
	var seq string
	var optimizedSeq string
	var err error
	for enzyme := range enzymesChan {
		// Remove spaces (from BtgZI)
		seq = strings.ReplaceAll(enzyme.Sequence, " ", "")

		// Add stop codons
		if string(seq[len(seq)-1]) != "*" {
			seq = seq + "*"
		}

		// Setup synthesis fixer funcs
		var functions []func(string, chan poly.DnaSuggestion, *sync.WaitGroup)
		functions = append(functions, poly.FindTypeIIS)

		// Optimize sequence
		optimizedSeq, err = poly.FixCds(":memory:", poly.Optimize(seq, finalTable), finalTable, functions)
		if err != nil {
			fmt.Println(err)
		}

		// Add stop codon, because for some reason it is removed when you optimize a sequence
		optimizedSeq = optimizedSeq + "TAA"

		for _, cutSite := range []string{"GAAGAC", "GGTCTC", "GCGATG", "CGTCTC", "GCTCTTC", "CACCTGC"} {
			if strings.Contains(optimizedSeq, cutSite) {
				fmt.Println(enzyme.Name + " contains " + cutSite)
			}
			if strings.Contains(poly.ReverseComplement(optimizedSeq), cutSite) {
				fmt.Println(enzyme.Name + " reverse complement contains " + cutSite)
			}
		}

		// Add BsaI cutting site
		optimizedSeq = "GGTCTCTA" + optimizedSeq + "TGAGCTTAGAGACC"

		// Add BtgZI/BbsI cutting site
		optimizedSeq = "GCGATGTTGAAGACCAGGAG" + optimizedSeq + "CGCTTAGTCTTCGCCATCGC"

		// Fragment if necessary. I know that none of these genes are THAT long, and Genparts go to 2000bp. An optimized
		// fragmenter will be required in the future for poly.
		fragments, _, err := fragment.Fragment(optimizedSeq, 770, 800, []string{"GGAG", "CGCT", "TTAG"})
		if err != nil {
			fmt.Println(err)
		}
		M13For := "GTAAAACGACGGCCAGT"
		M13Rev := "GTCATAGCTGTTTCCTG"
		switch {
		case len(fragments) == 1:
			enzymeFragments[enzyme.Name] = []string{M13For + fragments[0] + M13Rev}
		default:
			var newFragment string
			var newFragments []string
			for i, fragment := range fragments {
				switch {
				case i == 0:
					newFragment = M13For + fragment + "GGAGACCCGCCATCGC" + M13Rev
				case i == len(fragment)-1:
					newFragment = M13For + "GCGATGTTGGGTCTCA" + fragment + M13Rev
				default:
					newFragment = M13For + "GCGATGTTGGGTCTCA" + fragment + "GGAGACCCGCCATCGC" + M13Rev
				}
				newFragments = append(newFragments, newFragment)
			}
			enzymeFragments[enzyme.Name] = newFragments
		}
	}
	proteinFinish := time.Since(start)
	fmt.Println(fmt.Sprintf("Finished protein optimization and fixing in %fs", proteinFinish.Seconds()))
	start = time.Now()

	// The base plasmid that we are looking at is the `pPICZ(alpha) A` plasmid and the `pGAPZ(alpha) A` plasmid. These have the general features of the ColE1 origin with zeocin resistance, the AOX1 terminator, a specific expression promoter, and the alpha-factor secretion signal.
	// The TEF1 promoter driving expression of zeocin resistance has BsaI cut sites in the original versions of the plasmid, so I have switched these out with the pTEF1 promoters from the yeast toolkit (https://doi.org/10.1021/sb500366v FigS1). Small modifications are made to the intergenetic regions of this plasmid to remove BbsI and BsaI sites.
	// The construct is split into 4 parts:
	// - Promoter + secretion tag (https://doi.org/10.1021/acssynbio.6b00337 using alphaMF_no_EAEA for clean protein cleavage)
	// - Gene of interest
	// - Terminator + ZeoR
	// - ZeoR + Ori
	// The zeocin resistance marker is split into 2 different parts so only a successful cloning reaction can be selected for. The alphaMF_no_EAEA secretion tag is used so we can have a clean N-terminal tag, which is important for a lot of the proteins we wish to express.
	pGAP_alphaMF_no_EAEA := "GGAGTTTTTGTAGAAATGTCTTGGTGTCCTCGTCCAATCAGGTAGCCATCTCTGAAATATCTGGCTCCGTTGCAACTCCGAACGACCTGCTGGCAACGTAAAATTCTCCGGGGTAAAACTTAAATGTGGAGTAATGGAACCAGAAACGTCTCTTCCCTTCTCTCTCCTTCCACCGCCCGTTACCGTCCCTAGGAAATTTTACTCTGCTGGAGAGCTTCTTCTACGGCCCCCTTGCAGCAATGCTCTTCCCAGCATTACGTTGCGGGTAAAACGGAGGTCGTGTACCCGACCTAGCAGCCCAGGGATGGAAAAGTCCCGGCCGTCGCTGGCAATAATAGCGGGCGGACGCATGTCATGAGATTATTGGAAACCACCAGAATCGAATATAAAAGGCGAACACCTTTCCCAATTTTGGTTTCTCCTGACCCAAAGACTTTAAATTTAATTTATTTGTCCCTATTTCAATCAATTGAACAACTATTTCGAAACGATGAGATTTCCTTCAATTTTTACTGCTGTTTTATTCGCAGCATCCTCCGCATTAGCTGCTCCAGTCAACACTACAACAGAAGATGAAACGGCACAAATTCCGGCTGAAGCTGTCATCGGTTACTCAGATTTAGAAGGGGATTTCGATGTTGCTGTTTTGCCATTTTCCAACAGCACAAATAACGGGTTATTGTTTATAAATACTACTATTGCCAGCATTGCTGCTAAAGAAGAAGGGGTATCTCTCGAGAAAAGAATG"
	pAOX1_alphaMF_no_EAEA_1 := "GGAGGATCTAACATCCAAAGACGAAAGGTTGAATGAAACCTTTTTGCCATCCGACATCCACAGGTCCATTCTCACACATAAGTGCCAAACGCAACAGGAGGGGATACACTAGCAGCAGACCGTTGCAAACGCAGGACCTCCACTCCTCTTCTCCTCAACACCCACTTTTGCCATCGAAAAACCAGCCCAGTTATTGGGCTTGATTGGAGCTCGCTCATTCCAATTCCTTCTATTAGGCTACTAACACCATGACTTTATTAGCCTGTCTATCCTGGCCCCCCTGGCGAGGTTCATGTTTGTTTATTTCCGAATGCAACAAGCTCCGCATTACACCCGAACATCACTCCAGATGAGGGCTTTCTGAGTGTGGGGTCAAATAGTTTCATGTTCCCCAAATGGCCCAAAACTGACAGTTTAAACGCTGTCTTGGAACCTAATATGACAAAAGCGTGATCTCATCCAAGATGAACTAAGTTTGGTTCGTTGAAATGCTAACGGCCAGTTGGTCAAAAAGAAACTTCCAAAAGTCGGCATACCGTTTGTCTTGTTTGGTATTGATTGACGAATGCTCAAAAATAATCTCATTAATGCTTAG"
	pAOX1_alphaMF_no_EAEA_2 := "TTAGCGCAGTCTCTCTATCGCTTCTGAACCCCGGTGCACCTGTGCCGAAACGCAAATGGGGAAACACCCGCTTTTTGGATGATTATGCATTGTCTCCACATTGTATGCTTCCAAGATTCTGGTGGGAATACTGCTGATAGCCTAACGTTCATGATCAAAATTTAACTGTTCTAACCCCTACTTGACAGCAATATATAAACAGAAGGAAGCTGCCCTGTCTTAAACCTTTTTTTTTATCATCATTATTAGCTTACTTTCATAATTGCGACTGGTTCCAATTGACAAGCTTTTGATTTTAACGACTTTTAACGACAACTTGAGAAGATCAAAAAACAACTAATTATTCGAAACGATGAGATTTCCTTCAATTTTTACTGCTGTTTTATTCGCAGCATCCTCCGCATTAGCTGCTCCAGTCAACACTACAACAGAAGATGAAACGGCACAAATTCCGGCTGAAGCTGTCATCGGTTACTCAGATTTAGAAGGGGATTTCGATGTTGCTGTTTTGCCATTTTCCAACAGCACAAATAACGGGTTATTGTTTATAAATACTACTATTGCCAGCATTGCTGCTAAAGAAGAAGGGGTATCTCTCGAGAAAAGAATG"
	term_zeocinR := "GCTTGTCTTGCTAGATTCTAATCAAGAGGATGTCAGAATGCCATTTGCCTGAGAGATGCAGGCTTCATTTTTGATACTTTTTTATTTGTAACCTATATAGTATAGGATTTTTTTTGTCATTTTGTTTCTTCTCGTACGAGCTTGCTCCTGATCAGCCTATCTCGCAGCTGATGAATATCTTGTGGTAGGGGTTTGGGAAAATCATTCGAGTTTGATGTTTTTCTTGGTATTTCCCACTCCTCTTCAGAGTACAGAAGATTAAGTGAGACGATCCCCCACACACCATAGCTTCAAAATGTTTCTACTCCTTTTTTACTCTTCCAGATTTTCTCGGACTCCGCGCATCGCCGTACCACTTCAAAACACCCAAGCACAGCATACTAAATTTCCCCTCTTTCTTCCTCTAGGGTGTCGTTAATTACCCGTACTAAAGGTTTGGAAAAGAAAAAAGACACCGCCTCGTTTCTTTTTCTTCGTCGAAAAAGGCAATAAAAATTTTTATCACGTTTCTTTTTCTTGAAAATTTTT"
	term_zeocinR_2 := "TTTTTTGATTTTTTTCTCTTTCGATGACCTCCCATTGATATTTAAGTTAATAAACGGTCATCAATTTCTCAAGTTTCAGTTTCATTTTTCTTGTTCTATTACAACTTTTTTTACTTCTTGCTCATTAGAAAGAAAGCATAGCAATCTAATCTAAGGGCGGTGTTGACAATTAATCATCGGCATAGTATATCGGCATAGTATAATACGACAAGGTGAGGAACTAAACCATGGCCAAGTTGACCAGTGCCGTTCCAGTGCTCACCGCTAGAGACGTCGCCGGAGCTGTCGAGTTCTGGACCGACCGGCTCGGGTTCTCCCGGGACTTCGTGGAGGACGACTTCGCCGGTGTGGTCCGGGACGACGTGACCCTGTTCATCAGCGCT"
	zeocinR_colE1 := "CGCTGTCCAGGACCAGGTGGTGCCGGACAACACCCTGGCTTGGGTTTGGGTAAGAGGTCTGGACGAGCTGTACGCTGAGTGGTCGGAGGTCGTGTCCACGAACTTCAGAGACGCTTCCGGTCCAGCCATGACCGAGATCGGAGAGCAGCCATGGGGTAGAGAGTTCGCCCTGAGAGATCCAGCTGGTAACTGTGTTCACTTCGTGGCCGAGGAGCAGGACTGACACGTCCCCTTTTCCTTTGTCGATATCATGTAATTAGTTATGTCACGCTTACATTCACGCCCTCCCCCCACATCCGCTCTAACCGAAAAGGAAGGAGTTAGACAACCTGAAGTCTAGGTCCCTATTTATTTTTTTATAGTTATGTTAGTATTAAG"
	zeocinR_colE1_2 := "TAAGAACGTTATTTATATTTCAAATTTTTCTTTTTTTTCTGTACAGACGCGTGTACGCATGTAACATTATACTGAAAACCTTGCTTGAGAAGGTTTTGGGACGCTCGAAGGCTTTAATTTGCAAGCTGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCAATGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATCGGAGGGAG"
	ccdB := "GTAAAACGACGGCCAGTGGTCTCTAATGCGAGACCGGACTAAAAGCCAGATAACAGTATGCGTATTTGCGCGCTGATTTTTGCGGTATAAGAATATATACTGATATGTATACCCGAAGTATGTCAAAAAGAGGTATGCTATGAAGCAGCGTATTACAGTGACAGTTGACAGCGACAGCTATCAGTTGCTCAAGGCATATATGATGTCAATATCTCCGGTCTGGTAAGCACAACCATGCAGAATGAAGCCCGTCGTCTGCGTGCCGAACGCTGGAAAGCGGAAAATCAGGAAGGGATGGCTGAGGTCGCCCGGTTTATTGAAATGAACGGCTCTTTTGCTGACGAGAACAGGGGCTGGTGAAATGCAGTTTAAGGTTTACACCTATAAAAGAGAGAGCCGTTATCGTCTGTTTGTGGATGTACAGAGTGATATTATTGACACGCCCGGGCGACGGATGGTGATCCCCCTGGCCAGTGCACGTCTGCTGTCAGATAAAGTCTCCCGTGAACTTTACCCGGTGGTGCATATCGGGGATGAAAGCTGGCGCATGATGACCACCGATATGGCCAGTGTGCCGGTTTCCGTTATCGGGGAAGAAGTGGCTGATCTCAGCCACCGCGAAAATGACATCAAAAACGCCATTAACCTGATGTTCTGGGGAATATAAATGTCAGGCTCCCTTATACACAGCCGGTCTCTGCTTCGAGACCGTCATAGCTGTTTCCTG"

	// Add M13 primers + BbsI and BsaI cut sites
	pGAP_alphaMF_no_EAEA = "GTAAAACGACGGCCAGT" + "GCGATGTTGAAGACCAGGAGGTCTCT" + pGAP_alphaMF_no_EAEA + "AGAGACCCGCTTAGTCTTCGCCATCGC" + "GTCATAGCTGTTTCCTG"
	pAOX1_alphaMF_no_EAEA_1 = "GTAAAACGACGGCCAGT" + "GCGATGTTGAAGACCAGGAGGTCTCT" + pAOX1_alphaMF_no_EAEA_1 + "AGAGACCCGCTTAGTCTTCGCCATCGC" + "GTCATAGCTGTTTCCTG"
	pAOX1_alphaMF_no_EAEA_2 = "GTAAAACGACGGCCAGT" + "GCGATGTTGAAGACCAGGAGGTCTCT" + pAOX1_alphaMF_no_EAEA_2 + "AGAGACCCGCTTAGTCTTCGCCATCGC" + "GTCATAGCTGTTTCCTG"
	term_zeocinR = "GTAAAACGACGGCCAGT" + "GCGATGTTGAAGACCAGGAGGTCTCT" + term_zeocinR + "AGAGACCCGCTTAGTCTTCGCCATCGC" + "GTCATAGCTGTTTCCTG"
	term_zeocinR_2 = "GTAAAACGACGGCCAGT" + "GCGATGTTGAAGACCAGGAGGTCTCC" + term_zeocinR_2 + "AGAGACCCGCTTAGTCTTCGCCATCGC" + "GTCATAGCTGTTTCCTG"
	zeocinR_colE1 = "GTAAAACGACGGCCAGT" + "GCGATGTTGAAGACCAGGAGGTCTCT" + zeocinR_colE1 + "AGAGACCCGCTTAGTCTTCGCCATCGC" + "GTCATAGCTGTTTCCTG"
	zeocinR_colE1_2 = "GTAAAACGACGGCCAGT" + "GCGATGTTGAAGACCAGGAGGTCTCC" + zeocinR_colE1_2 + "AGAGACCCGCTTAGTCTTCGCCATCGC" + "GTCATAGCTGTTTCCTG"

	// Get pOpen_v3
	pOpen_v3 := poly.CloneSequence{Sequence: strings.ToUpper("TAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGGCCTACTATTAGCAACAACGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGAACCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACCTGCACCAGTCAGTAAAACGACGGCCAGTAGTCAAAAGCCTCCGACCGGAGGCTTTTGACTTGGTTCAGGTGGAGTGGGAGTAgtcttcGCcatcgCtACTAAAagccagataacagtatgcgtatttgcgcgctgatttttgcggtataagaatatatactgatatgtatacccgaagtatgtcaaaaagaggtatgctatgaagcagcgtattacagtgacagttgacagcgacagctatcagttgctcaaggcatatatgatgtcaatatctccggtctggtaagcacaaccatgcagaatgaagcccgtcgtctgcgtgccgaacgctggaaagcggaaaatcaggaagggatggctgaggtcgcccggtttattgaaatgaacggctcttttgctgacgagaacagggGCTGGTGAAATGCAGTTTAAGGTTTACACCTATAAAAGAGAGAGCCGTTATCGTCTGTTTGTGGATGTACAGAGTGATATTATTGACACGCCCGGGCGACGGATGGTGATCCCCCTGGCCAGTGCACGTCTGCTGTCAGATAAAGTCTCCCGTGAACTTTACCCGGTGGTGCATATCGGGGATGAAAGCTGGCGCATGATGACCACCGATATGGCCAGTGTGCCGGTCTCCGTTATCGGGGAAGAAGTGGCTGATCTCAGCCACCGCGAAAATGACATCAAAAACGCCATTAACCTGATGTTCTGGGGAATATAAATGTCAGGCTCCCTTATACACAGgcgatgttgaagaccaCGCTGAGGTGTCAATCGTCGGAGCCGCTGAGCAATAACTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTGCATGGTCATAGCTGTTTCCTGAGAGCTTGGCAGGTGATGACACACATTAACAAATTTCGTGAGGAGTCTCCAGAAGAATGCCATTAATTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGG"), Circular: true}

	// Simulate cloning of every gene into pOpen_v3
	//fragmentsToSimulate := append([]string{pGAP_alphaMF_no_EAEA, pAOX1_alphaMF_no_EAEA, term_zeocinR, zeocinR_colE1}, enzymeFrags...)
	for enzymeName, fragment := range enzymeFragments {
		var cloneSeqs []poly.CloneSequence
		for _, enzymeSeq := range fragment {
			cloneSeqs = append(cloneSeqs, poly.CloneSequence{Sequence: enzymeSeq, Circular: false})
		}
		simulationFragments := append(cloneSeqs, pOpen_v3)
		clonedFragments, err := poly.GoldenGate(simulationFragments, "BtgZI")
		if err != nil {
			fmt.Println(err)
		}
		if len(clonedFragments) != 1 {
			fmt.Println("!More than 1 cloned sequence simulated in basic genes")
			fmt.Println(enzymeName)
			fmt.Println(fragment)
		}
	}

	// Simulate cloning of backbone components into pOpen_v3
	for _, backbonePart := range []string{pGAP_alphaMF_no_EAEA, pAOX1_alphaMF_no_EAEA_1, pAOX1_alphaMF_no_EAEA_2, term_zeocinR, term_zeocinR_2, zeocinR_colE1, zeocinR_colE1_2} {
		clonedFragments, err := poly.GoldenGate([]poly.CloneSequence{poly.CloneSequence{Sequence: backbonePart, Circular: false}, pOpen_v3}, "BbsI")
		if err != nil {
			fmt.Println(err)
		}
		if len(clonedFragments) != 1 {
			fmt.Println("!More than 1 cloned sequence simulated in backbone parts")
			for _, clonedPlasmid := range clonedFragments {
				fmt.Println(clonedPlasmid)
			}
		}
	}

	// Simulate the cloning of every gene into a pichia plasmid
	var c strings.Builder
	backbone := []poly.CloneSequence{poly.CloneSequence{Sequence: term_zeocinR, Circular: false}, poly.CloneSequence{Sequence: term_zeocinR_2, Circular: false}, poly.CloneSequence{Sequence: zeocinR_colE1, Circular: false}, poly.CloneSequence{Sequence: zeocinR_colE1_2, Circular: false}}
	for enzymeName, enzymeSeqs := range enzymeFragments {
		var cloneSeqs []poly.CloneSequence
		for _, enzymeSeq := range enzymeSeqs {
			cloneSeqs = append(cloneSeqs, poly.CloneSequence{Sequence: enzymeSeq, Circular: false})
		}

		promoters := []string{"GAP", "AOX1"}
		for i, promoter := range promoters { //[]string{pGAP_alphaMF_no_EAEA, pAOX1_alphaMF_no_EAEA_1, pAOX1_alphaMF_no_EAEA_2} {
			var promoterCloneSeqs []poly.CloneSequence
			promoterCloneSeqs = append(backbone, cloneSeqs...)
			if promoter == "GAP" {
				promoterCloneSeqs = append(promoterCloneSeqs, poly.CloneSequence{Sequence: pGAP_alphaMF_no_EAEA, Circular: false})
			} else {
				promoterCloneSeqs = append(promoterCloneSeqs, poly.CloneSequence{Sequence: pAOX1_alphaMF_no_EAEA_1, Circular: false})
				promoterCloneSeqs = append(promoterCloneSeqs, poly.CloneSequence{Sequence: pAOX1_alphaMF_no_EAEA_2, Circular: false})
			}
			clonedFragments, err := poly.GoldenGate(promoterCloneSeqs, "BsaI")
			if err != nil {
				fmt.Println(err)
			}
			if len(clonedFragments) != 1 {
				fmt.Println("!More than 1 cloned sequence simulated")
				fmt.Println(enzymeName)
				fmt.Println(promoterCloneSeqs)
				for _, clonedPlasmid := range clonedFragments {
					fmt.Println(clonedPlasmid)
				}
			} else {
				fmt.Fprintf(&c, ">%s-%s\n%s\n", promoters[i], enzymeName, clonedFragments[0].Sequence)
			}
		}
	}
	e, err := os.Create("data/plasmids.fasta")
	defer e.Close()
	e.WriteString(c.String())

	cloningFinish := time.Since(start)
	fmt.Println(fmt.Sprintf("Finished cloning simulation in %fs", cloningFinish.Seconds()))
	start = time.Now()

	// Write fragments to a string
	var b strings.Builder
	var bases int
	for enzymeName, enzymeSeqs := range enzymeFragments {
		for i, enzymeSeq := range enzymeSeqs {
			bases = bases + len(enzymeSeq)
			fmt.Fprintf(&b, "%s-%d,%s\n", enzymeName, i, enzymeSeq)
		}
	}
	fmt.Fprintf(&b, "%s,%s\n", "term_zeocinR", term_zeocinR)
	fmt.Fprintf(&b, "%s,%s\n", "zeocinR_colE1", zeocinR_colE1)
	fmt.Fprintf(&b, "%s,%s\n", "zeocinR_colE1_2", zeocinR_colE1_2)
	fmt.Fprintf(&b, "%s,%s\n", "pGAP_alphaMF_no_EAEA", pGAP_alphaMF_no_EAEA)
	fmt.Fprintf(&b, "%s,%s\n", "pAOX1_alphaMF_no_EAEA_1", pAOX1_alphaMF_no_EAEA_1)
	fmt.Fprintf(&b, "%s,%s\n", "pAOX1_alphaMF_no_EAEA_2", pAOX1_alphaMF_no_EAEA_2)
	fmt.Fprintf(&b, "%s,%s\n", "term_zeocinR_2", term_zeocinR_2)
	fmt.Fprintf(&b, "%s,%s\n", "ccdB", ccdB)
	bases = bases + len(term_zeocinR)
	bases = bases + len(zeocinR_colE1)
	bases = bases + len(pGAP_alphaMF_no_EAEA)
	bases = bases + len(pAOX1_alphaMF_no_EAEA_1)
	bases = bases + len(pAOX1_alphaMF_no_EAEA_2)
	bases = bases + len(term_zeocinR_2)

	f, err := os.Create("data/output.csv")
	defer f.Close()
	f.WriteString(b.String())

	fmt.Println(fmt.Sprintf("Synthesizing %d bases", bases))
}
