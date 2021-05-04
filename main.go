package main

import (
	"fmt"
	"github.com/koeng101/poly"
	"os"
	"strings"
	"sync"
	"time"
)

func main() {
	// There are 4 chromosomes of Pichia. We're going to add them all together into a final codon table
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
		if len(optimizedSeq) > 1950 {
			enzymeFragments[enzyme.Name] = []string{"GTAAAACGACGGCCAGT" + optimizedSeq[:len(optimizedSeq)/2] + "GGAGACCCGCCATCGC" + "GTCATAGCTGTTTCCTG", "GTAAAACGACGGCCAGT" + "GCGATGTTGGGTCTCA" + optimizedSeq[(len(optimizedSeq)/2)-4:] + "GTCATAGCTGTTTCCTG"}
		} else {
			enzymeFragments[enzyme.Name] = []string{"GTAAAACGACGGCCAGT" + optimizedSeq + "GTCATAGCTGTTTCCTG"}
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
	pAOX1_alphaMF_no_EAEA := "GGAGGATCTAACATCCAAAGACGAAAGGTTGAATGAAACCTTTTTGCCATCCGACATCCACAGGTCCATTCTCACACATAAGTGCCAAACGCAACAGGAGGGGATACACTAGCAGCAGACCGTTGCAAACGCAGGACCTCCACTCCTCTTCTCCTCAACACCCACTTTTGCCATCGAAAAACCAGCCCAGTTATTGGGCTTGATTGGAGCTCGCTCATTCCAATTCCTTCTATTAGGCTACTAACACCATGACTTTATTAGCCTGTCTATCCTGGCCCCCCTGGCGAGGTTCATGTTTGTTTATTTCCGAATGCAACAAGCTCCGCATTACACCCGAACATCACTCCAGATGAGGGCTTTCTGAGTGTGGGGTCAAATAGTTTCATGTTCCCCAAATGGCCCAAAACTGACAGTTTAAACGCTGTCTTGGAACCTAATATGACAAAAGCGTGATCTCATCCAAGATGAACTAAGTTTGGTTCGTTGAAATGCTAACGGCCAGTTGGTCAAAAAGAAACTTCCAAAAGTCGGCATACCGTTTGTCTTGTTTGGTATTGATTGACGAATGCTCAAAAATAATCTCATTAATGCTTAGCGCAGTCTCTCTATCGCTTCTGAACCCCGGTGCACCTGTGCCGAAACGCAAATGGGGAAACACCCGCTTTTTGGATGATTATGCATTGTCTCCACATTGTATGCTTCCAAGATTCTGGTGGGAATACTGCTGATAGCCTAACGTTCATGATCAAAATTTAACTGTTCTAACCCCTACTTGACAGCAATATATAAACAGAAGGAAGCTGCCCTGTCTTAAACCTTTTTTTTTATCATCATTATTAGCTTACTTTCATAATTGCGACTGGTTCCAATTGACAAGCTTTTGATTTTAACGACTTTTAACGACAACTTGAGAAGATCAAAAAACAACTAATTATTCGAAACGATGAGATTTCCTTCAATTTTTACTGCTGTTTTATTCGCAGCATCCTCCGCATTAGCTGCTCCAGTCAACACTACAACAGAAGATGAAACGGCACAAATTCCGGCTGAAGCTGTCATCGGTTACTCAGATTTAGAAGGGGATTTCGATGTTGCTGTTTTGCCATTTTCCAACAGCACAAATAACGGGTTATTGTTTATAAATACTACTATTGCCAGCATTGCTGCTAAAGAAGAAGGGGTATCTCTCGAGAAAAGAATG"
	term_zeocinR := "GCTTGTCTTGCTAGATTCTAATCAAGAGGATGTCAGAATGCCATTTGCCTGAGAGATGCAGGCTTCATTTTTGATACTTTTTTATTTGTAACCTATATAGTATAGGATTTTTTTTGTCATTTTGTTTCTTCTCGTACGAGCTTGCTCCTGATCAGCCTATCTCGCAGCTGATGAATATCTTGTGGTAGGGGTTTGGGAAAATCATTCGAGTTTGATGTTTTTCTTGGTATTTCCCACTCCTCTTCAGAGTACAGAAGATTAAGTGAGACTTGCCAACAGGGAGTTCTTCAGAGACATGGAGGCTCAAAACGAAATTATTGACAGCCTAGACATCAATAGTCATACAACAGAAAGCGACCACCCAACTTTGGCTGATAATAGCGTATAAACAATGCATACTTTGTACGTTCAAAATACAATGCAGTAGATATATTTATGCATATTACATATAATACATATCACATAGGAAGCAACAGGCGCGTTGGACTTTTAATTTTCGAGGACCGCGAATCCTTACATCACACCCAATCCCCCACAAGTGATCCCCCACACACCATAGCTTCAAAATGTTTCTACTCCTTTTTTACTCTTCCAGATTTTCTCGGACTCCGCGCATCGCCGTACCACTTCAAAACACCCAAGCACAGCATACTAAATTTCCCCTCTTTCTTCCTCTAGGGTGTCGTTAATTACCCGTACTAAAGGTTTGGAAAAGAAAAAAGACACCGCCTCGTTTCTTTTTCTTCGTCGAAAAAGGCAATAAAAATTTTTATCACGTTTCTTTTTCTTGAAAATTTTTTTTTTTGATTTTTTTCTCTTTCGATGACCTCCCATTGATATTTAAGTTAATAAACGGTCATCAATTTCTCAAGTTTCAGTTTCATTTTTCTTGTTCTATTACAACTTTTTTTACTTCTTGCTCATTAGAAAGAAAGCATAGCAATCTAATCTAAGGGCGGTGTTGACAATTAATCATCGGCATAGTATATCGGCATAGTATAATACGACAAGGTGAGGAACTAAACCATGGCCAAGTTGACCAGTGCCGTTCCGGTGCTCACCGCGCGCGACGTCGCCGGAGCGGTCGAGTTCTGGACCGACCGGCTCGGGTTCTCCCGGGACTTCGTGGAGGACGACTTCGCCGGTGTGGTCCGGGACGACGTGACCCTGTTCATCAGCGCGCGCT"
	zeocinR_colE1 := "CGCTGTCCAGGACCAGGTGGTGCCGGACAACACCCTGGCCTGGGTGTGGGTGCGCGGCCTGGACGAGCTGTACGCCGAGTGGTCGGAGGTCGTGTCCACGAACTTCCGGGACGCCTCCGGGCCGGCCATGACCGAGATCGGCGAGCAGCCGTGGGGGCGGGAGTTCGCCCTGCGCGACCCGGCCGGCAACTGCGTGCACTTCGTGGCCGAGGAGCAGGACTGACACGTCCGACGGCGGCCCACGGGTCCCAGGCCTCGGAGATCCGTCCCCCTTTTCCTTTGTCGATATCATGTAATTAGTTATGTCACGCTTACATTCACGCCCTCCCCCCACATCCGCTCTAACCGAAAAGGAAGGAGTTAGACAACCTGAAGTCTAGGTCCCTATTTATTTTTTTATAGTTATGTTAGTATTAAGAACGTTATTTATATTTCAAATTTTTCTTTTTTTTCTGTACAGACGCGTGTACGCATGTAACATTATACTGAAAACCTTGCTTGAGAAGGTTTTGGGACGCTCGAAGGCTTTAATTTGCAAGCTGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCAATGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGCATGAGATCGGAG"

	// Add M13 primers + BbsI and BsaI cut sites
	pGAP_alphaMF_no_EAEA = "GTAAAACGACGGCCAGT" + "GCGATGTTGAAGACCAGGAGGTCTCT" + pGAP_alphaMF_no_EAEA + "AGAGACCCGCTTAGTCTTCGCCATCGC" + "GTCATAGCTGTTTCCTG"
	pAOX1_alphaMF_no_EAEA = "GTAAAACGACGGCCAGT" + "GCGATGTTGAAGACCAGGAGGTCTCT" + pAOX1_alphaMF_no_EAEA + "AGAGACCCGCTTAGTCTTCGCCATCGC" + "GTCATAGCTGTTTCCTG"
	term_zeocinR = "GTAAAACGACGGCCAGT" + "GCGATGTTGAAGACCAGGAGGTCTCT" + term_zeocinR + "AGAGACCCGCTTAGTCTTCGCCATCGC" + "GTCATAGCTGTTTCCTG"
	zeocinR_colE1 = "GTAAAACGACGGCCAGT" + "GCGATGTTGAAGACCAGGAGGTCTCT" + zeocinR_colE1 + "AGAGACCCGCTTAGTCTTCGCCATCGC" + "GTCATAGCTGTTTCCTG"

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
			for _, f := range fragment {
				fmt.Println(f)
			}
		}
	}

	// Simulate cloning of backbone components into pOpen_v3
	for _, backbonePart := range []string{pGAP_alphaMF_no_EAEA, pAOX1_alphaMF_no_EAEA, term_zeocinR, zeocinR_colE1} {
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
	backbone := []poly.CloneSequence{poly.CloneSequence{Sequence: term_zeocinR, Circular: false}, poly.CloneSequence{Sequence: zeocinR_colE1, Circular: false}}
	for _, enzymeSeqs := range enzymeFragments {
		var cloneSeqs []poly.CloneSequence
		for _, enzymeSeq := range enzymeSeqs {
			cloneSeqs = append(cloneSeqs, poly.CloneSequence{Sequence: enzymeSeq, Circular: false})
		}

		for _, promoter := range []string{pGAP_alphaMF_no_EAEA, pAOX1_alphaMF_no_EAEA} {
			var promoterCloneSeqs []poly.CloneSequence
			promoterCloneSeqs = append(backbone, cloneSeqs...)
			promoterCloneSeqs = append(promoterCloneSeqs, poly.CloneSequence{Sequence: promoter, Circular: false})
			clonedFragments, err := poly.GoldenGate(promoterCloneSeqs, "BsaI")
			if err != nil {
				fmt.Println(err)
			}
			if len(clonedFragments) != 1 {
				fmt.Println("!More than 1 cloned sequence simulated")
				fmt.Println(promoterCloneSeqs)
				for _, clonedPlasmid := range clonedFragments {
					fmt.Println(clonedPlasmid)
				}
			}
		}
	}
	cloningFinish := time.Since(start)
	fmt.Println(fmt.Sprintf("Finished cloning simulation in %fs", cloningFinish.Seconds()))
	start = time.Now()

	// Write fragments to a string
	var b strings.Builder
	for enzymeName, enzymeSeqs := range enzymeFragments {
		for i, enzymeSeq := range enzymeSeqs {
			fmt.Fprintf(&b, ">%s-%d\n%s\n", enzymeName, i, enzymeSeq)
		}
	}
	fmt.Fprintf(&b, ">%s\n%s\n", "term_zeocinR", term_zeocinR)
	fmt.Fprintf(&b, ">%s\n%s\n", "zeocinR_colE1", zeocinR_colE1)
	fmt.Fprintf(&b, ">%s\n%s\n", "pGAP_alphaMF_no_EAEA", pGAP_alphaMF_no_EAEA)
	fmt.Fprintf(&b, ">%s\n%s\n", "pAOX1_alphaMF_no_EAEA", pAOX1_alphaMF_no_EAEA)

	f, err := os.Create("data/output.fasta")
	defer f.Close()
	f.WriteString(b.String())
}