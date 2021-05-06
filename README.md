# Pichia toolkit

Pichia pastoris is a methylotrophic yeast that is popularly used in large-scale production of proteins. This repo contains designs and design software for a Pichia toolkit based on the pPICZ(alpha) and pGAPZ(alpha) plasmids. Our Pichia pastoris expression and secretion vectors can use the GAP constitutive promoter or AOX1 inducible promoter, have a zeocin resistance marker, and allow for MoClo-compatible proteins to be cloned in a single step. 

In total, this toolkit takes 19,831 base pairs to synthesize.

## Enzymes for expression

Sporenet Labs is building this toolkit for the specific purpose of lowering our enzyme costs. Pichia pastoris's secretion system should allow for simplified purification while increasing our end yield of enzyme, while creating useful parts for our customers to use. We are cloning the following enzymes for expression in Pichia pastoris:

- [Pfu-sso7d](https://patents.google.com/patent/US6627424B1/en) - the predecessor of Phusion. The patent finally expired on this enzyme. Phusion has a ~2.6x increase in fidelity over Pfu-sso7d, with [Pfu-sso7d having a 32x increase in fidelity over Taq](https://barricklab.org/twiki/bin/view/Lab/ProtocolsReagentsPfuSso7d). Polymerase is a large cost of DNA sequencing using Sporenet Lab's sequencing method, so we need to express this. We plan on using a [low cost purification method](https://pubmed.ncbi.nlm.nih.gov/31341574/).
- [CdnDI](https://patents.google.com/patent/US8748146B2/en) - a synthetic TypeIIS homing endonuclease. The patent for this enzyme is still active, but most interestingly they [cancelled claims](https://patents.google.com/patent/US20150024464A1/en) using different homing endonucleases for targeting. If we put this enzyme into production, we'll have to pay for licensing, but we may be able to engineer our own set of synthetic TypeIIS restriction enzymes. CdnDI will make for a good control.
- [p50-T4](https://patents.google.com/patent/US20120214208A1/en) - a better T4 ligase. This enzyme is under an abandoned patent, which is great. An expression plasmid is even available at [Addgene](https://www.addgene.org/87742/). Ligase is a large cost of DNA builds, so we need to express this.
- [T4-PNK](https://www.neb.com/products/m0201-t4-polynucleotide-kinase#Product%20Information) - T4 Polynucleotide Kinase. Phosphorylates DNA and is necessary for a step in Sporenet Lab's sequencing.
- [Taq](https://www.addgene.org/87742/) - Taq polymerase. I used the `pOpen_taq` plasmid + a His tag [since I know a his tag probably works](https://www.addgene.org/166944/). Taq is necessary for a step in Sporenet Lab's sequencing pipeline.
- [BsaI](http://rebase.neb.com/rebase/enz/Eco31I.html) - The most commonly used TypeIIS restriction enzyme. Technically, Eco31I. The folks at OpenBioeconomy lab found that Eco31I, an isoschizomer of BsaI, was unpatented, and they got the sequence for it. I'm attempting something that may not work with this enzyme and BtgZI - we know what overexpression of EcoRI [is toxic in yeast if not secreted](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC397259/) (likely due to diffusion into the nucleus), but will overexpression be toxic if we secrete the protein out of the cell? Only an experiment will tell.
- [BtgZI](https://patents.google.com/patent/US7029900) - an interesting TypeIIS restriction enzyme. It was patented by NEB, but patent has expired. This enzyme is an extremely interesting TypeIIS in that its recognition site is ~10bp from its cut site. If you overlap BtgZI and BsaI, you can even get an {BtgZI recognition}-{BsaI recognition}-{BsaI cut}-{BtgZI cut}, which we call "BtgZI skipping". We use a derivative of this method to be directly clone large genes either directly into expression vectors or into `pOpen_v3`, depending on the enzyme used. Normally, BtgZI is used in combination with BbsI for cloning small genes into `pOpen_v3`

## Cloning strategy 

We have two paths for cloning genes into Pichia pastoris expression vectors. We can either use BsaI and directly clone all fragments together into functional expression vectors, or we can use BtgZI and clone each individual gene and backbone component into `pOpen_v3`. From there, we can use BsaI to clone all fragments together into functional expression vectors. 

Later, customers will be able to order genes for expression in Pichia pastoris, and we will clone these genes into expression vectors using BsaI and prebuilt individual backbone components from `pOpen_v3`. While these parts do not take advantage of the larger Sporenet Labs modular ecosystem, they will work well enough for most customers and for our internal uses.

## Expression strategy

Both the GAP constitutive and AOX1 inducible promoter use the [alpha-no-EAEA](https://doi.org/10.1021/sb500366v) secretion tag to secrete proteins from Pichia pastoris. Following cleavage with Kex2, the full protein without the secretion tag is secreted into the media outside of the Pichia pastoris cells. After centrifugation, you can directly purify proteins from the resulting supernatant. 

