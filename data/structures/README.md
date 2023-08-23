# Visualizing SARS-CoV-2 Variant Fitness

This directory contains custom structures (not hosted on RCSB PDB) for visualizing the [fitness effects of mutations in every SARS-CoV-2 protein](https://github.com/jbloomlab/SARS2-mut-fitness) with [`dms-viz`](https://dms-viz.github.io/). Below are some notes on the structures we picked for each SARS-CoV-2 protein and proteins that we couldn't find reasonable structures for. 

The three stuctures included in this repository are: 
- [E](#e)
- [ORF6](#orf6)
- [ORF7b](#orf7b)

For each protein, there might be more applicable structures than the structures we chose. If you think this is the case, please raise an [issue](https://github.com/jbloomlab/SARS2-mut-fitness/issues) with your suggested alternative and a justification. 

## Structural Proteins

### S

Surface Glycoprotein. The structure is from from [this paper](https://doi.org/10.1016/j.cell.2020.02.058) from the Veesler lab that has Spike with one RBD in the up comfirmation.

Structure: [6VYB](https://www.rcsb.org/structure/6VYB)

### E

Envelope Protein. This is a structural protein that's part of the capsid and probably acts as a viroporin. The structure is likely pentameric. There aren't any solved structures of the whole thing. To make a structure for this protein, I took the computationally determined monomer from [this paper](https://zenodo.org/record/5521766) and aligned it to the NMR assmebly of the E protein from SARS-CoV-1 ([5X29](https://www.rcsb.org/structure/5X29)) to make a pentamer. I've saved the pentamer [locally](./E.pdb).

### M

Membrane protein. This forms a dimer. There is a Cryo-EM structure of the dimer from [this paper](https://doi.org/10.7554/eLife.81702): [8CTK](https://www.rcsb.org/structure/8CTK). Both chains (A, B) are the monomers.

### N

Nucleoprotein. This protein binds to the RNA in the virion. There is a fairly complete EM structure on the PDB from [this paper](https://doi.org/10.1093/micmic/ozac036): [8FD5](https://www.rcsb.org/structure/8FD5)

## Accessory Proteins

### ORF3a

A membrane protein that isn't a viroporin. It's probably involved in some endosomal pathway. There is [a great paper](https://elifesciences.org/articles/84477) characterizing the structure with Cryo-EM.

The best structure for this protein is here: [8EQJ](https://www.rcsb.org/structure/8EQJ)

### ORF6

This is a strange protein that isn't very large. It looks like one large helix. I've included a computationally modeled structure from [this paper](https://onlinelibrary.wiley.com/doi/10.1002/prot.26250). It seems like the protein binds to [STAT-1 and inhibits it's normal localization](https://www.nature.com/articles/s42003-022-03427-4). The structure is saved [locally](./ORF6.pdb). 

### ORF7a

This is an immunomodulatory protein. There is a crystal structure of the ectodomain of this protein from [this paper](https://doi.org/10.1016/j.isci.2021.102187): [7CI3](https://www.rcsb.org/structure/7CI3)

### ORF7b

This is a transmembrane protein. It's not super clear what it does. I don't know if there is really a structure either. I included a structure from RosettaTTAfold. I got the model from [here](https://www.ebi.ac.uk/interpro/entry/InterPro/IPR021532/). The structure is saved [locally](./ORF7b.pdb). 

### ORF8

It seems like this protein probably does something to the [ER stress response](https://journals.asm.org/doi/10.1128/jvi.00011-23). There is a crystal structure of this protein here: [7JX6](https://www.rcsb.org/structure/7JX6).

### ORF9b

This protein probably interacts with mitochondria in some way. There is a crystal structure of this protein  here: [6Z4U](https://www.rcsb.org/structure/6Z4U).

### ORF10

This protein is likely not even expressed.

## Non-structural Proteins

### nsp1

This protein seems to interfere with translation by binding to the ribosomal exit channel. There is a crystal structure of this protein that was characterized in [this paper](https://doi.org/10.1016/j.isci.2020.101903), and it can be found here: [7K3N](https://www.rcsb.org/structure/7k3n).

### nsp2

There is probably some role for this protein in modulating translation of viral proteins. There is full length structure of this protein determined with AlphaFold and Cyro-EM determined in [this paper](https://pubmed.ncbi.nlm.nih.gov/34013269/): [7MSX](https://www.rcsb.org/structure/7MSX)

### nsp3

Papin-like protease (PLpro). This is an active drug target for SARS-CoV-2. We might only really care about the PLpro domain of this protein. I don't think there is a full structure of this protein that includes all of the domains. This could be a decent structure for that, although there are a ton of other examples: [6WUU](https://www.rcsb.org/structure/6WUU).

### nsp4

This protein isn't well studied. It probably has some role in altering the behaviour of the ER. There isn't a strucutre. It's possible to model the structure with AlphaFold like they did in [this paper](https://www.nature.com/articles/s41392-022-00884-5#:~:text=Nsp2%20protein&text=SARS-CoV-2%20Nsp2%20comprises,4b). I didn't include this structure in the visualization. 

### nsp5

This is the SARS-CoV-2 Main Protease (3CLpro). There are hundreds of crystal structures of this homodimer. Here is the one I'm using right now: [6LU7](https://www.rcsb.org/structure/6lu7).

### nsp6

This protein probably has something to do with setting up replication centers for SARS-2 virions. There doesn't appear to be a structure for this protein. It might be possible to use AlphaFold like they did in [this paper](https://www.nature.com/articles/s41392-022-00884-5#:~:text=Nsp2%20protein&text=SARS-CoV-2%20Nsp2%20comprises,4b). I didn't include this structure in the visualization. 

### nsp7

Part of the RdRP polymerase complex. This is a co-factor of the polymerase complex. There is a structure of the whole complex. Nsp7 is specifically chain **C**: [6YYT](https://www.rcsb.org/structure/6YYT)

### nsp8

Part of the RdRP polymerase complex. This is a co-factor of the polymerase complex. There is a structure of the whole complex. Nsp8 are chains **B and D**: [6YYT](https://www.rcsb.org/structure/6YYT)

### nsp9

Is a DNA/RNA binding protein that forms part of the Replication/Transcription complex. It's a homodimer. There is a crystal structure here: [6WXD](https://www.rcsb.org/structure/6WXD).

### nsp10

Involved in methylation. This protein interacts with both nsp14 and nsp16. There is a structure in complex with both of these and a single crystal structure: [6ZPE](https://www.rcsb.org/structure/6ZPE)

### nsp12

Part of the RdRP polymerase complex. This is the main RdRp. There is a structure of the whole complex. Nsp12 is chain **A**: [6YYT](https://www.rcsb.org/structure/6YYT)

### nsp13

This protein has helicase activtiy. There is a structure of this protein here: [7NNG](https://www.rcsb.org/structure/7NNG). I think there are also some structures of this protein with other replicase proteins.

### nsp14

This protein associates with nsp10 and has some proofreading activity. There are multiple structures for this protein. Here is one in complex with nsp10 and RNA: [7N0C](https://www.rcsb.org/structure/7N0C).

### nsp15

This protein cleaves the 5′-polyuridine tracts in negative-strand RNA and prevents the activation of host pattern recognition receptor MDA5-mediated immune response. The structure is a homo-hexamer with D3 symmetry. Here is the structure: [6WXC](https://www.rcsb.org/structure/6wxc).

### nsp16

This protein is a 2'-O-methyltransferase that associates with nsp10. Here is a structure in complex with nsp10 where nsp16 is chain _A_: [6WVN](https://www.rcsb.org/structure/6WVN).

## Method

I relied on a few hand resources to pick a representative structure for each protein. [This feature from the PDB](https://www.rcsb.org/news/feature/5e74d55d2d410731e9944f52) is a really handy resource for the currently avaliable SARS-CoV-2 structures. This is also a very [helpful paper](https://www.nature.com/articles/s41392-022-00884-5#:~:text=Nsp2%20protein&text=SARS%2DCoV%2D2%20Nsp2%20comprises,4b) for finding structures. There are some structures (both experimentally and computationally determined) [here](https://zenodo.org/record/5521766).
