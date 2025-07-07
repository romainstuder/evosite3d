# Tutorial on Ancestral Sequence Reconstruction

This tutorial was part of a course on protein evolution done during ECCB 2014 in Strasbourg.

NB: If there are any question, comments or bugs, feel free to ask. ;) NB2: Revised in 2025

## Introduction

The slides of the introduction are available here: ancestral_sequence_reconstruction.pdf

In this practical, you will learn how to prepare files for CodeML, how to use it for reconstructing
ancestral sequences and how to compute the isoelectric point of protein sequences. There are a few
scripts you will need to download during the practical.

This practical will focus on the lysozyme, an enzyme (EC 3.2.1.17) that damages bacterial cell
walls. The Uniprot page is: <http://www.uniprot.org/uniprot/P61626> They evolved differently in
Primates:

## Part 1: Prepare alignment

The quality of the ancestral reconstruction will heavily depend on the quality of the alignment and
the tree topology (branch lengths are re-estimated during the reconstruction).

<https://www.ensembl.org/Homo_sapiens/Gene/Compara_Ortholog?db=core;g=ENSG00000090382;r=12:69348381-69354234>
Please download the sequence file:

```shell
cp Human_LYZ_orthologues.fa lysozyme_primates.seq
```

Make a multiple alignment, either with Mafft-L-INS-i or Clustal-Omega:

```shell
mafft-linsi lysozyme_primates.seq > lysozyme_primates.fasta
```

or

```shell
clustalo --in lysozyme_primates.seq --out lysozyme_primates.fasta
```

You can have a look at the alignment in Jalview:

```shell
jalview ./lysozyme_primates.fasta
```

The format of the resulting alignment is FASTA. However, most phylogenetic software use PHYLIP
format. So, you have to convert it into PHYLIP. Download the script "convert_fasta2phylip.py" and
execute it:

```shell
convert_fasta2phylip.py lysozyme_primates.fasta lysozyme_primates.phy
```

(If you are not familiar, have a look at the differences between the alignment in FASTA and PHYLIP
formats).

## Part 2: Prepare alignment

We can now generate a tree, either with PhyML (one of the most accurate tool) or FastTree (very fast
and pretty accurate):

```shell
phyml -i lysozyme_primates.phy -d aa -m WAG -c 8 -a e -b 0
mv lysozyme_primates.phy_phyml_tree.txt lysozyme_primates.tree
```

```shell
Option used:
-i = input file
-d aa: amino acid sequences
-m WAG: (substitution matrix). JTT works fine for most proteins, but other matrices (WAG, LG) can do slightly better.
-c 8: (numbers of categories for the gamma distribution)
-a e: (estimate alpha parameter for the gamma distribution)
-b 0: (we don't want boostrap, as this will cause trouble for further analyses in CodeML).
```

or run:

```shell
FastTree -nosupport lysozyme_primates.phy > lysozyme_primates.tree
```

Option used: `-nosupport`: (we don't want boostrap, as this will cause trouble for further analyses
in CodeML).

Finally, we could root the tree. Use `newick utilities` to root it by the group containing the
Marmoset sequence (Callithrix jacchus):

```shell
nw_reroot lysozyme_primates.tree ENSCATP00000036951 ENSMLEP00000013069
```

We can rerun and save the output as `lysozyme_primates_rooted.tree`:

```shell
nw_reroot lysozyme_primates.tree ENSCATP00000036951 ENSMLEP00000013069 > lysozyme_primates_rooted.tree
```

## Part 3: Run ancestral sequence reconstruction

The ancestral sequence reconstruction is done by CodeML, from the PAML package.

It is launched with

```shell
codeml control_file.ctl
```

You may have to copy the file "jones.dat" from the dat folder in the PAML package, or indicate its
location.

The control file contains many parameters:

```
      seqfile = lysozyme_primates.phy    * sequence data filename
     treefile = lysozyme_primates_root.tree   * tree structure file name
      outfile = lysozyme_primates.mlc    * main result file name

        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 2  * 0: concise; 1: detailed, 2: too much
      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 2  * 1:codons; 2:AAs; 3:codons-->AAs

        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
   aaRatefile = ./jones.dat  * only used for aa seqs with model=empirical(_F)

                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

        model = 2
                   * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                   * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
        Mgene = 0     * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
                   * AA: 0:rates, 1:separate

    fix_alpha = 0   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0.5 * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 1   * different alphas for genes
        ncatG = 4   * # of categories in dG of NSsites models


        getSE = 0  * 0: don\'t want them, 1: want S.E.s of estimates

 RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

   Small_Diff = .5e-6
    cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?
       method = 1  * Optimization method 0: simultaneous; 1: one branch a time


* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt.,
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt.,
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.


Explanation of some parameters:
runmode = 0 => We provide the tree.
clock = 0 => We don\'t set a molecular clock. We assume that the genes are evolving at different rate.
aaDist = 0 => We don\'t use the physicochemical properties of the amino acid.
aaRatefile = ./jones.dat => We use the JTT matrix. Other matrix could be used (WAG, etc...)
model = 2 => We use an empirical model (= substitutions matrix such as JTT).
fix_alpha = 0 => We estimated the alpha parameter of the gamma distribution.
alpha = 0.5 => We start the estimation from 0.5
RateAncestor = 1 => Force the estimation of ancestral states.
cleandata = 0 => Keep all ambigous data ("-", "X").
```

Please have a look at the output.

CodeML will also write into many files, but only two are of interest here:

- lysozyme_primates.mlc => Contains many information on evolutionary rates.
- rst => Contains ancestral states for sites and for nodes. Please have a look at both files.

In the rst file, there is also how the tree has been annotated, under the line "tree with node
labels for Rod Page's TreeView". You can copy this tree into a file (i.e.
"lysozyme_primates_annotated.tree"), open it with NJplot, and display the "bootstrap values". We can
see at which node correspond which ancestral sequence.

Questions:

- What is the estimated alpha parameter of the gamma curve?
- How many categories where used? And what are their frequencies?
- In Jalview, we observed at column 68 a mixture of E(Glu), Q(Gln) and R(Arg). But what was the most
  likely state of this position in the last common ancestor of all sequences? What are the
  probabilities of A(AlA), E(Glu), Q(Gln) and R(Arg)?
- What are the evolutionary events (amino acid substitutions) at the basis of the Hominoidea clade
  (Hylobates lar, Gorilla, Pan Paniscus, Homo sapiens)?

Now it is time to extract ancestral sequences and put them in a file. The rst file is quite
difficult to parse, hopefully, each ancestral sequence start by "node".

Use the script parse_rst.py to parse the rst:

```shell
parse_rst.py rst
```

It displays ancestral sequences in FASTA format. Let's put them in a file:

```shell
parse_rst.py rst > ancestral_sequences.fasta
```

## Part 4: Compute physico-chemical properties on ancestral sequences.

In Biopython, there is a function to compute the isoelectric point (pI):

analysed_protein = ProtParam.ProteinAnalysis(sequence) pI = analysed_protein.isoelectric_point()

The following script will compute the pI for all sequences in a FASTA file: compute_pI.py

By launching it, we can retrieve the pI for modern primate lysozymes:

```shell
compute_pI.py lysozyme_primates.fasta
```

And similarly for ancestral sequences:

```shell
compute_pI.py ancestral_sequences.fasta
```

## Part 5: Map properties on tree

We could easily map ancestral properties on the tree. The tree provided in rst contains nodes where
bootstrap information is. We just need to change the values of these nodes by the corresponding pI.
We will use the script `map_on_tree.py`:

```shell
map_on_tree.py ancestral_sequences.fasta lysozyme_primates_annotated.tree > lysozyme_primates_annotated_pI.tree
```

Have a look at both trees in a text editor.

You can load it in NJplot and see the different pI at the bootstrap place.

Alternatively, you can install FigTree: <http://tree.bio.ed.ac.uk/software/figtree/>

and load the tree.

In FigTree, set the following parameters: Appearance->Colour by:label Setup: Colours -> Scheme:
Colour gradient Tick gradient Line Weight 4
