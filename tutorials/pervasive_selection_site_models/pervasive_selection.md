# Detecting pervasive positive selection with site-models from CodeML / PAML

Disclaimer: Please don't hesitate to contact me if there is anything which is not working on your
computer, or any thing unclear, or even comments to improve it.
￼

## Theoretical principles:

When mutations are advantageous for the fitness, they are propagated at a higher rate in the
population. The selective pressure can be computed by the dN/dS ratio (ω). dS represents the
synonymous rate (keeping the same amino acid) and dN the non-synonymous rate (changing the amino
acid). In the absence of evolutionary pressure (genetic drift), the synonymous and non-synonymous
rates are supposed to be equal, so the dN/dS ratio is equal to 1. Under purifying selection, natural
selection prevents the replacement of amino acidS, so the dN will be lower than the dS, and
dN/dS < 1. And under positive selection, when mutations are advantageous for the fitness, they are
propagated at a higher rate in the population, so the replacement rate of amino acid is favoured by
selection, and dN/dS > 1.

We can distinguish two types of positive selection: pervasive positive selection and episodic
positive selection. The former implies that a site will be under continuous changes (i.e. adapting
to pathogens under arm-race), while the later implies that a site will change once and then be kept
in the clade (i.e. providing an advantage in a new environment). We can detect the later using the
branch-site model, for which I wrote the previous [tutorial on branch-site](..//positive_selection/branch_site.md).

To detect pervasive positive selection, we will use the site models from CodeML/PAML. Those models
allow the clustering of aligned columns (sites) in different groups, each group having a different
dN/dS value. There are many different sites models in CodeML, all assuming that the dN/dS ratio is
the same across branches, but different between sites.

Here are the different models we will use in this tutorial:

- M0: one unique dN/dS for all sites. This is the most basic model.
- M1a: assumes two categories of sites: sites with dN/dS<1 (negative selection) and sites with
  dN/dS =1 (neutral evolution).
- M2a: assumes three categories of sites: sites with dN/dS<1 (negative selection), sites with
  dN/dS=1 (neutral evolution) and sites with dN/dS>=1 (positive selection).
- M3: assumes multiple categories of selection, not necessarily positive selection.
- M7: assumes 10 categories following a beta-distribution of sites, all with different dN/dS <=1.
- M8: assumes 10 categories following a beta-distribution of sites, grouped, all with different
  dN/dS <=1, and an additional 11th category with dN/dS >=1 (positive selection allowed).
- M8a: assumes 10 categories following a beta-distribution of sites, grouped, all with dN/dS <=1,
  and an additional 11th category with dN/dS =1 (no positive selection allowed).

## Practical

We will focus on the major histocompatibility complex (MHC) protein, which detects peptides from
pathogens. As this gene is in the front line against invaders, it is submitted to strong selective
pressure to rapidly detect new antigenic peptides. Early work on positive selection was focused on
the MHC, so this is a very good example for this practical.

The uniprot code is HLA class II histocompatibility antigen, DQ beta 1 chain.
The Ensembl gene id is: ENSG00000179344.

### Download data from Ensembl: <http://www.ensembl.org>

Go on Ensembl website, and search for ENSG00000179344

#### 1) Download orthologues sequence

Comparative Genomics => Orthologues => Download orthologues

Then choose: Fasta, Unaligned sequences – CDS).
￼

#### 2) Download subtree

Comparative Genomics => Gene tree
Click on the blue node to select the following group: “Placental mammals ~100 MYA (Boreoeutheria)”
Gene Count 46

=> Export sub-tree Tree or Alignment
=> Format Newick, options "Full (web)" and "Final (merged) tree".

You have now two starting files:
- Sequences file: `Human_HLA_DQB1_orthologues.fa`
- Tree file: `HLA_DQB1_gene_tree.nh`

We are now going to process these files, in order to generate an alignment, visualise it, remove
spurious sequences and columns.

First, let’s rename these files to have shorter names

```shell
cp Human_HLA_DQB1_orthologues.fa HLA_DQB1.cds.fasta;
cp HLA_DQB1_gene_tree.nh HLA_DQB1.nh
```

Remove the species tag in the gene name

```shell
remove_ensembl_name_in_tree.py HLA_DQB1.nh > HLA_DQB1.tree
```

Extract gene names with Newick Utils (https://gensoft.pasteur.fr/docs/newick-utils/1.6/nwutils_tutorial.pdf)

```shell
nw_labels -I HLA_DQB1.tree > HLA_DQB1_names.txt
```

Extract CDS sequences that are in the tree, according to the extracted names
```shell
extract_sequences.py HLA_DQB1_names.txt HLA_DQB1.cds.fasta > HLA_DQB1_subset.cds.fasta
```

Translate CDS sequences to Amino Acid sequences using python script
```shell
translate_dna.py HLA_DQB1_subset.cds.fasta > HLA_DQB1_subset.aa.fasta
```

Make an alignment of these Amino Acid sequences
```shell
mafft-linsi HLA_DQB1_subset.aa.fasta > HLA_DQB1_subset.aa.mafft.fasta
```
Align CDS sequences by mapping them on the AA Alignment
```shell
realign_nuc_on_aa.py HLA_DQB1_subset.aa.mafft.fasta \
      HLA_DQB1_subset.cds.fasta \
      HLA_DQB1_subset.cds.mafft.fasta
```

With Jalview, load the alignment `Human_HLA_DQB1_subset.cds.mafft.fasta` and Visualise
it to see if there isn’t anything wrong

Move human sequence (`ENSP00000364080`) on top (use the key arrows). This is to inform CodeML to
use this sequence as reference.
Save the alignment as FASTA file again (CTRL-S)
#TODO: add script to automatise that.


Remove spurious sequences and columns with TrimAl
```shell
trimal -automated1 -in HLA_DQB1_subset.cds.mafft.fasta -resoverlap 0.75 -seqoverlap 85 -out HLA_DQB1_subset.cds.mafft.trimal.fasta -htmlout HLA_DQB1_subset.cds.mafft.trimal.html -colnumbering > HLA_DQB1_subset.cds.mafft.trimal.cols
```
Convert trimed sequences from FASTA to PHYLIP
```shell
convert_fasta2phylip.py HLA_DQB1_subset.cds.mafft.trimal.fasta HLA_DQB1_subset.cds.mafft.trimal.phy
```
Extract gene names to id.list
```shell
grep ">" HLA_DQB1_subset.cds.mafft.trimal.fasta | cut -c 2- > id.list
```
Using Newick Utilities, we load the id file to extract a pruned subtree from the starting tree
(contains only the taxa from the alignment)
```shell
id_list=$(cat id.list | xargs)
echo "$id_list"
eval "nw_prune -v HLA_DQB1.tree $id_list" > HLA_DQB1_subset.tree
```

Now we end up with two files:
- Alignment: `HLA_DQB1_subset.cds.mafft.trimal.phy`
- Tree: `HLA_DQB1_subset.tree`

## 2) Estimation of evolutionary values

This is the core of the tutorial. We will use codeml with three different control files (.ctl).
Each computation could take up to 30-60 minutes, depending of your CPU.

Compute many different site models: M0, M1a, M2a, M3 and M7.
Save the following commandS in HLA_DQB1_M0M1M2M3M7M8.ctl file

```
     seqfile = HLA_DQB1_subset.cds.mafft.trimal.phy  * sequence data file name
    treefile = HLA_DQB1_subset.tree                  * tree structure file name
     outfile = HLA_DQB1_M0M1M2M3M7M8.mlc             * main result file name

       noisy = 9   * 0,1,2,3,9: how much rubbish on the screen
     verbose = 1   * 1: detailed output, 0: concise output
     runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

     seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
   CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
       clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
      aaDist = 0   * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
       model = 0   * models for codons:
                   * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
     NSsites = 0 1 2 3 7 8 * 0:one w; 1:NearlyNeutral; 2:PositiveSelection;
                           * 3:discrete; 4:freqs; 5:gamma;6:2gamma;
                           * 7:beta;8:beta&w;9:beta&gamma;10:3normal
       icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
       Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all
   fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
       kappa = 2   * initial or fixed kappa
   fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
       omega = 1   * initial or fixed omega, for codons or codon-based AAs
       getSE = 0       * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0       * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
  Small_Diff = .45e-6  * Default value.
   cleandata = 0       * remove sites with ambiguity data (1:yes, 0:no)?
 fix_blength = 0       * 0: ignore, -1: random, 1: initial, 2: fixed
```

And execute it (this can take ~ 30-45 minutes):

```shell
codeml HLA_DQB1_M0M1M2M3M7M8.ctl
```

Important! Copy rst file to another name

```shell
cp rst HLA_DQB1_M0M1M2M3M7.rst.txt
```

One last thing is to compute site model M8a, which is the same as M8, except we fix the dN/dS to 1
(only negative selection and neutral evolution allowed). Save the following commandS in the file
"HLA_DQB1_M8a.ctl"

```
     seqfile = HLA_DQB1_subset.cds.mafft.trimal.phy   * sequence data file name
    treefile = HLA_DQB1_subset.tree                  * tree structure file name
     outfile = HLA_DQB1_M8a.mlc                      * main result file name

       noisy = 9   * 0,1,2,3,9: how much rubbish on the screen
     verbose = 1   * 1: detailed output, 0: concise output
     runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

     seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
   CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
       clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
      aaDist = 0   * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
       model = 0   * models for codons:
                   * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
     NSsites = 8   * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;
                   * 4:freqs; 5:gamma;6:2gamma;
                   * 7:beta;8:beta&w;9:beta&gamma;10:3normal
       icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
       Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all
   fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
       kappa = 2   * initial or fixed kappa
   fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate
       omega = 1   * initial or fixed omega, for codons or codon-based AAs
       getSE = 0       * 0: don\'t want them, 1: want S.E.s of estimates
RateAncestor = 0       * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
  Small_Diff = .45e-6  * Default value.
   cleandata = 0       * remove sites with ambiguity data (1:yes, 0:no)?
 fix_blength = 0       * 0: ignore, -1: random, 1: initial, 2: fixed
```

And execute (this can take ~ 5-10 minutes):
```shell
codeml HLA_DQB1_M8a.ctl
```

## 3) Results
### 3.1) Identification of positive selection

Have a look at the mlc files. If you want to retrieve the log-likelihood values:

```shell
grep "lnL" *.mlc
```

```
HLA_DQB1_M0M1M2M3M7M8.mlc:lnL(ntime: 32  np: 34):  -4838.327776      +0.000000
HLA_DQB1_M0M1M2M3M7M8.mlc:lnL(ntime: 32  np: 35):  -4711.107980      +0.000000
HLA_DQB1_M0M1M2M3M7M8.mlc:lnL(ntime: 32  np: 37):  -4692.732347      +0.000000
HLA_DQB1_M0M1M2M3M7M8.mlc:lnL(ntime: 32  np: 38):  -4692.078126      +0.000000
HLA_DQB1_M0M1M2M3M7M8.mlc:lnL(ntime: 32  np: 35):  -4718.466841      +0.000000
HLA_DQB1_M0M1M2M3M7M8.mlc:lnL(ntime: 32  np: 37):  -4690.545425      +0.000000
HLA_DQB1_M8a.mlc:         lnL(ntime: 32  np: 36):  -4706.268471      +0.000000
```

The order of lines being: M0, M1, M2, M3, M7, M8 and M8a. For each model, you directly get the number of parameters (np) and the log-likelihood value:

```
Model n.p. lnL
M0 34 -4838.3278
M1a 35 -4711.1080
M2a 37 -4692.7323
M3 38 -4692.0781
M7 35 -4718.4668
M8 37 -4690.5454
M8a 36 -4706.2685
```

Using these models, we can construct four likelihood-ratio tests (LRT), where three of them will
tell us if there is significant positive selection or not.

#### 1) M0-M3:
This one is an exception and will only tell us if there are different categories of
sites under different selective pressures. This test is not used to detect positive selection, and it is nearly always significant.

```
2x(L1-L0) = 2x[(-4692.0781) – (-4838.3278)] = 292.4993
d.f. = 38-34 = 4
=> 4.49405E-62
````

#### 2) M1a-M2a:
This test was the first site model developed to detect positive selection. We
contrast a model with 2 classes of sites against a model with 3 classes of sites.
Degree of freedom = 2.
```
2x(L1-L0) = 2x[(-4692.7323) – (-4711.1080)] = 36.7513
d.f. = 37-35 = 2
=> 1.04608E-08
```
The test is significant, so there is positive selection. This model is very conservative, and can
lack power under certain conditions.

#### 3) M7-M8:
This test also detects positive selection. We contrast a model with 10 classes of sites against a
model with 11 classes of sites.
Degree of freedom = 2.

```
2x(L1-L0) = 2x[(-4690.5454) – (-4838.3278)] = 55.842832
d.f. = 37-35 = 2
=> 7.47968E-13
```

The test is significant, so there is positive selection. However, this model can have problem power under certain conditions, and the following LRT is preferred.

#### 4) M8-M8a:
This is the latest test. We contrast a model with 11 classes of sites where positive
is not allowed (dN/dS=1) against a model with 11 classes of sites where positive is allowed
(dN/dS >=1).
Degree of freedom = 1.

```
2x(L1-L0) = 2x[(-4690.5454) – (-4706.2685)] = 31.446092
d.f. = 37-36 = 1
=> 1.48446E-07
```

This is the preferred test, combining power and robustness.

### 3.2) Identification of sites

As these tests are significant, we can move to the next step, which is the precise identification of
sites under positive selection. If the previous step was not significant, we should not move to this
stage.

In your mlc file, under the section of Model M2a and M8 (the only ones that allow positive
selection), you will find a section called “Bayes Empirical Bayes (BEB) analysis
(Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)”. This section contains the list that
have a BEB score [Pr(w>1)] higher than 50%. BEB values higher than 95% are indicated by *and sites
with values higher than 99% are indicated by **. Sometimes, there is no site detected, which means
there is probably a problem in your analysis or dataset if your test is significant. Sometimes, you
will find a lot of sites, which seems worrying, but it just means the average BEB (baseline) is
slightly above 50%. The most interesting sites are those with a BEB>95% (* or**).

```
Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)
Positively selected sites (*: P>95%; **: P>99%)
(amino acids refer to 1st sequence: ENSP00000364080)

            Pr(w>1)     post mean +- SE for w
     4 R      0.961*        2.858 +- 0.698
    36 F      0.671         2.150 +- 1.019
    53 L      0.999**       2.956 +- 0.610
    84 D      1.000**       2.957 +- 0.608
    97 G      0.649         2.098 +- 1.018
   112 V      0.969*        2.888 +- 0.695
   114 F      0.999**       2.955 +- 0.611
   115 R      0.743         2.367 +- 1.093
   116 G      1.000**       2.957 +- 0.609
   121 R      0.835         2.593 +- 0.993
   247 R      0.617         2.050 +- 1.110
```

The column correspondS to the position in the trimmed alignment (i.e. not related to the position in
the reference sequence, which is ENSP00000364080, the human). However, the amino acid correspondS
exactly to the reference sequence.

One of the site with the strongest BEB value is 84, with 1.000.  Its own dN/dS value is 2.957, with
a standard deviation of 0.608. As I said previously, the sites given by CodeML don’t correspond to
the human sequence. You can use the following script to extract the real position in the human
sequence:
```shell
get_position_cds_trimal.py HLA_DQB1_subset.cds.mafft.fasta HLA_DQB1_subset.cds.mafft.trimal.
cols "ENSP00000364080" 84
```

=> 84 89 D

=> This site 84 in CodeML (Trimal) correspondS to amino acid site 89 in the human sequence and codes
for an aspartic acid.

In total, we have six sites that are strongly interesting (BEB>95%): 4, 53, 84, 112, 114 and 116.
Let’s repeat the same for all these sites

```shell
site_list="4 53 84 112 114 116"
for site in $site_list; do ./get_position_cds_trimal.py HLA_DQB1_subset.cds.mafft.fasta HLA_DQB1_subset.cds.mafft.trimal.cols "ENSP00000364080" $site; done
```
=>
```
  4   8 R
 53  58 L
 84  89 D
112 117 V
114 119 F
116 121 G
```

We can do the same for sites with 50%<BEB<95%:
```shell
site_list="36 97 115 121 247"
for site in $site_list; do ./get_position_cds_trimal.py HLA_DQB1_subset.cds.mafft.fasta HLA_DQB1_subset.cds.mafft.trimal.cols "ENSP00000364080" $site; done
```
=>
```
 36  41 F
 97 102 G
115 120 R
121 126 R
247 252 R
```

We can also plot the dN/dS value per sites. At the bottom of the rst file
"HLA_DQB1_M0M1M2M3M7.rst.txt", extract the last part which looks like this and save as “beb.txt”

```
   1 K   0.02943 0.13133 0.20099 0.20092 0.16262 0.11574 0.07524 0.04540 0.02539 0.01292 0.00002 ( 3)  0.395 +-  0.198
   2 A   0.00278 0.02154 0.05527 0.09166 0.12189 0.14109 0.14752 0.14134 0.12323 0.09274 0.06093 ( 7)  0.737 +-  0.535
   3 L   0.00726 0.04637 0.09959 0.13896 0.15621 0.15361 0.13717 0.11292 0.08526 0.05671 0.00595 ( 5)  0.551 +-  0.270
   4 R   0.00000 0.00000 0.00000 0.00003 0.00019 0.00082 0.00248 0.00582 0.01104 0.01731 0.96230 (11)  2.860 +-  0.696
   5 I   0.44686 0.24476 0.14583 0.08007 0.04205 0.02141 0.01061 0.00510 0.00234 0.00099 0.00000 ( 1)  0.168 +-  0.149
   6 P   0.34022 0.22595 0.16343 0.10859 0.06876 0.04207 0.02495 0.01430 0.00780 0.00391 0.00001 ( 1)  0.221 +-  0.187
```

- 1st column = position in the trimmed alignment.
- 2nd column = amino acid from the reference sequence.
- 3rd to 13th column = BEB score for each class (10 neutral + 1 allowing positive selection).
- 14th = most likely class.
- 15th = estimated dN/dS value at this position.
- 16th = standard deviation for this dN/dS value.

You can see that position 4, the most likely class is 11th (BEB=0.96) with a dN/dS = 2.86+-0.70.

CodeML output uses a fixed delimitation. To parse it in R, we need to remove the space between
bracket and number

```shell
cat beb.txt | perl -pe "s/\( /\(/g" > tmp.txt; mv tmp.txt beb.txt
```

Here are some commands to use in R, to produce the following plot

```R
if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE)
   library(ggplot2)
}

df<-read.table("beb.txt", sep = "")

df$beb <- "No"
df$beb[df$V13 > 0.50] <- "Yes"

p <- ggplot(df, aes(V1, V15))
p + geom_point(aes(colour = factor(beb)))+
    geom_hline(yintercept = 1)+
    scale_color_manual(values = c("black", "red"))+
    labs(x = "Residue position")+
    labs(y = "Selective pressure [dN/dS]")+
    theme_bw()+
    theme(legend.position="none")+
    theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank())
ggsave("beb.png", height=3, width=4)
```

We can see that sites under positive selection represent only a small fraction, and most sites are
under strong purifying selection.

## 4) Visualisation of sites in 3D structure

We found that sites seems randomly distributed according to their residue position. It would be
interesting to see if they form a pattern in the 3D structure.
Go download the following pdb file 1uvq:
```shell
wget http://files.rcsb.org/download/1UVQ.pdb
```

Or download it from the webpage: <http://www.rcsb.org/pdb/explore.do?structureId=1UVQ>

Then load it in PyMOL:
```shell
pymol 1UVQ.pdb &
```

First, define the different molecules to their associated polypeptide chains
```
select HLA_DQA1, chain A
select HLA_DQB1, chain B
select peptide,  chain C
```
Highlight in cartoon/sticks and colour chains and peptide in different colours
```
hide everything
show cartoon, HLA_DQA1
show cartoon, HLA_DQB1
show sticks, peptide
colour white, HLA_DQA1
colour grey60, HLA_DQB1
colour green, peptide
```

To spice things up, the site numbering in many PDB files doesn’t correspond to the human sequence.
So we have to renumber it

```
alter HLA_DQB1, resi=int(resi)+32
```

Highlight sites with BEB>95%, and display as yellow spheres.

```
select sites_BEB95, HLA_DQB1 and resi 58+89+117+119+121
show spheres, sites_BEB95
colour yellow, sites_BEB95
```

Highlight sites with 50%<BEB<95%, and display them as yellow sticks.

```
select sites_BEB50, HLA_DQB1 and resi 41+102+120+126+252
show sticks, sites_BEB50
colour yellow, sites_BEB50
```

=> Most sites under positive selection (in yellow) are exactly in the binding site (in green),
facing the target peptide. This example with the MHC has been widely described in the literature
[Hugues AL et al. 1988]. To my disappointment, and unless there was some development since my
academic years, there are not so many examples that selective pressure, 3D structural and
experimental validation.


To have a nice finish as in the figure above, rotate the structure the way you want and type:

```
bg_color white
util.cnc
select none
set ray_opaque_background, 1
ray 2000
save 1UVQ_positive_sites.png
```

Et voila:
![Picture of MHC with positively selected site](1UVQ_positive_sites.png "1UVQ_positive_sites")
