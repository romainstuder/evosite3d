# Tutorial: estimating the stability effect of a mutation with FoldX

## Introduction:

Here is a brief tutorial on how to use FoldX to estimate the stability effect of a mutation in a 3D
structure. The stability (ΔG) of a protein is defined by the free energy, which is express in
kcal/mol. The lower it is, the more stable it is. G is difference of free energy between a wild-type
and mutant. A mutation that brings energy (ΔΔG>kcal/mol) will destabilise the structure, while a mutation that remove energy (ΔΔG<kcal/mol) will stabilise the structure. A common threshold is to say that a mutation has a significant effect if ΔΔG is >1 kcal/mol, which roughly corresponds to a single hydrogen bond.

A way to compute the free energy is to use molecular dynamics. Main problem: it can be very time-consuming.

FoldX uses an empirical method to estimate the stability effect of a mutation. The executable is available here: <http://foldx.crg.es/>

NB: I strongly encourage to read the manual (before or in parallel of this tutorial).

Foldx was used in many studies, i.e.:
Tokuriki N, Stricher F, Serrano L, Tawfik DS. How protein stability and new functions trade off. PLoS Comput Biol. 2008 Feb 29;4(2):e1000002 <http://dx.doi.org/10.1371/journal.pcbi.1000002>

Dasmeh P, Serohijos AW, Kepp KP, Shakhnovich EI. Positively selected sites in cetacean myoglobins contribute to protein stability. PLoS Comput Biol. 2013;9(3):e1002929. <http://dx.doi.org/10.1371/journal.pcbi.1002929>
And I personally used it in three of my studies:
Studer RA, Christin PA, Williams MA, Orengo CA. Stability-activity tradeoffs constrain the adaptive evolution of RubisCO. Proc Natl Acad Sci U S A. 2014 Feb 11;111(6):2223-8. <http://dx.doi.org/10.1073/pnas.1310811111>
Studer RA, Opperdoes FR, Nicolaes GA, Mulder AB, Mulder R. Understanding the functional difference between growth arrest-specific protein 6 and protein S: an evolutionary approach. Open Biol. 2014 Oct;4(10). pii: 140121. <http://dx.doi.org/10.1098/rsob.140121>

Rallapalli PM, Orengo CA, Studer RA, Perkins SJ. Positive selection during the evolution of the blood coagulation factors in the context of their disease-causing mutations. Mol Biol Evol. 2014 Nov;31(11):3040-56. <http://dx.doi.org/10.1093/molbev/msu248>

Example:

The structure is a bacterial cytochrome P450 (PDB:4TVF). You can download it PDB file (4TVF.pdb) from here: <http://www.rcsb.org/pdb/explore.do?structureId=4TVF>

We would to test the stability of mutatinf the leucine (L) at position 280 to an aspartic acid (D).

Here is the original structure, with Leu280 in green, and residues around 6Å in yellow:
￼

FoldX has different modes to run it, but I use the mode "runfile", which contains COMMANDS and OPTIONS in one single file.

FoldX works in two steps:

### 1) Repair the structure.

There are frequently problems in PDB structures, like steric clashes. FoldX will try to fix them and lower the global energy (ΔG). The command "RepairPDB" is better than the "Optimize" command. Here is the command file "foldx_repair.txt":
<TITLE>FOLDX_runscript;
<JOBSTART>#;
<PDBS>4TVF.pdb;
<BATCH>#;
<COMMANDS>FOLDX_commandfile;
<RepairPDB>#;
<END>#;
<OPTIONS>FOLDX_optionfile;
<Temperature>298;
<R>#;
<pH>7;
<IonStrength>0.050;
<water>-CRYSTAL;
<metal>-CRYSTAL;
<VdWDesign>2;
<OutPDB>true;
<pdb_hydrogens>false;
<END>#;
<JOBEND>#;
<ENDFILE>#;
We indicate which PDB file it needs to use, that we want to repair it (<RepairPDB>), that it will use water and metal bonds from the PDB file (<water>-CRYSTAL; <metal>-CRYSTAL;) and that we want a PDB as output (<OutPDB>true).

You have to run the command file "foldx_repair.txt" like this:
foldx3b6 -runfile foldx_repair.txt
This process is quite long (around 10 minutes). Here is the result (the original structure is now in white, while the repaired structure is in yellow/green):
￼

 We can see that some side chains have slightly move (in particular Phe16).

The starting free energy ΔG was 64.99 kcal/mol, and it was lowered to -48.15 kcal/mol, which is now stable (remember that a "+" sign means unstable while a "-" sign means stable).

Once it's finished, it will produce a file named "RepairPDB_4TVF.pdb", which you will use in the next step.

2) Perform the mutation

The mutation itself is perform by the BuildModel function. There are other methods, but the BuildModel is the most robust. You need also to specify the mutation in a separate file "individual_list.txt".

In the command file, you will see that is RepairPDB_4K33.pdb and not 4K33.pdb that is mutated. You will also notice "<numberOfRuns>3;". This is because some residues can have many rotamers and could have some convergence problems. You may to increase this values to 5 or 10, in case you are mutated long residues (i.e. Arginine) that have many rotamers.

Here the command file "foldx_build.txt":
<TITLE>FOLDX_runscript;
<JOBSTART>#;
<PDBS>RepairPDB_4TVF.pdb;
<BATCH>#;
<COMMANDS>FOLDX_commandfile;
<BuildModel>#,individual_list.txt;
<END>#;
<OPTIONS>FOLDX_optionfile;
<Temperature>298;
<R>#;
<pH>7;
<IonStrength>0.050;
<water>-CRYSTAL;
<metal>-CRYSTAL;
<VdWDesign>2;
<OutPDB>true;
<numberOfRuns>3;
<END>#;
<JOBEND>#;
<ENDFILE>#;
and the "individual_list.txt" (just one line):
LA280D;
It contains the starting amino acid (L), the chain (A), the position (280) and the amino acid you want at the end (D). One line correspond to one mutant. It means you can mutate many residues at the same per line (mutant) and also  produce different mutants by different numbers of lines.

You can run it by:
foldx3b6 -runfile foldx_build.txt
It is much faster this time (i.e. a few seconds) and will produce many files.

FoldX will first mutate the target residue (L) to itself (L) and move it as well as all neighbouring side chains multiple times. We can see that Leu280 (green) was rotated:
￼

=> This is will give the free energy of the wild-type (let's call it ΔGwt).

Then, it will mutate the target residue (L) to the desired mutant (D) and move it as well as all neighbouring side chains multiple times. We can see that Leu280 is mutated to Asp280 (see the two oxygen atoms in red):
￼

=> This is will give the free energy of the mutant (let's call it ΔGmut).

The difference in free energy (ΔΔG) is given by ΔGmut-ΔGwt.

In the file "Raw_BuildModel_RepairPDB_4TVF.fxout", you can retrieve the energy of the three runs for both WT and Mutant.

Run1:

* ΔGmut = RepairPDB_4TVF_1_0.pdb = -41.1377 kcal/mol
* ΔGwt = WT_RepairPDB_4TVF_1_0.pdb = -46.0464 kcal/mol
* => ΔΔG = ΔGmut-ΔGwt = (-41.1377)-(-46.0464) = +4.9087 kcal/mol

One file contains the average difference over all runs: "Average_BuildModel_RepairPDB_4K33.fxout".
You will notice that the difference in free energy ΔΔG is +4.84 kcal/mol (+- 0.06 kcal/mol).

=> It means the mutation L280D is highly destabilising (positive value, and much above 1.0 kcal/mol). Here is the final mutant:

PS: Another way to define the threshold is to use the SD deviation multiple times:

The reported accuracy of FoldX is 0.46 kcal/mol (i.e., the SD of the difference
between ΔΔGs calculated by FoldX and the experimental values). We can bin the ΔΔG values into seven categories:

1. highly stabilising (ΔΔG < −1.84 kcal/mol);
2. stabilising (−1.84 kcal/mol ≤ ΔΔG < −0.92 kcal/mol);
3. slightly stabilising (−0.92 kcal/mol ≤ ΔΔG < −0.46 kcal/mol);
4. neutral (−0.46 kcal/mol < ΔΔG ≤ +0.46 kcal/mol);
5. slightly destabilising (+0.46 kcal/mol < ΔΔG ≤ +0.92 kcal/mol);
6. destabilising (+0.92 kcal/mol < ΔΔG ≤ +1.84 kcal/mol);
7. highly destabilising (ΔΔG > +1.84 kcal/mol).
