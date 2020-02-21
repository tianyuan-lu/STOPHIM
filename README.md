# histone_modification_dynamics
Program for stochastic simulation of histone modification dynamics

We conceived a stochastic propagation model to simulate the dynamic process of deposition of H3K27 methylation marks by PRC2. In this study, we simulated a genomic region (length = 5,000 nucleosomes, corresponding to approximately 1 Mb of DNA) containing two PRC2 binding loci (at 2,000 and 3,000), in 20 cells. We introduced two actively transcribed regions marked with K36me3 and one intergenic region merged with K36me2.

## Bi-modal random walk of PRC2
In this conceptual model, PRC2 complexes are first recruited at specific genomic loci (nucleation sites) where they have a high binding affinity, and then randomly diffuse along the chromatin. We simulate the diffusion by a random walk process:

(1) We create a one-dimensional vector with continuous integers to represent coordinates of histones winding around the DNA sequence, and assign two narrow regions to harbor nucleation sites.

(2) We initiate eight PRC2 complexes at each of the two nucleation sites and record their coordinates, at time (t) = 0.

(3) For a PRC2 complex, we draw a random value: 
(3.1) with 95% probability, from a binomial distribution Binom(8, 0.5) to represent short-range diffusion due to spontaneous thermal motion; or
(3.2) with 5% probability, from a gamma distribution Gamma(0.5, 0.001) to represent possible long-range traveling of the complex after dissociating from the chromatin or across interacting genomic regions.

(4) We approximate this random value with the nearest (non-negative) integer and regard it as the distance the PRC2 complex travels within a time unit.

(5) We change the coordinate of the PRC2 complex at t = 1 if the travelling distance is non-zero, upstream or downstream with equal probability, based on the distance obtained in (4).

(6) For each PRC2 complex, we repeat (3) – (5) independently, so that we obtain the coordinates for all complexes at t = 1.

(7) We repeat (3) – (6) for a given period of time.

(8) We repeat (2) – (7) for a population of 20? cells, each cell independently yet having identical representation of the simulated chromatin and nucleation sites; We record the cumulative distribution of methylation marks until the relative abundance of different methylation marks reaches an equilibrium.
When a PRC2 complex travels beyond the border of our simulated region, we introduce a new complex initiating at a nucleation site to replace the one we lose track of, and independently simulate its trajectory together with other complexes.

## Catalytic features of PRC2
During the process described above, whenever a PRC2 complex encounters and resides on a histone, it methylates the unmethylated or lowly-methylated H3K27 to a more heavily-methylated status with differential catalytic activity. We simulate this dynamic process by conducting Bernoulli trials with pre-defined conversion rates at each time point at loci occupied by PRC2 complexes independently. We then iteratively record methylation status at each simulated locus in each cell. Generally, PRC2 more likely converts unmethylated H3K27 to mono-methylated H3K27, than mono-methylated H3K27 to di-methylated H3K27, than di-methylated H3K27 to tri-methylated H3K27. However, when an adjacent histone has tri-methylated H3K27, PRC2 will have largely enhanced catalytic capacity, especially for tri-methylating di-methylated H3K27, known as the allosteric activation process [reference]. The boosted capacity of an allosterically activated PRC2 will be maintained for a given short period of time after the complex dissociates from the histone. Meanwhile, antagonistic histone modifications, in particular H3K36me3 in actively transcribed regions and H3K36me2 in intergenic regions, decrease the catalytic activity of PRC2 complexes so that the conversion rates from lowly-methylated status to heavily-methylated status are reduced proportionally.

Cell division actively removes methylation marks
For each cell independently, we simulate events of cell division occurring periodically. During cell division, half of the methylation marks are randomly removed, replaced by unmethylated histones. Meanwhile, half of the PRC2 complexes are removed while the same number of complexes are initiated at the nucleation sites.

## Key parameters in PRC2 modeling 
First, PRC2 must be able to diffuse far enough from the nucleation sites to produce an approximately uniform distribution – outside of restricted regions – of all three K27me modifications. Second, within intergenic domains marked by K36me2, the antagonistic effect reduces the efficiency of deposition of K27me, with the strongest resistance to K27me3 and weakest to K27me1. As a result, the levels of K27me1 are slightly reduced and K27me3 greatly reduced in K36me2 domains as compared to adjacent intergenic regions. Seemingly paradoxically, the levels of K27me2, the intermediate mark are unaffected by the presence of K36me2. This is the result of relative differences in conversion efficiencies: although the levels of K27me1 are lower in those regions, once K27me1 is converted to K27me2 the conversion to K27me3 occurs very rarely, resulting in no change of K27me2 relative to the adjacent region. Third, within active genes, the deposition of K27me is hindered by the presence of K36me3 and/or histone turnover resulting from chromatin remodeling during transcription. Both of those phenomena have a strongest effect on K27me3 and the weakest on K27me1. As a result, in genes, K27me2 and K27me3 are depleted. Conversely, transcribed genes are the only genomic compartment where relatively high levels of K27me1 are retained; this can be explained not by preferential deposition of K27me1 in those regions, but the inability of PRC2 to convert the deposited K27me1 to higher methylation states. 

We tested whether plausible mechanisms for epigenome-wide re-patterning caused by K27M mutation can be realized by this model:

(1) We reduced the number of PRC2 complexes to examine whether this lack of enzyme could lead to observed phenomena.

(2) We explored whether systematically reducing the conversion rates and weakening allosteric activation could be the determinants.

(3) In the simulated scenario of harmful K27M mutants, we randomly introduced 5% K27M-mutant histones for each cell independently. Upon residing on a mutant histone which cannot be methylated, the catalytic capacity of PRC2 complexes will be permanently impaired. Namely, these PRC2 complexes become less likely to di-methylate mono-methylated H3K27, extremely less likely to tri-methylate di-methylated H3K27, and the allosteric activation effect will be inhibited. Moreover, the K27M-mutant histones are able to detain PRC2 complexes for a longer period of time than normal histones. As a consequence, the propagation of heavily-methylated marks is largely hindered.

## Additional features
More inhibitory mechanism induced by active transcription: Within actively transcribed regions, while stable deposition of H3K36me3 marks prevent deposition of di-methylated and tri-methylated H3K27 marks as the catalytic capacity of PRC2 decreases, methylation marks are also drastically removed through opening of double strands, exchange of nucleosomes and potential active de-methylation. Half of the methylated histones are substituted by unmethylated histones.
K27M mutants may induce higher possibility of dissociation from the chromatin.
Four figures will be generated: 

(1) Cumulative occupancy of PRC2 molecules

(2) Spatio-temporal distribution of H3K27 methylation marks

(3) Cumulative amount of H3K27 methylation marks

(4) Velocity of deposition of H3K27 methylation marks


Execute the following command for parameter details:
    
    Rscript histone_modification_dynamics.R -h

```ruby
usage: ./histone_modification_dynamics.R [-h] [-a ALPHA] [-b BETA]
                                         [--regionLength REGIONLENGTH]
                                         [--PRC2Peak1 PRC2PEAK1 [PRC2PEAK1 ...]]
                                         [--PRC2Peak2 PRC2PEAK2 [PRC2PEAK2 ...]]
                                         [--pme1 PME1] [--pme2 PME2]
                                         [--pme2adj PME2ADJ] [--pme3 PME3]
                                         [--pme3adj PME3ADJ]
                                         [--enhancedTime ENHANCEDTIME]
                                         [-N NPOP] [--Nenzyme1 NENZYME1]
                                         [--Nenzyme2 NENZYME2] [--Nstep NSTEP]
                                         [--period PERIOD]
                                         [--TSSTES TSSTES [TSSTES ...]]
                                         [--K36Me3 K36ME3 [K36ME3 ...]]
                                         [--K36Me2 K36ME2 [K36ME2 ...]]
                                         [--K36me3harm01 K36ME3HARM01]
                                         [--K36me3harm12 K36ME3HARM12]
                                         [--K36me3harm23 K36ME3HARM23]
                                         [--K36me2harm01 K36ME2HARM01]
                                         [--K36me2harm12 K36ME2HARM12]
                                         [--K36me2harm23 K36ME2HARM23]
                                         [-E EXPRESSION]
                                         [--periodExpression PERIODEXPRESSION]
                                         [--mutRate MUTRATE]
                                         [--K27Mharm01 K27MHARM01]
                                         [--K27Mharm12 K27MHARM12]
                                         [--K27Mharm23 K27MHARM23]
                                         [--K27MharmPermanent K27MHARMPERMANENT]
                                         [--sequestrationTime SEQUESTRATIONTIME]
                                         [--normalBackRate NORMALBACKRATE]
                                         [--K27MBackRate K27MBACKRATE]
                                         [--targetedRegion TARGETEDREGION [TARGETEDREGION ...]]
                                         [--shortRangeProb SHORTRANGEPROB]
                                         [--shortRangeLimit SHORTRANGELIMIT]
                                         [--shortRangeCentrality SHORTRANGECENTRALITY]
                                         [--K27MintroduceTime K27MINTRODUCETIME]
                                         [--equiPeriod EQUIPERIOD] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -a ALPHA, --alpha ALPHA
                        shape parameter of gamma distribution for simulating
                        PRC2 movement
  -b BETA, --beta BETA  rate parameter of gamma distribution for simulating
                        PRC2 movement; smaller beta allows longer-range
                        movement of PRC2
  --regionLength REGIONLENGTH
                        length of simulated genomic region
  --PRC2Peak1 PRC2PEAK1 [PRC2PEAK1 ...]
                        start/end coordinates indicating where a PRC2 binding
                        peak exists
  --PRC2Peak2 PRC2PEAK2 [PRC2PEAK2 ...]
                        start/end coordinates where a PRC2 binding peak exists
  --pme1 PME1           normal conversion rate from me0 to me1
  --pme2 PME2           normal conversion rate from me1 to me2
  --pme2adj PME2ADJ     increased conversion rate from me1 to me2 by
                        allosteric activation
  --pme3 PME3           normal conversion rate from me2 to me3
  --pme3adj PME3ADJ     increased conversion rate from me2 to me3 by
                        allosteric activation
  --enhancedTime ENHANCEDTIME
                        time for allosterically activated PRC2 to maintain its
                        enhanced catalytic ability; set to 0 to remove this
                        effect
  -N NPOP, --Npop NPOP  number of cells
  --Nenzyme1 NENZYME1   number of PRC2 molecules initially bound at peak1
  --Nenzyme2 NENZYME2   number of PRC2 molecules initially bound at peak2
  --Nstep NSTEP         total simulation time
  --period PERIOD       cell division interval
  --TSSTES TSSTES [TSSTES ...]
                        start/end coordinates for actively transcribed regions
  --K36Me3 K36ME3 [K36ME3 ...]
                        start/end coordinates for regions with K36me3 marks;
                        generally overlapping with actively transcribed
                        regions
  --K36Me2 K36ME2 [K36ME2 ...]
                        start/end coordinates for regions with K36me2 marks
  --K36me3harm01 K36ME3HARM01
                        K36me3 decreases me0 to me1 conversion rate by
                        K36me3harm01-fold
  --K36me3harm12 K36ME3HARM12
                        K36me3 decreases me1 to me2 conversion rate by
                        K36me3harm12-fold
  --K36me3harm23 K36ME3HARM23
                        K36me3 decreases me2 to me3 conversion rate by
                        K36me3harm23-fold
  --K36me2harm01 K36ME2HARM01
                        K36me2 decreases me0 to me1 conversion rate by
                        K36me2harm01-fold
  --K36me2harm12 K36ME2HARM12
                        K36me2 decreases me1 to me2 conversion rate by
                        K36me2harm12-fold
  --K36me2harm23 K36ME2HARM23
                        K36me2 decreases me2 to me3 conversion rate by
                        K36me2harm23-fold
  -E EXPRESSION, --Expression EXPRESSION
                        whether to include removal of methylated histones in
                        actively transcribed regions upon active gene
                        expression
  --periodExpression PERIODEXPRESSION
                        gene expression interval
  --mutRate MUTRATE     rate of K27M mutants; set to 0 to simulate a K27M-free
                        scenario
  --K27Mharm01 K27MHARM01
                        K27M decreases me0 to me1 conversion rate by
                        K27Mharm01-fold
  --K27Mharm12 K27MHARM12
                        K27M decreases me1 to me2 conversion rate by
                        K27Mharm12-fold
  --K27Mharm23 K27MHARM23
                        K27M decreases me2 to me3 conversion rate by
                        K27Mharm23-fold
  --K27MharmPermanent K27MHARMPERMANENT
                        whether K27M permanently damages the catalytic
                        capacity of PRC2 and prevents future allosteric
                        activation
  --sequestrationTime SEQUESTRATIONTIME
                        time to detain PRC2 on K27M; set to 0 to remove this
                        effect
  --normalBackRate NORMALBACKRATE
                        probability that a PRC2 automatically dissociates from
                        the genome and gets recruited back to the nucleation
                        site
  --K27MBackRate K27MBACKRATE
                        increased probability that a PRC2 returns to the
                        nucleation site when the molecule encounters K27M; set
                        to 0 to remove this effect
  --targetedRegion TARGETEDREGION [TARGETEDREGION ...]
                        start/end coordinates for the region used for
                        generating summary statistics; cropping out regions
                        close to the boundaries
  --shortRangeProb SHORTRANGEPROB
                        probability that a PRC2 takes small steps; set to 0 to
                        remove short-range movement
  --shortRangeLimit SHORTRANGELIMIT
                        maximum distance a PRC2 molecule can travel when it
                        takes small steps
  --shortRangeCentrality SHORTRANGECENTRALITY
                        centrality of the short-range movement; probability of
                        a PRC2 remains at its current locus is
                        shortRangeCentrality to the power of shortRangeLimit
  --K27MintroduceTime K27MINTRODUCETIME
                        time point when K27M is introduced; set to 0 to
                        introduce the mutants from the beginning of simulation
  --equiPeriod EQUIPERIOD
                        number of last time point(s) to be considered as
                        equilibrium stage; proportion of marks will be
                        reported by averaging across these time points. E.g.
                        set to 1 to report for only the last time point; set
                        to 100 to report average proportion across the last
                        100 time points
  -o OUTPUT, --output OUTPUT
                        prefix for output file names

```

Example script (equivalent to using default settings):

    Rscript histone_modification_dynamics.R -a 1.5 -b 0.001 --regionLength 5000 --PRC2Peak1 1490 1510 --PRC2Peak2 3190 3210 --pme1 0.95 --pme2 0.25 --pme2adj 0.5 --pme3 0.01 --pme3adj 0.15 --enhancedTime 2 --Npop 20 --Nenzyme1 8 --Nenzyme2 8 --Nstep 5000 --period 1000 --TSSTES 2100 2300 3500 3750 --K36Me3 2100 2300 3500 3750 --K36Me2 2301 2900 --K36me3harm01 5 --K36me3harm12 7.5 --K36me3harm23 15 --K36me2harm01 2 --K36me2harm12 3 --K36me2harm23 5 -E FALSE --periodExpression 100 --mutRate 0.05 --K27Mharm01 2 --K27Mharm12 2 --K27Mharm23 2 --K27MharmPermanent TRUE --sequestrationTime 5 --normalBackRate 0 --K27MBackRate 0 --targetedRegion 1001 4000 --shortRangeProb 0.95 --shortRangeLimit 8 --shortRangeCentrality 0.5 --K27MintroduceTime 0 --equiPeriod 1 --output ModificationPattern
