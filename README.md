# histone_modification_dynamics
Program for stochastic simulation of histone modification dynamics

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

    Rscript histone_modification_dynamics.R -a 1.2 -b 0.001 --regionLength 5000 --PRC2Peak1 1490 1510 --PRC2Peak2 3190 3210 --pme1 0.95 --pme2 0.2 --pme2adj 0.5 --pme3 0.01 --pme3adj 0.2 --enhancedTime 2 --Npop 5 --Nenzyme1 7 --Nenzyme2 7 --Nstep 5000 --period 1000 --TSSTES 2100 2400 3500 3700 --K36Me3 2100 2400 3500 3700 --K36Me2 2401 2700 --K36me3harm01 5 --K36me3harm12 7.5 --K36me3harm23 15 --K36me2harm01 2 --K36me2harm12 3 --K36me2harm23 5 -E FALSE --periodExpression 200 --mutRate 0.05 --K27Mharm01 1 --K27Mharm12 2.5 --K27Mharm23 5 --K27MharmPermanent FALSE --sequestrationTime 2 --normalBackRate 0 --K27MBackRate 0 --targetedRegion 1001 4000 --shortRangeProb 0.95 --shortRangeLimit 8 --shortRangeCentrality 0.5 --K27MintroduceTime 0 --equiPeriod 1 --output ModificationPattern
