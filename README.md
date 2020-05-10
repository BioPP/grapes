# grapes
(Greetings_here_is_yet_another) Rate of Adaptive Protein Evolution Software

Introduction
------------

**grapes** is a program estimating the of adaptive and non-adaptive amino-acid substitution rate from
polymorphism and divergence coding sequence data. It is written in C++ and based on the Bio++
libraries. This file is relevant to version 1.0.
**grapes** essentially re-implements Adam Eyre-Walker's DoFE program (Eyre-Walker and Keightley
2009) using a different optimization procedure and with a couple of additional features. The details of
the method are available from Galtier (2016). A similar approach and program have been developed by
Tataru et al. (2017). See Rousselle et al. (2018) for a comparison of the methods and clarification of
their relationships.

**multigrapes** is a program fitting a population genetic model to synonymous + non-synonymous Site
Frequency Spectrum (SFS) data, and is essentially a multi-SFS version of the program grapes. This file is relevant to version 1.1.
The main new feature of multi_grapes, compared to grapes, is the option of fitting model parameters
assumed to be shared by distinct data sets – specifically, the shape of the distribution of fitness effects
(DFE). This is the main focus of Galtier & Rousselle 2020, and the reason why this program was
developed. In addition, new models of the DFE are implemented in multi_grapes that are absent in
grapes.
These new option come with costs: (1) several of the models implemented in grapes are not available in
multi_grapes; (2) no estimate of the adaptive substitution rate is performed by multi_grapes.

Installation
------------

The distributed binaries for GNU/Linux are statically linked and ready to run.

Compilation
-----------

If, for any reason, a source compilation is needed, here are the dependencies :
* gcc >= 4.4
* libgsl >= 1.0
* Bio++. The source code is ready for the latest development version of the libraries. Stable version >= 2.0 (bpp-core, bpp-phyl, bpp-seq) can be used with little modification of the code.

Compilation can be achieved using **cmake**:
```bash
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local ..
make make
```
See the INSTALL file and cmake documentation for more detailed instructions.

Execution
---------

```bash
grapes -in input_file.dofe -out output_file.csv -model model_name [options]

multigrapes -in input_file.dofe -out output_file.csv -model model_name [options]
```

Input
-----
**grapes** needs information on the number of synonymous and non-synonymous substitutions between
focal and outgroup species, number and frequency of synonymous and non-synonymous SNPs (SFS),number of synonymous and non-synonymous sites in the divergence and polymorphism data sets.
These are passed via a DoFE file as described in the documentation of the DoFE program:
[http://www.lifesci.susx.ac.uk/home/Adam_Eyre-Walker/Website/Software.html][1].
**grapes** can handle both folded (as in DoFE) and unfolded SFS data. If SFS's are unfolded, this must be
indicated via an extra line in the input file containing: #unfolded . See examples at the bottom of this
file.

**multi_grapes** needs information on the number and frequency of synonymous and non-synonymous SNPs
(SFS) and number of synonymous and non-synonymous sites, for one or several data sets. These are
passed via a DoFE file as described in the documentation of the DoFE program:
[http://www.lifesci.susx.ac.uk/home/Adam_Eyre-Walker/Website/Software.html][1].
**multi_grapes** can handle both folded (as in DoFE) and unfolded SFS data. If SFS's are unfolded, this must
be indicated via an extra line in the input file containing: #unfolded. Multiple data sets should appear as
different lines of the DoFE file (see examples at the bottom of this file)

Analysis
--------
**grapes** will estimate the distribution of fitness effect of mutations (DFE), rate of adaptive evolution
(<img src="https://render.githubusercontent.com/render/math?math=\omega_a">), rate of non-adaptive evolution (<img src="https://render.githubusercontent.com/render/math?math=\omega_{na}">), and proportion of adaptive substitutions (a) by fitting a
population genetic model to SFS + divergence data in the maximum likelihood framework. grapes will
also perform more basic analyses, namely estimating a as <img src="https://render.githubusercontent.com/render/math?math=1-\frac{p_N / p_S}{d_N / d_S}">, referred to as Neutral
model, and using the corrected version of Fay, Wickoff and Wu (2002 Nature 415:1024), referred to as
FWW.

**multi_grapes** will estimate the distribution of fitness effect of mutations (DFE) by fitting a population
genetic model to SFS + divergence data in the maximum likelihood framework. Data sets will first be
analyzed separately (= with all parameters specific to each data set, as in grapes) then jointly (= with one parameter, the shape of the DFE, shared among data sets). 

Output
------
Basic output is written in the terminal. This corresponds to estimates of the adaptive and non-adaptive
rates, main model parameters and likelihoods. Detailed output is written in a **csv** file.

Grapes options
--------------

The most important option is the name of the assumed DFE model. Five distinct models are
implemented in addition to the Neutral model, namely `GammaZero` (=Gamma), `GammaExpo`, `DisplGamma`, `ScaledBeta`, and `FGMBesselK`.
See Galtier (2016) for details on what these models mean. One can either use one of the six models (e.g., `-model GammaExpo`), or do the six in a single
run by passing `-model all`. `GammaZero` with all default options is equivalent to the second approach
of Eyre-Walker and Keightley (2009), with nuisance parameters <img src="https://render.githubusercontent.com/render/math?math=r_i">'s.

### Additional options:

*`-nearly_neutral <float>`: defines the threshold of $Ne.s$ above which a mutation is considered
adaptive ($S_{adv}$ in Galtier 2016, default=5)
*`-FWW_threshold <float>`: minimal allele frequency in FWW alpha estimation (default=0.15)
*`-no_div_param`: implements so-called [-A] version in Galtier (2016); also called a DFE by Tataru et
al. (2017), and indicated with * in Rousselle et al. (2018); irrelevant if model = GammaZero;
(default=false)
*`-no_syn_orient_error`: force equal synonymous and non-synonymous mis-orientation rate(default=false)
*`-nb_rand_start <int>`: number of random starting values in model optimization (default=0);
setting positive values will slow down the program but decrease the probability of being trapped in
local optima.
*`-fixed_param <control_file_name>`: this option should be used if one does not want to
optimize every parameter, but rather have some parameters fixed to predefined values; parameter
names and predefined values are passed via a control file; see example at the bottom of this file
(default=none).

Multigrapes options
-------------------

The most important option is the name of the assumed DFE model. Two distinct models are implemented
in addition to the Neutral model, namely `GammaZero` (=Gamma) and `ReflectedGamma`. See Galtier
(2016) and Galtier & Rousselle (2020) for details on what these models mean. In addition, one can add a
class of lethal mutations using the `-p_lethal` option (see below).

### Additional options:

* `-fold` : fold the SFS (if unfolded)
* `-no_syn_orient_error` : force equal synonymous and non-synonymous mis-orientation rate (default=false)
* `-nb_rand_start <int>` : number of random starting values in model optimization (default=0); setting
positive values will slow down the program but decrease the probability of being trapped in local optima.
* `-no_separate` : do not perform separate analysis.
* `-no_shared` : do not perform shared analysis.
* `-shared_shape <int>,<int>,...`: fixed shared shape parameters (all the listed values will be fitted
successively); in the absence of this option, the shared shape parameter will be optimized.
* `-p_lethal <float>` : proportion of lethal mutations (this parameter cannot be optimized).

Developed by
------------

Nicolas Galtier
Institute of Evolutionary Sciences
CNRS – University Montpellier

Supported by
------------

Montpellier Bioinformatics & Biodiversity platform
[http://mbb.univ-montp2.fr]

References
----------

* Eyre-Walker A, Keightley PD. 2009 Estimating the rate of adaptive molecular evolution in the presence
of slightly deleterious mutations and population size change. Molecular Biology and Evolution.
26:2097–2108.
* Galtier N. 2016. Adaptive protein evolution in animals and the effective population size hypothesis.
PLoS Genetics 12:e1005774.
* Rousselle M., Mollion M., Nabholz B., Bataillon T., Galtier N. 2018. Overestimation of the adaptive
substitution rate in fluctuating populations. Biology Letters (in press)
* Tataru P, Mollion M, Glémin S, Bataillon T. 2017 Inference of distribution of fitness effects and
proportion of adaptive substitutions from polymorphism data. Genetics 207:1103–1119.
* Galtier N., Rousselle M. 2020. How much does Ne vary among species? biorxiv
[https://doi.org/10.1101/861849]

Example input files for Grapes:
-------------------------------

### folded DoFE file:

```
Allolobophora_chlorotica+Aporerctodea_icterica (4852 genes)
all_genes 12 254286 4360.47 1987.61 1339.62 1640.54 1087.89 581.969 92842.3 18049.3 9568.58 6722.38 6820.24 5353.69 2607.35 340000 10578 111000 24890
```
First line is a header/comment line.
Second line contains the data:
* entry 1 ( all_genes ): any string (dataset description)
* entry 2 ( 12 ): sample size (here, 6 diploid individuals)
* entry 3 ( 254286 ): number of non-synonymous sites, polymorphism data
* entry 4 → 3+n/2, where n equals sample size: non-synonymous SFS (number of non-synonymous
singletons, doubletons, etc...)
* entry 3+n/2+1 ( 92842.3 ): number of synonymous sites, polymorphism data
* entry 3+n/2+2 → 3+n/2+1+n/2, where n equals sample size: synonymous SFS (number of synonymous
singletons, doubletons, etc...)
* entry 3+n/2+1+n/2+1 ( 340000 ): number of non-synonymous sites, divergence data
* entry 3+n/2+1+n/2+2 ( 10578 ): number of non-synonymous substitutions, divergence data
* entry 3+n/2+1+n/2+3 ( 111000 ): number of synonymous sites, divergence data
* entry 3+n/2+1+n/2+4 ( 24890 ): number of synonymous substitutions, divergence data

### unfolded DoFE file:

```
Microtus_arvalis+Microtus_glareolus.fas (2943 genes)
#unfolded
all_genes 12 1.61188e+06 2008.19 590.033 165.33 109.824 85.956 73.6923 61.4066 63.5934 65.3077 104.637 119.604 361114 5046.4 2095.86 836.505 533.066 403.44 357.615 327.066 318.044 322.549 518.198 757.527 1.48228e+06 13089 432574 38940
```
Same as above with additional #unfolded line, and SFS's containing n-1 entries instead of n/2.                                 

Example control file:
---------------------

```
GammaExpo,negGshape,0.4
GammaExpo,negGmean,10000
```
Each line contains model name, parameter name and fixed parameter value.
Parameter names are listed in the table below:

Model | Parameter name | Parameter meaning
------|----------------|------------------
GammaZero | negGmean | Gamma distribution mean, negative effects
GammaZero | negGshape | Gamma distribution shape, negative effects
GammaExpo | negGmean | Gamma distribution mean, negative effects
GammaExpo | negGshape | Gamma distribution shape, negative effects
GammaExpo | posGmean | Exponential distribution mean, positive effects
GammaExpo | pos_prop | Proportion of beneficial mutations (positive s)
DisplGamma | negGmean | Un-displaced Gamma distribution mean
DisplGamma | negGshape | Un-displaced Gamma distribution shape
DisplGamma | s0 | Displacement
ScaledBeta | Ba | Beta distribution first shape parameter
ScaledBeta | Bb | Beta distribution second shape parameter
ScaledBeta | med_prop | Proportion of mutations with Ne.s above -25
FGMBesselK | BKm | m in equation 8, Lourenco et al 2011*
FGMBesselK | BKnsz2 | n/z 2 in equation 8, Lourenco et al 2011*
FGMBesselK | BKsigma | s in equation 8, Lourenco et al 2011*
FGMBesselK | BKscale | scaling parameter, FGM Bessel distribution

*Lourenco et al. 2011 Genetics 65:1559
