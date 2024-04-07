ods trace on;
ods trace off;

%let sasroot=C:\Users\jwq\Desktop\mph\1_k\version2code;

******Setting for run R function. Only run when MVTAPP=R***************************;
%let rroot=C:\Program Files\R\R-3.5.3;
options set=R_HOME="&rroot.";
option  SPOOL;
%macro RForQuantiles;
%include "&sasroot.\MVT_Quantiles.sas";**submit R function to calculate the quantiles in multivar normal distrubution;
%mend;
%macro RForProbabilities;
%include "&sasroot.\MVT_Probabilities.sas";**submit R function to calculate the quantiles in multivar normal distrubution;
%mend;
proc options option=rlang;run;
******End***************************************************************************;

******Setting for SAS IML code by Frank Bretz. Only run when MVTAPP=IML************;
%inc "&sasroot.\probmvt.sas";
******END***************************************************************************;

%inc "&sasroot.\mi1k.sas";

data aax;
 set ana818.HB818_ana818logous;
 if _n_/5=int( _n_/5);
 ID=_n_;
 if rannor(123)>0.5 then Sex="F";else sex="M";
 AGE=48+rannor(123)*18;
 put TRT ID SEX ESASTT AGE HB0 HB1 HB2 HB3 HB4 HB5 HB6 ;
run;

data sample;
input TRT $ ID SEX $ ESASTT $ AGE HB0 HB1 HB2 HB3 HB4 HB5 HB6;
cards;
BIW 10 M ESA-Naive 75 115 111 113 113 124 . .
BIW 20 M ESA-Naive 37 111.5 101 102 102 104 118 .
BIW 30 F ESA-Treated 36 117 118 120 120 114 . .
BIW 40 M ESA-Treated 34 122.5 . . . 123 104 119
BIW 50 F ESA-Treated 70 113 103 96 96 96 101 104
QW 60 M ESA-Naive 40 110.5 108 91 97 91 89 94
QW 70 F ESA-Naive 24 108 102 102 104 105 115 111
QW 80 M ESA-Treated 40 114.5 115 115 110 115 107 98
QW 90 M ESA-Treated 29 123 122 114 101 114 125 106
QW 100 F ESA-Treated 39 106.5 112 102 97 105 117 115
TIW 110 F ESA-Naive 61 118 116 117 115 120 99 101
TIW 120 F ESA-Naive 30 121 116 119 116 . 104 105
TIW 130 M ESA-Treated 53 120.5 113 103 107 104 105 102
TIW 140 M ESA-Treated 35 124 124 112 119 121 110 111
TIW 150 M ESA-Treated 72  113 105 111 112 99 99 .
TIW 160 M ESA-Treated 43 117 112 109 110 116 109 .
;
run;

PROC MI DATA= sample SEED=125636 NIMPUTE=10 out=sample_mcmc(rename=(HB0=BasewHB)); 
VAR HB0 HB1 HB2 HB3 HB4 HB5 HB6;
EM MAXITER=1000;
MCMC CHAIN=MULTIPLE IMPUTE=Full;
RUN;

***Sample for calculate LS mean difference, using R functions****;
%MI1K(IMPDAT=sample_mcmc,TRT=TRT,CONTROL=TIW,MEANTYPE=LSMEAN,RESVAR=HB6,CATVAR=SEX ESASTT,MODVAR=BasewHB AGE ESASTT SEX,
MVTAPP=R,ALPHA=0.05,LIMIT=BOTH,OUTDAT=res_sample);

***Sample for calculate LS mean difference, using SAS/IML code by Dr Bretz, test for Difference=-20****;
%MI1K(IMPDAT=sample_mcmc,TRT=TRT,CONTROL=TIW,MEANTYPE=LSMEAN,RESVAR=HB6,CATVAR=SEX ESASTT,MODVAR=BasewHB AGE ESASTT SEX,
MVTAPP=IML,ALPHA=0.05,TESTAT=-20,LIMIT=BOTH,OUTDAT=res_sample);

***Sample for calculate mean and LSM using SAS/IML code by Dr Bretz, 1-sided, upper limit****;
%MI1K(IMPDAT=sample_mcmc,TRT=TRT,CONTROL=TIW,MEANTYPE=BOTH,RESVAR=HB6,CATVAR=SEX ESASTT,MODVAR=BasewHB AGE ESASTT SEX,
MVTAPP=IML,ALPHA=0.025,LIMIT=LOWER,OUTDAT=res_sample);
