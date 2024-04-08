/* SAS Macro for the calculation of many-to-one comparison of arithmetic means or least-squared (LS) means 
   based on a multiple  imputed dataset created from PROC MI.

   Author : Wenqing Jiang
   Contact: 2111220022@bjmu.edu.cn; bjmujwq@163.com 
   Date   : 06 Sept 2023 - start version

   Arguments  : IMPDAT:   a data created by PROC MI with variable _IMPUTATION_ indicating the imputation ID.
                TRT    :  a numerical or categorical variable indicating the treatment allocation in the IMPDATA.
                CONTROL:  a number or a text string indicating the level of control's TRT in the IMPDATA.
	            MEANTYPE: the type of mean difference to calculate: 
	                      ARIMEAN = arithmetic  mean difference
	                      LSMEAN = LS mean difference;
                          BOTH = both arithmetic  and LS mean differences
                          The default string is BOTH.
	            RESVAR:   a numerical variable indicating the response value to be compared in the IMPDATA.
	            CATVAR:   categorical variable(s) to be adjusted in the general linear model(GLM) for LS mean difference in the IMPDATA. This argument is not needed  when MEANTYPE = ARIMEAN.
                MODVAR:   a statement of variable(s) except the TRT variable to build the GLM for LS means, that is, in PROC GLM, MODEL RESVAR=TRT MODVAR
                          However, please note this macro does not support any Nesting or Interaction Effects that contain TRT. Put TRT aside, users
                          can build the GLM with other nesting or interaction effects. This argument is not needed when MEANTYPE = ARIMEAN.
	            MVTAPP:   approach to solve  the probability and quantile of multi-variate t (MVT) distribution:
	                      IML = the SAS/IML code for MVT probability written by FRANK BRETZ and a SAS built-in root finding function for the quantile;
	                      R = call R package 'mvtnorm' written by ALAN GENZ. when specifying MVTAPP=R, user should install a right version of R and 
	                      'mvtnorm' package thus  SAS/IML can call R function. Default string is IML.
                ALPHA:    Overall significant level. Default value is 0.05
                TESTAT::  The test value for null hypothesis, that is, p-value is test for Diff = TESTAT. Default value is 0
                LIMIT:    Specify the limit(s) to calculate:
                          UPPER = 1-sided upper limit
                          LOWER = 1-sided lower limit
                          BOTH = 2-sided limits
                          Default string is BOTH.
	            OUTDAT:   a data to store the results, including the point estimation, interval (limit) estimation and their p-values of mean difference
                          and LS mean difference. When OUTDAT is skipped, the macro will print the results only  */



%MACRO MI1K(IMPDAT=,TRT=,CONTROL=,MEANTYPE=BOTH,RESVAR=,CATVAR=,MODVAR=,MVTAPP=IML,ALPHA=0.05,TESTAT=0,LIMIT=BOTH,OUTDAT=);
		data indata;set &impdat;obsnum=_n_;
        run;
		data null_;
        set sashelp.vcolumn end=eof; where libname="WORK" and memname="INDATA";if name=upcase("&trt.") or name = lowcase("&trt.") then call symputx("trttp",type);
		run;
		 *Create a standard treatment variable Standard_T;
		 %if &trttp^=char %then %do; data indata2;length Standard_T $50.; set indata; Standard_T=strip(put(&trt.,best.));run; %end;
		 %else %do; data indata2;length Standard_T $50.; set indata; Standard_T=&trt.;run; %end;

		ods output OneWayFreqs=OneWayFreqs;
		proc freq data=indata2;Table Standard_T ;run;

		data nonctrl_Standard_T;set OneWayFreqs;where Standard_T not in ("" "&CONTROL.") ;order_Standard_T=_n_;keep Standard_T order_Standard_T;run;
	    data ctrl_Standard_T;set OneWayFreqs;where Standard_T="&CONTROL." ;order_Standard_T=99;keep Standard_T order_Standard_T;run;
		data Levels_Standard_T;set ctrl_Standard_T nonctrl_Standard_T;run;*order of trt;

		proc sql noprint;
		  create table catord_trt as select a.obsnum, b.order_Standard_T from indata2 as a inner join Levels_Standard_T as b on a.Standard_T=b.Standard_T;*give order of TRT;
		  select strip(put(max(order_Standard_T),best.)) into: k from nonctrl_Standard_T; *treatment group number;
		  select strip(put(count(*),best.)) into: samplesize from indata2 where _Imputation_=1; *treatment group number;
		  select strip(put(max(_Imputation_),best.)) into: m from indata2; *imputation time;
		  select Standard_T into:trtname1 -:trtname&k from nonctrl_Standard_T; *name of treatment group;
		  select Standard_T into:ctrlname from ctrl_Standard_T; *name of treatment group;
		  create table indata3 as select a.*,b.order_Standard_T as Standard_TN from indata2 as a inner join catord_trt as b on a.obsnum=b.obsnum order by _Imputation_,Standard_TN;
		  select strip(put(count(*),best.)) into: n1 -: n&k  from indata3 where _Imputation_=1 and Standard_TN<99 group by Standard_TN; *samplesize of treatment groups;
		  select strip(put(count(*),best.)) into: n99 from indata3 where _Imputation_=1 and Standard_TN=99; *samplesize of treatment groups;
		quit;

		%if &MEANTYPE=ARIMEAN or &MEANTYPE=BOTH %then %do;
		ods output summary=TRTMEAN;
		proc means data=indata3;
		 by _Imputation_ Standard_TN;
		 var &RESVAR.;
		run;
 
		*Obtain MSE for arithimatic mean difference;
		ods output OverallANOVA=OANCOVA;
		proc anova data=indata3 PLOTS=NONE;
		by _Imputation_ ;
		class Standard_TN;
		model &RESVAR.=Standard_TN;
		run;
        %end;

		%if &MEANTYPE=LSMEAN or &MEANTYPE=BOTH %then %do;
		*Obtain Estimates, MSE, and Inverse of XPX for LS mean difference;
		ods output ParameterEstimates=PE OverallANOVA=OA InvXPX=InvXPX; *output XPX, MSE;
		proc glm data=indata3 plots=none;
		  by _Imputation_;
		  CLASS Standard_TN &CATVAR.;
		  MODEL &RESVAR. = Standard_TN &MODVAR./xpx solution INVERSE;
		quit;

		data invxpx_;
		 set invxpx;
		 if index(parameter,"Dummy") then do; paramid=input(substr(parameter,6),best.);end;
		 if 0<paramid<=&k.+1 then output;
		run;

		proc transpose data=invxpx_ out=trans_inv(drop=_NAME_);
		where  _imputation_=1;
		var parameter;
		run;

		data _null_;
		set trans_inv;
		combo = catx(',', of col:);
		call symputx("dummys",combo);
		run;
        %end;
		PROC IML;
		%if &MVTAPP=IML %then %do;
		load module = _all_;
		load _all_;
		%end;
		    L=j(&k,&k+1,0); L[1:&k,1:&k]=i(&k); L[1:&k,&k+1]=repeat(-1,&k,1);
			prob=1-&ALPHA.;
			%if &LIMIT=UPPER or &LIMIT=LOWER %then %do;
			prob=1-2*&ALPHA.;
			%end;
			START CALMAT(Q,U) GLOBAL(QBA,UBA,BM,TM, RM,V,R,PROB,TIV,RESMAT);
				 SUMQ=J(&K,1,0); 
		         SUMU=J(&K,&K,0);
				 DO i=1 TO &m;
				   SUMQ=SUMQ+Q[,i]; 
		           SUMU=SUMU+U[,(i-1)*&k+1 :i*&k];
				 END;
				 QBA=SUMQ/&m; 
		         UBA=SUMU/&m;
		         SUMBM=J(&K,&K,0);
				 DO i=1 TO &m;
				   deltaQ=Q[,i]-QBA; SUMBM=SUMBM+ deltaQ*deltaQ`;
				 END;
				 BM=SUMBM/(&m-1);
				 TM=UBA+ (1+1/&m)*BM;
				 RM=((1+1/&m)/&k) * trace(BM*inv(UBA));
				 V=floor((&m-1)*(1+1/RM)**2);
				 R=inv(sqrt(diag(TM)))*Tm*inv(sqrt(diag(TM)));
				 AII=((vecdiag(TM))##0.5);
				 TESTC=J(&K,1,&TESTAT.);
		         TIV=abs(QBA-TESTC)/AII;
				 Print QBA TM TIV;
				 RESMAT=j(&k,5,0);


					%if &MVTAPP=R %then %do;%RForQuantiles;
     					call ExportMatrixToR(r, "r"); 
					call ExportMatrixToR(v, "v");
					call ExportMatrixToR(prob, "prob");
                                        %end;
                                        %else %if &MVTAPP=IML %then %do;
						Start Func(x)  global(R,v,prob);
								  DELTA = J(1,&k,0);
								  LOWER = J(1,&k,-x);
								  UPPER = J(1,&k,x);
								  INFIN = J(1,&k,2);
						          COVAR = R;
								  MAXPTS = 2000*&k*&k*&k;
								  ABSEPS = .001;
								  RELEPS = 0; 
								  RUN MVN_DIST(&k, v, DELTA, LOWER, UPPER, INFIN, COVAR, MAXPTS, ABSEPS, RELEPS,  ERROR, VALUE, NEVALS, INFORM );
								  y=VALUE-prob;	 
						Return(y);
						Finish;

					g = froot("Func", {0 5});
					%end;
					%if &MVTAPP=R %then %do; call ImportMatrixFromR(t2data, "t2data");  %end; 
					%do i=1 %to &k;
						TI=j(&k,1,TIV[&i,1]);
							%if &MVTAPP=R %then %do;
                                                        %RForProbabilities;
						        call ExportMatrixToR(ti, "ti"); 
	                                                %end;
							%else %if &MVTAPP=IML %then %do;
							TII=TIV[&i,1];
							DELTA = J(1,&k,0);
							LOWER = J(1,&k,-TII);
							UPPER = J(1,&k,TII);
							INFIN = J(1,&k,2);
							COVAR = R;
							MAXPTS = 2000*&k*&k*&k;
							ABSEPS = .001;
							RELEPS = 0;
                            RUN MVN_DIST(&k, v, DELTA, LOWER, UPPER, INFIN, COVAR, MAXPTS, ABSEPS, RELEPS,  ERROR, VALUE, NEVALS, INFORM );
                            t2pdata=VALUE;
							%end;
						%if &MVTAPP=R %then %do; call ImportMatrixFromR(t2pdata, "t2pdata"); %end;
						RESMAT[&i,1]=&i;
						RESMAT[&i,5]=1-t2pdata;
                    Print TI t2pdata;
		            %end;

				  RESMAT[,2]=QBA;
				  %if &MVTAPP=R %then %do;
								  RESMAT[,3]=QBA - t2data[1,1]*AII;
								  RESMAT[,4]=QBA + t2data[1,1]*AII;
				  %end;
				  %else %if &MVTAPP=IML %then %do;
								  RESMAT[,3]=QBA - g*AII;
								  RESMAT[,4]=QBA + g*AII;
				  %end;
			FINISH;

		**LS Means Case;
			%if &MEANTYPE=LSMEAN or &MEANTYPE=BOTH %then %do;
		    use PE; read all var {Estimate} where(parameter ? "Standard") into ESTTRT ; close PE;
			use OA; read all var {MS} where(Source ? "Error")  into MSE; close OA;
			use invxpx_; read all var {&dummys.} where(Parameter?"Dummy")  into INVXPX; close invxpx_;
			QLS=j(&k,&m,.);
			ULS=j(&k,&m*&k,.);
			do i=1 to &m; 
		        BETA_=ESTTRT[(i-1)*(&k+1)+1: i*(&k+1),];
				InvXPX_=InvXPX[(i-1)*(&k+1)+1: i*(&k+1), 1:&k+1 ];
			    MSE_=MSE[i,];
			    QLS[,i]=L * BETA_;
		        ULS[, (i-1)*&k+1 :i*&k ]=MSE_*L*InvXPX_*L`;
			end;
			RUN CALMAT(QLS,ULS);
		    Create RES_LS from RESMAT;append from RESMAT;
            %end;

		***ARIMEAN CASE;
		  %if &MEANTYPE=ARIMEAN or &MEANTYPE=BOTH %then %do;
		   use TRTMEAN ; read all var {&RESVAR._Mean} into TRTMEAN; 
		                 read all var {&RESVAR._StdDev} into TRTSTD;     close TRTMEAN;
		   use OANCOVA; read all var {MS} where(Source ? "Error")  into MSE; close OANCOVA;
		   QAM=j(&k,&m,0);
		   UAM=j(&k,&m*&k,0);
		   UAMEQV=j(&k,&m*&k,0);
		   do i=1 to &m; 
		        BETA_=TRTMEAN[(i-1)*(&k+1)+1: i*(&k+1),];
			    QAM[,i]=L * BETA_;
				DG=j(&k,1,0);DGEQV=j(&k,1,0);
				  %do j=1 %to &k;
		             DG[&j,1]=TRTSTD[(i-1)*(&k+1)+&j,1]##2/&&n&j;
					 DGEQV[&j,1]=MSE[i,1]/&&n&j;
				  %end;
		        UAM[, (i-1)*&k+1 :i*&k ]=j(&k,&k,TRTSTD[i*(&k+1),1]**2/&n99.)+DIAG(DG);
				UAMEQV[, (i-1)*&k+1 :i*&k ]=j(&k,&k,MSE[i,1]/&n99.)+DIAG(DGEQV);
			end;
			RUN CALMAT(QAM,UAM);Create RES_AM from RESMAT;append from RESMAT;
			RUN CALMAT(QAM,UAMEQV);Create RES_AMEQV from RESMAT;append from RESMAT;
		Quit;
        %end;

		%if &LIMIT=UPPER or &LIMIT=UPPER %then %do;
			prob=1-2*&ALPHA.;
			%end;

		Data RES &OUTDAT;
		 Retain Type Comparison col2 
              %if &LIMIT^=UPPER %then %do; col3 %end;
			  %if &LIMIT^=LOWER %then %do; col4 %end;
              col5;
		 length Type $100 Comparison $100;
		 Set   %if &MEANTYPE=ARIMEAN or &MEANTYPE=BOTH %then %do; RES_AM(in=a) RES_AMEQV(in=b) %end;
               %if &MEANTYPE=LSMEAN or &MEANTYPE=BOTH %then %do; RES_LS(in=c) %end;;
         %if &MEANTYPE=ARIMEAN or &MEANTYPE=BOTH %then %do;
         if a then Type="Arithimatic Mean (Heterogeneity Assumption)";
		 if b then Type="Arithimatic Mean (Homogeneity Assumption)";
		 %end;
         %if &MEANTYPE=LSMEAN or &MEANTYPE=BOTH %then %do; 
		 if c then Type="Least Squares Means";
		 %end;
		  %do i=1 %to &k;
		   if col1=&i then Comparison="&&trtname&i vs. &ctrlname";
		  %end;
		  Rename col2=Estimate 
			%if &LIMIT^=UPPER %then %do;  col3=Lower  %end;
		     %if &LIMIT^=LOWER %then %do;  col4=Upper  %end; 
          col5=pValue; 
		  Keep Type Comparison col2 
              %if &LIMIT^=UPPER %then %do; col3 %end;
			  %if &LIMIT^=LOWER %then %do; col4 %end;
         col5;
		run;

		Proc print data=RES noobs;run;
%MEND MI1K;
