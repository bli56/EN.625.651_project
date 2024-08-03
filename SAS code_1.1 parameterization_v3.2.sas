DM LOG 'CLEAR';
DM LST 'CLEAR';

/*** Data preparation ***/
/* Find subjects with definite El Escorial Diagnosis */
proc import out=proact_ELESC
	datafile= "D:\Baini_files\ZM\Studies\CB03-154-201\Sample size\Data\PRO-ACT database\PROACT_ELESCORIAL.csv"
	dbms=csv
	replace;
	getnames=YES;
quit;

data proact_ELESC;
	keep subject_id;
	set proact_ELESC;
	if el_escorial = 'Definite';
proc sort;
	by subject_id;
quit;


/* Find subjects with non-missing age values */
proc import out=proact_DM
	datafile= "D:\Baini_files\ZM\Studies\CB03-154-201\Sample size\Data\PRO-ACT database\PROACT_DEMOGRAPHICS.csv"
	dbms=csv
	replace;
	getnames=YES;
quit;

data proact_DM;
	keep subject_id Age;
	set proact_DM;
	if Age ^=.;
proc sort;
	by subject_id;
quit;

/* Find subjects with SVC>60% at the screening visit */
/* ATTENTION: this condition is not usedv for filtering because no active subjects were SVC>60% at the screening visit  */
/*
proc import out=proact_SVC
	datafile= "D:\Baini_files\ZM\Studies\CB03-154-201\Sample size\Data\PRO-ACT database\PROACT_SVC.csv"
	dbms=csv
	replace;
	getnames=YES;
quit;

data proact_SVC;
	keep subject_id pct_of_Normal_Trial_1 slow_vital_Capacity_Delta;
	set proact_SVC;
	if pct_of_Normal_Trial_1 ^=. and slow_vital_Capacity_Delta =< 0 and slow_vital_Capacity_Delta ^=.;
proc sort;
	by subject_id slow_vital_Capacity_Delta;
quit;

data proact_SVC;
	set proact_SVC;
	by subject_id slow_vital_Capacity_Delta;
	if last.subject_id;
run;
*/

/* Retrieve alocation info */
proc import out=proact_treat
	datafile= "D:\Baini_files\ZM\Studies\CB03-154-201\Sample size\Data\PRO-ACT database\PROACT_TREATMENT.csv"
	dbms=csv
	replace;
	getnames=YES;
quit;

data ls;
	keep subject_id Study_Arm;
	set proact_treat;
run;


/* Retrieve ALSFRS-R data  */
proc import out=proact_alsfrs
	datafile= "D:\Baini_files\ZM\Studies\CB03-154-201\Sample size\Data\PRO-ACT database\PROACT_ALSFRS.csv"
	dbms=csv
	replace;
	getnames=YES;
quit;

data proact_alsfrs;
	keep subject_id ALSFRS_R_Total ALSFRS_Delta;
	set proact_alsfrs;
proc sort;
	by subject_id ALSFRS_Delta;
quit;


/* Retrieve ALS history info */
proc import out=PROACT_ALSH
	datafile= "D:\Baini_files\ZM\Studies\CB03-154-201\Sample size\Data\PRO-ACT database\PROACT_ALSHISTORY.csv"
	dbms=csv
	replace;
	getnames=YES;
quit;

data PROACT_ALSH;
	keep subject_id Onset_Bulbar Site_of_Onset Onset_Delta;
	set PROACT_ALSH;
	if Site_of_Onset = '' then Onset_Bulbar =.;  /* This derived info is no used for now */
		else if Site_of_Onset ^= '' then do;
			if Site_of_Onset = 'Onset: Bulbar' then Onset_Bulbar = 1;
			else Onset_Bulbar = 0;
		end;
	if Onset_Delta ^=.;
proc sort;
	by subject_id Onset_Delta;
quit;

data Proact_alsh2;
	set Proact_alsh;
	by subject_id Onset_Delta;
	if first.subject_id;
run;

/* Combine */
data alsfrs0;
	merge proact_alsfrs(in=a) ls(in=b);
	by subject_id;
	if a and b;
run;


data alsfrs1;
	merge alsfrs0(in=a) proact_ELESC(in=b);
	by subject_id;
	if a and b;
run; 

data alsfrs2;
	merge alsfrs1(in=a) proact_DM(in=b);
	by subject_id;
	if a and b;
run; 

/*
data alsfrs3;
	merge alsfrs2(in=a) proact_SVC(in=b);
	by subject_id;
	if a and b;
run;
*/

data alsfrs3;
	merge alsfrs2(in=a) PROACT_ALSH(in=b);
	by subject_id;
	if a and b;
run;

data alsfrs4;
	set alsfrs3;
	if ALSFRS_R_Total ^=. and Study_Arm ^= '' and ALSFRS_Delta ^=. and Onset_Delta ^=. and age ^=.;
	if abs(Onset_Delta)>. and abs(Onset_Delta) < 540;  /* Subjects who are < 540 days from symptom onset  */
	if ALSFRS_Delta >=0 and ALSFRS_Delta =< 334;
	if ALSFRS_Delta = 0 then do; week = 0; month = 0; weekflag = 0; end;  /* baseline. NOTE: In protocol the baseline is day 1  */
		else if ALSFRS_Delta >= 4  and ALSFRS_Delta =< 10  then do; week = 2; month = 0.5; weekflag = 2; end;
		else if ALSFRS_Delta >= 16  and ALSFRS_Delta =< 26  then do; week = 4; month = 1; weekflag = 2; end;
		else if ALSFRS_Delta >= 44  and ALSFRS_Delta =< 54  then do; week = 8; month = 2; weekflag = 2; end;
		else if ALSFRS_Delta >= 72  and ALSFRS_Delta =< 82  then do; week = 12; month = 3; weekflag = 2; end;
		else if ALSFRS_Delta >= 100  and ALSFRS_Delta =< 110  then do; week = 16; month = 4; weekflag = 2; end;
		else if ALSFRS_Delta >= 128  and ALSFRS_Delta =< 138  then do; week = 20; month = 5; weekflag = 2; end;
		else if ALSFRS_Delta >= 156  and ALSFRS_Delta =< 166  then do; week = 24; month = 6; weekflag = 2; end;

		else if ALSFRS_Delta >= 184  and ALSFRS_Delta =< 194  then do; week = 28; month = 7; weekflag = 2; end;
		else if ALSFRS_Delta >= 212  and ALSFRS_Delta =< 222  then do; week = 32; month = 8; weekflag = 2; end;
		else if ALSFRS_Delta >= 240  and ALSFRS_Delta =< 250  then do; week = 36; month = 9; weekflag = 2; end;
		else if ALSFRS_Delta >= 268  and ALSFRS_Delta =< 278  then do; week = 40; month = 10; weekflag = 2; end;
		else if ALSFRS_Delta >= 296  and ALSFRS_Delta =< 306  then do; week = 44; month = 11; weekflag = 2; end;
		else if ALSFRS_Delta >= 324  and ALSFRS_Delta =< 334  then do; week = 48; month = 12; weekflag = 2; end;
	if week ^=.;
proc sort;
	by subject_id ALSFRS_Delta;
quit;

/* Delete the subjects that having only one record */
proc sql;
	create table eff_data1 as select * from alsfrs4
	group by subject_id
	having count(subject_id)>1;
quit;

/* Calculate del-FS */
data del_FS;
	keep subject_id del_FS;
	set eff_data1;
	del_FS = ALSFRS_R_Total/abs(Onset_Delta);
	if ALSFRS_Delta = 0;
run;

data eff_data2;
	merge eff_data1(in=a) del_FS(in=b);
	by subject_id;
	if a and b;
proc sort;
	by subject_id ALSFRS_Delta;
quit;


/* Check availability of data */
%macro freq_chk(name=);
	proc sort data=&name out=chk_subject nodupkey;
		by subject_id;
	quit;
	title1 "Check subjects";
	proc freq data=chk_subject;
		table study_Arm;
	quit;
	proc delete data = chk_subject;
	quit;

	proc sort data=&name out=chk_obs;
		by subject_id alsfrs_delta;
	quit;
	title1 "Check observations";
	proc freq data=chk_obs;
		table study_Arm;
	quit;
	proc delete data = chk_obs;
	quit;
	title1;
%mend;

/*
%freq_chk(name=Alsfrs0);
%freq_chk(name=Alsfrs1);
%freq_chk(name=Alsfrs2);
%freq_chk(name=Alsfrs3);
%freq_chk(name=Alsfrs4);
*/
%freq_chk(name=eff_data2);


/* delete datasets of no use */
%macro deletdat(library=, name=);
	proc delete data = &library..&name;
	quit;
%mend;

%deletdat(library=work, name=Proact_alsfrs);
%deletdat(library=work, name=Proact_alsh);
%deletdat(library=work, name=Proact_alsh2);
%deletdat(library=work, name=Proact_dm);
%deletdat(library=work, name=Proact_elesc);
/* %deletdat(library=work, name=Proact_svc); */
%deletdat(library=work, name=Proact_treat);

%deletdat(library=work, name=Alsfrs0);
%deletdat(library=work, name=Alsfrs1);
%deletdat(library=work, name=Alsfrs2);
%deletdat(library=work, name=Alsfrs3);
%deletdat(library=work, name=Alsfrs4);
%deletdat(library=work, name=eff_data1);
%deletdat(library=work, name=Del_fs);
%deletdat(library=work, name=ls);


/**************************************************************************/
/***** Extract info to assist the parameterization of data simulation *****/
/**************************************************************************/

proc import out=proact_die
	datafile= "D:\Baini_files\ZM\Studies\CB03-154-201\Sample size\Data\PRO-ACT database\PROACT_DEATHDATA.csv"
	dbms=csv
	replace;
	getnames=YES;
quit;

proc sort data=proact_die;
	by subject_id;
quit;

data Eff_die;
	merge Eff_data2(in=a) Proact_die(in=b);
	by subject_id;
	if a;
run;

proc sort data=Eff_die out=Eff_die2 (keep=subject_id Study_Arm Death_Days) nodupkey;
	by subject_id;
quit;


/*** 1) Without considering the length of follow-up time and censoring time window (e.g., foltime=anc.cens=beta0.cens= 1000 for simple.surv.sim() ) ***/
proc means data=Eff_die2 noprint;
	var Death_Days;
	output out=stat_die mean= STD= STDERR = mode= P5= P25= P50= P75= P95= /autoname;
quit;
/* Mean(Death_Days) = 332.36363636, so 'beta0.ev' of simple.surv.sim() = log(332.36363636) = 5.80623 */
/* SD(Death_Days) = 89.944364117, SE(Death_Days) = 15.657304072 */

/*** 2) Apply the censoring time window: 0~(335-1) ***/
data Eff_die3;
	set Eff_die2;
	if .<death_days and death_days =< 334 then do;
		deathD = death_days;
		deathFL = 1;
		CNSR = 0;
	end;
	else do;
		deathD = 334;
		deathFL = 0;
		CNSR = 1;
	end;
proc sort;
	by Study_Arm deathFL subject_id;
quit;


data Eff_die3b;
	set Eff_die3;
	if deathFL = 1;
run;

proc means data=Eff_die3b noprint;
	var deathD;
	output out=stat_die2 mean= STD= STDERR = mode= P5= P25= P50= P75= P95= /autoname;
quit;
/* Mean(deathD) = 260.4375, so 'beta0.ev' of simple.surv.sim() = log(260.4375) = 5.562363 */
/* SD(deathD) = 68.324684412, SE(deathD) = 17.081171103 */


/* Pooled freq */
title1;
proc freq data=Eff_die3;
	table deathFL /out=freq_die;
quit;

proc freq data=Eff_die3;
	table deathFL /out=freq_die2;
	by Study_Arm;
quit;
/* In filtered PRAO-ACT data, considering the censoring time window [0, 330+5-1], death occurred in 1 subject (1/22=4.545455%) who received active drug, and in 15 subjects (15/149=10.06711%) who received placebo. So, the overall death rate is 16/171 = 9.356725% */



/*** Shared-baseline mixed model #3 (Model used by Paganoni et al. 2020 (also about AMX-0035): "A shared-baseline, mixed-effects regression model was used, with no added covariates.") ***/

%let time_var = month;
/* Specify effect size */
%let effect_size = 0.3;

/*       t --> week or month  */
/*     Z_i --> Study_Arm      */
/*   Age_i --> age            */
/* DelFS_i --> del_FS         */

title1;
proc mixed data=eff_data2 method=reml alpha=0.05;
	class Study_Arm(ref='Placebo') subject_id;  /*** Note! Adding time variable at class statement (just like proc mixed for RCBD design) would lead to failure in converge ***/
	model ALSFRS_R_Total = &time_var  study_Arm*&time_var / s cl solution residual outp=outp;
	random intercept &time_var  / subject=subject_id type=un solution;
	store store_mixed;
	contrast 'Contrast: Arm(1) - Arm(0)' study_Arm*&time_var 1 -1;
	estimate 'Estimate: Arm(1) - Arm(0)' study_Arm*&time_var 1 -1;
	ods output  tests3=t3
				SolutionF=Solu_Fix
				SolutionR=Solu_Rand
				contrasts=contrasts
				estimates=estimates;
quit;

/* Fixed effect "Intercept" is 37.5393 with SE = 0.4376 */
/* Fixed effect "month" is -1.2342 with SE = 0.06387 */
/* Fixed effect "month*Study_Arm" is 0.1992 with SE = 0.1951 */

proc plm restore = store_mixed;
	estimate 'Estimate: Arm(1) - Arm(0)' study_Arm*&time_var 1 -1 /testvalue=0.3 upper cl alpha=0.05;
quit;


/*** Calculate the mean and SD for random opponents ***/
proc sort data=Solu_Rand;
	by effect subject_id;
quit;

proc means data=Solu_Rand noprint;
	by effect;
	var estimate;
	output out=dis_rand n=n mean= STD= STDERR= min= max= /autoname;
quit;

/* Randm effect "Intercept"：mean=0, SD = 5.2072 */
/* Randm effect "Month"：mean=0, SD = 0.6952 */


/*** Calculate the mean and SD for residual ***/
proc means data=outp noprint;
	var resid;
	output out=stat_resid mean= STD= STDERR= mode= P5= P25= P50= P75= P95= /autoname;
quit;

proc datasets lib=work memtype=data noprint;
	modify stat_resid; 
	attrib _all_ label='';
quit;

/* resid_mean = 0, and resid_SD = 1.6801010817, so individual-level random error per each visit e ~ rnorm(1, mean=0, sd=1.6801010817) */

/*** Power analysis (old method): Obtaining subject-specific parameter estimates ***/
data pw_F;
	set t3 contrasts;
	alpha = 0.05; /* F-test is one-sided */
	Noncen = NumDF*FValue;
	FCrit = finv(1-Alpha, NumDF, DenDF, 0);
	Power = 1 - probf(FCrit, NumDF, DenDF, Noncen);
title1 "Power: F-test, contrast";
proc print;
quit;
title1;

data pw_t;
	set estimates;
	do alpha = 0.1;
		do effectsize = 0.3; /* Attention! this is already greater than the mean estimate of slope. We shall determine the effect size after the interim analysis. */
			tcrit = tinv(1-alpha/2, df); /* Returns the two-sided quantile from t-distribution */
			noncen = effectsize/stderr;
			power = sum(1, -cdf("t", tcrit,df,noncen),cdf("t",-tcrit,df,noncen));
			output;
		end;
	end;
run;

title1 "Power: Estimate (t-test)";
proc print data= pw_t;
quit;
title1 "";


/*** Power analysis (new method): Obtaining subject-specific parameter estimates ***/
/*** Obtain subject-specific parameter estimates ***/
proc sql noprint;
	select Estimate into : fix_incpt
		from Solu_fix
		where Effect = 'Intercept';
	select Estimate into : fix_act
		from Solu_fix
		where Effect = "&time_var.*Study_Arm" and Study_Arm = 'Active';
	select Estimate into : fix_time
		from Solu_fix
		where Effect = "&time_var";

quit;

%put &fix_incpt;
%put &fix_act;
%put &fix_time (per &time_var);

proc sort data=eff_data2 (keep=subject_id Study_Arm) out=ls_subj nodupkey;
	by subject_id Study_Arm;
quit;

proc sort data=solu_rand;
	by subject_id;
quit;

data solu_rand2;
	merge solu_rand(in=a) ls_subj(in=b);
	by subject_id;
	if a;
proc sort data=solu_rand2;
	by subject_id Effect;
quit;

data solu_comb;
	set solu_rand2;
	if Effect = 'Intercept' then  Estim_comb = Estimate + &fix_incpt;
		else if Effect = "&time_var" and Study_Arm = 'Active' then Estim_comb = Estimate + &fix_time + &fix_act;
		else if Effect = "&time_var" and Study_Arm = 'Placebo' then Estim_comb = Estimate + &fix_time;
run;


/* Obtain pooled CI for the COMBINED effect */
proc sort data=solu_comb;
	by effect subject_id;
quit;

proc means data=solu_comb noprint;
	var Estim_comb;
	by effect;
	output out=stat_slop_comb n=n mean=mean std=std lclm=lclm uclm=uclm;
quit;

/* Combined intercept: mean=37.5393, SD=5.3266292574 */
/* Combined slope：mean=-1.20857193, SD=0.5914529297 */


/* Obtain CI for the COMBINED effect by arm */
proc sort data=solu_comb;
	by Study_Arm effect subject_id;
quit;

proc means data=solu_comb noprint;
	var Estim_comb;
	by study_arm effect;
	output out=stat_slop_comb2 n=n mean=mean std=std  lclm=lclm uclm=uclm;
quit;

/* Active: */
/* Combined intercept: mean=38.443803358, SD=3.6198544592 */
/* Combined slope：mean=-0.999668892, SD=0.3885615451 */

/* Placebo: */
/* Combined intercept: mean=37.405749169, SD=5.5309814647 */
/* Combined slope：mean=-1.239416674, SD=0.6106951407 */


/* Obtain pooled CI for the RANDOM effect */
proc sort data=solu_comb;
	by effect subject_id;
quit;
proc means data=solu_comb noprint;
	var estimate;
	by effect;
	output out=stat_slop_rand n=n mean=mean std=std lclm=lclm uclm=uclm;
quit;

/* Random intercept: mean=0, SD=5.3266 */
/* Random slope：mean=0, SD=0.5861 */


/* Obtain CI for the RANDOM effect by arm */
proc sort data=solu_comb;
	by study_arm effect subject_id;
quit;
proc means data=solu_comb noprint;
	var estimate;
	by study_arm effect;
	output out=stat_slop_rand2 n=n mean=mean std=std lclm=lclm uclm=uclm;
quit;

/* Active: */
/* Random intercept: mean=0.9045, SD=3.6199 */
/* Random slope：mean=0.03533, SD=0.3886 */

/* Placebo: */
/* Random intercept: mean=-0.1336, SD=5.5310 */
/* Random slope：mean=-0.00522, SD=0.6107 */


/* Calculate power (sub-method 1) */
proc sql noprint;
	
	/* Pooled SD */
	select std into : slope_time_SD_pooled
		from Stat_slop_comb
		where Effect = "&time_var";
	select std into : intcpt_SD_pooled
		from Stat_slop_comb
		where Effect = "Intercept";


	/* Stata for each group */
	select n into : n_act
		from Stat_slop_comb2
		where Study_Arm = 'Active' and Effect = "&time_var";
	select n into : n_pla
		from Stat_slop_comb2
		where Study_Arm = 'Placebo' and Effect = "&time_var";

	select mean into : slope_time_mean_act
		from Stat_slop_comb2
		where Study_Arm = 'Active' and Effect = "&time_var";
	select mean into : intcpt_mean_pla
		from Stat_slop_comb2
		where Study_Arm = 'Placebo' and Effect = "&time_var";

quit;

%put &slope_time_SD_pooled;
%put &intcpt_SD_pooled;

%put &n_act;
%put &n_pla;

%put &slope_time_mean_act;
%put &intcpt_mean_pla;



/* slope_time_SD = 0.591453, intcpt_SD = 5.326629, n_act = 22, n_pla = 149 */

proc power;
	twosamplemeans
		test = diff
		meandiff = &effect_size.
		STDDEV = &slope_time_SD_pooled   /* Only pooled SD is allowed */
		alpha = 0.05
		SIDES = U
		power =.
		npergroup = &n_act &n_pla;
quit;


/* Calculate power (sub-method 2); and obtain CI for difference in means  */
data solu_comb2;
	set solu_comb;
	if effect="&time_var";
run;

ods output  Statistics = t_stat
			ConfLimits = t_CI
			TTests = t_test
			Equality = t_equal;
proc ttest data=solu_comb2 alpha=0.05 sides=2;
	class Study_Arm;
	var Estim_comb;
quit;

proc sql noprint;
	select Mean into : diff_mean
		from T_ci
		where Class = 'Diff (1-2)';

	select StdDev into : diff_SD
		from T_ci
		where Class = 'Diff (1-2)';

	select DF into : ttest_DF
		from T_test
		where Method = 'Pooled';
quit;

%put &diff_mean;
%put &diff_SD;
%put &ttest_DF;

/* diff_mean = 0.2397 */
/* diff_SD = 0.5877 */
/* ttest_DF = 169 */



/*** Prepare the dataset for Bayesian power analysis on overall survival ***/
data Eff_die3;


















































































































































































































































































































































































