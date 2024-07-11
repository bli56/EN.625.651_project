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

data surv1;
	merge Eff_data2(in=a) Proact_die(in=b);
	by subject_id;
	if a;
run;

/*** Apply the censoring time window: 0~(335-1) ***/
data Surv1;
	set Surv1;
	if .<death_days and death_days =< 334 then do;
		deathD = death_days;
		deathFL = 1;
	end;
	else do;
		deathD = .;
		deathFL = 0;
	end;
proc sort;
	by subject_id week;
quit;


data Surv2;
	set Surv1;
	by subject_id week;
	if last.subject_id;
proc sort;
	by subject_id week;
quit;


/*** Find age at baseline  ***/
data BL_age;
	keep subject_id age;
	set Surv1;
	by subject_id week;
	rename age = BL_age;
	if first.subject_id;
proc sort;
	by subject_id;
quit;


data Surv3;
	merge Surv2(in=a) BL_age(in=b);
	by subject_id;
	if a and b;
proc sort;
	by subject_id week;
quit;


/**************************/
/*** Classical Analysis ***/
/**************************/

proc phreg data = Surv3
			plots=survival;
	class Study_Arm ;
	model week*deathFL(0) = Study_Arm;
	hazardratio Study_Arm;
quit;


proc phreg data = Surv3
			plots=survival;
	class Study_Arm ;
	model week*deathFL(0) = Study_Arm BL_age del_FS;
	hazardratio Study_Arm;
quit;


/*************************/
/*** Bayesian Analysis ***/
/*************************/

ods trace on;
ods output 
	CensoredSummary = CensoredSummary    /* Summary of the Number of Event and Censored Value */
	ParameterEstimates = ParameterEstimates    /* Maximum Likelihood Estimates */
	FitStatistics = FitStatistics   /* Fit Statistics */
	PostSumInt = PostSumInt   /* Posterior Summaries and Intervals */
	Gelman = Gelman   /* Gelman-Rubin Diagnostics */
	Geweke = Geweke   /* Geweke Diagnostics */ 
	ESS = ESS    /* Effective Sample Sizes */
	HazardRatios = HazardRatios;   /* Hazard Ratios */
proc phreg data = Surv3;
	class Study_Arm ;
	model week*deathFL(0) = Study_Arm;
	bayes seed=1 diagnostic=all plots=all;
	hazardratio Study_Arm;
quit;
ods trace off;






proc phreg data = Surv3;
	class Study_Arm ;
	model week*deathFL(0) = Study_Arm BL_age del_FS;
	bayes seed=1 diagnostic=all plots=all;
	hazardratio Study_Arm;
quit;







































































































































































































































































































































































































































































































































































































































































