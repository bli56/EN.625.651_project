DM LOG 'CLEAR';
DM LST 'CLEAR';


proc printto log="D:\Baini_files\ZM\Studies\CB03-154-201\Sample size\Log\SASlog.txt";
quit;

/*** Read-in the simulated data ***/
proc import out=simu_data
	datafile= "D:\Baini_files\ZM\Studies\CB03-154-201\Sample size\Data_simulated\data_simu.csv"
	dbms=csv
	replace;
	getnames=YES;
quit;

data simu_data;
	drop VAR1;
	set simu_data;
run;

proc sort data=simu_data nodupkey out=datalength (keep=data_id);
	by data_id;
quit;

/* Count the total number of datasets */
%let dsid=%sysfunc(open(datalength));
%let N_data=%sysfunc(attrn(&dsid, nlobs));
%let dsid=%sysfunc(close(&dsid));
%put The total number of datasets is &N_data;

/*** Set parameters ***/
%let visit = visit_month;
%let N_act = 98;
%let N_pla = 49;
%let time_var = month;
%let effectsize = 0.3;


%macro batch_run(max_data_id=, N_act=, N_pla=, time_var=);
	%do i = 1 %to &max_data_id.;

		/*** Data preparation ***/
		data subset;
			set simu_data;
			if data_id = &i.;
		proc sort;
			by data_id trt subject_id visit_day;
		quit;


		proc sort data=subset out=lst (keep=subject_id trt) nodupkey;
			by subject_id trt;
		quit;

		data lst_act lst_pla;
			set lst;
			if trt = 0 then output lst_pla;
			else if trt = 1 then output lst_act;
		run;

		/* Select the specified number of subjects for each group */
		data lst_act_pick;
			set lst_act (obs=&N_act.);
		run;

		data lst_pla_pick;
			set lst_pla (obs=&N_pla.);
		run;

		data lst_pick;
			set lst_act_pick lst_pla_pick;
		proc sort;
			by subject_id;
		quit;

		proc delete data = work.lst;
		quit;

		proc delete data = work.lst_act;
		quit;

		proc delete data = work.lst_pla;
		quit;

		proc delete data = work.lst_act_pick;
		quit;

		proc delete data = work.lst_pla_pick;
		quit;

		proc sort data=subset;
			by data_id subject_id trt visit_day;
		quit;

		proc sort data=lst_pick;
			by subject_id;
		quit;

		data data_pick_&i.;
			merge subset (in=a) lst_pick(in=b);
			by subject_id;
			if b;
		run;

		data data_pick_&i.;
			length Study_Arm $12.;
			set data_pick_&i.;
			if trt = 0 then Study_Arm = 'Placebo';
				else if trt = 1 then Study_Arm = 'Active';
			rename  visit_month = month
					visit_week = week
					visit_day = day;
		run;
		
		/*** Shared-baseline mixed model #3 (Model used by Paganoni et al. 2020 (also about AMX-0035): "A shared-baseline, mixed-effects regression model was used, with no added covariates.") ***/
		/*       t --> week or month  */
		/*     Z_i --> Study_Arm      */
		/*   Age_i --> age            */
		/* DelFS_i --> del_FS         */

		title1;
		proc mixed data=data_pick_&i. method=reml alpha=0.05 maxiter=1000 maxfunc=5000;
			class Study_Arm(ref='Placebo') subject_id;
			model ALSFRS_R_Total = &time_var study_Arm*&time_var / cl solution residual outp=outp;
			random intercept &time_var  / subject=subject_id type=un solution;
			contrast 'Contrast: Arm(1) - Arm(0)' study_Arm*&time_var 1 -1;
			estimate 'Estimate: Arm(1) - Arm(0)' study_Arm*&time_var 1 -1;
			ods output  tests3=t3
						SolutionF=Solu_Fix
						SolutionR=Solu_Rand
						contrasts=contrasts
						estimates=estimates;
		quit;

		/* Calculate APPROXIMATE SD for effect coefficient */
		data Solu_Fix;
			set Solu_Fix;
			if Study_Arm ^= 'Placebo' then do;
				SD_appr = StdErr * sqrt(DF+1);
			end;
		run;

		/*** Power analysis (Method 1): Not obtain subject-specific parameter estimates ***/
		data pw1f_&i.;
			data_id = &i.;
			set t3 contrasts;
			alpha = 0.05; /* F-test is one-sided */
			Noncen = NumDF*FValue;
			FCrit = finv(1-Alpha, NumDF, DenDF, 0);
			Power = 1 - probf(FCrit, NumDF, DenDF, Noncen);  /*** ERROR message when estimating power for Effect = "Month" ***/
		run;

		data pw1t_&i.;
			data_id = &i.;
			set estimates;
			do alpha = 0.1;
				do effectsize = &effectsize.;
					tcrit = tinv(1-alpha/2, df); /* Returns the two-sided quantile from t-distribution */
					noncen = effectsize/stderr;
					power = sum(1, -cdf("t", tcrit,df,noncen),cdf("t",-tcrit,df,noncen));
					output;
				end;
			end;
		run;

		/*** Power analysis (Method 2): Obtain subject-specific parameter estimates ***/
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

		proc sort data=data_pick_&i. (keep=subject_id Study_Arm) out=ls_subj nodupkey;
			by subject_id Study_Arm;
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

		/* Obtain pooled CI for the combined slope */
		proc sort data=solu_comb;
			by effect subject_id;
		quit;
		proc means data=solu_comb noprint;
			var estim_comb;
			by effect;
			output out=stat_slop_comb_&i. n=n mean=mean std=std stderr=stderr lclm=lclm uclm=uclm;
		quit;

		/* Obtain CI for the combined slope by arm */
		proc sort data=solu_comb;
			by Study_Arm effect subject_id;
		quit;
		proc means data=solu_comb noprint;
			var Estim_comb;
			by study_arm effect;
			output out=stat_slop_comb2_&i. n=n mean=mean std=std stderr=stderr lclm=lclm uclm=uclm;
		quit;

		/* Calculate power (sub-method 2.1) */
	
		proc sql noprint;
			
			/* Pooled SD */
			select std into : slope_time_SD_pooled
				from stat_slop_comb_&i.
				where Effect = "&time_var";
			select std into : intcpt_SD_pooled
				from stat_slop_comb_&i.
				where Effect = "Intercept";

			/* Stata for each group */
			select n into : n_act
				from stat_slop_comb2_&i.
				where Study_Arm = 'Active' and Effect = "&time_var";
			select n into : n_pla
				from stat_slop_comb2_&i.
				where Study_Arm = 'Placebo' and Effect = "&time_var";

			select mean into : slope_time_mean_act
				from stat_slop_comb2_&i.
				where Study_Arm = 'Active' and Effect = "&time_var";
			select mean into : intcpt_mean_pla
				from stat_slop_comb2_&i.
				where Study_Arm = 'Placebo' and Effect = "&time_var";

		quit;

		%put &slope_time_SD_pooled;
		%put &intcpt_SD_pooled;

		%put &n_act;
		%put &n_pla;

		%put &slope_time_mean_act;
		%put &intcpt_mean_pla;

		ods output  FixedElements=fe21_&i. output=pw21_&i.;
		proc power;
			twosamplemeans
				test = diff
				meandiff = &effectsize.
				STDDEV = &slope_time_SD_pooled.   /* Only pooled SD is allowed */
				alpha = 0.05
				SIDES = U
				power =.
				npergroup = &n_act. &n_pla.;
		quit;

		data fe21_&i.;
			data_id = &i.;
			set fe21_&i.;
		run;

		data pw21_&i.;
			data_id = &i.;
			set pw21_&i.;
		run;

		/* Calculate power (sub-method 2.2); and obtain CI for difference in means  */
		data solu_comb2;
			set solu_comb;
			if effect="&time_var";
		run;

		ods output  Statistics = tsta22
					ConfLimits = tCI22
					TTests = tt22
					Equality = teq22;
		proc ttest data=solu_comb2 alpha=0.05 sides=2;
			class Study_Arm;
			var Estim_comb;
		quit;

		proc sql noprint;
			select Mean into : diff_mean
				from tCI22
				where Class = 'Diff (1-2)';

			select StdDev into : diff_SD
				from tCI22
				where Class = 'Diff (1-2)';

			select DF into : ttest_DF
				from tt22
				where Method = 'Pooled';
		quit;

		%put &diff_mean;
		%put &diff_SD;
		%put &ttest_DF;

		data pw22t_&i.;
			data_id = &i.;
			do alpha = 0.1;  /*NOTE: this alpha is different from the alpha used in PROC TTEST for CI calculation */
				do effectsize = &effectsize.;
					tcrit = tinv(1-alpha/2, &ttest_DF.); /* Returns the two-sided quantile from t-distribution */
					noncen = effectsize/&diff_SD.;
					power = sum(1, -cdf("t", tcrit, &ttest_DF., noncen),cdf("t", -tcrit, &ttest_DF., noncen));
					output;
				end;
			end;
		run;

	%end;
%mend;

%batch_run(max_data_id=&N_data., N_act=&N_act., N_pla=&N_pla., time_var=&time_var.); 


/**************************************************/
/*** Stat summary of picked datasets: mortality ***/
/**************************************************/

/*** Collect picked simu datasets ***/
data all_picked;
	set data_pick_:;
run;

proc datasets lib=work nolist nowarn;
	delete data_pick_:;
quit;

proc sort data=all_picked;
	by data_id trt subject_id day;
quit;


proc sort data=all_picked nodupkey out=mort1 (keep=data_id subject_id trt event_status);
	by data_id subject_id trt;
quit;

data mort1;
	usubjid+1;
	set mort1;
run;

/*** 1 Pooled mortality ***/
/* Overal mortality */
proc freq data=mort1 noprint;
	table event_status/out=freq_mort1;
quit;

/* Mortality by trt */
proc sort data=mort1;
	by trt data_id subject_id;
quit;
proc freq data=mort1 noprint;
	table event_status/out=freq_mort1_1;
	by trt;
quit;

/*** 2 Mortality by data_id ***/
/* Overal mortality by data_id */
proc sort data=mort1;
	by data_id trt subject_id;
quit;

proc freq data=mort1 noprint;
	table event_status/out=freq_mort2;
	by data_id;
quit;

data freq_mort2b;
	set freq_mort2;
	if event_status = 1;
run;

proc means data=freq_mort2b noprint;
	var percent;
	output out=stat_mort2 n= median= mean= std= stderr= min= max= /autoname;
quit;

proc datasets lib=work memtype=data noprint;
	modify stat_mort2; 
	attrib _all_ label='';
quit;


/* Mortality by data_id trt */
proc sort data=mort1;
	by data_id trt subject_id;
quit;

proc freq data=mort1 noprint;
	table event_status/out=freq_mort2_1;  
	by data_id trt;
quit;

data freq_mort2_1b;
	set freq_mort2_1;
	if event_status = 1;
proc sort;
	by trt;
quit;

proc means data=freq_mort2_1b noprint;
	var percent;
	by trt;
	output out=stat_mort2_1 n= median= mean= std= stderr= min= max= /autoname;
quit;

proc datasets lib=work memtype=data noprint;  /* Some plaebo groups have no death, so they are not shown in the dataset */
	modify stat_mort2_1; 
	attrib _all_ label='';
quit;


/*******************************************************/
/*** Stat summary of picked datasets: ALSFRS_R_Total ***/
/*******************************************************/

%let visit = month;

/*** 1 Pooled ALSFRS_R_Total by trt visit ***/
proc sort data=all_picked;
	by trt  &visit.;
quit;

proc means data=all_picked noprint;
	var ALSFRS_R_Total;
	by trt &visit.;
	output out=stat_als1 n= median= mean= std= stderr= min= max= /autoname;
quit;

proc datasets lib=work memtype=data noprint;
	modify stat_als1; 
	attrib _all_ label='';
quit;

/*** 2 ALSFRS_R_Total by data_id trt visit ***/
proc sort data=all_picked;
	by data_id trt  &visit.;
quit;

proc means data=all_picked noprint;
	var ALSFRS_R_Total;
	by data_id trt &visit.;
	output out=stat_als1_2 n= median= mean= std= stderr= min= max= /autoname;
quit;


/*********************************/
/*** Stat summary of modelling ***/
/*********************************/

/*** Collect power analysis results ***/
data all_Pw1f;
	set Pw1f_:;
run;
data all_Pw1t;
	set Pw1t_:;
run;
proc datasets lib=work nolist nowarn;
	delete Pw1f_:;
quit;
proc datasets lib=work nolist nowarn;
	delete Pw1t_:;
quit;


data all_fe21;
	set fe21_:;
run;
data all_pw21;
	set pw21_:;
run;
proc datasets lib=work nolist nowarn;
	delete fe21_:;
quit;
proc datasets lib=work nolist nowarn;
	delete pw21_:;
quit;


data all_pw22t;
	set pw22t_:;
run;
proc datasets lib=work nolist nowarn;
	delete pw22t_:;
quit;


/*** Summarize power analysis results ***/
/* Method 1 */
data All_pw1f;
	set All_pw1f;
	if Effect = '' then Effect = 'Contrast';
proc sort;
	by Effect;
quit;

proc means data=All_pw1f noprint;
	var Power;
	by Effect;
	output out= All_pw1f_stat n=n mean=mean median=median std=std stderr=stderr min=min max=max;
quit;

proc export data=All_pw1f_stat
    outfile="D:\Baini_files\ZM\Studies\CB03-154-201\Sample size\Output\All_pw1f_stat.xlsx"
    dbms=xlsx
    replace;
    sheet="1";
quit;


proc means data=All_pw1t noprint; /* Good */
	var Power;
	output out= All_pw1t_stat n=n mean=mean median=median std=std stderr=stderr min=min max=max;
quit;

proc export data=All_pw1t_stat
    outfile="D:\Baini_files\ZM\Studies\CB03-154-201\Sample size\Output\All_pw1t_stat.xlsx"
    dbms=xlsx
    replace;
    sheet="1";
quit;



/* Method 2.1 */
proc sort data=All_pw21;
	by Index;
quit;
proc means data=All_pw21 noprint;
	by Index;
	var Power;
	output out= All_pw21_stat n=n mean=mean median=median std=std stderr=stderr min=min max=max;
quit;


proc export data=All_pw21_stat
    outfile="D:\Baini_files\ZM\Studies\CB03-154-201\Sample size\Output\All_pw21_stat.xlsx"
    dbms=xlsx
    replace;
    sheet="1";
quit;



/* Method 2.2 */
proc means data=All_pw22t noprint;
	var power;
	output out= ALL_pw22t_stat n=n mean=mean median=median std=std stderr=stderr min=min max=max;
quit;


proc export data=ALL_pw22t_stat
    outfile="D:\Baini_files\ZM\Studies\CB03-154-201\Sample size\Output\ALL_pw22t_stat.xlsx"
    dbms=xlsx
    replace;
    sheet="1";
quit;





