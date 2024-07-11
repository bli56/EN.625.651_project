DM LOG 'CLEAR';
DM LST 'CLEAR';

/**************************************/
/***** Read-in the simulated data *****/
/**************************************/
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


/*******************************************/
/*** Stat summary of datasets: mortality ***/
/*******************************************/

proc sort data=simu_data nodupkey out=mort1 (keep=data_id subject_id trt event_status);
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


/************************************************/
/*** Stat summary of datasets: ALSFRS_R_Total ***/
/************************************************/

%let visit = visit_month;

/*** 1 Pooled ALSFRS_R_Total by trt visit ***/
proc sort data=simu_data;
	by trt  &visit.;
quit;

proc means data=simu_data noprint;
	var ALSFRS_R_Total;
	by trt &visit.;
	output out=stat_als1 n= median= mean= std= stderr= min= max= /autoname;
quit;

proc datasets lib=work memtype=data noprint;
	modify stat_als1; 
	attrib _all_ label='';
quit;

/*** 2 ALSFRS_R_Total by data_id trt visit ***/
proc sort data=simu_data;
	by data_id trt  &visit.;
quit;

proc means data=simu_data noprint;
	var ALSFRS_R_Total;
	by data_id trt &visit.;
	output out=stat_als1_2 n= median= mean= std= stderr= min= max= /autoname;
quit;

































