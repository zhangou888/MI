/*****************************************************************************/
/***             Topic: Measurement Invariance Tests in SAS                ***/
/***                     Psychometric Conference 2013                      ***/
/***                     Pearson Clinical Assessment                       ***/
/***                     Presenter: Ou Zhang                               ***/
/***                     Date: 2013-11-02                                  ***/
/*****************************************************************************/

/*This example data from the book chapter by Thompson and Green (2006, p.139,
Table 5.2, Dataset 2) were borrowed, which contain six measured variables aiming to
assess preschool children academic (V1-V3) and social school readiness (V4-V6).
Preschool children were divided into two groups: Group 1—day-care and Group 2—
home-care. We use the two matrices of means and covariances for Group 1 and Group 2,
and the sample sizes are 250 and 150 for Group 1 and Group 2, respectively.*/

/*Creating two separate datasets in the type of covariance matrix in the SAS system*/
data group1(type=cov); 
/*** (type=cov) is necessary to tell SAS this is a covariance matrix 
     instead of raw data collected in rows. ***/
infile datalines missover;
input _NAME_ $ _TYPE_ $ V1-V6;
datalines;
. MEAN 49.14 82.60 104.95 78.58 54.95 119.91
V1 COV 154.54
V2 COV 44.75 90.23
V3 COV 40.98 22.77 78.76
V4 COV 41.35 2.83 7.92 220.12
V5 COV 23.28 9.12 0.75 61.46 159.88
V6 COV 47.08 22.88 16.05 125.08 84.31 332.26
;
data group2(type=cov);
infile datalines missover;
input _NAME_ $ _TYPE_ $ V1-V6;
datalines;
. MEAN 55.01 81.22 97.01 72.67 46.57 124.28
V1 COV 124.93
V2 COV 52.19 80.67
V3 COV 64.45 42.85 83.56
V4 COV 59.95 33.10 38.34 290.65
V5 COV 32.32 16.09 18.29 124.71 169.17
V6 COV 87.15 39.70 51.82 174.20 108.39 355.22
;
run;

/*******************************************************************************/
/*Step 0—Evaluainge same factor model in each of the groups and combined groups*/
/*******************************************************************************/
*** Group 1 ***;
proc calis data=group1 cov method=ml outstat=stat1 outfit=fit1 nobs=250 maxiter=1000 maxfunc=1000;
/* List all the goodness-of-fit indices you want to report */
fitindex NOINDEXTYPE on(only)=[CHISQ DF PROBCHI SRMSR CFI RMSEA];
 var V1-V6;
 lineqs
	v1 = b1 * f1 + e1,
	v2 = b2 * f1 + e2,
	v3 = b3 * f1 + e3,
	v4 = b4 * f2 + e4,
	v5 = b5 * f2 + e5,
	v6 = b6 * f2 + e6;	
 variance
	f1 = 1.0,         /* factor variance is set as 1.0 */
	f2 = 1.0,
	e1-e6 = ve1-ve6;  /* error variances are saved to ve1-ve6 */
 cov
	F1 F2 = d12;      /* factor covariance */ 
run;

/* Group 2*/
proc calis data=group2 cov method=ml outstat=stat2 outfit=fit2 nobs=150 maxiter=1000 maxfunc=1000;
/* List all the goodness-of-fit indices you want to report */
fitindex NOINDEXTYPE on(only)=[CHISQ DF PROBCHI SRMSR CFI RMSEA];
 var V1-V6;
 lineqs
	v1 = b1 * f1 + e1,
	v2 = b2 * f1 + e2,
	v3 = b3 * f1 + e3,
	v4 = b4 * f2 + e4,
	v5 = b5 * f2 + e5,
	v6 = b6 * f2 + e6;	
 variance
	f1 = 1.0,        /* factor variance is set as 1.0 */
	f2 = 1.0,
	e1-e6 = ve1-ve6;  /* error variances are saved to ve1-ve6 */
 cov
	f1 f2 = d12;   /* factor covariance */ 
run;

*** Combined groups for co/var ***;
proc calis cov method=ml outstat=stat3 outfit=fit3 maxiter=1000 maxfunc=1000 pall;
 fitindex NOINDEXTYPE on(only)=[CHISQ DF PROBCHI SRMSR CFI RMSEA];
 group 1 / data=group1 nobs=250;
 group 2 / data=group2 nobs=150;
 model 1 / group=1;
 var V1-V6;
 lineqs
	v1 = b1 * f1 + e1,
	v2 = b2 * f1 + e2,
	v3 = b3 * f1 + e3,
	v4 = b4 * f2 + e4,
	v5 = b5 * f2 + e5,
	v6 = b6 * f2 + e6;	
 variance
	f1 = 1.0,         /* factor variance is set as 1.0 */
	f2 = 1.0,
	e1-e6 = ve1-ve6;  /* error variances are saved to ve1-ve6 */
 cov
	f1 f2 = d12;      /* factor covariance */ 
 model 2 / group=2;
 refmodel 1 / AllNewParms; /* reference group is group 1 */ 
run;

/*******************************************************************************/
/*      Step 1—Evaluating Configural Invariance Model                          */
/*******************************************************************************/
proc calis cov method=ml outstat=stat4 outfit=fit4 pall maxiter=1000 maxfunc=1000;
fitindex NOINDEXTYPE on(only)=[CHISQ DF PROBCHI SRMSR CFI RMSEA];
 group 1 / data=group1 nobs=250;
 group 2 / data=group2 nobs=150;
 model 1 / group=1;
 lineqs
	v1 = a1 * intercept + 1  * f1 + e1,   /* 1 of the factor loading for each factor is set as 1.0 */
	v2 = a2 * intercept + b2 * f1 + e2,
	v3 = a3 * intercept + b3 * f1 + e3,
	v4 = a4 * intercept + 1  * f2 + e4,   /* 1 of the factor loading for each factor is set as 1.0 */
	v5 = a5 * intercept + b5 * f2 + e5,
	v6 = a6 * intercept + b6 * f2 + e6;	
 variance
	f1  = vf1,        /* factor variance freely estimated */
	f2  = vf2,
	e1-e6 = ve1-ve6;
 cov
	f1 f2 = d12;      /* factor covariance */ 
 mean
 	f1 = 0,          /* factor mean fix to 0 */  
    f2 = 0;  
 model 2 / group=2;
 refmodel 1;         /* reference group is group 1 */
 lineqs
	v1 = a7  * intercept + 1   * f1 + e1,   /* 1 of the factor loading is set as 1.0 */
	v2 = a8  * intercept + b8  * f1 + e2,
	v3 = a9  * intercept + b9  * f1 + e3,
	v4 = a10 * intercept + 1   * f2 + e4,    /* 1 of the factor loading is set as 1.0 */
	v5 = a11 * intercept + b11 * f2 + e5,
	v6 = a12 * intercept + b12 * f2 + e6;	
 variance
	f1  = nvf1,         
	f2  = nvf2,
	e1-e6  = nve1-nve6;  
 cov
	f1 f2 = nd12;   /* factor covariance */ 
 mean
 	f1 = 0,         /* F1 and F2 have own variance but mean fixed to 0 */ 
    f2 = 0; 
run;

/*******************************************************************************/
/*      Step 2—Evaluating between-group equivalence of factor loadings         */
/*                          Metric Invariance Model                            */
/*******************************************************************************/
proc calis cov method=ml outstat=stat5 outfit=fit5 pall maxiter=1000 maxfunc=1000;
fitindex NOINDEXTYPE on(only)=[CHISQ DF PROBCHI SRMSR CFI RMSEA];
 group 1 / data=group1 nobs=250;
 group 2 / data=group2 nobs=150;
 model 1 / group=1;
 lineqs
	v1 = a1 * intercept + b1 * f1 + e1,
	v2 = a2 * intercept + b2 * f1 + e2,
	v3 = a3 * intercept + b3 * f1 + e3,
	v4 = a4 * intercept + b4 * f2 + e4,
	v5 = a5 * intercept + b5 * f2 + e5,
	v6 = a6 * intercept + b6 * f2 + e6;	
 variance
	f1 = 1.0,         /* factor variance is set as 1.0 */
	f2 = 1.0,
	e1-e6 = ve1-ve6;  /* error variances are saved to ve1-ve6 */
 cov
	f1 f2 = d12;      /* factor covariance */
 mean
 	f1 = 0,           /* factor mean fix to 0 */  
    f2 = 0; 
 model 2 / group=2;
  refmodel 1;
  lineqs
	v1 = a21 * intercept + b1 * f1 + e1,
	v2 = a22 * intercept + b2 * f1 + e2,
	v3 = a23 * intercept + b3 * f1 + e3,
	v4 = a24 * intercept + b4 * f2 + e4,
	v5 = a25 * intercept + b5 * f2 + e5,
	v6 = a26 * intercept + b6 * f2 + e6;	
 variance
	f1 = vf1,           /* factor variance is freely estimated */
	f2 = vf2,
	e1-e6 = nve1-nve6;  /* error variances are saved to nve1-nve6 */
 cov
	f1 f2 = nd12;       /* factor covariance */ 
 mean
 	f1 = 0,             /* factor mean fix to 0 */  
    f2 = 0;
run;

/*******************************************************************************/
/*      Step 3—Evaluating between-group equivalence of factor Intercepts       */
/*                       Scalar Invariance Model                               */
/*******************************************************************************/
proc calis cov method=ml outstat=stat6 outfit=fit6 maxiter=1000 maxfunc=1000 pall;
fitindex NOINDEXTYPE on(only)=[CHISQ DF PROBCHI SRMSR CFI RMSEA];
 group 1 / data=group1 nobs=250;
 group 2 / data=group2 nobs=150;
 model 1 / group=1;
 lineqs
	v1 = a1 * intercept + b1 * f1 + e1,
	v2 = a2 * intercept + b2 * f1 + e2,
	v3 = a3 * intercept + b3 * f1 + e3,
	v4 = a4 * intercept + b4 * f2 + e4,
	v5 = a5 * intercept + b5 * f2 + e5,
	v6 = a6 * intercept + b6 * f2 + e6;	
 variance
	f1 = 1.0,          /* factor variance is set as 1.0 */
	f2 = 1.0,
	e1-e6 = ve1-ve6;  /* error variances are saved to ve1-ve6 */
 cov
	f1 f2 = d12;      /* factor covariance */ 
 mean
 	f1 = 0,
    f2 = 0;
 model 2 / group=2;
  refmodel 1;
  lineqs
	v1 = a1 * intercept + b1 * f1 + e1,
	v2 = a2 * intercept + b2 * f1 + e2,
	v3 = a3 * intercept + b3 * f1 + e3,
	v4 = a4 * intercept + b4 * f2 + e4,
	v5 = a5 * intercept + b5 * f2 + e5,
	v6 = a6 * intercept + b6 * f2 + e6;	
 variance
	f1 = vf1,         
	f2 = vf2,
	e1-e6 = nve1-nve6;  /* error variances are saved to nve1-nve6 */
 cov
	f1 f2 = nd12;      /* factor covariance */ 
  mean
 	f1 = m_f1,
    f2 = m_f2;
run;


/*******************************************************************************/
/*      Step 4—Evaluating between-group equivalence of residual variance       */
/*                       Residual Invariance Model                               */
/*******************************************************************************/
proc calis cov method=ml outstat=stat7 outfit=fit7 maxiter=1000 maxfunc=1000 pall;
fitindex NOINDEXTYPE on(only)=[CHISQ DF PROBCHI SRMSR CFI RMSEA];
 group 1 / data=group1 nobs=250;
 group 2 / data=group2 nobs=150;
 model 1 / group=1;
 lineqs
	v1 = a1 * intercept + b1 * f1 + e1,
	v2 = a2 * intercept + b2 * f1 + e2,
	v3 = a3 * intercept + b3 * f1 + e3,
	v4 = a4 * intercept + b4 * f2 + e4,
	v5 = a5 * intercept + b5 * f2 + e5,
	v6 = a6 * intercept + b6 * f2 + e6;	
 variance
	f1 = 1.0,          /* factor variance is set as 1.0 */
	f2 = 1.0,
	e1-e6 = ve1-ve6;  /* error variances are saved to ve1-ve6 */
 cov
	f1 f2 = d12;      /* factor covariance */ 
 mean
 	f1 = 0,
    f2 = 0;
 model 2 / group=2;
  refmodel 1;
  lineqs
	v1 = a1 * intercept + b1 * f1 + e1,
	v2 = a2 * intercept + b2 * f1 + e2,
	v3 = a3 * intercept + b3 * f1 + e3,
	v4 = a4 * intercept + b4 * f2 + e4,
	v5 = a5 * intercept + b5 * f2 + e5,
	v6 = a6 * intercept + b6 * f2 + e6;	
 variance
	f1 = vf1,         
	f2 = vf2,
	e1-e6 = ve1-ve6;  /* error variances are saved to ve1-ve6 */
 cov
	f1 f2 = nd12;      /* factor covariance */ 
  mean
 	f1 = m_f1,
    f2 = m_f2;
run;

/* Model fit Indices */
%macro fit(dat,var); 
data &dat&var;set &dat&var;
keep FitIndex FitValue;
run;
proc transpose data=&dat&var
			   out=&dat&var.1;
run;

data &dat&var.1;set &dat&var.1;
	keep COL11 COL12 COL13 COL17 COL21 COL22 COL23 COL24 COL25 COL29;
	rename COL11 = chi;
    rename COL12 = df;
	rename COL13 = p;
    rename COL17 = SRMSR;
    rename COL21 = RMSEA;
	rename COL22 = RMSEA_lo;
	rename COL23 = RMSEA_hi;
	rename COL24 = Pclose;
	rename COL25 = AIC;
    rename COL29 = CFI;
run;

data &dat&var.1;set &dat&var.1;
	RMSEA90 = compress('['||round(RMSEA_lo,.001)||','||round(RMSEA_hi,.001)||']');
	drop RMSEA_lo RMSEA_hi;
run;

data &dat&var.1;
	retain chi df p CFI RMSEA RMSEA90 SRMSR AIC Pclose;
    set &dat&var.1;
	label RMSEA90  ="RMSEA 90%CI"; 
run;
%mend;

%fit(fit,3);
%fit(fit,4);
%fit(fit,5);
%fit(fit,6);
%fit(fit,7);

data fit;set fit31 fit41 fit51 fit61 fit71;run;

proc export data=fit
		    outfile="C:\temp\MI_fit.xls" 
		    dbms=excel2000 replace; 
		    sheet="model fit"; 
run;









