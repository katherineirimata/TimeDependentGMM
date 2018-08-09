/*****************************************************************/
/**	SAS Code: %GMM Macro   									    **/
/** Programmer: Katherine Cai									**/
/** Description: Performs GMM regression, using method by		**/
/**		Lalonde, Wilson, and Yin (2014)  		 				**/
/** Notes: Requires macro %MVINTEGRATION to be run first		**/
/**		See SAS® Macro for Generalized Method of Moments 		**/
/**		Estimation for Longitudinal Data with Time-Dependent    **/
/**		Covariates, Cai and Wilson (2016) for reference         **/
/*****************************************************************/


%macro GMM(ds=, file=, reflib=, timeVar=, outVar=, predVar=, idVar=, alpha=, distr=);
LIBNAME DS &ds.;  

proc sql;
select max(&timeVar.) into :timepts from DS.&file.;
quit;

DATA mydata; SET DS.&file.; 
%let n = &timepts.;
  ARRAY dummys {*} 3.  t1 - t&n.;
 
  DO i=2 TO &n.;			      
    dummys(i) = 0;
  END;
  dummys( &timeVar.  ) = 1;		
  drop i t1;
  %let timeDummyVars = t2 - t&n.;
RUN;

TITLE "Pooled Regression by &timeVar.";
PROC SORT DATA=mydata OUT=mydatasorted; 
BY &timeVar.; RUN;

%if &distr.=normal %then
 	%do;
		proc reg data=mydatasorted NOPRINT;
		BY &timeVar.;
		MODEL &outVar. = &predVar.;
		OUTPUT OUT=outpool3 PREDICTED=mu;
		RUN;

		DATA outpool3; SET outpool3;
		wt = 1/mu; 
		rsdraw = &outVar.-mu; 
%end;
%else %if &distr.=bin %then
      %do;
		PROC logistic DATA=mydatasorted NOPRINT; 
		BY &timeVar.;
		MODEL &outVar. (event='1') = &predVar. / aggregate scale=none;
		OUTPUT OUT=outpool3 P=mu;
		RUN;

		DATA outpool3; SET outpool3;
		wt = mu*(1-mu);
		rsdraw = &outVar.-mu; 
%end;
%else %do;
	%put ERROR must be normal or binomial distributions; 
	%return;
%end;

PROC SORT DATA=outpool3 OUT=outpool3 ;
  BY &idVar. &timeVar.; RUN;

QUIT;
proc iml;
libname reflib &reflib.;
RESET STORAGE=reflib.MVIntegration;
LOAD;

use outpool3;
read all VARIABLES {&predVar.} into Zmat; 
read all var {wt} into wt;
read all var {rsdraw} into rsd;
read all VARIABLES {&idVar. &timeVar.};
close outpool3;

corrind = 1; 

N=ncol(unique(&idVar.));
T = max(&timeVar.); 
Np = ncol(Zmat);
alpha = &alpha.;

print N; *N=number of unique subjects;
print T; *T=number of observations/measurements per subject;

cutoff2 = 0.05; cutoff3 = 0.01;
level1 = 25000;level2 = 250000;level3 = 2500000; 
ABSEPS = 0.000001;
RELEPS = 0;


start rho(a,rsd) global(N,T);
abm = j(N,2*T,.);
abm[,1:T] = shape(rsd,N);
abm[,T+1:2*T] = shape(a,N);
corr = corr(abm);  
rho = corr[1:T,T+1:2*T]; 
return(rho);
finish rho;

start stddev(a,rsd) global(N,T);
bm = shape(rsd,N); 
bdev = bm-j(N,1,1)*bm[:,]; 
bdev2 = bdev#bdev;      
am = shape(a,N);   
adev = am-j(N,1,1)*am[:,];  
adev2 = adev#adev;      
stddev = sqrt( (1/N)*t(bdev2)*adev2 );
return(stddev);
finish stddev;

*Corrected standardization;
start stdzn(x) global(N,T);
xrows = shape(x,N);  
y = xrows - xrows[:,]; 
vcv = (1/(N-1))*t(y)*y; 
v = diag(vcv);
sinv = sqrt(inv(v));
x2 = y*sinv;  
x2  = shape(x2,N*T,1);
return(x2);
finish stdzn;

pvec = j(Np*T*T,1,.);sevec = j(Np*T*T,1,.); 
r4out = j(T,T*Np,.); se4out = j(T,T*Np,.); z4out = j(T,T*Np,.); p4out = j(T,T*Np,.); 

y = rsd;
y_std = stdzn(y);
DO i=1 TO Np;
x = wt#Zmat[,i]; 
x_std = stdzn(x);
r = rho(x_std,y_std);
se = stddev(x_std,y_std);
z = sqrt(N)*(r/se);
p = 2*(1-cdf('normal',abs(z))); 

r4out[,T*(i-1)+1:T*i] = r;
se4out[,T*(i-1)+1:T*i] = se;
z4out[,T*(i-1)+1:T*i] = z;
p4out[,T*(i-1)+1:T*i] = p;

DO j = 1 TO T;
p[j,j] = 1; 
END;
pvec[T*T*(i-1)+1:T*T*i,1] = shape(p,T*T,1); 
sevec[T*T*(i-1)+1:T*T*i,1] = shape(se,T*T,1); 
END;

pse = j(Np*T*T,2,.); pse[,1] = pvec; pse[,2] = sevec; 
call sort(pse,{1}); 

stop = 0;
rnking = 1;

DO WHILE (stop<1);

	pmin = pse[rnking,1];

	se4test = pse[rnking:Np*T*T,2]; 

	L = Np*T*T - rnking +1;


	PRINT rnking;
	PRINT pmin;

	if pmin = 0 then DO; pact = 0 ; print 'Integration Skipped'; print pact; END;
	if pmin >=0.5 then DO; pact = 1; print 'Integration Skipped'; print pact; END;

	if pmin>0 & pmin <0.5 then DO;

		 LOWER = probit(pmin*0.5)*J(1,L,1); 
		 UPPER = probit(1-pmin*0.5)*J(1,L,1);
		 INFIN = J(1,L,-1);
		 Corr1 = I(L); Corr2 = j(L,L,1); if corrind = 1 then Corr = Corr1; else Corr = Coor2;
		 COVAR = diag(se4test)*Corr*diag(se4test); *Variance matrix for test;

				MAXPTS = level3;
				RUN MVN_DIST( L, LOWER, UPPER, INFIN, COVAR, MAXPTS, ABSEPS, RELEPS,  ERROR, VALUE, NEVALS, INFORM );
			pact = 1- VALUE;
			print 'Level 3 Integration';	print pact;

	END; 

	 if pact >= alpha | L=1 then stop = 1;
	 if pact < alpha & L>1 then DO;	rnking = rnking +1; END;

END; 


print 'final pvalue';
print pact, L;


*Correlation Test - by Multiple Test;
TypeVec = (pvec >= pmin*j(Np*T*T,1,1) );
Type = shape(TypeVec, Np, T*T);

TypeMtx = j(Np+T,T*T,.);
Tint = I(T);
TypeMtx[1,]= shape(Tint,1,T*T);
TypeMtx[2:1+Np,] = Type;
DO i = 2 to T;					
	Tt = j(T,T,0);Tt[,i] = Tint[,i];
	TypeMtx[Np+i,] = shape(Tt,1,T*T);
END;

CREATE TypeMtx from TypeMtx;
APPEND from TypeMtx;
CLOSE TypeMtx;

print 'Each row of TypeMtx is the type vector for each of the predictors (by multiple test)';
print 'The intercept admitted T valid equations';
print TypeMtx;

*Correlation Test - by Individual Tests;
TypeVec2 = (pvec >= alpha*j(Np*T*T,1,1) );
Type2 = shape(TypeVec2, Np, T*T);
TypeMtx2 = TypeMtx;
TypeMtx2[2:1+Np,] = Type2;
CREATE TypeMtx2 from TypeMtx2;
APPEND from TypeMtx2;
CLOSE TypeMtx2;

print 'Each row of TypeMtx2 is the type vector for each of the predictors (by individual tests)';
print 'The intercept admitted T valid equations';

print TypeMtx2;


Type4out = j(T,T*(Np+T),.);
DO i=1 to Np+T;
Type4out[,T*(i-1)+1:T*i] = shape(TypeMtx[i,],T,T);
END;


print '========================================================================================';
print 'Note: r4out, se4out, z4out, p4out consider only non-intercept, non-time, dummy variables';
print '========================================================================================';
print r4out;
print se4out;
print z4out;
print p4out;
print 'Type4out is the TypeMtx rearranged (for Excel output); each TxT block is the type vector for a predictor';
print 'Correlation results by multiple test, the intercept admitted only T valid equations';
print Type4out; 

x = wt; 
x_std = stdzn(x);
r = rho(x_std,y_std);
se = stddev(x_std,y_std);
z = sqrt(N)*(r/se);
p = 2*(1-cdf('normal',abs(z)));
print r, se, z, p;
T_wt = p>0.05;
print T_wt;
T_wt = shape(T_wt,1); 
print T_wt;
TypeMtx3 = j(Np+T,T*T,.);
TypeMtx3[1,] = T_wt;
TypeMtx3[2:Np+1,] = Type2;
DO i = 2 to T;
	T_wtT = j(T,T,0);
	T_wtT[,i] = T_wt[,i];
	T_wtT = shape(T_wtT,1); 
	print "T_twT #"; print i;
	print T_wtT;
	TypeMtx3[Np+i,] = T_wtT;
END;

CREATE TypeMtx3 from TypeMtx3;
APPEND from TypeMtx3;
CLOSE TypeMtx3;
print 'Each row of TypeMtx3 is the type vector for each of the predictors (by individual tests)';
print 'The intercept and time dummy variables may admit more than T valid equations';
print TypeMtx3;
Quit;

PROC SORT DATA=mydata;
BY &idVar. &timeVar.; RUN;

data mydata_time;
set mydata(keep = &timeDummyVars.);
run;

PROC genmod DATA=mydata descend;
      CLASS &idVar. &timeVar.;
      MODEL &outVar.= &predVar. &timeDummyVars. / DIST=&distr. ;
      REPEATED SUBJECT=&idVar. /WITHIN=&timeVar. CORR=indep CORRW;
	  OUTPUT OUT=GEEout XBETA=xb RESRAW = rraw;
	  ods output GEEEmpPEst=betaGEE;
RUN;
ods output close;

PROC IML;   
USE mydata_time;
READ all VAR _ALL_ INTO Zmat2[colname = Timenames];
CLOSE mydata_time;
USE mydata;
READ all VAR {&predVar.} INTO Zmat1[colname = Prednames];
READ all VAR {&outVar.} INTO yvec;
read all VARIABLES {&idVar. &timeVar.};
CLOSE mydata;
Zmat = Zmat1 || Zmat2;
PRINT 'Method: 2SGMM YWL';

N=ncol(unique(&idVar.)); *N=number of observations;
Pn=ncol(Zmat)+1; *Pn=number of covariates/parameters to estimate; 
						
Tn=max(&timeVar.); *Tn=number of observations/measurements per subject;

nr = nrow(Zmat);     
nc = ncol(Zmat);   
int = j(nr,1,1);
Xmat =j(nr,nc+1,.); Xmat[,1]=int; Xmat[,2:nc+1]=Zmat;

USE betaGEE;
READ all VAR {Estimate} INTO beta0;
CLOSE betaGEE;
beta0 = t(beta0);

USE TypeMtx;
READ all INTO TypeMtx;
CLOSE TypeMtx;
print 'Each row of TypeMtx is the type vector for each of the predictors (by multiple test)' ;
print TypeMtx;

neq = j(Pn,1,0);
DO p =1 TO Pn;
    neq[p] = ncol(loc(TypeMtx[p,]^=0));
END;

nloc = j(1,Pn+1,0);
DO p =1 TO Pn;
  nloc[p+1] = sum(neq[1:p]);
END;
nreg = sum(neq);

Wn = I(nreg); 
S = j(nreg,nreg,0);

START TSGMM(beta) global(Pn,Tn,N,Xmat,yvec,nreg,TypeMtx,nloc,Wn,S);      
Gn = j(nreg,1,0);                
S = j(nreg,nreg,0);
eq = j(nreg,N,0);              

DO i = 1 TO N;  
  x = Xmat[(i-1)*Tn+1:i*Tn,]; 
  y = yvec[(i-1)*Tn+1:i*Tn];  
  if "&distr." = "bin" then 
  	mu = exp(x*t(beta)) / ( 1+exp(x*t(beta)) );
  else if "&distr." = "normal" then
 	mu = x*t(beta);
  Rsd = y - mu;               
  DO p = 1 TO Pn;
  	if "&distr." = "bin" then
    	D = x[,p]#mu#(1- mu);
	else if "&distr." = "normal" then
		D= x[,p]#(mu##(-1));
    Eqmtx = Rsd*t(D);
    eq[nloc[p]+1:nloc[p+1],i] = Eqmtx[loc(TypeMtx[p,]^=0)];   
  END;
  S = S + eq[,i]*t(eq[,i]); 
END;  
Gn = eq[,:];               
f = t(Gn)*Wn*Gn; *Objective function to be minimized; 
RETURN(f);
FINISH TSGMM;

tc = {2000 2000}; optn = {0 2}; *optn[1]=0 specifies minimization;
  CALL NLPNRA(rc, xres,"TSGMM", beta0,optn, , tc); *Nonlinear optimization using Newton Raphson;
  beta0 = xres; 
  Wn = ginv(S/N);  

  CALL NLPNRA(rc, xres,"TSGMM", beta0,optn, , tc);
  beta = xres;
  Wn = ginv(S/N);

DG = j(nreg,Pn,.); 
DO k = 1 TO Pn;
  DGi = j(nreg,N,0); 
  DO i = 1 TO N; 
    x = Xmat[(i-1)*Tn+1:i*Tn,]; 
    y = yvec[(i-1)*Tn+1:i*Tn];   
	if "&distr." = "bin" then 
    	mu = exp(x*t(beta)) / ( 1+exp(x*t(beta)) );
	else if "&distr." = "normal" then
		mu = x*t(beta);
    Rsd = y - mu;  
	if "&distr." = "bin" then do;
    	Dk =  x[,k]#mu#(1- mu);
		Dkz =  x[,k]#(1- 2*mu);
		end;
	else if "&distr." = "normal" then do;
		Dk =  x[,k]#(mu##(-1));
		Dkz =  x[,k]#(-mu##(-2));
		end;
    DO p = 1 TO Pn;
	  if "&distr." = "bin" then 
      	Dp = x[,p]#mu#(1- mu);
	  else if "&distr." = "normal" then
		Dp = x[,p]#(mu##(-1));
      Dkzp = Dkz#Dp;        
	  DGmtx = Rsd*t(Dkzp)-Dk*t(Dp);
      DGi[nloc[p]+1:nloc[p+1],i] = DGmtx[loc(TypeMtx[p,]^=0)];   
    END;
  END;
  DG[,k]= DGi[,:]; 
END;   

AsymVar = (1/N)*ginv(t(DG)*Wn*DG); 
AVvec = vecdiag(AsymVar);
StdDev = sqrt(AVvec);

zvalue = t(beta)/StdDev;
pvalue = 2*(1-cdf('normal',abs(zvalue)));

Outmtx = j(Pn,4,.);
Outtitle={'Estimate'  'StdDev'  'Zvalue'  'Pvalue'};
Varnames_int = {'Intercept'};
Varnames=t(Varnames_int || Prednames || Timenames);

Outmtx[,1]=t(beta);
Outmtx[,2]=StdDev;
Outmtx[,3]=zvalue;
Outmtx[,4]=pvalue;
PRINT Outmtx[colname=Outtitle rowname=Varnames];

resvec = yvec - Xmat*t(beta);

CREATE resdata FROM resvec[colname={"res"}];
APPEND FROM resvec;
CLOSE resdata;


betavec = shape(beta,1);
print betavec;
CREATE betaGMM var {betavec};;
APPEND;
CLOSE betaGMM;


QUIT;

%mend GMM;
