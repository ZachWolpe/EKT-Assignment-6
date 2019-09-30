
* __________________________ Assignment 6: Question 2 __________________________;
data cdata;
set '/folders/myfolders/sasuser.v94/EKT 720/Assignment 6/cdata.sas7bdat';
run;

proc print data=cdata (obs=10);

proc sgplot data=cdata;
	scatter y=y x=x;
run;



proc iml;
use cdata;
read all into xy;
n = nrow(xy);
y_orig = xy[,2];
x_orig = J(n,1,1)||xy[,1];
y = y_orig;
x = x_orig;
iters = 1000;

start reg;
	k=ncol(x);
	bh=inv(x`*x)*x`*y;
	sse=ssq(y-x*bh);
	cssy=ssq(y-(sum(y)/n));
	r2=(cssy-sse)/cssy;
	F = (r2/(k-1))/((1-r2)/(n-k));
	if pp=1 then do;
		print bh r2 f;
	end;
finish reg;




* _______ Bootstrap Pairs _______;

start bootstrap_pairs;
	in = sample(1:n, n, 'replace')`;
	x = x_orig[in,];
	y = y_orig[in,];
	call reg;
finish bootstrap_pairs;



do i=1 to iters;
	call bootstrap_pairs;
	res = res // (bh` || r2);
end;

create pair_betas from res[colname={'b0 b1 r2'}];
append from res;


* _______ Bootstrap Errors _______;

start bootstrap_errors;	
	in = sample(1:n, n, 'replace')`;
	y = y_orig + e[in,];
	call reg;
finish bootstrap_errors;

x = x_orig;
y = y_orig;
call reg;
e = y - x*bh;

do i=1 to iters;
	call bootstrap_errors;
	res_e = res_e // (bh` || r2);
end;

create error_betas from res_e[colname={'b0 b1 r2'}];
append from res_e;



* _______ Results _______;

proc sgplot data=error_betas;
	histogram b0 / binwidth=10;
	title 'Error Resampling';
	title2 'Beta 0';
run;

proc sgplot data=error_betas;
	histogram b1 / binwidth=.1;
	title 'Error Resampling';
	title2 'Beta 1';
run;
	
proc sgplot data=error_betas;
	histogram r2 / binwidth=0.001;
	title 'Error Resampling';
	title2 'R2: Coefficient of Determination';
run;


proc sgplot data=pair_betas;
	histogram b0 / binwidth=10;
	title 'Error Resampling';
	title2 'Beta 0';
run;

proc sgplot data=pair_betas;
	histogram b1 / binwidth=.1;
	title 'Error Resampling';
	title2 'Beta 1';
run;
	
proc sgplot data=pair_betas;
	histogram r2 / binwidth=.001;
	title 'Error Resampling';
	title2 'R2: Coefficient of Determination';
run;
	
	
	



















