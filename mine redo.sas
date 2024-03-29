
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
pp=0;

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

upper_r2 = mean(res[,3]) + 1.96*std(res[,3]);
lower_r2 = mean(res[,3]) - 1.96*std(res[,3]);
res = res || J(nrow(res), 1, lower_r2) || J(nrow(res), 1, upper_r2);


create pair_betas from res[colname={'b0' 'b1' 'r2' 'lower_r2' 'upper_r2'}];
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

upper_r2 = mean(res_e[,3]) + 1.96*std(res_e[,3]);
lower_r2 = mean(res_e[,3]) - 1.96*std(res_e[,3]);
res_e = res_e || J(nrow(res_e), 1, lower_r2) || J(nrow(res_e), 1, upper_r2);


create error_betas from res_e[colname={'b0' 'b1' 'r2' 'lower_r2' 'upper_r2'}];
append from res_e;
quit;

proc print data=error_betas (obs=10);
proc print data=pair_betas (obs=10);


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
	refline  lower_r2 / axis=x ;
	refline  upper_r2 / axis=x ;
	title 'Error Resampling';
	title2 'R2: Coefficient of Determination';
run;


proc sgplot data=pair_betas;
	histogram b0 / binwidth=10;
	title 'Pair Resampling';
	title2 'Beta 0';
run;

proc sgplot data=pair_betas;
	histogram b1 / binwidth=.1;
	title 'Pair Resampling';
	title2 'Beta 1';
run;
	
proc sgplot data=pair_betas;
	histogram r2 / binwidth=.001;
	refline  lower_r2 / axis=x ;
	refline  upper_r2 / axis=x ;
	title 'Pair Resampling';
	title2 'R2: Coefficient of Determination';
run;








* __________________________ Assignment 6: Question 3 __________________________;
	
data cdata;
set '/folders/myfolders/sasuser.v94/EKT 720/Assignment 6/cdata.sas7bdat';
run;

proc sgplot data=cdata;
	scatter x=x y=y;
run;




proc iml;
use cdata;
read all into xy;
n = nrow(xy);
x_orig = J(n,1,1)||xy[,1];
y_orig = xy[,2];


start reg;
	n=nrow(x);
	k=ncol(x);
	bh=inv(x`*x)*x`*y;
	sse=ssq(y-x*bh);
	cssy=ssq(y-(sum(y)/n));
	r2=(cssy-sse)/cssy;
finish reg;



x1 = 175;
x2 = 255;
xa = (x_orig[,2]-x1)#(x_orig[,2]>x1);
xb = (x_orig[,2]-x2)#(x_orig[,2]>x2);
x = x_orig || xa || xb;

bh=inv(x_orig`*x_orig)*x_orig`*y_orig;
e = y_orig - x_orig*bh;

do i=1 to 1000;
	in = sample(1:n, n, 'replace')`;
	y = y_orig + e[in,];
	call reg;
	results = results // (bh` || r2);
end;

nm = {'b0' 'b1' 'b2' 'b3' 'r2'};

do i=1 to ncol(results);
	low = mean(results[,i]) - 1.96*std(results[,i]);
	high = mean(results[,i]) + 1.96*std(results[,i]);
	
	print 'confidence intervals';
	print (nm[i]) low high;
	* plot;
	call Histogram(results[,i]);
end;
quit;










* __________________________ Assignment 6: Question 4 __________________________;


* __________________________   OBSERVATIONAL SEARCH   __________________________;

data cdata;
set '/folders/myfolders/sasuser.v94/EKT 720/Assignment 6/cdata.sas7bdat';
run;




proc iml;
use cdata;
read all into xy;
n = nrow(xy);
call sort(xy, {1});
x_orig = J(n,1,1)||xy[,1];
y_orig = xy[,2];
y = y_orig;


start reg;
	if inv(x`*x) ^= 0 then; do;
		k=ncol(x);
		bh=inv(x`*x)*x`*y;
		yh = x*bh;
		sse=ssq(y-x*bh);
		mse = sse/(n-k);
		cssy=ssq(y-(sum(y)/n));
		r2=(cssy-sse)/cssy;
		F = (r2/(k-1))/((1-r2)/(n-k));
	end;
	else; do;
		mse = 90009;
		r2 = 90009;
	end;
finish reg;



start def_x(x,x1,x2);
	return x || (x[,2]-x1)#(x[,2]>x1) || (x[,2]-x2)#(x[,2]>x2);
finish def_x;


* define min number of observations per segment;
m = 5;

do i=m to nrow(xy)-2*m;
	do j=(i+m) to nrow(xy)-m;
		x1 = x_orig[i,2];
		x2 = x_orig[j,2];
		x = def_x(x_orig,x1,x2);		

		call reg;
		results = results // (i || j || mse || r2);
	end;
end;


print "best breakpoints using MSE";
in = results[,3][>:<];
b1_mse = x_orig[results[in,1],2];
b2_mse = x_orig[results[in,2],2];
print in b1_mse b2_mse;


print "best breakpoints using R2";
in = results[,4][<:>];
b1_r2 = x_orig[results[in,1],2];
b2_r2 = x_orig[results[in,2],2];
print in b1_r2 b2_r2;


* create datasets with optimal breakpoints;
x = def_x(x_orig, b1_mse, b2_mse);
call reg;
d1 = xy || yh;

 
x = def_x(x_orig, b1_r2, b2_r2);
call reg;
d2 = xy || yh;

create d1 from d1[colname={'x' 'y' 'yh'}]; append from d1;
create d2 from d2[colname={'x' 'y' 'yh'}]; append from d2;
quit;

proc print data=d1 (obs=5);
proc print data=d2 (obs=5);

 

proc sgplot data=d1;
	scatter y=y x=x;
	scatter y=yh x=x;
run;










