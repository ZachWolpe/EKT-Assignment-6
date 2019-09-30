
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


start bootstrap_pairs;
	in = sample(1:n, n, 'replace')`;
	x = x_orig[in,];
	y= y_orig[in,];
	call reg;
finish bootstrap_pairs;

 
start bootstrap_errors;	
	in = sample(1:n, n, 'replace')`;
	y = y_orig + e[in,];
	call reg;
finish bootstrap_errors;



do i=1 to 10000;
	call bootstrap_pairs;
	betas = betas // bh`;
end;

create pair_betas from betas[colname={'b0 b1'}];
append from betas;


x = x_orig;
y = y_orig;
call reg;
e = y - bh;
do i=1 to 10000;
	call bootstrap_errors;
	betas_e = betas_e // bh`;
end;

create error_betas from betas_e[colname={'b0 b1'}];
append from betas_e;








e = y - bh;















