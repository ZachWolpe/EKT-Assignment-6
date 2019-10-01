
data a;
set '/folders/myfolders/sasuser.v94/EKT 720/Assignment 6/cdata.sas7bdat';
run;





* ________________________ Bootstrap Pairs Class ________________________;
options ls=72 nodate pageno=1 ;


data a (keep = x y);
do i=1 to 100;
	x = ranuni(567)*200 + 300;
	u = rannor(4576)*800;
	y = -15 + 25*x + u;
	output;
end;


proc reg data=a;
	model y=x / clb;
run;


proc iml;
title 'Bootstrap X,Y pairs';
use a;
read all into xy;
n=nrow(xy);
x=J(n,1,1)||xy[,1];
y=xy[,2];


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


*bootstrap routine ;
start bs1;
	u = sample(1:n,n);
	bs = in[u,];
	if pp=1 then do;
		print u bs;
	end;
finish bs1;


pp=1;
print 'Model';
call reg;

pp=0;
it=1000;


do i=1 to it;
	in=xy;
	call bs1;
	x = J(n,1,1) || bs[,1];
	y = bs[,2];
	
	call reg;
	bh_sum = bh_sum // (bh` || r2);
end;


nm = {'bh1' 'bh2' 'r2'};
print bh_sum[colname=nm];

nm={"bh1" "bh2" "r2"};
create b from bh_sum[colname=nm];
	append from bh_sum;


alpha=0.05;

p= alpha/2 || 1-alpha/2;
print p ;


call qnt1(CI, bh_sum, p);
print ((1-alpha)*100) "% confidence intervals", CI ;
quit ;



proc univariate data=b normal;
	var bh1 bh2 r2;
	histogram bh1 bh2 r2;
run;











* ________________________ Bootstrap Errors Class ________________________;
options ls=72 nodate pageno=1 ; 

data a (keep=x y);
do i=1 to 100;
	x = ranuni(567)*200 + 300;
	u = rannor(4567)*800;
	y = -15 + 25*x + u;
	output;
end;

run;

proc reg data=a;
	model y=x;
run;


proc iml;
print 'Bootstrap Errors';

use a;
read all into xy;
n = nrow(xy);
x = J(n,1,1)||xy[,1];
y = xy[,2];

* regression module;
start reg;
	k = ncol(x);
	bh = (inv(x`*x)*x`*y);
	sse = ssq(y - x*bh);
	cssy = ssq(y - (sum(y)/n));
	r2 = (cssy-sse)/cssy;
	F = (r2/(k-1))/((1-r2)/(n-k));
	if pp=1 then do; print bh r2 f; end;
finish reg;

*bootstrap routine ;
start bs1;
	u = sample(1:n,n);
	bs = in[u,];
	if pp=1 then do; print u bs; end;
finish bs1;

pp=1;
print 'Model';
call reg;
yh = x*bh;
e = y - yh;

print e;

pp=0;
it=1000;

do i=1 to it;
	in = e;
	call bs;
	x = J(n,1,1)||xy[,1];
	y = yh+bs[,1];
	*print x y;
	*print "BS" i;
	call reg;
	bh_sum = bh_sum // (bh` || r2);
end;

print bh_sum;

bm = {'bh1' 'bh2' 'r2'};
create b from bh_sum[colname=nm];
	append from bh_sum;
	
alpha=0.05l
p = alpha/2 || 1-alpha/2;
print p;

call qnt1(CI, bh_sum, p);
print ((1-alpha)*100) '% confidence intervals', CI;
quit;



proc univariate data=b normal;
	var bh1 bh2 r2;
	histogram bh1 bh2 r2;
run;

