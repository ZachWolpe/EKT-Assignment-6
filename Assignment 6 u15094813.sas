options ls=72 nodate pageno=1 ; 
************** ***************  ASSIGNMENT 6  *************** **************;

data cdata;
set '/folders/myfolders/sasuser.v94/EKT 720/Assignment 6/cdata.sas7bdat';
run;
proc print data=cdata (obs=10);
run;

proc sgplot data=cdata;
	scatter x=x y=y;
run;


************** Question 1: Logic - Monte Carlo Integration **************;

proc iml;
do i=1 to 1000;
da = da // rand('integer', 1,10); * inclusive;
end;
print da;

start function_x(x);
	do i=1 to nrow(x);
		fx = 10 + 6.32972*x[i] - 1.72728*(x[i]**2) + 0.2017*(x[i]**3) 
		- 0.00996*(x[i]**4) + 0.00017*(x[i]**5);
	end;
	return fx;	
finish function_x;

start MC_integration(fx,n,lower_limit, upper_limit);
	/* Monte Carlo Integration Sampler
		Parameters:
			fx: function to integrate
			n: number of iterations
			lower_limit: lower itegration limit
			upper_limit: upper itegration limit
	*/
	do i=1 to n;
		x = rand('integer', lower_limit, upper_limit);
		y = function_x(x);
		if 
	end;
finish MC_integration;

* Question 1.1;
x_1 = do(3, 8, 0.01)`;
y_1 = function_x(x_1);

* define sampling space;
buffer = 2.5
upper = max(y_1) + buffer;
lower = min(y_1) - buffer;
left 

print (max(y_1));


xy_1 = x_1 || y_1;
cn = {'x' 'y'};
create xy_1 from xy_1[colname=cn];
	append from xy_1;
run;

proc sgplot data=xy_1;
	scatter y=y x=x;
run;



