options ls=72 nodate pageno=1 ; 

************** ***********  ASSIGNMENT 6  *********** **************;
test add;

data cdata;
set '/folders/myfolders/sasuser.v94/EKT 720/Assignment 6/cdata.sas7bdat';
run;
proc print data=cdata (obs=10);
run;

proc sgplot data=cdata;
	scatter x=x y=y;
run;


proc iml;
start function_x(x);
	do i=1 to nrow(x);
		y = y // 
		10 + 6.32972*x[i] - 1.72728*(x[i]**2) + 0.2017*(x[i]**3) 
		- 0.00996*(x[i]**4) + 0.00017*(x[i]**5);
	end;
	return y;	
finish function_x;

start sampler(x, y);
	
finish sampler;

* Question 1.1;
x_1 = do(3, 8, 0.1)`;
y_1 = function_x(x);





xy = x || y;
cn = {'x' 'y'};
create xy from xy[colname=cn];
	append from xy;
run;

proc sgplot data=xy;
	scatter y=y x=x;
run;



