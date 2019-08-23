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

start function_x(x);
	fx = 10 + 6.32972*x - 1.72728*(x**2) + 0.2017*(x**3) 
		- 0.00996*(x**4) + 0.00017*(x**5);
	return fx;	
finish function_x;



start MC_integration(function_x, n, x_lower_limit, x_upper_limit);

	t = do(3,8,.1)`;
	do i=1 to nrow(t);
		yt = yt // function_x(t[i]);
	end;

	approx_height_u = round(max(yt))+3;
	approx_height_l = round(min(yt))-3;
	
	results = J(n,3,999);
	do i=1 to n;
		x = rand('integer', round(x_lower_limit), round(x_upper_limit));
		y = rand('integer', approx_height_l, approx_height_u);
		y_bound = function_x(x);
		results[i,1] = x; 
		results[i,2] = y;
		if y > y_bound then do; 
			results[i,3] = 0; 
			end;
		else do; 
			results[i,3] = 1; 
			end;
	end;
	print results;
	return results;
finish MC_integration;


*********** Question 1.1 ***********;
q11 = MC_integration(function_x, 10, 3, 8);

cm = {'x' 'y' 'above'};
create q11 from q11[colname=cm];
	append from q11;
run;

print q11;





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



