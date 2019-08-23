options ls=72 nodate pageno=1 ; 
* - - - - - - - - - - - - - - - - - ASSIGNMENT 6 - - - - - - - - - - - - - - - - - - - - - ;    


************** ******* Question 1: Logic - Monte Carlo Integration ******* **************;
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
	
	total_area = (approx_height_u - approx_height_l) * 
		(round(x_upper_limit) - round(x_lower_limit));
 	area_under_curve = mean(results[,3])*total_area;
	return area_under_curve;
finish MC_integration;


n = 10000;
* Question 1.3;
q11 = MC_integration(function_x, n, 3, 8);
print 'Area Under Curve: ' (q11);

* Question 1.2; 
q12 = MC_integration(function_x, n, 1, 10);
print 'Area Under Curve: ' (q12);

* Question 1.3;
q13 = MC_integration(function_x, n, 0, 20);
print 'Area Under Curve: ' (q13);






************** ********** Question 2: Bootstrap Regression ********** **************;
data cdata;
set '/folders/myfolders/sasuser.v94/EKT 720/Assignment 6/cdata.sas7bdat';
run;
proc print data=cdata (obs=10);
run;

proc sgplot data=cdata;
	scatter x=x y=y;
run;



















