options ls=72 nodate pageno=1 ; 

************** ***********  ASSIGNMENT 6  *********** **************;


data cdata;
set '/folders/myfolders/sasuser.v94/EKT 720/Assignment 5/cdata.sas7bdat';
run;
proc print data=cdata (obs=10);
run;
