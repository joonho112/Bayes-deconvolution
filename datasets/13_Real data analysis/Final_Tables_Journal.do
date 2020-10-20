/*****
*** This program is set up to produce the final tables for the Long-Term Subsidios paper
*****/

# delimit;
clear all;
set more off;
*set mem 2g;
set matsize 11000;
capture log close;
global strformat "%8.3f";


/*****
**** Final data set
*****/

use "Final_Data_CCT_Long_Barrera-Osorio_et_al.dta";


*Declaration of baseline variables;

global baseline_vars "
   i.s_teneviv_fill i.s_utilities_fill i.s_durables_fill s_infraest_hh_fill 
   s_age_sorteo_fill s_sexo_fill s_yrs_fill i.grade_fill s_single_fill s_edadhead_fill
   s_yrshead_fill s_tpersona_fill s_num18_fill i.s_estrato_fill s_puntaje_fill s_ingtotal_fill

   s_teneviv_miss s_utilities_miss s_durables_miss s_infraest_hh_miss
   s_age_sorteo_miss s_sexo_miss s_yrs_miss grade_miss s_single_miss s_edadhead_miss
   s_yrshead_miss s_tpersona_miss s_num18_miss s_estrato_miss s_puntaje_miss s_ingtotal_miss
   ";

/****
*** Table 1: Match Rates
****/

gen row_num = _n;
gen str30 col0 = "";
forvalues i = 1(1)10 {;
   gen str20 col`i' = "";
   };

local row = 0;

foreach var in enrolled_2006 enrolled_2007 enrolled_2008 icfes spadies spadies_longterm {;
   local row = `row' + 1;
   replace col0 = "`var'" if row_num == `row';
   summarize `var' if suba == 0;
   replace col1 = string(r(mean), "$strformat")  if row_num == `row';
   summarize `var' if suba == 1;
   replace col2 = string(r(mean), "$strformat")  if row_num == `row';
   summarize `var' if suba < 2;
   replace col3 = string(r(mean), "$strformat")  if row_num == `row';

   };

  *Grades 9-11;

local row = `row' + 2;
replace col0 = "icfes" if row_num == `row';
summarize icfes if grade > 8 & suba == 0;
replace col1 = string(r(mean), "$strformat")  if row_num == `row';
summarize icfes if grade > 8 & suba == 1;
replace col2 = string(r(mean), "$strformat")  if row_num == `row';
summarize icfes if grade > 8 & suba < 2;
replace col3 = string(r(mean), "$strformat")  if row_num == `row';

local row = `row' + 1;
replace col0 = "spadies" if row_num == `row';
summarize spadies if grade > 8 & suba == 0;
replace col1 = string(r(mean), "$strformat")  if row_num == `row';
summarize spadies if grade > 8 & suba == 1;
replace col2 = string(r(mean), "$strformat")  if row_num == `row';
summarize spadies if grade > 8 & suba < 2;
replace col3 = string(r(mean), "$strformat")  if row_num == `row';

local row = `row' + 1;
replace col0 = "spadies_longterm" if row_num == `row';
summarize spadies_longterm if grade > 8 & suba == 0;
replace col1 = string(r(mean), "$strformat")  if row_num == `row';
summarize spadies_longterm if grade > 8 & suba == 1;
replace col2 = string(r(mean), "$strformat")  if row_num == `row';
summarize spadies_longterm if grade > 8 & suba < 2;
replace col3 = string(r(mean), "$strformat")  if row_num == `row';


*Grades 6-8;
local row = `row' + 2;
replace col0 = "icfes" if row_num == `row';
summarize icfes if grade < 9 & suba == 0;
replace col1 = string(r(mean), "$strformat")  if row_num == `row';

local row = `row' + 1;
replace col0 = "spadies" if row_num == `row';
summarize spadies if grade < 9 & suba == 0;
replace col1 = string(r(mean), "$strformat")  if row_num == `row';

local row = `row' + 1;
replace col0 = "spadies_longterm" if row_num == `row';
summarize spadies_longterm if grade < 9 & suba == 0;
replace col1 = string(r(mean), "$strformat")  if row_num == `row';
  
sort row_num;
outsheet col0 col1-col4 using "Table_01.xls" if _n <= `row'+1, replace comma;
drop col* row_num;

/****
*** Table 2: Balance in availability of matching info;
****/

xi: reg any_ID T1_treat T2_treat if suba == 0, cluster(school_code);
sum any_ID if T1_treat == 0 & T2_treat == 0 & suba == 0;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2 using "Table_02_Panel_A.xls", keep(T1_treat T2_treat) addstat("Mean", `ymean', "T1=T2p", `pval') nocons se replace  bdec(3) sdec(3) adec(2);

xi: reg full_last_name T1_treat T2_treat if suba == 0, cluster(school_code);
sum full_last_name if T1_treat == 0 & T2_treat == 0 & suba == 0;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2 using "Table_02_Panel_A.xls", keep(T1_treat T2_treat) addstat("Mean", `ymean', "T1=T2p", `pval') nocons se append  bdec(3) sdec(3) adec(2);

xi: reg any_ID T3_treat if suba == 1, cluster(school_code);
su any_ID if T3_treat == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_02_Panel_B.xls", keep(T3_treat) addstat("Mean", `ymean') nocons se  replace  bdec(3) sdec(3)  adec(2);

xi: reg full_last_name T3_treat if suba == 1, cluster(school_code);
su full_last_name if T3_treat == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_02_Panel_B.xls", keep(T3_treat) addstat("Mean", `ymean') nocons se  append  bdec(3) sdec(3)  adec(2);

/****
*** Table 3: Secondary Education outcomes (On-Time Enrollment and Taking the ICFES exam);
*** For format purposes, we are dividing the table in two panels, columns 1 to 8, panel A; rest of columns panel B;  
****/

/****
*** On-Time Enrollment;
****/

*SC: Socio-Econ Controls and FE's;
xi: areg on_time T1_treat T2_treat $baseline_vars if suba == 0, cluster(school_code) absorb(school_code); 
su on_time if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2 using "Table_03.xls",  keep(T1_treat T2_treat) se nocons replace addstat("Mean", `ymean', "T1=T2p", `pval') bdec(3) sdec(3)  adec(2);

*Suba: Socio-Econ Controls and FE's;
xi: areg on_time T3_treat $baseline_vars if suba == 1, cluster(school_code) absorb(school_code); 
su on_time if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_03.xls", keep(T3_treat) se nocons append addstat("Mean", `ymean') bdec(3) sdec(3)  adec(2);

*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg on_time T1_treat T2_treat T3_treat $baseline_vars suba if grade > 8, cluster(school_code) absorb(school_code); 
su on_time if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval1 = r(p);
test T1_treat == T3_treat;
local pval2 = r(p);
test T2_treat == T3_treat;
local pval3 = r(p);
outreg2 using "Table_03.xls", keep(T1_treat T2_treat T3_treat) se nocons append 
             addstat("Mean", `ymean', "T1=T2p", `pval1', "T1=T3p", `pval2', "T2=T3p", `pval3') bdec(3) sdec(3)  adec(2);

   
*Grades 6-8: Socio-Econ Controls and FE's;
xi: areg on_time T1_treat T2_treat $baseline_vars if suba == 0 & grade < 9, cluster(school_code) absorb(school_code); 
su on_time if treatment == 0 & suba == 0 & grade < 9;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2 using "Table_03.xls", keep(T1_treat T2_treat ) se nocons append addstat("Mean", `ymean', "T1=T2p", `pval') bdec(3) sdec(3)  adec(2);


/****
*** Heldback;
****/


*SC: Socio-Econ Controls and FE's;
xi: areg heldback T1_treat T2_treat $baseline_vars if suba == 0, cluster(school_code) absorb(school_code); 
su heldback if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2 using "Table_03.xls",  keep(T1_treat T2_treat) se nocons append addstat("Mean", `ymean', "T1=T2p", `pval') bdec(3) sdec(3)  adec(2);

*Suba: Socio-Econ Controls and FE's;
xi: areg heldback T3_treat $baseline_vars if suba == 1, cluster(school_code) absorb(school_code); 
su heldback if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_03.xls", keep(T3_treat) se nocons append addstat("Mean", `ymean') bdec(3) sdec(3)  adec(2);

/****
*** Dropout;
****/

*SC: Socio-Econ Controls and FE's;
xi: areg dropout T1_treat T2_treat $baseline_vars if suba == 0, cluster(school_code) absorb(school_code); 
su dropout if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2 using "Table_03.xls",  keep(T1_treat T2_treat) se nocons append addstat("Mean", `ymean', "T1=T2p", `pval') bdec(3) sdec(3)  adec(2);

*Suba: Socio-Econ Controls and FE's;
xi: areg dropout T3_treat $baseline_vars if suba == 1, cluster(school_code) absorb(school_code); 
su dropout if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_03.xls", keep(T3_treat) se nocons append addstat("Mean", `ymean') bdec(3) sdec(3)  adec(2);

/****
*** Icfes;
****/

*SC: Socio-Econ Controls and FE's;
xi: areg icfes T1_treat T2_treat $baseline_vars if suba == 0, cluster(school_code) absorb(school_code); 
su icfes if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2  using "Table_03.xls", keep(T1_treat T2_treat) se nocons append addstat("Mean", `ymean', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);


*Suba: Socio-Econ Controls and FE's;
xi: areg icfes T3_treat $baseline_vars if suba == 1, cluster(school_code) absorb(school_code); 
su icfes if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_03.xls", keep(T3_treat) se nocons append addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);


*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg icfes T1_treat T2_treat T3_treat $baseline_vars suba if grade > 8, cluster(school_code) absorb(school_code); 
su icfes if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval1 = r(p);
test T1_treat == T3_treat;
local pval2 = r(p);
test T2_treat == T3_treat;
local pval3 = r(p);
outreg2 using "Table_03.xls", keep(T1_treat T2_treat T3_treat) se append nocons addstat("Mean", `ymean', 
    "T1=T2p", `pval1', "T1=T3p", `pval2', "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);

   
*Grades 6-8: Socio-Econ Controls and FE's;
xi: areg icfes T1_treat T2_treat $baseline_vars if suba == 0 & grade < 9, cluster(school_code) absorb(school_code); 
su icfes if treatment == 0 & suba == 0 & grade < 9;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2  using "Table_03.xls", keep(T1_treat T2_treat) se nocons append addstat("Mean", `ymean', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

/****
*** Table 4: SPADIES Medium-term;
*** For format purposes, we are dividing the table in two panels, columns 1 to 4, panel A; rest of columns panel B;
****/

*SC: Socio-Econ Controls and FE's;
xi: areg spadies T1_treat T2_treat $baseline_vars if suba == 0, cluster(school_code) absorb(school_code); 
su spadies if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2  using "Table_04.xls", keep(T1_treat T2_treat) se nocons replace addstat("Mean", `ymean', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

*Suba: Socio-Econ Controls and FE's;
xi: areg spadies T3_treat $baseline_vars if suba == 1, cluster(school_code) absorb(school_code); 
su spadies if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_04.xls", keep(T3_treat) se nocons append addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);

*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg spadies T1_treat T2_treat T3_treat $baseline_vars suba if grade > 8, cluster(school_code) absorb(school_code); 
su spadies if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval1 = r(p);
test T1_treat == T3_treat;
local pval2 = r(p);
test T2_treat == T3_treat;
local pval3 = r(p);
outreg2 using "Table_04.xls", keep(T1_treat T2_treat T3_treat) se append nocons addstat("Mean", `ymean', 
   "T1=T2p", `pval1',  "T1=T3p", `pval2', "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);

*Grades 6-8: Socio-Econ Controls and FE's;
xi: areg spadies T1_treat T2_treat $baseline_vars if suba == 0 & grade < 9, cluster(school_code) absorb(school_code); 
su spadies if treatment == 0 & suba == 0 & grade < 9;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2  using "Table_04.xls", keep(T1_treat T2_treat) se nocons append addstat("Mean", `ymean', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

*Grades 9-11: effects by type of institution; 

foreach var in university vocational not_reported {;

	xi: areg `var' T1_treat T2_treat $baseline_vars if suba == 0 & grade > 8, cluster(school_code) absorb(school_code); 
	su `var' if treatment == 0 & suba == 0 & grade > 8;
	local ymean = r(mean);
	test T1_treat == T2_treat;
	local pval = r(p);
	outreg2 using "Table_04.xls", keep(T1_treat T2_treat) se nocons append addstat("Mean", `ymean', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);
	};

foreach var in university vocational not_reported {;

	xi: areg `var' T3_treat $baseline_vars if suba == 1 & grade > 8, cluster(school_code) absorb(school_code); 
	su `var' if treatment == 0 & suba == 1  & grade > 8;
	local ymean = r(mean);
	outreg2 using "Table_04.xls", keep(T3_treat) se nocons append addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);
	};

	
/****
*** Table 5: SPADIES, long-term effects; 
*** For format purposes, we are dividing the table in three panels, columns 1 to 4, panel A; columns 5 to 8, panel B; rest of columns panel C;
****/

/****
*** Enrollment; 
****/

*SC: Socio-Econ Controls and FE's;
xi: areg spadies_longterm T1_treat T2_treat $baseline_vars if suba == 0, cluster(school_code) absorb(school_code); 
su spadies_longterm if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2  using "Table_05.xls", keep(T1_treat T2_treat) se nocons replace addstat("Mean", `ymean', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

*Suba: Socio-Econ Controls and FE's;
xi: areg spadies_longterm T3_treat $baseline_vars if suba == 1, cluster(school_code) absorb(school_code); 
su spadies_longterm if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_05.xls", keep(T3_treat) se nocons append addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);

*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg spadies_longterm T1_treat T2_treat T3_treat $baseline_vars suba if grade > 8, cluster(school_code) absorb(school_code); 
su spadies_longterm if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval1 = r(p);
test T1_treat == T3_treat;
local pval2 = r(p);
test T2_treat == T3_treat;
local pval3 = r(p);
outreg2 using "Table_05.xls", keep(T1_treat T2_treat T3_treat) se append nocons addstat("Mean", `ymean', 
   "T1=T2p", `pval1', "T1=T3p", `pval2', "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);

*Grades 6-8: Socio-Econ Controls and FE's;
xi: areg spadies_longterm T1_treat T2_treat $baseline_vars if suba == 0 & grade < 9, cluster(school_code) absorb(school_code); 
su spadies_longterm if treatment == 0 & suba == 0 & grade < 9;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2  using "Table_05.xls", keep(T1_treat T2_treat) se nocons append addstat("Mean", `ymean', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

/****
*** On-time enrollment, unconditional; 
****/

*SC: Socio-Econ Controls and FE's;
xi: areg on_time_uncond_longterm T1_treat T2_treat $baseline_vars if suba == 0, cluster(school_code) absorb(school_code); 
su on_time_uncond_longterm if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2  using "Table_05.xls", keep(T1_treat T2_treat) se nocons append addstat("Mean", `ymean', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

*Suba: Socio-Econ Controls and FE's;
xi: areg on_time_uncond_longterm T3_treat $baseline_vars if suba == 1, cluster(school_code) absorb(school_code); 
su on_time_uncond_longterm if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_05.xls", keep(T3_treat) se nocons append addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);

*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg on_time_uncond_longterm T1_treat T2_treat T3_treat $baseline_vars suba if grade > 8, cluster(school_code) absorb(school_code); 
su on_time_uncond_longterm if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval1 = r(p);
test T1_treat == T3_treat;
local pval2 = r(p);
test T2_treat == T3_treat;
local pval3 = r(p);
outreg2 using "Table_05.xls", keep(T1_treat T2_treat T3_treat) se append nocons addstat("Mean", `ymean', 
   "T1=T2p", `pval1', "T1=T3p", `pval2', "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);

*Grades 6-8: Socio-Econ Controls and FE's;
xi: areg on_time_uncond_longterm T1_treat T2_treat $baseline_vars if suba == 0 & grade < 9, cluster(school_code) absorb(school_code); 
su on_time_uncond_longterm if treatment == 0 & suba == 0 & grade < 9;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2  using "Table_05.xls", keep(T1_treat T2_treat) se nocons append addstat("Mean", `ymean', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

/****
*** Tertiary graduation, unconditional; 
****/

*SC: Socio-Econ Controls and FE's;
xi: areg grad_uncond_longterm T1_treat T2_treat $baseline_vars if suba == 0, cluster(school_code) absorb(school_code); 
su grad_uncond_longterm if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2  using "Table_05.xls", keep(T1_treat T2_treat) se nocons append addstat("Mean", `ymean', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

*Suba: Socio-Econ Controls and FE's;
xi: areg grad_uncond_longterm T3_treat $baseline_vars if suba == 1, cluster(school_code) absorb(school_code); 
su grad_uncond_longterm if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_05.xls", keep(T3_treat) se nocons append addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);

*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg grad_uncond_longterm T1_treat T2_treat T3_treat $baseline_vars suba if grade > 8, cluster(school_code) absorb(school_code); 
su grad_uncond_longterm if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval1 = r(p);
test T1_treat == T3_treat;
local pval2 = r(p);
test T2_treat == T3_treat;
local pval3 = r(p);
outreg2 using "Table_05.xls", keep(T1_treat T2_treat T3_treat) se append nocons addstat("Mean", `ymean', 
   "T1=T2p", `pval1', "T1=T3p", `pval2',  "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);

*Grades 6-8: Socio-Econ Controls and FE's;
xi: areg grad_uncond_longterm T1_treat T2_treat $baseline_vars if suba == 0 & grade < 9, cluster(school_code) absorb(school_code); 
su grad_uncond_longterm if treatment == 0 & suba == 0 & grade < 9;
local ymean = r(mean);
test T1_treat == T2_treat;
local pval = r(p);
outreg2  using "Table_05.xls", keep(T1_treat T2_treat) se nocons append addstat("Mean", `ymean', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

/*****
*** Table 6: Joint test of of effects, Original estimation paper;
*****/

*Basic and Savings individually, medium term;
xi: reg on_time T1_treat T2_treat $baseline_vars i.school_code if suba == 0;
estimates store on_time;
xi: reg icfes T1_treat T2_treat $baseline_vars i.school_code if suba == 0;
estimates store icfes;
xi: reg spadies T1_treat T2_treat $baseline_vars i.school_code if suba == 0; 
estimates store spadies;

suest on_time icfes spadies, cluster(school_code);

test [on_time_mean]T1_treat [icfes_mean]T1_treat [spadies_mean]T1_treat;

local t1chi = r(chi2);
local t1pval = r(p);

outreg2 using "Table_06.xls", keep (T1_treat) se nocons replace addstat("T1 MT", `t1chi', "T1p MT", `t1pval')  bdec(3) sdec(3)  adec(3); 

test [on_time_mean]T2_treat [icfes_mean]T2_treat [spadies_mean]T2_treat;

local t2chi = r(chi2);
local t2pval = r(p);

outreg2 using "Table_06.xls", keep (T2_treat) se nocons append addstat("T2 MT", `t2chi', "T2p MT", `t2pval')  bdec(3) sdec(3)  adec(3); 

estimates clear;

*Basic and Savings individually, long term only;
xi: reg spadies_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0;
estimates store spadieslt;
xi: reg on_time_uncond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0;
estimates store ontimelt;
xi: reg grad_uncond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0; 
estimates store gradlt;

suest spadieslt ontimelt gradlt, cluster(school_code);

test [spadieslt_mean]T1_treat [ontimelt_mean]T1_treat [gradlt_mean]T1_treat;

local t1chi = r(chi2);
local t1pval = r(p);

outreg2 using "Table_06.xls", keep (T1_treat) se nocons append addstat("T1 MT", `t1chi', "T1p MT", `t1pval')  bdec(3) sdec(3)  adec(3); 

test [spadieslt_mean]T2_treat [ontimelt_mean]T2_treat [gradlt_mean]T2_treat;

local t2chi = r(chi2);
local t2pval = r(p);

outreg2 using "Table_06.xls", keep (T1_treat) se nocons append addstat("T2 MT", `t2chi', "T2p MT", `t2pval')  bdec(3) sdec(3)  adec(3); 

estimates clear;

*Basic and Savings individually, both medium and long term;

xi: reg on_time T1_treat T2_treat $baseline_vars i.school_code if suba == 0;
estimates store on_time;
xi: reg icfes T1_treat T2_treat $baseline_vars i.school_code if suba == 0;
estimates store icfes;
xi: reg spadies_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0;
estimates store spadieslt;
xi: reg on_time_uncond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0;
estimates store ontimelt;
xi: reg grad_uncond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0; 
estimates store gradlt;

suest on_time icfes spadieslt ontimelt gradlt, cluster(school_code);

test [on_time_mean]T1_treat [icfes_mean]T1_treat [spadieslt_mean]T1_treat [ontimelt_mean]T1_treat [gradlt_mean]T1_treat;

local t1chi = r(chi2);
local t1pval = r(p);

outreg2 using "Table_06.xls", keep (T1_treat) se nocons append addstat("T1 MT", `t1chi', "T1p MT", `t1pval')  bdec(3) sdec(3)  adec(3); 

test [on_time_mean]T2_treat [icfes_mean]T2_treat [spadieslt_mean]T2_treat [ontimelt_mean]T2_treat [gradlt_mean]T2_treat;

local t2chi = r(chi2);
local t2pval = r(p);

outreg2 using "Table_06.xls", keep (T2_treat) se nocons append addstat("T2 MT", `t2chi', "T2p MT", `t2pval')  bdec(3) sdec(3)  adec(3); 

estimates clear;


*Tertiary, medium term;
xi: reg on_time T3_treat $baseline_vars i.school_code if suba == 1;
estimates store on_time;
xi: reg icfes T3_treat $baseline_vars i.school_code if suba == 1;
estimates store icfes;
xi: reg spadies T3_treat $baseline_vars i.school_code if suba == 1; 
estimates store spadies;

suest on_time icfes spadies, cluster(school_code);

test [on_time_mean]T3_treat [icfes_mean]T3_treat [spadies_mean]T3_treat;

local t3chi = r(chi2);
local t3pval = r(p);

outreg2 using "Table_06.xls", keep (T3_treat) se nocons append addstat("T3 MT", `t3chi', "T3p MT", `t3pval')  bdec(3) sdec(3)  adec(3); 

estimates clear;

*Tertiary, long term only;

xi: reg spadies_longterm T3_treat $baseline_vars i.school_code if suba == 1; 
estimates store spadieslt;
xi: reg on_time_uncond_longterm T3_treat $baseline_vars i.school_code if suba == 1; 
estimates store ontimelt;
xi: reg grad_uncond_longterm T3_treat $baseline_vars i.school_code if suba == 1; 
estimates store gradlt;

suest spadieslt ontimelt gradlt, cluster(school_code);

test [spadieslt_mean]T3_treat [ontimelt_mean]T3_treat [gradlt_mean]T3_treat;

local t3chi = r(chi2);
local t3pval = r(p);

outreg2 using "Table_06.xls", keep (T3_treat) se nocons append addstat("T3 MT", `t3chi', "T3p MT", `t3pval')  bdec(3) sdec(3)  adec(3); 

estimates clear;

*Tertiary, both medium and long term;
xi: reg on_time T3_treat $baseline_vars i.school_code if suba == 1;
estimates store on_time;
xi: reg icfes T3_treat $baseline_vars i.school_code if suba == 1;
estimates store icfes;
xi: reg spadies_longterm T3_treat $baseline_vars i.school_code if suba == 1; 
estimates store spadieslt;
xi: reg on_time_uncond_longterm T3_treat $baseline_vars i.school_code if suba == 1; 
estimates store ontimelt;
xi: reg grad_uncond_longterm T3_treat $baseline_vars i.school_code if suba == 1; 
estimates store gradlt;

suest on_time icfes spadieslt ontimelt gradlt, cluster(school_code);

test [on_time_mean]T3_treat [icfes_mean]T3_treat [spadieslt_mean]T3_treat [ontimelt_mean]T3_treat [gradlt_mean]T3_treat;

local t3chi = r(chi2);
local t3pval = r(p);

outreg2 using "Table_06.xls", keep (T3_treat) se nocons append addstat("T3 MT", `t3chi', "T3p MT", `t3pval')  bdec(3) sdec(3)  adec(3); 

estimates clear;

*Basic versus Savings, medium term -- all grades spadies w/ interact;
estimates clear;
xi: reg on_time T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0;
estimates store on_time;
xi: reg icfes T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0;
estimates store icfes;
xi: reg spadies T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0; 
estimates store spadies;

suest on_time icfes spadies, cluster(school_code);
test [on_time_mean]T1_treat [on_time_mean]T1_treat_upper
	[icfes_mean]T1_treat [icfes_mean]T1_treat_upper
	[spadies_mean]T1_treat [spadies_mean]T1_treat_upper;

local t12chi = r(chi2);
local t12pval = r(p);

outreg2 using "Table_06.xls", keep (T1_treat) se nocons append addstat("T1T2 MT", `t12chi', "T1T2p MT", `t12pval')  bdec(3) sdec(3)  adec(3); 

	
*Test for difference in effects by grade (in paper but not table);
test [on_time_mean]T1_treat_upper [icfes_mean]T1_treat_upper [spadies_mean]T1_treat_upper;

*Basic versus Savings, only long term -- all grades spadies w/ interact;
estimates clear;   
xi: reg spadies_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0; 
estimates store spadieslt;
xi: reg on_time_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0; 
estimates store ontimelt;
xi: reg grad_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0; 
estimates store gradlt;

suest spadieslt ontimelt gradlt, cluster(school_code);

test [spadieslt_mean]T1_treat [spadieslt_mean]T1_treat_upper
   	[ontimelt_mean]T1_treat [ontimelt_mean]T1_treat_upper
	[gradlt_mean]T1_treat [gradlt_mean]T1_treat_upper;

local t12chi = r(chi2);
local t12pval = r(p);

outreg2 using "Table_06.xls", keep (T3_treat) se nocons append addstat("T1T2 LT", `t12chi', "T1T2p LT", `t12pval')  bdec(3) sdec(3)  adec(3); 


*Basic versus Savings, both medium and long term -- all grades spadies w/ interact;
estimates clear;   
xi: reg on_time T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0;
estimates store on_time;
xi: reg icfes T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0;
estimates store icfes;
xi: reg spadies_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0; 
estimates store spadieslt;
xi: reg on_time_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0; 
estimates store ontimelt;
xi: reg grad_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0; 
estimates store gradlt;

suest on_time icfes spadieslt ontimelt gradlt, cluster(school_code);

test [on_time_mean]T1_treat [on_time_mean]T1_treat_upper
	[icfes_mean]T1_treat [icfes_mean]T1_treat_upper
	[spadieslt_mean]T1_treat [spadieslt_mean]T1_treat_upper
   	[ontimelt_mean]T1_treat [ontimelt_mean]T1_treat_upper
	[gradlt_mean]T1_treat [gradlt_mean]T1_treat_upper;

local t12chi = r(chi2);
local t12pval = r(p);

outreg2 using "Table_06.xls", keep (T3_treat) se nocons append addstat("T1T2 LT", `t12chi', "T1T2p LT", `t12pval')  bdec(3) sdec(3)  adec(3); 


*Test for difference in effects by grade (in paper but not table);
*Footnote 22;
test [on_time_mean]T1_treat_upper [icfes_mean]T1_treat_upper [spadieslt_mean]T1_treat_upper [ontimelt_mean]T1_treat_upper [gradlt_mean]T1_treat_upper;

estimates clear;

*Basic versus Savings, medium term -- upper secondary;
xi: reg on_time T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0;
estimates store on_time;
xi: reg icfes T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0;
estimates store icfes;
xi: reg spadies T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0; 
estimates store spadies;

suest on_time icfes spadies, cluster(school_code);

test [on_time_mean]T1_treat
	[icfes_mean]T1_treat
	[spadies_mean]T1_treat;

local t12chi = r(chi2);
local t12pval = r(p);

outreg2 using "Table_06.xls", keep (T3_treat) se nocons append addstat("T1T2 Up MT", `t12chi', "T1T2p Up MT", `t12pval')  bdec(3) sdec(3)  adec(3); 
	
	
estimates clear; 

*Basic versus Savings, only long term -- upper secondary;
xi: reg spadies_longterm T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0; 
estimates store spadieslt;
xi: reg on_time_uncond_longterm T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0; 
estimates store ontimelt;
xi: reg grad_uncond_longterm T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0; 
estimates store gradlt;

suest spadieslt ontimelt gradlt, cluster(school_code);

test [spadieslt_mean]T1_treat
	[ontimelt_mean]T1_treat
	[gradlt_mean]T1_treat;

local t12chi = r(chi2);
local t12pval = r(p);

outreg2 using "Table_06.xls", keep (T3_treat) se nocons append addstat("T1T2 Up LT", `t12chi', "T1T2p Up LT", `t12pval')  bdec(3) sdec(3)  adec(3); 


estimates clear;   

*Basic versus Savings, both medium and long term -- upper secondary;
xi: reg on_time T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0;
estimates store on_time;
xi: reg icfes T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0;
estimates store icfes;
xi: reg spadies_longterm T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0; 
estimates store spadieslt;
xi: reg on_time_uncond_longterm T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0; 
estimates store ontimelt;
xi: reg grad_uncond_longterm T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0; 
estimates store gradlt;

suest on_time icfes spadieslt ontimelt gradlt, cluster(school_code);

test [on_time_mean]T1_treat
	[icfes_mean]T1_treat
	[spadieslt_mean]T1_treat
	[ontimelt_mean]T1_treat
	[gradlt_mean]T1_treat;
local t12chi = r(chi2);
local t12pval = r(p);

outreg2 using "Table_06.xls", keep (T3_treat) se nocons append addstat("T1T2 Up LT", `t12chi', "T1T2p Up LT", `t12pval')  bdec(3) sdec(3)  adec(3); 

estimates clear;   

*Basic versus Savings, medium term -- lower secondary;
xi: reg on_time T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0;
estimates store on_time;
xi: reg icfes T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0;
estimates store icfes;
xi: reg spadies T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0; 
estimates store spadies;

suest on_time icfes spadies, cluster(school_code);

test [on_time_mean]T1_treat
	[icfes_mean]T1_treat
	[spadies_mean]T1_treat;
local t12chi = r(chi2);
local t12pval = r(p);

outreg2 using "Table_06.xls", keep (T3_treat) se nocons append addstat("T1T2 Low MT", `t12chi', "T1T2p Low MT", `t12pval')  bdec(3) sdec(3)  adec(3); 

estimates clear;   

*Basic versus Savings, only long term -- lower secondary;
xi: reg spadies_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0; 
estimates store spadieslt;
xi: reg on_time_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0; 
estimates store ontimelt;
xi: reg grad_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0; 
estimates store gradlt;

suest spadieslt ontimelt gradlt, cluster(school_code);

test [spadieslt_mean]T1_treat
	[ontimelt_mean]T1_treat
	[gradlt_mean]T1_treat;

local t12chi = r(chi2);
local t12pval = r(p);

outreg2 using "Table_06.xls", keep (T3_treat) se nocons append addstat("T1T2 Low LT", `t12chi', "T1T2p Low LT", `t12pval')  bdec(3) sdec(3)  adec(3); 
		
estimates clear;   

*Basic versus Savings, both medium and long term -- lower secondary;
xi: reg on_time T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0;
estimates store on_time;
xi: reg icfes T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0;
estimates store icfes;
xi: reg spadies_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0; 
estimates store spadieslt;
xi: reg on_time_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0; 
estimates store ontimelt;
xi: reg grad_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0; 
estimates store gradlt;

suest on_time icfes spadieslt ontimelt gradlt, cluster(school_code);

test [on_time_mean]T1_treat
	[icfes_mean]T1_treat
	[spadieslt_mean]T1_treat
	[ontimelt_mean]T1_treat
	[gradlt_mean]T1_treat;

local t12chi = r(chi2);
local t12pval = r(p);

outreg2 using "Table_06.xls", keep (T3_treat) se nocons append addstat("T1T2 Low LT", `t12chi', "T1T2p Low LT", `t12pval')  bdec(3) sdec(3)  adec(2); 
		
estimates clear;   


*Secondary versus terciary: difference between treatments, Medium term;

xi: reg on_time treatment T1_treat T2_treat $baseline_vars suba i.school_code if grade > 8;  
estimates store on_time;
xi: reg icfes treatment T1_treat T2_treat $baseline_vars suba i.school_code if grade > 8; 
estimates store icfes;
xi: reg spadies treatment T1_treat T2_treat $baseline_vars suba i.school_code if grade > 8; 
estimates store spadies;

suest on_time icfes spadies, cluster(school_code);

test [on_time_mean]T2_treat
	[icfes_mean]T2_treat
	[spadies_mean]T2_treat;

local t23chi = r(chi2);
local t23pval = r(p);

outreg2 using "Table_06.xls", keep (treatment) se nocons append addstat("T2T3 MT", `t23chi', "T2T3p MT", `t23pval')  bdec(3) sdec(3)  adec(3); 


*Secondary versus terciary only long term;

xi: reg spadies_longterm treatment T1_treat T2_treat $baseline_vars suba i.school_code if grade > 8 ;  
estimates store spadieslt;
xi: reg on_time_uncond_longterm treatment T1_treat T2_treat $baseline_vars suba i.school_code if grade > 8 ; 
estimates store ontimelt;
xi: reg grad_uncond_longterm treatment T1_treat T2_treat $baseline_vars suba i.school_code if grade > 8 ; 
estimates store gradlt;

suest spadieslt ontimelt gradlt, cluster(school_code);

test [spadieslt_mean]T2_treat
	 [ontimelt_mean]T2_treat
	[gradlt_mean]T2_treat;

local t23chi = r(chi2);
local t23pval = r(p);

outreg2 using "Table_06.xls", keep (treatment) se nocons append addstat("T2T3 MT", `t23chi', "T2T3p MT", `t23pval')  bdec(3) sdec(3)  adec(3); 

*Secondary versus terciary both medium and long term;
xi: reg on_time treatment T1_treat T2_treat $baseline_vars suba i.school_code if grade > 8 ;  
estimates store on_time;
xi: reg icfes treatment T1_treat T2_treat $baseline_vars suba i.school_code if grade > 8; 
estimates store icfes;
xi: reg spadies_longterm treatment T1_treat T2_treat $baseline_vars suba i.school_code if grade > 8 ;  
estimates store spadieslt;
xi: reg on_time_uncond_longterm treatment T1_treat T2_treat $baseline_vars suba i.school_code if grade > 8 ; 
estimates store ontimelt;
xi: reg grad_uncond_longterm treatment T1_treat T2_treat $baseline_vars suba i.school_code if grade > 8 ; 
estimates store gradlt;

suest on_time icfes spadieslt ontimelt gradlt, cluster(school_code);

test  [on_time_mean]T2_treat
     [icfes_mean]T2_treat
     [spadieslt_mean]T2_treat
	 [ontimelt_mean]T2_treat
	[gradlt_mean]T2_treat;

local t23chi = r(chi2);
local t23pval = r(p);

outreg2 using "Table_06.xls", keep (treatment) se nocons append addstat("T2T3 MT", `t23chi', "T2T3p MT", `t23pval')  bdec(3) sdec(3)  adec(3); 


/*****
*** APPENDIXES;
*****/


/*****
*A. Internal Validity; 
*****/

/*****
*B. Heterogenous effects ; 
*** For format purposes, we are dividing the two tables of heterogeniery effects in two panels, columns 1 to 6, panel A; rest of columns panel B;
*****/


/*****
*** Heterogenous effects by sisben score;
*The other way: Standarization of income;
sum s_ingtotal_fill;
replace s_ingtotal_fill=(s_ingtotal_fill-r(mean))/r(sd); 

gen T1_inco=T1_treat*s_ingtotal_fill;
gen T2_inco=T2_treat*s_ingtotal_fill;
gen T3_inco=T3_treat*s_ingtotal_fill;

*Results are very similar;
*****/


*Generation interaction treatment with sisben score;

gen T1_sisb=T1_treat*s_puntaje_fill;
gen T2_sisb=T2_treat*s_puntaje_fill;
gen T3_sisb=T3_treat*s_puntaje_fill;

*****;
*OUTCOME;
*A. on time enrollemtne;
*****;

*SC: Socio-Econ Controls and FE's;
xi: areg on_time T1_treat T2_treat T1_sisb T2_sisb s_puntaje_fill $baseline_vars if suba == 0, cluster(school_code) a(school_code); 
su on_time if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_sisb == T2_sisb;
local fstat = r(F);
local pval = r(p);
outreg2 using "Table_Appen_B1.xls", keep(T1_treat T2_treat T1_sisb T2_sisb s_puntaje_fill) se nocons excel replace addstat("Mean", `ymean', "T1_sisb=T2_sisb", `fstat', "T1_sisb=T2p_sisb", `pval') bdec(3) sdec(3)  adec(2);


*Suba: Socio-Econ Controls and FE's;
xi: areg on_time T3_treat T3_sisb s_puntaje_fill $baseline_vars if suba == 1, cluster(school_code) a(school_code); 
su on_time if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_Appen_B1.xls", keep(T3_treat T3_sisb s_puntaje_fill) se nocons excel append addstat("Mean", `ymean') bdec(3) sdec(3)  adec(2);


*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg on_time T1_treat T2_treat T3_treat T1_sisb T2_sisb T3_sisb s_puntaje_fill $baseline_vars suba if grade > 8, cluster(school_code) a(school_code); 
su on_time if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_sisb == T2_sisb;
local fstat1 = r(F);
local pval1 = r(p);
test T1_sisb == T3_sisb;
local fstat2 = r(F);
local pval2 = r(p);
test T2_sisb == T3_sisb;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_B1.xls", keep(T1_treat T2_treat T3_treat T1_sisb T2_sisb T3_sisb s_puntaje_fill) se nocons excel append addstat("Mean", `ymean', 
   "T1_sisb=T2_sisb", `fstat1', "T1_sisb=T2p_sisb", `pval1',  "T1_sisb=T3_sisb", `fstat2', "T1_sisb=T3p_sisb", `pval2',  "T2=_sexo=T3_sisb", `fstat3', "T2_sisb=T3p_sisb", `pval3') bdec(3) sdec(3)  adec(2);

*****;
*OUTCOME;
*B. Icfes;
*****;

*SC: Socio-Econ Controls and FE's;
xi: areg icfes T1_treat T2_treat T1_sisb T2_sisb s_puntaje_fill $baseline_vars if suba == 0, cluster(school_code) a(school_code); 
su icfes if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_sisb == T2_sisb;
local fstat = r(F);
local pval = r(p);
outreg2 using "Table_Appen_B1.xls", keep(T1_treat T2_treat T1_sisb T2_sisb s_puntaje_fill) se nocons excel append addstat("Mean", `ymean', "T1_sisb=T2_sisb", `fstat', "T1_sisb=T2p_sisb", `pval') bdec(3) sdec(3)  adec(2);

*Suba: Socio-Econ Controls and FE's;
xi: areg icfes T3_treat T3_sisb s_puntaje_fill $baseline_vars if suba == 1, cluster(school_code) a(school_code); 
su icfes if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_Appen_B1.xls", keep(T3_treat T3_sisb s_puntaje_fill) se nocons excel append addstat("Mean", `ymean') bdec(3) sdec(3)  adec(2);

*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg icfes T1_treat T2_treat T3_treat T1_sisb T2_sisb T3_sisb s_puntaje_fill $baseline_vars suba if grade > 8, cluster(school_code) a(school_code); 
su icfes if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_sisb == T2_sisb;
local fstat1 = r(F);
local pval1 = r(p);
test T1_sisb == T3_sisb;
local fstat2 = r(F);
local pval2 = r(p);
test T2_sisb == T3_sisb;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_B1.xls", keep(T1_treat T2_treat T3_treat T1_sisb T2_sisb T3_sisb s_puntaje_fill) se nocons excel append addstat("Mean", `ymean', 
   "T1_sisb=T2_sisb", `fstat1', "T1_sisb=T2p_sisb", `pval1',  "T1_sisb=T3_sisb", `fstat2', "T1_sisb=T3p_sisb", `pval2',  "T2=_sexo=T3_sisb", `fstat3', "T2_sisb=T3p_sisb", `pval3') bdec(3) sdec(3)  adec(2);


*****;
*OUTCOME;
*C. SPADIES;
*****;

*SC: Socio-Econ Controls and FE's;
xi: areg spadies T1_treat T2_treat T1_sisb T2_sisb s_puntaje_fill $baseline_vars if suba == 0, cluster(school_code) a(school_code); 
su spadies if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_sisb == T2_sisb;
local fstat = r(F);
local pval = r(p);
outreg2 using "Table_Appen_B1.xls", keep(T1_treat T2_treat T1_sisb T2_sisb s_puntaje_fill) se nocons excel append addstat("Mean", `ymean', "T1_sisb=T2_sisb", `fstat', "T1_sisb=T2p_sisb", `pval') bdec(3) sdec(3)  adec(2);


*Suba: Socio-Econ Controls and FE's;
xi: areg spadies T3_treat T3_sisb s_puntaje_fill $baseline_vars if suba == 1, cluster(school_code) a(school_code); 
su spadies if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_Appen_B1.xls", keep(T3_treat T3_sisb s_puntaje_fill) se nocons excel append addstat("Mean", `ymean') bdec(3) sdec(3)  adec(2);

*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg spadies T1_treat T2_treat T3_treat T1_sisb T2_sisb T3_sisb s_puntaje_fill $baseline_vars suba if grade > 8, cluster(school_code) a(school_code); 
su spadies if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_sisb == T2_sisb;
local fstat1 = r(F);
local pval1 = r(p);
test T1_sisb == T3_sisb;
local fstat2 = r(F);
local pval2 = r(p);
test T2_sisb == T3_sisb;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_B1.xls", keep(T1_treat T2_treat T3_treat T1_sisb T2_sisb T3_sisb s_puntaje_fill) se nocons excel append addstat("Mean", `ymean', 
   "T1_sisb=T2_sisb", `fstat1', "T1_sisb=T2p_sisb", `pval1',  "T1_sisb=T3_sisb", `fstat2', "T1_sisb=T3p_sisb", `pval2',  "T2=_sexo=T3_sisb", `fstat3', "T2_sisb=T3p_sisb", `pval3') bdec(3) sdec(3)  adec(2);

*****;
*OUTCOME;
*D. SPADIES long term;
*****;

*SC: Socio-Econ Controls and FE's;
xi: areg spadies_longterm T1_treat T2_treat T1_sisb T2_sisb s_puntaje_fill $baseline_vars if suba == 0, cluster(school_code) a(school_code); 
su spadies_longterm if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_sisb == T2_sisb;
local fstat = r(F);
local pval = r(p);
outreg2 using "Table_Appen_B1.xls", keep(T1_treat T2_treat T1_sisb T2_sisb s_puntaje_fill) se nocons excel append addstat("Mean", `ymean', "T1_sisb=T2_sisb", `fstat', "T1_sisb=T2p_sisb", `pval') bdec(3) sdec(3)  adec(2);


*Suba: Socio-Econ Controls and FE's;
xi: areg spadies_longterm T3_treat T3_sisb s_puntaje_fill $baseline_vars if suba == 1, cluster(school_code) a(school_code); 
su spadies_longterm if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_Appen_B1.xls", keep(T3_treat T3_sisb s_puntaje_fill) se nocons excel append addstat("Mean", `ymean') bdec(3) sdec(3)  adec(2);

*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg spadies_longterm T1_treat T2_treat T3_treat T1_sisb T2_sisb T3_sisb s_puntaje_fill $baseline_vars suba if grade > 8, cluster(school_code) a(school_code); 
su spadies_longterm if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_sisb == T2_sisb;
local fstat1 = r(F);
local pval1 = r(p);
test T1_sisb == T3_sisb;
local fstat2 = r(F);
local pval2 = r(p);
test T2_sisb == T3_sisb;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_B1.xls", keep(T1_treat T2_treat T3_treat T1_sisb T2_sisb T3_sisb s_puntaje_fill) se nocons excel append addstat("Mean", `ymean', 
   "T1_sisb=T2_sisb", `fstat1', "T1_sisb=T2p_sisb", `pval1',  "T1_sisb=T3_sisb", `fstat2', "T1_sisb=T3p_sisb", `pval2',  "T2=_sexo=T3_sisb", `fstat3', "T2_sisb=T3p_sisb", `pval3') bdec(3) sdec(3)  adec(2);   


 *** Heterogenous effects by gender;
   
*Generation interaction treatment with gender (cleaner than xi reg);

gen T1_sexo=T1_treat*s_sexo_fill;
gen T2_sexo=T2_treat*s_sexo_fill;
gen T3_sexo=T3_treat*s_sexo_fill;

*****;
*OUTCOME;
*A. on time enrollemtne;
*****;


*SC: Socio-Econ Controls and FE's;
xi: areg on_time T1_treat T2_treat T1_sexo T2_sexo s_sexo_fill $baseline_vars if suba == 0, cluster(school_code) a(school_code); 
su on_time if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_sexo == T2_sexo;
local fstat = r(F);
local pval = r(p);
outreg2 using "Table_Appen_B2.xls", keep(T1_treat T2_treat T3_treat T1_sexo T2_sexo T3_sexo s_sexo_fill) se replace excel nocons addstat("Mean", `ymean', "T1=T2", `fstat', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);


*Suba: Socio-Econ Controls and FE's;
xi: areg on_time T3_treat T3_sexo s_sexo_fill $baseline_vars if suba == 1, cluster(school_code) a(school_code); 
su on_time if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_Appen_B2.xls", keep(T3_treat T3_sexo s_sexo_fill) se excel append nocons addstat("Mean", `ymean') bdec(3) sdec(3)  adec(2);


*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg on_time T1_treat T2_treat T3_treat T1_sexo T2_sexo T3_sexo s_sexo_fill $baseline_vars suba if grade > 8, cluster(school_code) a(school_code); 
su on_time if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_sexo == T2_sexo;
local fstat1 = r(F);
local pval1 = r(p);
test T1_sexo == T3_sexo;
local fstat2 = r(F);
local pval2 = r(p);
test T2_sexo == T3_sexo;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_B2.xls", keep(T1_treat T2_treat T3_treat T1_sexo T2_sexo T3_sexo s_sexo_fill) se excel nocons append addstat("Mean", `ymean', 
   "T1_sexo=T2_sexo", `fstat1', "T1_sexo=T2p_sexo", `pval1',  "T1_sexo=T3_sexo", `fstat2', "T1_sexo=T3p_sexo", `pval2',  "T2=_sexo=T3_sexo", `fstat3', "T2_sexo=T3p_sexo", `pval3');

   
*****;
*OUTCOME;
*B. Icfes;
*****;

*SC: Socio-Econ Controls and FE's;
xi: areg icfes T1_treat T2_treat T1_sexo T2_sexo s_sexo_fill $baseline_vars if suba == 0, cluster(school_code) a(school_code); 
su icfes if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_sexo == T2_sexo;
local fstat = r(F);
local pval = r(p);
outreg2 using "Table_Appen_B2.xls", keep(T1_treat T2_treat T3_treat T1_sexo T2_sexo T3_sexo s_sexo_fill) se append excel nocons addstat("Mean", `ymean', "T1=T2", `fstat', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);


*Suba: Socio-Econ Controls and FE's;
xi: areg icfes T3_treat T3_sexo s_sexo_fill $baseline_vars if suba == 1, cluster(school_code) a(school_code); 
su icfes if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_Appen_B2.xls", keep(T3_treat T3_sexo s_sexo_fill) se excel append nocons addstat("Mean", `ymean') bdec(3) sdec(3)  adec(2);

*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg icfes T1_treat T2_treat T3_treat T1_sexo T2_sexo T3_sexo s_sexo_fill $baseline_vars suba if grade > 8, cluster(school_code) a(school_code); 
su icfes if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_sexo == T2_sexo;
local fstat1 = r(F);
local pval1 = r(p);
test T1_sexo == T3_sexo;
local fstat2 = r(F);
local pval2 = r(p);
test T2_sexo == T3_sexo;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_B2.xls", keep(T1_treat T2_treat T3_treat T1_sexo T2_sexo T3_sexo s_sexo_fill) se excel nocons append addstat("Mean", `ymean', 
   "T1_sexo=T2_sexo", `fstat1', "T1_sexo=T2p_sexo", `pval1',  "T1_sexo=T3_sexo", `fstat2', "T1_sexo=T3p_sexo", `pval2',  "T2=_sexo=T3_sexo", `fstat3', "T2_sexo=T3p_sexo", `pval3');


*****;
*OUTCOME;
*C. SPADIES;
*****;

*SC: Socio-Econ Controls and FE's;
xi: areg spadies T1_treat T2_treat T1_sexo T2_sexo s_sexo_fill $baseline_vars if suba == 0, cluster(school_code) a(school_code); 
su spadies if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_sexo == T2_sexo;
local fstat = r(F);
local pval = r(p);
outreg2 using "Table_Appen_B2.xls", keep(T1_treat T2_treat T3_treat T1_sexo T2_sexo T3_sexo s_sexo_fill) se append excel nocons addstat("Mean", `ymean', "T1=T2", `fstat', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);


*Suba: Socio-Econ Controls and FE's;
xi: areg spadies T3_treat T3_sexo s_sexo_fill $baseline_vars if suba == 1, cluster(school_code) a(school_code); 
su spadies if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_Appen_B2.xls", keep(T3_treat T3_sexo s_sexo_fill) se excel append nocons addstat("Mean", `ymean') bdec(3) sdec(3)  adec(2);

*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg spadies T1_treat T2_treat T3_treat T1_sexo T2_sexo T3_sexo s_sexo_fill $baseline_vars suba if grade > 8, cluster(school_code) a(school_code); 
su spadies if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_sexo == T2_sexo;
local fstat1 = r(F);
local pval1 = r(p);
test T1_sexo == T3_sexo;
local fstat2 = r(F);
local pval2 = r(p);
test T2_sexo == T3_sexo;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_B2.xls", keep(T1_treat T2_treat T3_treat T1_sexo T2_sexo T3_sexo s_sexo_fill) se excel nocons append addstat("Mean", `ymean', 
   "T1_sexo=T2_sexo", `fstat1', "T1_sexo=T2p_sexo", `pval1',  "T1_sexo=T3_sexo", `fstat2', "T1_sexo=T3p_sexo", `pval2',  "T2=_sexo=T3_sexo", `fstat3', "T2_sexo=T3p_sexo", `pval3') bdec(3) sdec(3)  adec(2);

*****;
*OUTCOME;
*D. SPADIES long term;
*****;

*SC: Socio-Econ Controls and FE's;
xi: areg spadies_longterm T1_treat T2_treat T1_sexo T2_sexo s_sexo_fill $baseline_vars if suba == 0, cluster(school_code) a(school_code); 
su spadies_longterm if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_sexo == T2_sexo;
local fstat = r(F);
local pval = r(p);
outreg2 using "Table_Appen_B2.xls", keep(T1_treat T2_treat T3_treat T1_sexo T2_sexo T3_sexo s_sexo_fill) se append excel nocons addstat("Mean", `ymean', "T1=T2", `fstat', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);


*Suba: Socio-Econ Controls and FE's;
xi: areg spadies_longterm T3_treat T3_sexo s_sexo_fill $baseline_vars if suba == 1, cluster(school_code) a(school_code); 
su spadies_longterm if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_Appen_B2.xls", keep(T3_treat T3_sexo s_sexo_fill) se excel append nocons addstat("Mean", `ymean') bdec(3) sdec(3)  adec(2);

*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg spadies_longterm T1_treat T2_treat T3_treat T1_sexo T2_sexo T3_sexo s_sexo_fill $baseline_vars suba if grade > 8, cluster(school_code) a(school_code); 
su spadies_longterm if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_sexo == T2_sexo;
local fstat1 = r(F);
local pval1 = r(p);
test T1_sexo == T3_sexo;
local fstat2 = r(F);
local pval2 = r(p);
test T2_sexo == T3_sexo;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_B2.xls", keep(T1_treat T2_treat T3_treat T1_sexo T2_sexo T3_sexo s_sexo_fill) se excel nocons append addstat("Mean", `ymean', 
   "T1_sexo=T2_sexo", `fstat1', "T1_sexo=T2p_sexo", `pval1',  "T1_sexo=T3_sexo", `fstat2', "T1_sexo=T3p_sexo", `pval2',  "T2=_sexo=T3_sexo", `fstat3', "T2_sexo=T3p_sexo", `pval3') bdec(3) sdec(3)  adec(2);
   


/*****
*** Appen C1: Randomization Inference;
*****/

*Run at the end of the file. Higly demanding in computational time (2000 iteractions);  
  
/*****
*** Appen C2: Medium-term tertiary enrollment conditional on reaching grade 9 ;
*****/
  
*Var for grades 6-8 at baseline who arrive to grade 9: reach_grade_9;
*Var that capture all people who arrive to grade 9. ;

gen all_reach_grade_9=1 if reach_grade_9==1 | grade>8; 

*Regression of spadies ("short term") for all studants who arrived to grade (all), only for san cristobal;

xi: areg spadies T1_treat T2_treat $baseline_vars if suba == 0, cluster(school_code) absorb(school_code);
su spadies if suba == 0 & treatment == 0;
local ymean = r(mean);
outreg2 using "Table_Appen_C2.xls", keep (T1_treat T2_treat) se nocons replace addstat("Mean control", `ymean')  bdec(3) sdec(3)  adec(2); 

xi: areg spadies T1_treat T2_treat $baseline_vars if suba == 0 & all_reach_grade_9==1, cluster(school_code) absorb(school_code);
su spadies if suba == 0 & treatment == 0 & all_reach_grade_9==1;
local ymean = r(mean);
outreg2 using "Table_Appen_C2.xls", keep (T1_treat T2_treat) se nocons append addstat("Mean control", `ymean')  bdec(3) sdec(3)  adec(2); 


/*****
* Append C3: Tertiary enrollment, medium term, Lee-bound estimation;
*****/

set seed 2023;

*Short term spadies; 

*Regression, Treatment 1;
reg spadies T1_treat if suba == 0 & T2_treat~=1 & grade>8, cluster(school_code); 
outreg2  using "Table_Appen_C3.xls", keep(T1_treat) se nocons replace bdec(3) sdec(3);

*Bound, Treatment 1;
leebounds spadies T1_treat if T2_treat~=1 & suba == 0 & grade>8 , select(icfes) tight(s_sexo_fill s_estrato_fill) vce(boot);
outreg2  using "Table_Appen_C3.xls", se nocons append  bdec(3) sdec(3);

*Regression, Treatment 2;
reg spadies T2_treat if suba == 0 & T1_treat~=1 & grade>8, cluster(school_code); 
outreg2  using "Table_Appen_C3.xls", keep(T2_treat) se nocons append bdec(3) sdec(3);

*Bound, Treatment 2;
leebounds spadies T2_treat if T1_treat~=1 & suba == 0 & grade>8, select(icfes) tight(s_sexo_fill s_estrato_fill) vce(boot);
outreg2  using "Table_Appen_C3.xls", se nocons append  bdec(3) sdec(3);

*Regression, Treatment 3;
reg spadies T3_treat if suba == 1 & grade>8, cluster(school_code); 
outreg2  using "Table_Appen_C3.xls", keep(T3_treat) se nocons append bdec(3) sdec(3);

*Bound, Treatment 3;
leebounds spadies T3_treat if suba == 1 & grade>8 , select(icfes) tight(s_sexo_fill s_estrato_fill) vce(boot);
outreg2  using "Table_Appen_C3.xls", se nocons append  bdec(3) sdec(3);

*Long term spadies;
*Regression, Treatment 1;
reg spadies_longterm T1_treat if suba == 0 & T2_treat~=1 & grade>8, cluster(school_code); 
outreg2  using "Table_Appen_C3.xls", keep(T1_treat) se nocons append bdec(3) sdec(3);

*Bound, Treatment 1;
leebounds spadies_longterm T1_treat if T2_treat~=1 & suba == 0  & grade>8, select(icfes) tight(s_sexo_fill s_estrato_fill) vce(boot);
outreg2  using "Table_Appen_C3.xls", se nocons append  bdec(3) sdec(3);

*Regression, Treatment 2;
reg spadies_longterm T2_treat if suba == 0 & T1_treat~=1 & grade>8, cluster(school_code); 
outreg2  using "Table_Appen_C3.xls", keep(T2_treat) se nocons append bdec(3) sdec(3);

*Bound, Treatment 2;
leebounds spadies_longterm T2_treat if T1_treat~=1 & suba == 0 & grade>8, select(icfes) tight(s_sexo_fill s_estrato_fill) vce(boot);
outreg2  using "Table_Appen_C3.xls",  se nocons append  bdec(3) sdec(3);

*Regression, Treatment 3;
reg spadies_longterm T3_treat if suba == 1 & grade>8, cluster(school_code); 
outreg2  using "Table_Appen_C3.xls", keep(T3_treat) se nocons append bdec(3) sdec(3);

*Bound, Treatment 3;
leebounds spadies_longterm T3_treat if suba == 1 & grade>8, select(icfes) tight(s_sexo_fill s_estrato_fill) vce(boot);
outreg2  using "Table_Appen_C3.xls", se nocons append  bdec(3) sdec(3);


/*****
* Append C4: Long term tertiary: on-time and graduaction conditional on ICFES and on enrollment;
*** For format purposes, we are dividing the table in two panels, columns 1 to 6, panel A; rest of columns panel B;
*****/


/****
*** Ontime Estimation CONDITIONAL on ICFES;
****/

*SC: Socio-Econ Controls and FE's;
xi: reg on_time_uncond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0 & icfes==1, cluster(school_code); 
su on_time_uncond_longterm if treatment == 0 & suba == 0 & icfes==1;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat = r(F);
local pval = r(p);
outreg2 using "Table_Appen_C4.xls", replace keep(T1_treat T2_treat) se nocons addstat("Mean", `ymean', "T1=T2", `fstat', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

*Suba: Socio-Econ Controls and FE's;
xi: reg on_time_uncond_longterm T3_treat $baseline_vars i.school_code if suba == 1 & icfes==1, cluster(school_code); 
su on_time_uncond_longterm if treatment == 0 & suba == 1 & icfes==1;
local ymean = r(mean);
outreg2 using "Table_Appen_C4.xls", append keep(T3_treat) se nocons addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);


*Grades 9-11: Socio-Econ Controls and FE's;
xi: reg on_time_uncond_longterm T1_treat T2_treat T3_treat $baseline_vars suba i.school_code if grade > 8 & icfes==1, cluster(school_code); 
su on_time_uncond_longterm if treatment == 0 & grade > 8 & icfes==1;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat1 = r(F);
local pval1 = r(p);
test T1_treat == T3_treat;
local fstat2 = r(F);
local pval2 = r(p);
test T2_treat == T3_treat;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_C4.xls", append keep(T1_treat T2_treat T3_treat) se nocons addstat("Mean", `ymean', 
   "T1=T2", `fstat1', "T1=T2p", `pval1',  "T1=T3", `fstat2', "T1=T3p", `pval2',  "T2=T3", `fstat3', "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);

/****
*** Graduation Estimation CONDITIONAL on ICFES;
****/

*SC: Socio-Econ Controls and FE's;
xi: reg grad_uncond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0 & icfes==1, cluster(school_code); 
su grad_uncond_longterm if treatment == 0 & suba == 0 & icfes==1;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat = r(F);
local pval = r(p);
outreg2 using "Table_Appen_C4.xls", append keep(T1_treat T2_treat) se nocons addstat("Mean", `ymean', "T1=T2", `fstat', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);


*Suba: Socio-Econ Controls and FE's;
xi: reg grad_uncond_longterm T3_treat $baseline_vars i.school_code if suba == 1 & icfes==1, cluster(school_code); 
su grad_uncond_longterm if treatment == 0 & suba == 1 & icfes==1;
local ymean = r(mean);
outreg2 using "Table_Appen_C4.xls", append keep(T3_treat) se nocons addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);


*Grades 9-11: Socio-Econ Controls and FE's;
xi: reg grad_uncond_longterm T1_treat T2_treat T3_treat $baseline_vars suba i.school_code if grade > 8 & icfes==1, cluster(school_code); 
su grad_uncond_longterm if treatment == 0 & grade > 8 & icfes==1;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat1 = r(F);
local pval1 = r(p);
test T1_treat == T3_treat;
local fstat2 = r(F);
local pval2 = r(p);
test T2_treat == T3_treat;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_C4.xls", append keep(T1_treat T2_treat T3_treat) se nocons addstat("Mean", `ymean', 
   "T1=T2", `fstat1', "T1=T2p", `pval1',  "T1=T3", `fstat2', "T1=T3p", `pval2',  "T2=T3", `fstat3', "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);

/****
*** On time enrollment Variable created CONDITIONAL on teriary enrollment;
****/

*SC: Socio-Econ Controls and FE's;
xi: areg on_time_longterm T1_treat T2_treat $baseline_vars if suba == 0, cluster(school_code) absorb(school_code); 
su on_time_longterm if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat = r(F);
local pval = r(p);
outreg2 using "Table_Appen_C4.xls", append keep(T1_treat T2_treat) se nocons addstat("Mean", `ymean', "T1=T2", `fstat', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

*Suba: Socio-Econ Controls and FE's;
xi: areg on_time_longterm T3_treat $baseline_vars if suba == 1, cluster(school_code) absorb(school_code); 
su on_time_longterm if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_Appen_C4.xls", append keep(T3_treat) se nocons addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);


*Grades 9-11: Socio-Econ Controls and FE's;
xi: areg on_time_longterm T1_treat T2_treat T3_treat $baseline_vars suba if grade > 8, cluster(school_code) absorb(school_code); 
su on_time_longterm if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat1 = r(F);
local pval1 = r(p);
test T1_treat == T3_treat;
local fstat2 = r(F);
local pval2 = r(p);
test T2_treat == T3_treat;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_C4.xls", append keep(T1_treat T2_treat T3_treat) se nocons addstat("Mean", `ymean', 
   "T1=T2", `fstat1', "T1=T2p", `pval1',  "T1=T3", `fstat2', "T1=T3p", `pval2',  "T2=T3", `fstat3', "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);

   
/****
*** Graduation Variable created CONDITIONAL on enrollment;
****/

*SC: Socio-Econ Controls and FE's;
xi: reg grad_cond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code); 
su grad_cond_longterm if treatment == 0 & suba == 0;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat = r(F);
local pval = r(p);
outreg2 using "Table_Appen_C4.xls", append keep(T1_treat T2_treat) se nocons addstat("Mean", `ymean', "T1=T2", `fstat', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

*Suba: Socio-Econ Controls and FE's;
xi: reg grad_cond_longterm T3_treat $baseline_vars i.school_code if suba == 1, cluster(school_code); 
su grad_cond_longterm if treatment == 0 & suba == 1;
local ymean = r(mean);
outreg2 using "Table_Appen_C4.xls", append keep(T3_treat) se nocons addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);


*Grades 9-11: Socio-Econ Controls and FE's;
xi: reg grad_cond_longterm T1_treat T2_treat T3_treat $baseline_vars suba i.school_code if grade > 8, cluster(school_code); 
su grad_cond_longterm if treatment == 0 & grade > 8;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat1 = r(F);
local pval1 = r(p);
test T1_treat == T3_treat;
local fstat2 = r(F);
local pval2 = r(p);
test T2_treat == T3_treat;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_C4.xls", append keep(T1_treat T2_treat T3_treat) se nocons addstat("Mean", `ymean', 
   "T1=T2", `fstat1', "T1=T2p", `pval1',  "T1=T3", `fstat2', "T1=T3p", `pval2',  "T2=T3", `fstat3', "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);

   
/*****
*** Append C5: exclusion of sisben 1. Four outcomes: on_time, icfes, spadies and spadies_longterm;
*** For format purposes, we are dividing the table in two panels, columns 1 to 6, panel A; rest of columns panel B;
*****/

*On-time: SC: Socio-Econ Controls and FE's;
xi: areg on_time T1_treat T2_treat $baseline_vars if suba == 0 & r_be_sisb_punt>=11, cluster(school_code) absorb(school_code);
su on_time if treatment == 0 & suba == 0 & r_be_sisb_punt>=11;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat = r(F);
local pval = r(p);
outreg2  using "Table_Appen_C5.xls", keep(T1_treat T2_treat) se nocons replace addstat("Mean", `ymean', "T1=T2", `fstat', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

*On-time: Suba: Socio-Econ Controls and FE's and conditional on icfes;
xi: areg on_time T3_treat $baseline_vars if suba == 1 & r_be_sisb_punt>=11, cluster(school_code) absorb(school_code); 
su on_time if treatment == 0 & suba == 1 & r_be_sisb_punt>=11;
local ymean = r(mean);
outreg2 using "Table_Appen_C5.xls", keep(T3_treat) se nocons append addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);

*On-time, Grades 9-11: Socio-Econ Controls and FE's;
xi: areg on_time T1_treat T2_treat T3_treat $baseline_vars suba if grade > 8 & r_be_sisb_punt>=11, cluster(school_code) absorb(school_code); 
su on_time if treatment == 0 & grade > 8 & r_be_sisb_punt>=11;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat1 = r(F);
local pval1 = r(p);
test T1_treat == T3_treat;
local fstat2 = r(F);
local pval2 = r(p);
test T2_treat == T3_treat;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_C5.xls", keep(T1_treat T2_treat T3_treat) se append nocons addstat("Mean", `ymean', 
   "T1=T2", `fstat1', "T1=T2p", `pval1',  "T1=T3", `fstat2', "T1=T3p", `pval2',  "T2=T3", `fstat3', "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);


*Icfes: SC: Socio-Econ Controls and FE's;
xi: areg icfes T1_treat T2_treat $baseline_vars if suba == 0 & r_be_sisb_punt>=11, cluster(school_code) absorb(school_code);
su icfes if treatment == 0 & suba == 0 & r_be_sisb_punt>=11;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat = r(F);
local pval = r(p);
outreg2  using "Table_Appen_C5.xls", keep(T1_treat T2_treat) se nocons append addstat("Mean", `ymean', "T1=T2", `fstat', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

*Icfes: Suba: Socio-Econ Controls and FE's and conditional on icfes;
xi: areg icfes T3_treat $baseline_vars if suba == 1 & r_be_sisb_punt>=11, cluster(school_code) absorb(school_code); 
su icfes if treatment == 0 & suba == 1 & r_be_sisb_punt>=11;
local ymean = r(mean);
outreg2 using "Table_Appen_C5.xls", keep(T3_treat) se nocons append addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);

*Icfes, Grades 9-11: Socio-Econ Controls and FE's;
xi: areg icfes T1_treat T2_treat T3_treat $baseline_vars suba if grade > 8 & r_be_sisb_punt>=11, cluster(school_code) absorb(school_code); 
su icfes if treatment == 0 & grade > 8 & r_be_sisb_punt>=11;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat1 = r(F);
local pval1 = r(p);
test T1_treat == T3_treat;
local fstat2 = r(F);
local pval2 = r(p);
test T2_treat == T3_treat;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_C5.xls", keep(T1_treat T2_treat T3_treat) se append nocons addstat("Mean", `ymean', 
   "T1=T2", `fstat1', "T1=T2p", `pval1',  "T1=T3", `fstat2', "T1=T3p", `pval2',  "T2=T3", `fstat3', "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);

*Spadies: SC: Socio-Econ Controls and FE's;
xi: areg spadies T1_treat T2_treat $baseline_vars if suba == 0 & r_be_sisb_punt>=11, cluster(school_code) absorb(school_code);
su spadies if treatment == 0 & suba == 0 & r_be_sisb_punt>=11;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat = r(F);
local pval = r(p);
outreg2  using "Table_Appen_C5.xls", keep(T1_treat T2_treat) se nocons append addstat("Mean", `ymean', "T1=T2", `fstat', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

*Spadies: Suba: Socio-Econ Controls and FE's and conditional on icfes;
xi: areg spadies T3_treat $baseline_vars if suba == 1 & r_be_sisb_punt>=11, cluster(school_code) absorb(school_code); 
su spadies if treatment == 0 & suba == 1 & r_be_sisb_punt>=11;
local ymean = r(mean);
outreg2 using "Table_Appen_C5.xls", keep(T3_treat) se nocons append addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);

*Spadies, Grades 9-11: Socio-Econ Controls and FE's;
xi: areg spadies T1_treat T2_treat T3_treat $baseline_vars suba if grade > 8 & r_be_sisb_punt>=11, cluster(school_code) absorb(school_code); 
su spadies if treatment == 0 & grade > 8 & r_be_sisb_punt>=11;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat1 = r(F);
local pval1 = r(p);
test T1_treat == T3_treat;
local fstat2 = r(F);
local pval2 = r(p);
test T2_treat == T3_treat;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_C5.xls", keep(T1_treat T2_treat T3_treat) se append nocons addstat("Mean", `ymean', 
   "T1=T2", `fstat1', "T1=T2p", `pval1',  "T1=T3", `fstat2', "T1=T3p", `pval2',  "T2=T3", `fstat3', "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);

   
*Spadies LT: SC: Socio-Econ Controls and FE's;
xi: areg spadies_longterm T1_treat T2_treat $baseline_vars if suba == 0 & r_be_sisb_punt>=11, cluster(school_code) absorb(school_code);
su spadies_longterm if treatment == 0 & suba == 0 & r_be_sisb_punt>=11;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat = r(F);
local pval = r(p);
outreg2  using "Table_Appen_C5.xls", keep(T1_treat T2_treat) se nocons append addstat("Mean", `ymean', "T1=T2", `fstat', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

*Spadies LT: Suba: Socio-Econ Controls and FE's and conditional on icfes;
xi: areg spadies_longterm T3_treat $baseline_vars if suba == 1 & r_be_sisb_punt>=11, cluster(school_code) absorb(school_code); 
su spadies_longterm if treatment == 0 & suba == 1 & r_be_sisb_punt>=11;
local ymean = r(mean);
outreg2 using "Table_Appen_C5.xls", keep(T3_treat) se nocons append addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);

*Spadies LT, Grades 9-11: Socio-Econ Controls and FE's;
xi: areg spadies_longterm T1_treat T2_treat T3_treat $baseline_vars suba if grade > 8 & r_be_sisb_punt>=11, cluster(school_code) absorb(school_code); 
su spadies_longterm if treatment == 0 & grade > 8 & r_be_sisb_punt>=11;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat1 = r(F);
local pval1 = r(p);
test T1_treat == T3_treat;
local fstat2 = r(F);
local pval2 = r(p);
test T2_treat == T3_treat;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_C5.xls", keep(T1_treat T2_treat T3_treat) se append nocons addstat("Mean", `ymean', 
   "T1=T2", `fstat1', "T1=T2p", `pval1',  "T1=T3", `fstat2', "T1=T3p", `pval2',  "T2=T3", `fstat3', "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);

   
   
/*****
* Append C6: Result on ICFES and on SPADIES, restricting sample with enrollment information;
*****/
   
xi: areg icfes T1_treat T2_treat $baseline_vars  if suba == 0 & m_enrolled ~= ., cluster(school_code) absorb(school_code);  
su icfes if treatment == 0 & suba == 0 & m_enrolled ~= .;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat = r(F);
local pval = r(p);
outreg2 using "Table_Appen_C6.xls", replace keep(T1_treat T2_treat) se nocons addstat("Mean", `ymean', "T1=T2", `fstat', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

xi: areg icfes T3_treat $baseline_vars if suba == 1 & m_enrolled ~= ., cluster(school_code) absorb(school_code); 
su icfes if treatment == 0 & suba == 1 & m_enrolled ~= .;
local ymean = r(mean);
outreg2 using "Table_Appen_C6.xls", append keep(T3_treat) se nocons addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);

xi: areg icfes T1_treat T2_treat T3_treat $baseline_vars suba if grade > 8 & m_enrolled ~= ., cluster(school_code) absorb(school_code); 
su icfes if treatment == 0 & grade > 8 & m_enrolled ~= .;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat1 = r(F);
local pval1 = r(p);
test T1_treat == T3_treat;
local fstat2 = r(F);
local pval2 = r(p);
test T2_treat == T3_treat;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_C6.xls", append keep(T1_treat T2_treat T3_treat) se nocons addstat("Mean", `ymean', 
   "T1=T2", `fstat1', "T1=T2p", `pval1',  "T1=T3", `fstat2', "T1=T3p", `pval2',  "T2=T3", `fstat3', "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);

xi: areg spadies T1_treat T2_treat $baseline_vars  if suba == 0 & m_enrolled ~= ., cluster(school_code) absorb(school_code);  
su spadies if treatment == 0 & suba == 0 & m_enrolled ~= .;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat = r(F);
local pval = r(p);
outreg2 using "Table_Appen_C6.xls", append keep(T1_treat T2_treat) se nocons addstat("Mean", `ymean', "T1=T2", `fstat', "T1=T2p", `pval')  bdec(3) sdec(3)  adec(2);

xi: areg spadies T3_treat $baseline_vars if suba == 1 & m_enrolled ~= ., cluster(school_code) absorb(school_code); 
su spadies if treatment == 0 & suba == 1 & m_enrolled ~= .;
local ymean = r(mean);
outreg2 using "Table_Appen_C6.xls", append keep(T3_treat) se nocons addstat("Mean", `ymean')  bdec(3) sdec(3)  adec(2);

xi: areg spadies T1_treat T2_treat T3_treat $baseline_vars suba if grade > 8 & m_enrolled ~= ., cluster(school_code) absorb(school_code); 
su spadies if treatment == 0 & grade > 8 & m_enrolled ~= .;
local ymean = r(mean);
test T1_treat == T2_treat;
local fstat1 = r(F);
local pval1 = r(p);
test T1_treat == T3_treat;
local fstat2 = r(F);
local pval2 = r(p);
test T2_treat == T3_treat;
local fstat3 = r(F);
local pval3 = r(p);
outreg2 using "Table_Appen_C6.xls", append keep(T1_treat T2_treat T3_treat) se nocons addstat("Mean", `ymean', 
   "T1=T2", `fstat1', "T1=T2p", `pval1',  "T1=T3", `fstat2', "T1=T3p", `pval2',  "T2=T3", `fstat3', "T2=T3p", `pval3')  bdec(3) sdec(3)  adec(2);

/*
   
/*****
*** Appen C1: Randomization Inference;
*** Aprox. time of running: 5 days;
*** Set seed at an arbitrary number for replication;
*** Generataion of another dataset with the results of the estimation; 
*****/

keep T1_treat T2_treat T3_treat on_time icfes spadies spadies_longterm on_time_uncond_longterm grad_uncond_longterm school_code 
suba T1_treat_upper control control_upper T1_treat_lower control_lower treatment

   s_teneviv_fill s_utilities_fill s_durables_fill s_infraest_hh_fill 
   s_age_sorteo_fill s_sexo_fill s_yrs_fill grade_fill s_single_fill s_edadhead_fill
   s_yrshead_fill s_tpersona_fill s_num18_fill s_estrato_fill s_puntaje_fill s_ingtotal_fill
   s_teneviv_miss s_utilities_miss s_durables_miss s_infraest_hh_miss
   s_age_sorteo_miss s_sexo_miss s_yrs_miss grade_miss s_single_miss s_edadhead_miss
   s_yrshead_miss s_tpersona_miss s_num18_miss s_estrato_miss s_puntaje_miss s_ingtotal_miss
;

global baseline_vars "
   i.s_teneviv_fill i.s_utilities_fill i.s_durables_fill s_infraest_hh_fill 
   s_age_sorteo_fill s_sexo_fill s_yrs_fill i.grade_fill s_single_fill s_edadhead_fill
   s_yrshead_fill s_tpersona_fill s_num18_fill i.s_estrato_fill s_puntaje_fill s_ingtotal_fill

   s_teneviv_miss s_utilities_miss s_durables_miss s_infraest_hh_miss
   s_age_sorteo_miss s_sexo_miss s_yrs_miss grade_miss s_single_miss s_edadhead_miss
   s_yrshead_miss s_tpersona_miss s_num18_miss s_estrato_miss s_puntaje_miss s_ingtotal_miss
   ";

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

* Medium term; 

set seed 2049; 
randcmd ((T1_treat) xi: reg on_time T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T1_treat) xi: reg icfes T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T1_treat) xi: reg spadies T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code)), treatvars(T1_treat) reps(2000); 

matrix T1_joint = e(ROmni);
svmat double T1_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;
		   
set seed 2049;
randcmd ((T2_treat) xi: reg on_time T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T2_treat) xi: reg icfes T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T2_treat) xi: reg spadies T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code)), treatvars(T2_treat) reps(2000); 

matrix T2_joint = e(ROmni);
svmat double T2_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

set seed 2049;
randcmd ((T3_treat) xi: reg on_time T3_treat $baseline_vars i.school_code if suba == 1, cluster(school_code)) 
        ((T3_treat) xi: reg icfes T3_treat $baseline_vars i.school_code if suba == 1, cluster(school_code)) 
		((T3_treat) xi: reg spadies T3_treat $baseline_vars i.school_code if suba == 1, cluster(school_code)), treatvars(T3_treat) reps(2000);
		
matrix T3_joint = e(ROmni);
svmat double T3_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

* Long term only; 

set seed 2049;
randcmd ((T1_treat) xi: reg spadies_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T1_treat) xi: reg on_time_uncond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T1_treat) xi: reg grad_uncond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code)), 
		      treatvars(T1_treat) reps(2000); 

matrix T1_lt_joint = e(ROmni);
svmat double T1_lt_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

set seed 2049;
randcmd ((T2_treat) xi: reg spadies_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T2_treat) xi: reg on_time_uncond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T2_treat) xi: reg grad_uncond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code)), 
		      treatvars(T2_treat) reps(2000); 

matrix T2_lt_joint = e(ROmni);
svmat double T2_lt_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

set seed 2049;
randcmd ((T3_treat) xi: reg spadies_longterm T3_treat $baseline_vars i.school_code if suba == 1, cluster(school_code) ) 
           ((T3_treat) xi: reg on_time_uncond_longterm T3_treat $baseline_vars i.school_code if suba == 1, cluster(school_code)) 
		   ((T3_treat) xi: reg grad_uncond_longterm T3_treat $baseline_vars i.school_code if suba == 1, cluster(school_code)), 
		      treatvars(T3_treat) reps(2000); 

matrix T3_lt_joint = e(ROmni);
svmat double T3_lt_joint; 

* Long term and medium term joinly;

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

set seed 2049;
randcmd ((T1_treat) xi: reg on_time T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T1_treat) xi: reg icfes T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T1_treat) xi: reg spadies_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code))
		   ((T1_treat) xi: reg on_time_uncond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code))
		   ((T1_treat) xi: reg grad_uncond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code)), 
		      treatvars(T1_treat) reps(2000); 

matrix T1_both_joint = e(ROmni);
svmat double T1_both_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

set seed 2049;
randcmd ((T2_treat) xi: reg on_time T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T2_treat) xi: reg icfes T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T2_treat) xi: reg spadies_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code))
		   ((T2_treat) xi: reg on_time_uncond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code))
		   ((T2_treat) xi: reg grad_uncond_longterm T1_treat T2_treat $baseline_vars i.school_code if suba == 0, cluster(school_code)), 
		       treatvars(T2_treat) reps(2000); 

matrix T2_both_joint = e(ROmni);
svmat double T2_both_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

set seed 2049;
randcmd ((T3_treat) xi: reg on_time T3_treat $baseline_vars i.school_code if suba == 1, cluster(school_code) ) 
           ((T3_treat) xi: reg icfes T3_treat $baseline_vars i.school_code if suba == 1, cluster(school_code)) 
		   ((T3_treat) xi: reg spadies_longterm T3_treat $baseline_vars i.school_code if suba == 1, cluster(school_code))
		   ((T3_treat) xi: reg on_time_uncond_longterm T3_treat $baseline_vars i.school_code if suba == 1, cluster(school_code))
		   ((T3_treat) xi: reg grad_uncond_longterm T3_treat $baseline_vars i.school_code if suba == 1, cluster(school_code)), 
				treatvars(T3_treat) reps(2000); 

matrix T3_both_joint = e(ROmni);
svmat double T3_both_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;


/*****
*** Test of H0: T2_treat==T1_treat, all grades;
*****/

*Medium term;


set seed 2049;
randcmd ((T1_treat T1_treat_upper) xi: reg on_time T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T1_treat T1_treat_upper) xi: reg icfes T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T1_treat T1_treat_upper) xi: reg spadies T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)), 
		     treatvars(T1_treat T1_treat_upper) reps(2000); 

matrix T1_T2_joint = e(ROmni);
svmat double T1_T2_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

*Long term;

set seed 2049;
randcmd ((T1_treat T1_treat_upper) xi: reg spadies_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T1_treat T1_treat_upper) xi: reg on_time_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T1_treat T1_treat_upper) xi: reg grad_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)), 
		     treatvars(T1_treat T1_treat_upper) reps(2000); 

matrix T1_T2_lt_joint = e(ROmni);
svmat double T1_T2_lt_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

*Both term;

set seed 2049;
randcmd ((T1_treat T1_treat_upper) xi: reg on_time T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code))
		((T1_treat T1_treat_upper) xi: reg icfes T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		((T1_treat T1_treat_upper) xi: reg spadies_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
        ((T1_treat T1_treat_upper) xi: reg on_time_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		((T1_treat T1_treat_upper) xi: reg grad_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)), 
		 treatvars(T1_treat T1_treat_upper) reps(2000); 

matrix T1_T2_both_joint = e(ROmni);
svmat double T1_T2_both_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

/*****
*** Test of H0: T2_treat==T1_treat, upper;
*****/


*Basic versus Savings, medium term -- upper secondary;

set seed 2049;
randcmd ((T1_treat) xi: reg on_time T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T1_treat) xi: reg icfes T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T1_treat) xi: reg spadies T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0, cluster(school_code)), 
		     treatvars(T1_treat) reps(2000); 

matrix T12_upper_joint = e(ROmni);
svmat double T12_upper_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

*Basic versus Savings, only long term -- upper secondary;

set seed 2049;
randcmd ((T1_treat) xi: reg spadies_longterm T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T1_treat) xi: reg on_time_uncond_longterm T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T1_treat) xi: reg grad_uncond_longterm T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0, cluster(school_code)), 
		     treatvars(T1_treat) reps(2000); 

matrix T12_upper_lt_joint = e(ROmni);
svmat double T12_upper_lt_joint; 


save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

*Basic versus Savings, both medium and long term -- upper secondary;

set seed 2049;
randcmd  ((T1_treat) xi: reg on_time T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T1_treat) xi: reg icfes T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T1_treat) xi: reg spadies_longterm T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T1_treat) xi: reg on_time_uncond_longterm T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T1_treat) xi: reg grad_uncond_longterm T1_treat T1_treat_lower control control_lower $baseline_vars i.school_code if suba == 0, cluster(school_code)), 
		     treatvars(T1_treat) reps(2000); 

matrix T12_upper_both_joint = e(ROmni);
svmat double T12_upper_both_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

/*****
*** Test of H0: T2_treat==T1_treat, lower;
*****/

*Basic versus Savings, medium term -- lower secondary;

set seed 2049;
randcmd ((T1_treat) xi: reg on_time T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T1_treat) xi: reg icfes T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T1_treat) xi: reg spadies T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)), 
		     treatvars(T1_treat) reps(2000); 

matrix T12_lower_joint = e(ROmni);
svmat double T12_lower_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;


*Basic versus Savings, only long term -- lower secondary;

set seed 2049;
randcmd ((T1_treat) xi: reg spadies_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T1_treat) xi: reg on_time_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T1_treat) xi: reg grad_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)), 
		     treatvars(T1_treat) reps(2000); 

matrix T12_lower_lt_joint = e(ROmni);
svmat double T12_lower_lt_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;


*Basic versus Savings, both medium and long term -- lower secondary;

set seed 2049;
randcmd  ((T1_treat) xi: reg on_time T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T1_treat) xi: reg icfes T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T1_treat) xi: reg spadies_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code) ) 
           ((T1_treat) xi: reg on_time_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)) 
		   ((T1_treat) xi: reg grad_uncond_longterm T1_treat T1_treat_upper control control_upper $baseline_vars i.school_code if suba == 0, cluster(school_code)), 
		     treatvars(T1_treat) reps(2000); 

matrix T12_lower_both_joint = e(ROmni);
svmat double T12_lower_both_joint; 



save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;


/*****
*** Test of H0: Difference T3_treat==T2_treat;
*****/

*Medium term;

set seed 2049;
randcmd ((T2_treat) xi: areg on_time treatment T1_treat T2_treat $baseline_vars suba if grade_fill > 8, cluster(school_code) absorb(school_code)) 
        ((T2_treat) xi: areg icfes treatment T1_treat T2_treat $baseline_vars suba if grade_fill > 8, cluster(school_code) absorb(school_code))
        ((T2_treat)xi: areg spadies treatment T1_treat T2_treat $baseline_vars suba if grade_fill > 8, cluster(school_code) absorb(school_code)), 
		 treatvars(T2_treat) reps(2000); 

matrix T2_T3_joint = e(ROmni);
svmat double T2_T3_joint; 		 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

*Long term; 

set seed 2049;
randcmd ((T2_treat) xi: areg spadies_longterm treatment T1_treat T2_treat $baseline_vars suba if grade_fill > 8, cluster(school_code) absorb(school_code)) 
        ((T2_treat) xi: areg on_time_uncond_longterm treatment T1_treat T2_treat $baseline_vars suba if grade_fill > 8, cluster(school_code) absorb(school_code))
        ((T2_treat)xi: areg grad_uncond_longterm treatment T1_treat T2_treat $baseline_vars suba if grade_fill > 8, cluster(school_code) absorb(school_code)), 
		 treatvars(T2_treat) reps(2000);

matrix T2_T3_lt_joint = e(ROmni);
svmat double T2_T3_lt_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;
 		 
*Both terms;

set seed 2049;
randcmd ((T2_treat) xi: areg on_time treatment T1_treat T2_treat $baseline_vars suba if grade_fill > 8, cluster(school_code) absorb(school_code)) 
        ((T2_treat) xi: areg icfes treatment T1_treat T2_treat $baseline_vars suba if grade_fill > 8, cluster(school_code) absorb(school_code))
        ((T2_treat) xi: areg spadies_longterm treatment T1_treat T2_treat $baseline_vars suba if grade_fill > 8, cluster(school_code) absorb(school_code)) 
        ((T2_treat) xi: areg on_time_uncond_longterm treatment T1_treat T2_treat $baseline_vars suba if grade_fill > 8, cluster(school_code) absorb(school_code))
        ((T2_treat)xi: areg grad_uncond_longterm treatment T1_treat T2_treat $baseline_vars suba if grade_fill > 8, cluster(school_code) absorb(school_code)), 
		 treatvars(T2_treat) reps(2000); 

matrix T2_T3_both_joint = e(ROmni);
svmat double T2_T3_both_joint; 

save "Final_Data_CCT_Long_Barrera-Osorio_et_al_plus_RI.dta", replace;

*Generate Table Variables;

foreach var of varlist T1_joint* T2_joint* T3_joint* T1_lt_joint* T2_lt_joint* T3_lt_joint*  
					   T1_both_joint* T2_both_joint* T3_both_joint* T1_T2_joint* T1_T2_lt_joint* 
					   T1_T2_both_joint* T12_upper_joint* T12_upper_lt_joint* 
					   T12_upper_both_joint* T12_lower_joint* T12_lower_lt_joint* 
					   T12_lower_both_joint* T2_T3_joint* T2_T3_lt_joint* T2_T3_both_joint*  { ;
egen `var'_b=mean(`var') ;
replace `var' = `var'_b ;
drop `var'_b ;
} ; 

   gen row_num = _n ;
   gen str30 col0 = "" ;
   forvalues i = 1(1)6 {  ;
      gen str20 col`i' = ""  ;
      }  ;
 
   *Values of the Omnius p test;

replace col1 = string(T1_joint1) if row_num == 1 ; 
replace col1 = string(T1_joint2) if row_num == 2 ;
replace col1 = string(T1_joint3) if row_num == 3 ;

replace col2 = string(T1_lt_joint1) if row_num == 1 ; 
replace col2 = string(T1_lt_joint2) if row_num == 2 ;
replace col2 = string(T1_lt_joint3) if row_num == 3 ;

replace col3 = string(T1_both_joint1) if row_num == 1 ; 
replace col3 = string(T1_both_joint2) if row_num == 2 ;
replace col3 = string(T1_both_joint3) if row_num == 3 ;
   
replace col1 = string(T2_joint1) if row_num == 5 ; 
replace col1 = string(T2_joint2) if row_num == 6 ;
replace col1 = string(T2_joint3) if row_num == 7 ;

replace col2 = string(T2_lt_joint1) if row_num == 5 ; 
replace col2 = string(T2_lt_joint2) if row_num == 6 ;
replace col2 = string(T2_lt_joint3) if row_num == 7 ;

replace col3 = string(T2_both_joint1) if row_num == 5 ; 
replace col3 = string(T2_both_joint2) if row_num == 6 ;
replace col3 = string(T2_both_joint3) if row_num == 7 ;
   
replace col1 = string(T3_joint1) if row_num == 9 ; 
replace col1 = string(T3_joint2) if row_num == 10 ;
replace col1 = string(T3_joint3) if row_num == 11 ;
   
replace col2 = string(T3_lt_joint1) if row_num == 9 ; 
replace col2 = string(T3_lt_joint2) if row_num == 10 ;
replace col2 = string(T3_lt_joint3) if row_num == 11 ;

replace col3 = string(T3_both_joint1) if row_num == 9 ; 
replace col3 = string(T3_both_joint2) if row_num == 10 ;
replace col3 = string(T3_both_joint3) if row_num == 11 ;

replace col1 = string(T1_T2_joint1) if row_num == 13 ; 
replace col1 = string(T1_T2_joint2) if row_num == 14 ;
replace col1 = string(T1_T2_joint3) if row_num == 15 ;

replace col2 = string(T1_T2_lt_joint1) if row_num == 13 ; 
replace col2 = string(T1_T2_lt_joint2) if row_num == 14 ;
replace col2 = string(T1_T2_lt_joint3) if row_num == 15 ;

replace col3 = string(T1_T2_both_joint1) if row_num == 13 ; 
replace col3 = string(T1_T2_both_joint2) if row_num == 14 ;
replace col3 = string(T1_T2_both_joint3) if row_num == 15 ;

replace col1 = string(T12_upper_joint1) if row_num == 17 ; 
replace col1 = string(T12_upper_joint2) if row_num == 18 ;
replace col1 = string(T12_upper_joint3) if row_num == 19 ;

replace col2 = string(T12_upper_lt_joint1) if row_num == 17 ; 
replace col2 = string(T12_upper_lt_joint2) if row_num == 18 ;
replace col2 = string(T12_upper_lt_joint3) if row_num == 19 ;

replace col3 = string(T12_upper_both_joint1) if row_num == 17 ; 
replace col3 = string(T12_upper_both_joint2) if row_num == 18 ;
replace col3 = string(T12_upper_both_joint3) if row_num == 19 ;

replace col1 = string(T12_lower_joint1) if row_num == 21 ; 
replace col1 = string(T12_lower_joint2) if row_num == 22 ;
replace col1 = string(T12_lower_joint3) if row_num == 23 ;

replace col2 = string(T12_lower_lt_joint1) if row_num == 21 ; 
replace col2 = string(T12_lower_lt_joint2) if row_num == 22 ;
replace col2 = string(T12_lower_lt_joint3) if row_num == 23 ;

replace col3 = string(T12_lower_both_joint1) if row_num == 21 ; 
replace col3 = string(T12_lower_both_joint2) if row_num == 22 ;
replace col3 = string(T12_lower_both_joint3) if row_num == 23 ;

replace col1 = string(T2_T3_joint1) if row_num == 25 ; 
replace col1 = string(T2_T3_joint2) if row_num == 26 ;
replace col1 = string(T2_T3_joint3) if row_num == 27 ;

replace col2 = string(T2_T3_lt_joint1) if row_num == 25 ; 
replace col2 = string(T2_T3_lt_joint2) if row_num == 26 ;
replace col2 = string(T2_T3_lt_joint3) if row_num == 27 ;

replace col3 = string(T2_T3_both_joint1) if row_num == 25 ; 
replace col3 = string(T2_T3_both_joint2) if row_num == 26 ;
replace col3 = string(T2_T3_both_joint3) if row_num == 27 ;

   sort row_num  ;
   outsheet col0 col1-col4 using "Table_Appen_C1.xls" if _n <= 30, replace  ;
   drop col* row_num  ;
 
*/

   
log close;


