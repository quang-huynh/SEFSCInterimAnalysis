//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##
//##  SEDAR 28 Assessment: Spanish mackerel, May 2012
//##
//##  NMFS, Beaufort Lab, Sustainable Fisheries Branch
//##
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>


//Inputs selectivity for discards (as known) and estimates different time dep. F's for discards & landings 
//Uses interpolated landings instead of interpolated F's
//Uses an input value for historical F rather than setting it = F(1950)
//Inputs historical selectivity as well.


DATA_SECTION

!!cout << "Starting Spanish Mackerel Assessment Model" << endl;

// Starting and ending year of the model (year data starts)
init_int styr;
init_int endyr;
//Starting year to estimate recruitment deviation from S-R curve
init_int styr_rec_dev;

//Total number of ages
init_int nages;

!!cout << "Number of ages" << nages << endl;

// Vector of ages for age bins
init_vector agebins(1,nages);
!!cout << "Age bins" << agebins << endl;
//number assessment years
//int styrR;
number nyrs;
number nyrs_rec;
//this section MUST BE INDENTED!!!
 LOCAL_CALCS
   nyrs=endyr-styr+1.;
   nyrs_rec=endyr-styr_rec_dev+1.;
 END_CALCS

//Total number of length bins for each matrix
init_int nlenbins;

// Vector of lengths for length bins (cm)(midpoint)
init_vector lenbins(1,nlenbins);

//discard mortality constants
//init_number set_Dmort_HL;
//init_number set_Dmort_GN;
init_number set_Dmort_MRFSS;

//initial proportion female
init_number set_prop_f_a0;
!!cout << "proportion female" << set_prop_f_a0 << endl;
//Total number of iterations for spr calcs
init_int n_iter_spr;
!!cout << "spr iterationsn" << n_iter_spr << endl;
//Total number of iterations for msy calcs
init_int n_iter_msy;
!!cout << "MSY iterations" << n_iter_msy << endl;
//starting index of ages for exploitation rate: if model has age-0s, ages of E are (value-1) to oldest    
init_int set_E_age_st; 
!!cout << "Exploitation age" << set_E_age_st << endl;
//bias correction (set to 1.0 for no bias correction or 0.0 to compute from rec variance)
init_number set_BiasCor;
!!cout << "Bias correction" << set_BiasCor << endl;
// Von Bert parameters for males and females, respectively
//init_number set_Linf_m;
//init_number set_K_m;
//!!cout<<"Linf males"<< set_Linf_m<<endl;
//init_number set_t0_m;
//init_number set_Linf_f;
//!!cout<<"Linf females"<< set_Linf_f<<endl;
//init_number set_K_f;
//init_number set_t0_f;
//!!cout<<"t0 females"<< set_t0_f<<endl;
//CV of length at age
init_number set_len_cv;
 
//length(mm)-weight(whole weight in g) relationship: W=aL^b
init_number wgtpar_a;
init_number wgtpar_b;
init_number wgtpar_am;
init_number wgtpar_bm;
init_number wgtpar_af;
init_number wgtpar_bf;

//weight-weight relationship:whole weight to gutted weight -- gutted=a*whole
//init_number wgtpar_whole2gutted  NOT CURRENTLY USED
init_vector meanlen_fishery(1,nages);
init_vector meanlen_m(1,nages);
init_vector meanlen_f(1,nages);
!!cout<<"mean length fishery"<<meanlen_fishery<<endl;
!!cout<<"mean length males"<<meanlen_m<<endl;
!!cout<<"mean length females"<<meanlen_f<<endl;

//Female maturity and proportion female at age
init_vector maturity_f_obs(1,nages);  
!!cout << "maturity" << maturity_f_obs(1,nages) << endl;          //total maturity of females
init_vector prop_f_obs(1,nages);                //proportion female at age
//#############################################################################
//##############################Commercial Hand Lines fishery ########################
//CPUE
//FL trip ticket
init_int styr_FL_HL_cpue;
//  !!cout << "start year of FL HL" << styr_FL_HL_cpue << endl;                                             
init_int endyr_FL_HL_cpue;                                            
init_vector obs_FL_HL_cpue(styr_FL_HL_cpue,endyr_FL_HL_cpue);//Observed CPUE
init_vector FL_HL_cpue_cv(styr_FL_HL_cpue,endyr_FL_HL_cpue); //CV of cpue
//  !!cout << "FL HL cpue cvs " << FL_HL_cpue_cv << endl;
// Landings (1000 lb whole weight)
init_int styr_HL_L;
init_int endyr_HL_L;
init_vector obs_HL_L(styr_HL_L,endyr_HL_L);
init_vector HL_L_cv(styr_HL_L,endyr_HL_L);    //vector of CV of landings by year
//  !!cout << "start year of HL Landings" <<styr_HL_L << endl;
// Discards (1000s)
//init_int styr_HL_D;
//init_int endyr_HL_D;
//init_vector obs_HL_released(styr_HL_D,endyr_HL_D); //vector of observed releases by year, multiplied by discard mortality for fitting  
//init_vector HL_D_cv(styr_HL_D,endyr_HL_D);    //vector of CV of discards by year
// Length Compositions (1cm bins)
//init_int nyr_HL_lenc;
//init_ivector yrs_HL_lenc(1,nyr_HL_lenc);
//init_vector nsamp_HL_lenc(1,nyr_HL_lenc);
//init_matrix obs_HL_lenc(1,nyr_HL_lenc,1,nlenbins);
  !!cout << "start year of FL HL" <<styr_FL_HL_cpue << endl;
// Age Compositions 
init_int nyr_HL_agec;
init_ivector yrs_HL_agec(1,nyr_HL_agec);
init_vector nsamp_HL_agec(1,nyr_HL_agec);
init_matrix obs_HL_agec(1,nyr_HL_agec,1,nages);
//  !!cout << "Handline data" << obs_HL_agec << endl;
//#############################################################################
//################################Poundnet fishery ########################################
// Landings (1000 lb whole weight)
init_int styr_PN_L;
init_int endyr_PN_L;
init_vector obs_PN_L(styr_PN_L,endyr_PN_L);
init_vector PN_L_cv(styr_PN_L,endyr_PN_L);
// Length Compositions (1cm bins)
//init_int nyr_PN_lenc;
//init_ivector yrs_PN_lenc(1,nyr_PN_lenc);
//init_vector nsamp_PN_lenc(1,nyr_PN_lenc);
//init_matrix obs_PN_lenc(1,nyr_PN_lenc,1,nlenbins);
// Age compositions
init_int nyr_PN_agec;
init_ivector yrs_PN_agec(1,nyr_PN_agec);
init_vector nsamp_PN_agec(1,nyr_PN_agec);
init_matrix obs_PN_agec(1,nyr_PN_agec,1,nages);
//  !!cout << "Poundnet data" <<nyr_PN_agec << endl;
//#############################################################################
//###################Commercial Gillnet fishery #########################
// Landings - (1000 lb whole weight)
init_int styr_GN_L;
init_int endyr_GN_L;
init_vector obs_GN_L(styr_GN_L,endyr_GN_L); //vector of observed landings by year 
init_vector GN_L_cv(styr_GN_L,endyr_GN_L);    //vector of CV of landings by year

// Length Compositions (1cm bins)
//init_int nyr_GN_lenc;
//init_ivector yrs_GN_lenc(1,nyr_GN_lenc);
//init_vector nsamp_GN_lenc(1,nyr_GN_lenc);
//init_matrix obs_GN_lenc(1,nyr_GN_lenc,1,nlenbins);
// Age compositions
init_int nyr_GN_agec;
init_ivector yrs_GN_agec(1,nyr_GN_agec);
init_vector nsamp_GN_agec(1,nyr_GN_agec);
init_matrix obs_GN_agec(1,nyr_GN_agec,1,nages);
  !!cout << "Gillnet data" <<obs_GN_agec << endl;
//#############################################################################
//################################Castnet fishery ########################################
// Landings (1000 lb whole weight)
init_int styr_CN_L;
init_int endyr_CN_L;
init_vector obs_CN_L(styr_CN_L,endyr_CN_L);
init_vector CN_L_cv(styr_CN_L,endyr_CN_L);
// Length Compositions (1cm bins)
//init_int nyr_CN_lenc;
//init_ivector yrs_CN_lenc(1,nyr_CN_lenc);
//init_vector nsamp_CN_lenc(1,nyr_CN_lenc);
//init_matrix obs_CN_lenc(1,nyr_CN_lenc,1,nlenbins);
// Age compositions
init_int nyr_CN_agec;
init_ivector yrs_CN_agec(1,nyr_CN_agec);
init_vector nsamp_CN_agec(1,nyr_CN_agec);
init_matrix obs_CN_agec(1,nyr_CN_agec,1,nages);
//  !!cout << "Castnet data" <<obs_CN_agec << endl;
//#############################################################################
//############################MRFSS landings #################################
//CPUE
init_int styr_MRFSS_cpue;
init_int endyr_MRFSS_cpue;
init_vector obs_MRFSS_cpue(styr_MRFSS_cpue,endyr_MRFSS_cpue);//Observed CPUE
init_vector MRFSS_cpue_cv(styr_MRFSS_cpue,endyr_MRFSS_cpue); //CV of cpue
//  !!cout << "MRFSS cpue" <<MRFSS_cpue_cv << endl;
// Landings (1000 of fish)
init_int styr_MRFSS_L;
init_int endyr_MRFSS_L;
init_vector obs_MRFSS_L(styr_MRFSS_L,endyr_MRFSS_L);
init_vector MRFSS_L_cv(styr_MRFSS_L,endyr_MRFSS_L);
  !!cout << "MRFSS Landings cv" <<MRFSS_L_cv << endl;
// Discards (1000s)
init_int styr_MRFSS_D;
init_int endyr_MRFSS_D;
init_vector obs_MRFSS_released(styr_MRFSS_D,endyr_MRFSS_D); //vector of observed releases by year, multiplied by discard mortality for fitting 
init_vector MRFSS_D_cv(styr_MRFSS_D,endyr_MRFSS_D);    //vector of CV of discards by year
// Length Compositions (1cm bins)
//init_int styr_MRFSS_lenc;
//init_int endyr_MRFSS_lenc;
//init_int nyr_MRFSS_lenc;
//init_vector nsamp_MRFSS_lenc(styr_MRFSS_lenc,endyr_MRFSS_lenc);
//init_matrix obs_MRFSS_lenc(styr_MRFSS_lenc,endyr_MRFSS_lenc,1,nlenbins);
// !!cout << "MRFSS length comps" <<obs_MRFSS_lenc << endl;
// Age Compositions 
init_int styr_MRFSS_agec;
init_int endyr_MRFSS_agec;
init_int nyr_MRFSS_agec;
init_ivector yrs_MRFSS_agec(1,nyr_MRFSS_agec);
init_vector nsamp_MRFSS_agec(1,nyr_MRFSS_agec);
init_matrix obs_MRFSS_agec(1,nyr_MRFSS_agec,1,nages);
//  !!cout << "MRFSS" <<obs_MRFSS_agec << endl; 
	!!cout<<"nyrs mrfss"<< nyr_MRFSS_agec<<endl;
//#############################################################################
//############################Shrimp Bycatch #################################
// Bycatch (1000s of fish)
init_int styr_shrimp_B;
init_int endyr_shrimp_B;
init_vector obs_shrimp_B(styr_shrimp_B,endyr_shrimp_B);
init_vector shrimp_B_cv(styr_shrimp_B,endyr_shrimp_B);
//init_vector obs_shrimp_B_effort(styr_shrimp_B,endyr_shrimp_B);
//init_vector shrimp_B_cv(styr_shrimp_B,endyr_shrimp_B);
//init_int styr_shrimp_B_fit;  //start year for computing geometric mean of shrimp bycatch for fitting
//init_int endyr_shrimp_B_fit; //end year for computing geometric mean of shrimp bycatch for fitting

//   !!cout << "Shrimp bycatch" <<shrimp_B_cv << endl;
//#############################################################################
//############################SEAMAP#################################
//YOY CPUE
init_int styr_SMAP_YOY_cpue;
init_int endyr_SMAP_YOY_cpue;
init_vector obs_SMAP_YOY_cpue(styr_SMAP_YOY_cpue,endyr_SMAP_YOY_cpue);//Observed CPUE
init_vector SMAP_YOY_cpue_cv(styr_SMAP_YOY_cpue,endyr_SMAP_YOY_cpue); //CV of cpue
//1-yr-old CPUE
//init_int styr_SMAP_1YR_cpue;
//init_int endyr_SMAP_1YR_cpue;
//init_vector obs_SMAP_1YR_cpue(styr_SMAP_1YR_cpue,endyr_SMAP_1YR_cpue);//Observed CPUE
//init_vector SMAP_1YR_cpue_cv(styr_SMAP_1YR_cpue,endyr_SMAP_1YR_cpue); //CV of cpue
//   !!cout << "SEAMAP 1-yr" <<obs_SMAP_1YR_cpue << endl;
//#############################################################################
//##################Parameter values and initial guesses #################################
//--weights for likelihood components-------------------------------------------------------------------------------
init_number set_w_L;
init_number set_w_D;
//init_number set_w_shrimp;
//init_number set_w_lc;
//init_number set_w_lc_HL;
//init_number set_w_lc_PN;
//init_number set_w_lc_GN;
//init_number set_w_lc_CN;
//init_number set_w_lc_MRFSS;
//init_number set_w_ac;
init_number set_w_ac_HL;
init_number set_w_ac_PN;
init_number set_w_ac_GN;
init_number set_w_ac_CN;
init_number set_w_ac_MRFSS;
init_number set_w_I_FL_HL;
init_number set_w_I_MRFSS;
init_number set_w_I_SMAP_YOY;
//init_number set_w_I_SMAP_1YR;
init_number set_w_R;
init_number set_w_R_init;
init_number set_w_R_end;
init_number set_w_F;
init_number set_w_B1dB0;         // weight on B1/B0
init_number set_w_fullF;         //penalty for any fullF>5
init_number set_w_cvlen_dev;         //penalty on cv deviations at age
init_number set_w_cvlen_diff;       //penalty on first difference of cv deviations at age

//Initial guesses or fixed values
init_number set_steep;
init_number set_steep_se;      //SE of recruitment steepness
//init_number set_M;
init_vector set_M(1,nages); //age-dependent: used in model
init_number set_M_constant; //age-independent: used only for MSST
init_number set_rec_sigma;     //recruitment standard deviation in log space
init_number set_rec_sigma_se;  //SE of recruitment standard deviation in log space

//--index catchability------------------------------------------------------------------------------------------------------------
 
init_number set_logq_FL_HL;         //catchability coefficient (log) for FL hand lines           
init_number set_logq_MRFSS;   
init_number set_logq_SMAP_YOY;       
//init_number set_logq_SMAP_1YR;           
//init_number set_logq_shrimp_B;
//--F's--------------------------------
init_number set_F_hist;
  !!cout << "historic F" <<set_F_hist << endl;
init_number set_log_avg_F_HL;  //hand lines
init_vector set_log_F_dev_HL(styr_HL_L,endyr_HL_L);  
init_number set_log_avg_F_PN;  //pound nets
init_vector set_log_F_dev_PN(styr_PN_L,endyr_PN_L);  
init_number set_log_avg_F_GN;  //gill nets
  !!cout << "AVG F GN" <<set_log_avg_F_GN << endl;
init_vector set_log_F_dev_GN(styr_GN_L,endyr_GN_L);  
init_number set_log_avg_F_CN;  //cast nets
init_vector set_log_F_dev_CN(styr_CN_L,endyr_CN_L);
//  !!cout << "F devs CN" <<set_log_F_dev_CN << endl;
init_number set_log_avg_F_MRFSS;  //mrfss
init_vector set_log_F_dev_MRFSS(styr_MRFSS_L,endyr_MRFSS_L);
init_number set_log_avg_F_MRFSS_D;  //mrfss discards
  !!cout << "AVG F MRFSS_D" <<set_log_avg_F_MRFSS_D << endl;
init_vector set_log_F_dev_MRFSS_D(styr_MRFSS_D,endyr_MRFSS_D);
  !!cout << "F devs MRFSS _D" <<set_log_F_dev_MRFSS_D << endl;
init_number set_log_avg_F_shrimp;  //shrimp
init_vector set_log_F_dev_shrimp(styr_shrimp_B,endyr_shrimp_B);

//Set some more initial guesses of estimated parameters
init_number set_log_R0;
init_number set_R_autocorr;
  !!cout << "R autocorr" <<set_R_autocorr << endl;

//Initial guesses of estimated selectivity parameters
init_number set_selpar_L50_HL_keep;
  !!cout << "parm L50 for HL" <<set_selpar_L50_HL_keep << endl;
init_number set_selpar_slope_HL;
	!!cout << "parm slope for HL"<<set_selpar_slope_HL<<endl;
init_number set_selpar_L501_PN;  
	!!cout<<"PN a501"<<set_selpar_L501_PN<<endl;
init_number set_selpar_slope1_PN;
	!!cout<<"PN slope1"<<set_selpar_slope1_PN<<endl;
init_number set_selpar_L502_PN;
	!!cout<<"PN a502"<<set_selpar_L502_PN<<endl;
init_number set_selpar_slope2_PN;
	!!cout<<"PN slope2"<<set_selpar_slope2_PN<<endl;
init_number set_selpar_L50_GN_keep;
!!cout<<"GN a50"<<set_selpar_L50_GN_keep<<endl;
init_number set_selpar_slope_GN;
init_number set_selpar_L501_CN;
init_number set_selpar_slope1_CN;
init_number set_selpar_L502_CN;
init_number set_selpar_slope2_CN;
//init_number set_selpar_sigma_CN;
//init_number set_selpar_afull_CN;
init_number set_selpar_L501_MRFSS_keep;
init_number set_selpar_slope1_MRFSS_keep; 
init_number set_selpar_L502_MRFSS_keep;
init_number set_selpar_slope2_MRFSS_keep;   
//init_number set_selpar_sigma_MRFSS;
//init_number set_selpar_afull_MRFSS;
init_vector set_sel_historical_F(1,nages);
   !!cout << "historical F selex" <<set_sel_historical_F << endl;
init_vector set_sel_MRFSS_D_F(1,nages);  //selectivity vectors for female, male discards input directly
init_vector set_sel_MRFSS_D_M(1,nages);
init_vector set_sel_shrimp_B(1,nages); //shrimp bycatch selectivity input directly
  !!cout << "shrimp bycatch selex" <<set_sel_shrimp_B << endl;

init_number set_L50_diff; //difference between males and females

//aging error matrix
init_matrix set_age_error_matrix(1,nages,1,nages);

//historic recreational landings multiplier
init_number L_rec_multiply;

//threshold sample sizes for length comps 
//init_number minSS_HL_lenc;
//init_number minSS_PN_lenc;
//init_number minSS_GN_lenc;
//init_number minSS_CN_lenc;
//init_number minSS_MRFSS_lenc;
 //cout << "min SS MRFSS" <<minSS_MRFSS_lenc<< endl;
//threshold sample sizes for age comps
init_number minSS_HL_agec;
init_number minSS_PN_agec;
init_number minSS_GN_agec;
init_number minSS_CN_agec;
init_number minSS_MRFSS_agec;
	!!cout<<"minSS for MRFSS"<<minSS_MRFSS_agec<<endl;

// #######Indices for year(iyear), age(iage),length(ilen) ###############
int iyear;
int iyear2;
int iage;
int ilen;
int E_age_st;   //starting age for exploitation rate: (value-1) to oldest 
int ff;
 
init_number end_of_data_file;
//this section MUST BE INDENTED!!!
 LOCAL_CALCS
   if(end_of_data_file!=999)
   {
       cout << "*** WARNING: Data File NOT READ CORRECTLY ****" << endl;
       cout << "" <<endl;
       exit(0);
     
   }
   else
   {
    cout << "Data File read correctly" << endl;
   } 
 END_CALCS 
   
PARAMETER_SECTION 
//number Linf_m; 
//number K_m; 
//number t0_m;  
//number Linf_f;
//number K_f; 
//number t0_f; 
//init_bounded_number Linf_m(300,1000,-3); 
//init_bounded_number K_m(0.15,1.0,-3); 
//init_bounded_number t0_m(-5,0,-3); 
//init_bounded_number Linf_f(300,1000,-3); 
//init_bounded_number K_f(0.15,1.0,-3); 
//init_bounded_number t0_f(-5,0,-3);   


//vector meanlen_fishery(1,nages);
//vector meanlen_m(1,nages);
//vector meanlen_f(1,nages);

vector wgt_g_f(1,nages);        //whole wgt in g - females 
vector wgt_kg_f(1,nages);       //whole wgt in kg 
vector wgt_f(1,nages);          //whole wgt in mt 
vector wgt_klb_f(1,nages); //whole wgt in 1000 lb 
vector wgt_g_land(1,nages);        //whole wgt in g - females 
vector wgt_kg_land(1,nages);       //whole wgt in kg 
vector wgt_land(1,nages);          //whole wgt in mt 
vector wgt_klb_land(1,nages); //whole wgt in 1000 lb 
vector wgt_g_m(1,nages);        //whole wgt in g - males 
vector wgt_kg_m(1,nages);       //whole wgt in kg 
vector wgt_m(1,nages); //whole wgt in mt 
vector wgt_klb_m(1,nages);      //whole wgt in 1000 lb 
//vector meanlen_m(1,nages);        //mean length at age -males 
//vector meanlen_f(1,nages);        //mean length at age -females 
number sqrt2pi; 
number g2mt;                    //conversion of grams to metric tons 
number g2kg; //conversion of grams to kg   
number mt2klb;  //conversion of metric tons to 1000 lb
  
number nlenbins2;
  matrix lenprob_m(1,nages,1,nlenbins);           //distn of size at age (age-length key, 1cm bins) - males
  matrix lenprob_f(1,nages,1,nlenbins);           //distn of size at age (age-length key, 1cm bins) - females
  matrix lenprob_m2(1,nages,1,nlenbins+20);           //distn of size at age (age-length key, 1cm bins) - males
  matrix lenprob_f2(1,nages,1,nlenbins+20);           //distn of size at age (age-length key, 1cm bins) - females
  vector lenbins2(1,nlenbins+20);

  init_bounded_number log_len_cv(-5,-0.3,-4);  //KWS add a fourth phase
  //init_bounded_number log_len_cv(-4.6,-0.7,2) //cv expressed in log-space, bounds correspond to 0.01, 0.5
  //init_bounded_dev_vector log_len_cv_dev(1,nages,-2,2,3)
  vector len_cv(1,nages);

//----Predicted length and age compositions
  number prop_f_a0; //proportion female at age 0

  //matrix pred_HL_lenc(1,nyr_HL_lenc,1,nlenbins);
  //matrix pred_PN_lenc(1,nyr_PN_lenc,1,nlenbins);
  //matrix pred_GN_lenc(1,nyr_GN_lenc,1,nlenbins);
  //matrix pred_CN_lenc(1,nyr_CN_lenc,1,nlenbins);
  //matrix pred_MRFSS_lenc(styr_MRFSS_lenc,endyr_MRFSS_lenc,1,nlenbins);
   
  matrix pred_HL_agec(1,nyr_HL_agec,1,nages);
  matrix pred_PN_agec(1,nyr_PN_agec,1,nages);
  matrix pred_GN_agec(1,nyr_GN_agec,1,nages);
  matrix pred_CN_agec(1,nyr_CN_agec,1,nages);
  matrix pred_MRFSS_agec(1,nyr_MRFSS_agec,1,nages);
  
  //nsamp_X_allyr vectors used only for R output of comps with nonconsecutive yrs
  //vector nsamp_HL_lenc_allyr(styr,endyr);
  vector nsamp_HL_agec_allyr(styr,endyr);
  //vector nsamp_PN_lenc_allyr(styr,endyr);
  vector nsamp_PN_agec_allyr(styr,endyr);
  //vector nsamp_GN_lenc_allyr(styr,endyr);
  vector nsamp_GN_agec_allyr(styr,endyr);
  //vector nsamp_CN_lenc_allyr(styr,endyr);
  vector nsamp_CN_agec_allyr(styr,endyr);
  //vector nsamp_MRFSS_lenc_allyr(styr,endyr);
  vector nsamp_MRFSS_agec_allyr(styr,endyr);
   
//----Aging error
  matrix age_error_matrix(1,nages,1,nages);

//-----Population-----------------------------------------------------------------------------------
  matrix N_F(styr,endyr+1,1,nages);           //Population numbers for females by year and age at start of yr
  matrix N_M(styr,endyr+1,1,nages);           //Population numbers for males by year and age at start of yr
  matrix N_F_mdyr(styr,endyr+1,1,nages);        //Population numbers by year and age at mdpt of yr: used for comps and SSB
  matrix N_M_mdyr(styr,endyr+1,1,nages);
  matrix B(styr,endyr+1,1,nages);             //Biomass by year and age - sexes combined
  vector totB(styr,endyr+1);                //Total biomass by year
  vector SSB(styr,endyr);       //Spawning biomass by year
  number SSB_extra;
  vector rec(styr,endyr+1);       //Recruits by year
  matrix prop_f(styr,endyr+1,1,nages);        //Proportion female by year and age
  vector prop_f_F0(1,nages);                //proportion of females in unexploited pop
  vector maturity_f(1,nages);               //Proportion of female mature at age
  vector reprod(1,nages);                 //recruitment calcs

//---Stock-Recruit Function (Beverton-Holt, steepness parameterization)----------
  init_bounded_number log_R0(5,100,1);       //log(virgin Recruitment)
  //number log_R0;
  number R0;
  init_bounded_number steep(0.25,0.9,-3);     //steepness
  //number steep;  //uncomment to fix steepness, comment line directly above
  init_bounded_dev_vector log_rec_dev(styr_rec_dev,endyr,-5,5,3); //log recruitment deviations
  //vector log_dev_N_rec(styr_rec_dev,endyr);
  vector log_dev_R(styr,endyr+1);           //used in output. equals zero except for yrs in log_dev_N_rec
  number var_rec_dev; 
  number sigma_rec_dev;                      //variance of log recruitment deviations. 
  init_bounded_number rec_sigma(0.5,1.2,-4);  //sd recruitment residuals
  number rec_sigma_sq;                        //square of rec_sigma      
  number rec_logL_add;                        //additive term in -logL term   
                                                 //Estimate from yrs with unconstrainted S-R(XXXX-XXXX)
  number BiasCor;                           //Bias correction in equilibrium recruits
  init_bounded_number R_autocorr(0,1.0,-3);  //autocorrelation in SR  //KWS turned off
  //number R_autocorr;
  number R_autocorr_sd;
  number steep_sd;                 //steepness for stdev report
  number S0;                                //equal to spr_F0*R0 = virgin SSB - for males = 0 
  number B0;                              //equal to bpr_F0*R0 = virgin Biomass (combined sex)
  number R1;                                //Recruits in styr

  number S1S0;                     //SSB(styr) / virgin SSB
  number popstatus;                //SSB(endyr) / virgin SSB

//---Selectivity-------------------------------------------------------------------------

  number L50_diff; //shift in selectivity b/w males and females (same for all gears)
  vector sel_historical_F(1,nages);
    
   //Commercial handline
  vector sel_HL_keep_M(1,nages);  //time invariant - fish that are kept - males
  vector sel_HL_keep_F(1,nages);  //time invariant - fish that are kept - females
  init_bounded_number selpar_slope_HL(0.1,10,3); 
  init_bounded_number selpar_L50_HL_keep(0.5,8,2);
    
  //Poundnets
  vector sel_PN_M(1,nages);  //time invariant - males
  vector sel_PN_F(1,nages);  //females
  init_bounded_number selpar_slope1_PN(0.1,10,3); 
  init_bounded_number selpar_L501_PN(0.4,8,2);
  init_bounded_number selpar_slope2_PN(0.1,10,3); 
  init_bounded_number selpar_L502_PN(-2.0,8,2);
  
  //Commercial gillnet
  vector sel_GN_keep_M(1,nages);  //pre-1995 - fish that are kept - males
  vector sel_GN_keep_F(1,nages);  //pre-1995 - fish that are kept - females
  init_bounded_number selpar_slope_GN(0.5,10.0,3); 
  init_bounded_number selpar_L50_GN_keep(1.0,8,2);
  
  //Castnets
  vector sel_CN_M(1,nages);  //time invariant - males
  vector sel_CN_F(1,nages);  //females
  init_bounded_number selpar_slope1_CN(0.5,10.0,3); 
  init_bounded_number selpar_L501_CN(2.0,8,2);  
  init_bounded_number selpar_slope2_CN(0.5,10.0,3); 
  init_bounded_number selpar_L502_CN(2.0,10.0,2);      
  //init_bounded_number selpar_sigma_CN(0.1,10000.0,-3);
  //number selpar_afull_CN;
  
  //MRFSS
  vector sel_MRFSS_keep_M(1,nages);  //time invariant - fish that are kept - males
  vector sel_MRFSS_keep_F(1,nages);  //time invariant - fish that are kept - females
  vector sel_MRFSS_D_M(1,nages);  //selectivity for discards - males
  vector sel_MRFSS_D_F(1,nages);  //selectivity for discards - females    
  init_bounded_number selpar_slope1_MRFSS_keep(0.1,15,3); 
  init_bounded_number selpar_L501_MRFSS_keep(0.5,8,2);
	init_bounded_number selpar_slope2_MRFSS_keep(0.1,15,3);     
	init_bounded_number selpar_L502_MRFSS_keep(-2.0,8,2);   
  //init_bounded_number selpar_sigma_MRFSS(0.1,10000.0,3);
  //number selpar_afull_MRFSS;
 
  //shrimp
  vector sel_shrimp(1,nages);
 
  //effort-weighted, recent selectivities
  vector sel_wgted_L_F(1,nages);  //toward landings,females
  vector sel_wgted_D_F(1,nages);  //toward discards  
  vector sel_wgted_tot_F(1,nages);//toward Z, landings plus dead discards
  number max_sel_wgted_tot_F;
  vector sel_wgted_L_M(1,nages);  //toward landings,males
  vector sel_wgted_D_M(1,nages);  //toward discards  
  vector sel_wgted_tot_M(1,nages);//toward Z, landings plus dead discards
  number max_sel_wgted_tot_M;

//-------CPUE Predictions--------------------------------
  vector pred_FL_HL_cpue(styr_FL_HL_cpue,endyr_FL_HL_cpue);             //predicted FL handline index (pounds/trip)
  matrix N_FL_HL(styr_FL_HL_cpue,endyr_FL_HL_cpue,1,nages);             //used to compute above index
  vector pred_MRFSS_cpue(styr_MRFSS_cpue,endyr_MRFSS_cpue);             //predicted MRFSS index (number per angler-trip)
  matrix N_MRFSS(styr_MRFSS_cpue,endyr_MRFSS_cpue,1,nages);             //used to compute above index
  vector pred_SMAP_YOY_cpue(styr_SMAP_YOY_cpue,endyr_SMAP_YOY_cpue);             //predicted SMAP_YOY index (number per tow)
  //vector pred_SMAP_1YR_cpue(styr_SMAP_1YR_cpue,endyr_SMAP_1YR_cpue);             //predicted SMAP_1YR index (number per tow)  

//---Catchability (CPUE q's)----------------------------------------------------------
  
  init_bounded_number log_q_FL_HL(-20,-5,1);
  init_bounded_number log_q_MRFSS(-40,-5,1);  //KWs adjusted lower bound
  init_bounded_number log_q_SMAP_YOY(-40,-5,1);
  //init_bounded_number log_q_SMAP_1YR(-40,-5,1);
  //init_bounded_number log_q_shrimp_B(-30,-5,1);

  //init_bounded_number q_rate(-0.1,0.1,-3);
  number q_rate;
  
//---Landings Bias------------------------------------------------------------------
  //init_bounded_number L_early_bias(0.1,10.0,3);  
  //number L_early_bias;
  //number L_commHAL_bias;

//---C is landings in (numbers), L is landings in wgt (mt)--------------------------------------------------

  matrix C_HL_M(styr,endyr,1,nages);               //landings (numbers) at age
  matrix C_HL_F(styr,endyr,1,nages);               //landings (numbers) at age
  matrix L_HL(styr,endyr,1,nages);                 //landings (mt) at age
  matrix L_HL_klb(styr,endyr,1,nages);             //landings (1000 lb whole weight) at age
  vector pred_HL_L(styr_HL_L,endyr_HL_L);          //yearly landings (klb) summed over ages 
  
  matrix C_PN_M(styr,endyr,1,nages);               //landings (numbers) at age
  matrix C_PN_F(styr,endyr,1,nages);               //landings (numbers) at age
  matrix L_PN(styr,endyr,1,nages);                 //landings (mt) at age
  matrix L_PN_klb(styr,endyr,1,nages);             //landings (1000 lb whole weight) at age    
  vector pred_PN_L(styr_PN_L,endyr_PN_L);          //yearly landings (klb) summed over ages

  matrix C_GN_M(styr,endyr,1,nages);               //landings (numbers) at age
  matrix C_GN_F(styr,endyr,1,nages);               //landings (numbers) at age
  matrix L_GN(styr,endyr,1,nages);                 //landings (mt) at age
  matrix L_GN_klb(styr,endyr,1,nages);             //landings (1000 lb whole weight) at age
  vector pred_GN_L(styr_GN_L,endyr_GN_L);          //yearly landings (klb) summed over ages 
  
  matrix C_total(styr,endyr,1,nages);              //catch in number 
  matrix L_total(styr,endyr,1,nages);              //landings in klb 
  vector L_total_yr(styr,endyr);                   //total landings (klb) by yr summed over ages
  matrix D_total(styr,endyr,1,nages);              //discards in klb 
  vector D_total_yr(styr,endyr);                   //total discards (klb) by yr summed over ages
  matrix B_total(styr,endyr,1,nages);              //bycatch in klb 
  vector B_total_yr(styr,endyr);                   //total bycatch (klb) by yr summed over ages

  matrix C_CN_M(styr,endyr,1,nages);               //landings (numbers) at age
  matrix C_CN_F(styr,endyr,1,nages);               //landings (numbers) at age
  matrix L_CN(styr,endyr,1,nages);                 //landings (mt) at age
  matrix L_CN_klb(styr,endyr,1,nages);             //landings (1000 lb whole weight) at age    
  vector pred_CN_L(styr_CN_L,endyr_CN_L);          //yearly landings (klb) summed over ages 

  matrix C_MRFSS_M(styr,endyr,1,nages);            //landings (numbers) at age
  matrix C_MRFSS_F(styr,endyr,1,nages);            //landings (numbers) at age
  matrix L_MRFSS(styr,endyr,1,nages);              //landings (mt) at age
  matrix L_MRFSS_klb(styr,endyr,1,nages);          //landings (1000 lb whole weight) at age    
  vector pred_MRFSS_L(styr_MRFSS_L,endyr_MRFSS_L); //yearly landings (numbers in 1000's) summed over ages
  
//---Discards (number dead fish) --------------------------------------------------
  
 
  matrix C_MRFSS_D_M(styr_MRFSS_D,endyr_MRFSS_D,1,nages);     //discards (numbers) at age
  matrix C_MRFSS_D_F(styr_MRFSS_D,endyr_MRFSS_D,1,nages);     //discards (numbers) at age
  matrix L_MRFSS_D(styr,endyr,1,nages);                       //discards in mt whole weight 
  vector pred_MRFSS_D(styr_MRFSS_D,endyr_MRFSS_D);            //yearly discards summed over ages
  vector obs_MRFSS_D(styr_MRFSS_D,endyr_MRFSS_D);             //observed releases multiplied by discard mortality

//---Bycatch  ------------------------------------------------------------------
  matrix C_shrimp_M(styr,endyr,1,nages);
  matrix C_shrimp_F(styr,endyr,1,nages);
  matrix L_shrimp(styr,endyr,1,nages);        //bycatch in mt whole weight 
  vector pred_shrimp_B(styr,endyr);
 // matrix C_shrimp_M(styr_shrimp_B,endyr_shrimp_B,1,nages);
 // matrix C_shrimp_F(styr_shrimp_B,endyr_shrimp_B,1,nages);
 // matrix L_shrimp(styr,endyr,1,nages);        //bycatch in mt whole weight 
 // vector pred_shrimp_B(styr_shrimp_B,endyr_shrimp_B);
  
  //vector log_obs_shrimp_B(styr_shrimp_B_fit,endyr_shrimp_B_fit);
  //vector log_pred_shrimp_B(styr_shrimp_B_fit,endyr_shrimp_B_fit);
  //number obs_shrimp_geomean;
  //number pred_shrimp_geomean;
//---MSY calcs----------------------------------------------------------------------------

  number F_HL_prop;       //proportion of F_full attributable to hand lines, last three yrs
  number F_PN_prop;
  number F_GN_prop;
  number F_CN_prop; 
  number F_MRFSS_prop; 
  number F_MRFSS_D_prop;
  number F_shrimp_prop;
  number F_temp_sum;      //sum of geom mean full Fs in last yrs, used to compute F_fishery_prop

  number SSB_msy_out;           //SSB at msy
  number F_msy_out;             //F at msy
  number msy_out;               //max sustainable yield
  number B_msy_out;             //total biomass at MSY 
  number E_msy_out;             //exploitation rate at MSY (ages E_age_st plus)  
  number R_msy_out;             //equilibrium recruitment at F=Fmsy
  number D_msy_out;             //equilibrium dead discards at F=Fmsy
  number spr_msy_out;           //spr at F=Fmsy

  vector N_age_msy_F(1,nages);         //numbers at age for MSY calculations: beginning of yr - Females
  vector N_age_msy_mdyr_F(1,nages);    //numbers at age for MSY calculations: mdpt of yr  
  vector C_age_msy_F(1,nages);         //catch at age for MSY calculations
  vector Z_age_msy_F(1,nages);         //total mortality at age for MSY calculations
  vector D_age_msy_F(1,nages);         //discard mortality (dead discards) at age for MSY calculations
  vector F_L_age_msy_F(1,nages);       //fishing mortality (landings, not discards) at age for MSY calculations
  vector F_D_age_msy_F(1,nages); 
  vector N_age_msy_M(1,nages);         //numbers at age for MSY calculations: beginning of yr - Males
  vector N_age_msy_mdyr_M(1,nages);    //numbers at age for MSY calculations: mdpt of yr  
  vector C_age_msy_M(1,nages);         //catch at age for MSY calculations
  vector Z_age_msy_M(1,nages);         //total mortality at age for MSY calculations
  vector D_age_msy_M(1,nages);         //discard mortality (dead discards) at age for MSY calculations
  vector F_L_age_msy_M(1,nages);       //fishing mortality (landings, not discards) at age for MSY calculations
  vector F_D_age_msy_M(1,nages);   
  vector F_msy(1,n_iter_msy);        //values of full F to be used in per-recruit and equilibrium calculations
  vector spr_msy(1,n_iter_msy);      //reproductive capacity-per-recruit values corresponding to F values in F_msy
  vector R_eq(1,n_iter_msy);     //equilibrium recruitment values corresponding to F values in F_msy
  vector L_eq(1,n_iter_msy);     //equilibrium landings(mt) values corresponding to F values in F_msy
  vector SSB_eq(1,n_iter_msy);   //equilibrium reproductive capacity values corresponding to F values in F_msy
  vector B_eq(1,n_iter_msy);     //equilibrium biomass values corresponding to F values in F_msy
  vector E_eq(1,n_iter_msy);     //equilibrium exploitation rates corresponding to F values in F_msy
  vector D_eq(1,n_iter_msy);     //equilibrium discards (1000s) corresponding to F values in F_msy

  vector FdF_msy(styr,endyr);
  vector EdE_msy(styr,endyr);  
  vector SdSSB_msy(styr,endyr);
  number SdSSB_msy_end;
  number FdF_msy_end;
  number EdE_msy_end;

//--------Mortality------------------------------------------------------------------
  vector M(1,nages);                         //age-dependent natural mortality
  number M_constant;                         //age-indpendent: used only for MSST
  matrix F_M(styr,endyr,1,nages);            //males
  matrix F_F(styr,endyr,1,nages);            //females  

  number F_hist;                        //historical F input into assessment model
  vector F_F_hist(1,nages);                           //historical F used to get stable age dist in 1st year
  vector Z_hist(1,nages);
  vector fullF(styr,endyr);                   //Fishing mortality rate by year
  vector E(styr,endyr);                       //Exploitation rate by year    
  vector fullF_sd(styr,endyr);
  vector E_sd(styr,endyr);
  matrix Z_M(styr,endyr,1,nages);             //males
  matrix Z_F(styr,endyr,1,nages);             //females

  init_bounded_number log_avg_F_HL(-10,1,1);
  init_bounded_dev_vector log_F_dev_HL(styr_HL_L,endyr_HL_L,-20,5,2);
  matrix F_HL_M(styr,endyr,1,nages);
  matrix F_HL_F(styr,endyr,1,nages);
  vector F_HL_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  //number log_F_init_commDV;
  
  init_bounded_number log_avg_F_PN(-10,0,1);
  init_bounded_dev_vector log_F_dev_PN(styr_PN_L,endyr_PN_L,-20,5,2);
  matrix F_PN_M(styr,endyr,1,nages);
  matrix F_PN_F(styr,endyr,1,nages);
  vector F_PN_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality

  init_bounded_number log_avg_F_GN(-10,0,1);
  init_bounded_dev_vector log_F_dev_GN(styr_GN_L,endyr_GN_L,-20,5,2);
  matrix F_GN_M(styr,endyr,1,nages);
  matrix F_GN_F(styr,endyr,1,nages);
  vector F_GN_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  
  init_bounded_number log_avg_F_CN(-10,0,1);
  init_bounded_dev_vector log_F_dev_CN(styr_CN_L,endyr_CN_L,-20,5,2);
  matrix F_CN_M(styr,endyr,1,nages);
  matrix F_CN_F(styr,endyr,1,nages);
  vector F_CN_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality

  init_bounded_number log_avg_F_MRFSS(-10,0,1);
  init_bounded_dev_vector log_F_dev_MRFSS(styr_MRFSS_L,endyr_MRFSS_L,-20,5,2);
  matrix F_MRFSS_M(styr,endyr,1,nages);
  matrix F_MRFSS_F(styr,endyr,1,nages);
  vector F_MRFSS_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number MRFSS_ratio;  //ratio of F_MRFSS/F_GN for saltwater report years to apply to missing years;
                       //straight interpolation used from last angler report to first 'real' year of landings
  
//--Discard mortality stuff------------------------------------------------------------------------------
  
  init_bounded_number log_avg_F_MRFSS_D(-15,0,1);
  init_bounded_dev_vector log_F_dev_MRFSS_D(styr_MRFSS_D,endyr_MRFSS_D,-20,5,2);
  matrix F_MRFSS_D_M(styr,endyr,1,nages);
  matrix F_MRFSS_D_F(styr,endyr,1,nages);
  vector F_MRFSS_D_out(styr,endyr); //used for intermediate calculations in fMRFSS get_mortality
  number F_MRFSS_D_ratio;  //ratio of average discard F to fishery F, for projection discards back before data
   
  number Dmort_MRFSS;
  
//--------Bycatch -----------------------------;
   init_bounded_number log_avg_F_shrimp(-15,0,1);
  init_bounded_dev_vector log_F_dev_shrimp(styr_shrimp_B,endyr_shrimp_B,-20,6,2);
  //vector log_obs_shrimp_B_effort(styr_shrimp_B,endyr_shrimp_B);
  //matrix F_shrimp(styr,endyr,1,nages);
  //vector F_shrimp_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  //vector log_F_shrimp_out(endyr-2,endyr); //used for getting weighted sel
  matrix F_shrimp(styr,endyr,1,nages);
  vector F_shrimp_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  
////---Per-recruit stuff----------------------------------------------------------------------------------
  vector N_age_spr_F(1,nages);         //numbers at age for SPR calculations: beginning of year, females only
  vector N_age_spr_mdyr_F(1,nages);    //numbers at age for SPR calculations: midyear, females only  
  vector C_age_spr_F(1,nages);         //catch at age for SPR calculations, females only
  vector Z_age_spr_F(1,nages);         //total mortality at age for SPR calculations, females only
  vector F_L_age_spr_F(1,nages);      //fishing mortality (landings, not discards) at age for SPR calculations, females 
  vector N_age_spr_M(1,nages);         //numbers at age for SPR calculations: beginning of year, males only
  vector N_age_spr_mdyr_M(1,nages);    //numbers at age for SPR calculations: midyear, males only  
  vector C_age_spr_M(1,nages);         //catch at age for SPR calculations, males only
  vector Z_age_spr_M(1,nages);         //total mortality at age for SPR calculations, males only
  vector F_L_age_spr_M(1,nages);      //fishing mortality (landings, not discards) at age for SPR calculations, males 
  vector spr_static(styr,endyr);     //vector of static SPR values by year
  vector F_spr(1,n_iter_spr);        //values of full F to be used in per-recruit and equilibrium calculations
  vector spr_spr(1,n_iter_spr);      //reporductive capacity-per-recruit values corresponding to F values in F_spr
  vector L_spr(1,n_iter_spr);        //landings(mt)-per-recruit values corresponding to F values in F_spr
  vector E_spr(1,n_iter_spr);        //exploitation rate values corresponding to F values in F_spr

  vector N_spr_F0(1,nages);          //Used to compute spr at F=0: midpt of year
  vector N_bpr_F0(1,nages);          //Used to compute bpr at F=0: start of year
  number spr_F0;                     //Spawning biomass per recruit at F=0
  vector N_spr_F_init(1,nages);          //Used to compute spr at F=0: start of year 
  vector N_spr_F_init_mdyr(1,nages);          //Used to compute spr at F=0: midpt of year 
  number spr_F_init;                 //Spawning biomass per recruit at F=F_F_hist

  number bpr_F0;                     //Biomass per recruit at F=0 

//-------Objective function components-----------------------------------------------------------------------------
  number w_L;
  number w_D;
  //  number w_shrimp;
  //number w_lc;
  //number w_lc_HL;
  //number w_lc_PN;
  //number w_lc_GN;
  //number w_lc_CN;
  //number w_lc_MRFSS;
  //number w_ac;  
  number w_ac_HL;
  number w_ac_PN;
  number w_ac_GN;
  number w_ac_CN;
  number w_ac_MRFSS;
  number w_I_FL_HL;
  number w_I_MRFSS;
 // number w_I_SMAP_1YR;
  number w_I_SMAP_YOY;
  number w_R;
  number w_R_init;
  number w_R_end;
  number w_F;
  number w_B1dB0;
  number w_fullF;
  number w_cvlen_dev;
  number w_cvlen_diff;

  number f_FL_HL_cpue;     
  number f_MRFSS_cpue;
  number f_SMAP_YOY_cpue;
  //number f_SMAP_1YR_cpue;
    
  number f_HL_L;
  number f_PN_L;
  number f_GN_L;
  number f_CN_L;
  number f_MRFSS_L;
   
  number f_MRFSS_D;
  
  number f_shrimp_B;
   
  //number f_HL_lenc;
  //number f_PN_lenc;
  //number f_GN_lenc;
  //number f_CN_lenc;
  //number f_MRFSS_lenc;

  number f_HL_agec;
  number f_PN_agec;
  number f_GN_agec;
  number f_CN_agec;
  number f_MRFSS_agec;
  
  number f_N_dev;               //weight on recruitment deviations to fit S-R curve
  number f_N_dev_early;         //extra weight against deviations before styr
  number f_N_dev_end;         //extra constraint on last 3 years of recruitment variability
  number f_Fend_constraint;     //penalty for F deviation in last 5 years
  number f_B1dB0_constraint;    //penalty to fix B(styr)/K  
  number f_fullF_constraint;    //penalty for fullF>5
  number f_cvlen_dev_constraint; //deviation penalty on cv's of length at age
  number f_cvlen_diff_constraint;//first diff penalty on cv's of length at age
	number f_priors;								//prior information on parameters
  
  objective_function_value fval;
  number fval_unwgt;

//--Dummy arrays for output convenience  --------------------------
  vector xdum(styr,endyr);
  vector xdum2(styr,endyr+1);
//--Other dummy variables ----
  number sel_diff_dum;
  number zero_dum;
  number dzero_dum; 
  number huge_number;
 // init_number x_dum; //used only during model development. can be removed.
   
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
INITIALIZATION_SECTION


//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
GLOBALS_SECTION
  #include "admodel.h"          // Include AD class definitions   KWS
  #include "admb2r.cpp"    // Include S-compatible output functions (needs preceding)  KWS

//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
RUNTIME_SECTION
 maximum_function_evaluations 1000, 3000, 10000, 20000;
 //maximum_function_evaluations 1, 1, 1;
 convergence_criteria 1e-4, 1e-4, 1e-4;
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PRELIMINARY_CALCS_SECTION
// Set values of fixed parameters or set initial guess of estimated parameters
  Dmort_MRFSS=set_Dmort_MRFSS;

  obs_MRFSS_D=Dmort_MRFSS*obs_MRFSS_released;
  
  obs_MRFSS_L(styr_MRFSS_L,1979)=L_rec_multiply*obs_MRFSS_L(styr_MRFSS_L,1979);
  obs_MRFSS_D(styr_MRFSS_L,1979)=L_rec_multiply*obs_MRFSS_D(styr_MRFSS_L,1979);

  E_age_st=set_E_age_st;   //E computed over ages E_age_st +   [E_age_st-1+ if model starts with age 0]

  //Linf_f=set_Linf_f;
  //Linf_m=set_Linf_m;
  //K_f=set_K_f;
  //K_m=set_K_m;
  //t0_f=set_t0_f;
  //t0_m=set_t0_m;
 
  nlenbins2=nlenbins+20;  //lenbins,nlenbins for plus group; lenbins2,nlenbins2 assume truncation at 90cm
  lenbins2(1)=lenbins(1);
  for(ilen=2;ilen<=nlenbins2;ilen++){
    lenbins2(ilen)=lenbins2(ilen-1)+1;
  }
  prop_f_a0=set_prop_f_a0;
  //selpar_L50_D_2=t0-log(1.0-304.8/Linf)/K; //age at size limit: 304.8 = limit in mm

  M=set_M; 
  M_constant=set_M_constant;
  steep=set_steep;
  //log_dev_N_rec=0.0;
  R_autocorr=set_R_autocorr;
  rec_sigma=set_rec_sigma;
  
  log_q_FL_HL=set_logq_FL_HL;
  log_q_MRFSS=set_logq_MRFSS;
  log_q_SMAP_YOY=set_logq_SMAP_YOY;
  //log_q_SMAP_1YR=set_logq_SMAP_1YR;
	//log_q_shrimp_B=set_logq_shrimp_B;
  //q_rate=set_q_rate;

  //L_early_bias=set_L_early_bias;
  //L_commHAL_bias=1.0;

  w_L=set_w_L;          
  w_D=set_w_D; 
  //w_shrimp=set_w_shrimp;          
  //w_lc_=set_w_lc;         
 //w_lc_HL = set_w_lc_HL;
 //w_lc_PN = set_w_lc_PN;
 //w_lc_GN = set_w_lc_GN;
 //w_lc_CN = set_w_lc_CN;
 //w_lc_MRFSS = set_w_lc_MRFSS;  
  //w_ac=set_w_ac;         
	w_ac_HL = set_w_ac_HL;
  w_ac_PN = set_w_ac_PN;
  w_ac_GN = set_w_ac_GN;
  w_ac_CN = set_w_ac_CN;
  w_ac_MRFSS = set_w_ac_MRFSS;
  w_I_FL_HL=set_w_I_FL_HL;  
  w_I_MRFSS=set_w_I_MRFSS;  
  w_I_SMAP_YOY=set_w_I_SMAP_YOY;  
  //w_I_SMAP_1YR=set_w_I_SMAP_1YR;  

  w_R=set_w_R;          
  w_R_init=set_w_R_init;     
  w_R_end=set_w_R_end;      
  w_F=set_w_F;          
  w_B1dB0=set_w_B1dB0;      
  w_fullF=set_w_fullF;      
  w_cvlen_dev=set_w_cvlen_dev;  
  w_cvlen_diff=set_w_cvlen_diff; 

  F_hist=set_F_hist;
  log_avg_F_HL=set_log_avg_F_HL;
  log_F_dev_HL=set_log_F_dev_HL;
  log_avg_F_PN=set_log_avg_F_PN;
  log_F_dev_PN=set_log_F_dev_PN; 
  log_avg_F_GN=set_log_avg_F_GN;
  log_F_dev_GN=set_log_F_dev_GN;
  log_avg_F_CN=set_log_avg_F_CN;
  log_F_dev_CN=set_log_F_dev_CN;
  log_avg_F_MRFSS=set_log_avg_F_MRFSS;
  log_F_dev_MRFSS=set_log_F_dev_MRFSS;   
  
  log_avg_F_MRFSS_D=set_log_avg_F_MRFSS_D;  
  log_F_dev_MRFSS_D=set_log_F_dev_MRFSS_D;  
   
  log_avg_F_shrimp=set_log_avg_F_shrimp;  
  log_F_dev_shrimp=set_log_F_dev_shrimp;
  //log_obs_shrimp_B_effort=log(obs_shrimp_B_effort);
  //log_obs_shrimp_B=log(obs_shrimp_B(styr_shrimp_B_fit,endyr_shrimp_B_fit));
  //obs_shrimp_geomean=mfexp(sum(log_obs_shrimp_B)/(endyr_shrimp_B_fit-styr_shrimp_B_fit+1.0));
    
   
  log_len_cv=log(set_len_cv);
  log_R0=set_log_R0;

  L50_diff=set_L50_diff; //shift in selectivity b/w males and females (same for all gears)
  sel_historical_F=set_sel_historical_F;
  
  selpar_slope_HL=set_selpar_slope_HL; 
  selpar_L50_HL_keep=set_selpar_L50_HL_keep;
  
  selpar_slope1_PN=set_selpar_slope1_PN; 
  selpar_L501_PN=set_selpar_L501_PN;   
  selpar_slope2_PN=set_selpar_slope2_PN;   
  selpar_L502_PN=set_selpar_L501_PN;      
                                                    
  selpar_slope_GN=set_selpar_slope_GN; 
  selpar_L50_GN_keep=set_selpar_L50_GN_keep;
  
  selpar_slope1_CN=set_selpar_slope1_CN; 
  selpar_L501_CN=set_selpar_L501_CN;     
  selpar_slope2_CN=set_selpar_slope2_CN; 
  selpar_L502_CN=set_selpar_L502_CN;     
  
  selpar_slope1_MRFSS_keep=set_selpar_slope1_MRFSS_keep; 
  selpar_L501_MRFSS_keep=set_selpar_L501_MRFSS_keep;
  selpar_slope2_MRFSS_keep=set_selpar_slope2_MRFSS_keep;
  selpar_L502_MRFSS_keep=set_selpar_L502_MRFSS_keep;  
  //selpar_sigma_MRFSS=set_selpar_sigma_MRFSS;
  //selpar_afull_MRFSS=set_selpar_afull_MRFSS
  
  sel_MRFSS_D_M=set_sel_MRFSS_D_M;  //selectivity for discards - males
  sel_MRFSS_D_F=set_sel_MRFSS_D_F;  //selectivity for discards - females 
  sel_shrimp=set_sel_shrimp_B;  //selectivity for bycatch (same for both sexes)
  
 sqrt2pi=sqrt(2.*3.14159265);
 g2mt=0.000001;         //conversion of grams to metric tons
 g2kg=0.001;         //conversion of grams to kg 
 mt2klb=2.20462;   //converstion of metric tons to 1000 lb 
 //df=0.001; //difference for msy derivative approximations
 zero_dum=0.0;
 
 //additive constant to prevent division by zero
 dzero_dum=0.00001;

 SSB_msy_out=0.0;
 
 maturity_f=maturity_f_obs;
 
////Fill in maturity matrix for calculations for styr to styr
//    for(iyear=styr; iyear<=styr-1; iyear++)
//      {
//         maturity_f(iyear)=maturity_f_obs;
//         maturity_m(iyear)=maturity_m_obs;
//         prop_m(iyear)=prop_m_obs(styr);
//         prop_f(iyear)=1.0-prop_m_obs(styr);
//      }
//    for (iyear=styr;iyear<=endyr;iyear++)
//      {
//         maturity_f(iyear)=maturity_f_obs;
//         maturity_m(iyear)=maturity_m_obs;
//         prop_m(iyear)=prop_m_obs(iyear);
//         prop_f(iyear)=1.0-prop_m_obs(iyear);
//      }

//Fill in sample sizes of comps sampled in nonconsec yrs. 
//Used only for output in R object   
      //nsamp_HL_lenc_allyr=missing;
      //nsamp_CN_lenc_allyr=missing;
      
      nsamp_HL_agec_allyr=missing;
      nsamp_PN_agec_allyr=missing;
      nsamp_CN_agec_allyr=missing;
      nsamp_GN_agec_allyr=missing;
      nsamp_MRFSS_agec_allyr=missing;
      
     // for (iyear=1; iyear<=nyr_HL_lenc; iyear++)
     //    {
     //      nsamp_HL_lenc_allyr(yrs_HL_lenc(iyear))=nsamp_HL_lenc(iyear);
     //    }   
     //for (iyear=1; iyear<=nyr_CN_lenc; iyear++)
      //   {
      //     nsamp_CN_lenc_allyr(yrs_CN_lenc(iyear))=nsamp_CN_lenc(iyear);
      //   } 
      
      for (iyear=1; iyear<=nyr_HL_agec; iyear++)
         {
           nsamp_HL_agec_allyr(yrs_HL_agec(iyear))=nsamp_HL_agec(iyear);
         }  
      for (iyear=1; iyear<=nyr_PN_agec; iyear++)
         {
           nsamp_PN_agec_allyr(yrs_PN_agec(iyear))=nsamp_PN_agec(iyear);
         }  
      for (iyear=1; iyear<=nyr_CN_agec; iyear++)
         {
           nsamp_CN_agec_allyr(yrs_CN_agec(iyear))=nsamp_CN_agec(iyear);
         } 
         //cout<<"nsamp vector CN"<<nsamp_CN_agec_allyr<<endl; 
			for (iyear=1; iyear<=nyr_GN_agec; iyear++)
         {
           nsamp_GN_agec_allyr(yrs_GN_agec(iyear))=nsamp_GN_agec(iyear);
         }  
        // cout<<"nsamp vector GN"<<nsamp_GN_agec_allyr<<endl;
      for (iyear=1; iyear<=nyr_MRFSS_agec; iyear++)
         {
           nsamp_MRFSS_agec_allyr(yrs_MRFSS_agec(iyear))=nsamp_MRFSS_agec(iyear);
         }  
     //  cout<<"nsamp vector MRFSS"<<nsamp_MRFSS_agec_allyr<<endl;   
  age_error_matrix=set_age_error_matrix;

//fill in reproductive caculations with zero's
  reprod.initialize();
  prop_f.initialize();
  prop_f_F0.initialize();

//fill in F's, Catch matrices, and log rec dev with zero's
  F_HL_M.initialize();
  F_HL_F.initialize();
  C_HL_M.initialize();
  C_HL_F.initialize();
  F_PN_M.initialize();
  F_PN_F.initialize();
  C_PN_F.initialize();
  C_PN_M.initialize();  
  F_GN_M.initialize();
  F_GN_F.initialize();
  C_GN_M.initialize();
  C_GN_F.initialize();
  F_CN_M.initialize();
  F_CN_F.initialize();
  C_CN_M.initialize();
  C_CN_F.initialize();
  F_MRFSS_M.initialize();
  F_MRFSS_F.initialize();
  C_MRFSS_M.initialize();
  C_MRFSS_F.initialize();
  F_shrimp.initialize();
  C_shrimp_M.initialize();
  C_shrimp_F.initialize();
  
  C_MRFSS_D_M.initialize();
  C_MRFSS_D_F.initialize();  
  
  L_MRFSS_D.initialize();        //discards in 1000lb whole weight 
  L_shrimp.initialize();        //discards in 1000lb whole weight 
  
  F_MRFSS_D_M.initialize();
  F_MRFSS_D_F.initialize();
      
  log_dev_R.initialize();

//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
TOP_OF_MAIN_SECTION
  arrmblsize=20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(1600);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);
  

//>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PROCEDURE_SECTION
	R0=mfexp(log_R0);
 //		cout<<"start"<<endl;
	get_length_and_weight_at_age();
  //	cout << "got length and weight transitions" <<endl;
	get_selectivity();
//  	cout << "got selectivity" << endl;
	get_length_at_age_dist();
  //	cout<< "got predicted length at age distribution"<<endl;
	get_spr_F0();
  //	cout << "got F0 spr" << endl;
	get_mortality();
  //	cout << "got mortalities" << endl;
	get_bias_corr();
  //	cout<< "got recruitment bias correction" << endl;
	get_numbers_at_age();
  //	cout << "got numbers at age" << endl;  
	get_landings_numbers();
  //	cout << "got catch at age" << endl;
	get_landings_wgt();
  //	cout << "got landings" << endl;
	get_dead_discards();
  //	cout << "got discards" << endl;
	//get_length_comps();
  //	cout<< "got length comps"<< endl;
	get_age_comps();
  //	cout<< "got age comps"<< endl;
	get_sel_weighted_current();
  //	cout<<"got sel weighted"<<endl;
	get_indices();
  //	cout << "got indices" << endl;
 	evaluate_objective_function();
  //	cout << "objective function calculations complete" << endl;


FUNCTION get_length_and_weight_at_age
  //compute mean length (mm) and weight (whole and gutted) at age
    //meanlen_m=Linf_m*(1.0-mfexp(-K_m*(agebins-t0_m)));   //length in mm - males
    //cout<<"mean length at age--males"<<meanlen_m<<endl;
    
    //meanlen_f=Linf_f*(1.0-mfexp(-K_f*(agebins-t0_f)));   //length in mm - females
    //cout<<"mean length at age--females"<<meanlen_f<<endl;
    
    wgt_kg_land=wgtpar_a*pow(meanlen_fishery,wgtpar_b); //wgt in kg in the fishery
    //wgt_kg_land=g2kg*wgt_g_land;										//wgt in kilos in the fishery
    wgt_land=0.001*wgt_kg_land;												//mt of whole wgt:
    wgt_klb_land=mt2klb*wgt_land;										//1000lb of whole wgt
    
    wgt_kg_m=wgtpar_am*pow(meanlen_m,wgtpar_bm);        //wgt in kg 
    //wgt_kg_m=g2kg*wgt_g_m;                           //wgt in kilos 
    wgt_m=0.001*wgt_kg_m;     													//mt of whole wgt: g2mt converts g to mt
    wgt_klb_m=mt2klb*wgt_m;                          //1000 lb of whole wgt
    
    wgt_kg_f=wgtpar_af*pow(meanlen_f,wgtpar_bf);        //wgt in grams
    //wgt_kg_f=g2kg*wgt_g_f;                           //wgt in kilos 
    wgt_f=0.001*wgt_kg_f;     												//mt of whole wgt: g2mt converts g to mt
    wgt_klb_f=mt2klb*wgt_f;   											//1000s lb

FUNCTION get_selectivity
// ------- compute landings and discard selectivities 
	//	for (iyear=styr; iyear<=endyr; iyear++)
		//{
			//	sel_HL_keep_M(iyear)=logistic(agebins, selpar_L50_HL, selpar_slope_HL);
				//sel_HL_keep_M(iyear)=logistic(agebins, selpar_L50_HL, selpar_slope_HL);
				//sel_HL_keep_M(iyear)=logistic(agebins, selpar_L50_HL, selpar_slope_HL);
				
				
  for (iage=1; iage<=nages; iage++)
  {   
  	 sel_HL_keep_M(iage)=1./(1.+mfexp(-1.*selpar_slope_HL*(double(agebins(iage))-(selpar_L50_HL_keep-L50_diff)))); 
     sel_HL_keep_F(iage)=1./(1.+mfexp(-1.*selpar_slope_HL*(double(agebins(iage))-selpar_L50_HL_keep))); 
     //sel_PN_M(iage)=1./(1.+mfexp(-1.*selpar_slope_PN*(double(agebins(iage))-selpar_L50_PN-L50_diff))); 
     //sel_PN_F(iage)=1./(1.+mfexp(-1.*selpar_slope_PN*(double(agebins(iage))-selpar_L50_PN)));     
     sel_PN_M(iage)=(1./(1.+mfexp(-1.*selpar_slope1_PN*(double(agebins(iage))-
                       (selpar_L501_PN-L50_diff)))))*(1-(1./(1.+mfexp(-1.*selpar_slope2_PN*
                       (double(agebins(iage))-((selpar_L501_PN-L50_diff)+selpar_L502_PN)))))); //double logistic 
     sel_PN_F(iage)=(1./(1.+mfexp(-1.*selpar_slope1_PN*(double(agebins(iage))-
                       selpar_L501_PN))))*(1-(1./(1.+mfexp(-1.*selpar_slope2_PN*
                       (double(agebins(iage))-(selpar_L501_PN+selpar_L502_PN)))))); //double logistic 
     sel_GN_keep_M(iage)=1./(1.+mfexp(-1.*selpar_slope_GN*(double(agebins(iage))-(selpar_L50_GN_keep-L50_diff)))); 
     sel_GN_keep_F(iage)=1./(1.+mfexp(-1.*selpar_slope_GN*(double(agebins(iage))-selpar_L50_GN_keep))); 
     //sel_GN_keep_M2(iage)=(1./(1.+mfexp(-1.*selpar_slope_GN2*(double(agebins(iage))-
     //                  selpar_L50_GN2))))*(1-(1./(1.+mfexp(-1.*selpar_slope2_GN2*
     //                  (double(agebins(iage))-(selpar_L50_GN2+selpar_L502_GN2)))))); //double logistic 
     
           
     //sel_CN_M(iage)=1./(1.+mfexp(-1.*selpar_slope_CN*(double(agebins(iage))-selpar_L50_CN-L50_diff))); 
     //sel_CN_F(iage)=1./(1.+mfexp(-1.*selpar_slope_CN*(double(agebins(iage))-selpar_L50_CN)));    
   } 
   //sel_CN_M=logistic_exponential(agebins, selpar_L50_CN, selpar_slope_CN, selpar_sigma_CN, selpar_afull_CN); 
   //sel_CN_F=logistic_exponential(agebins, selpar_L50_CN, selpar_slope_CN, selpar_sigma_CN, selpar_afull_CN); 
   	//--Implement Butterworth's Logistic exponential selex to avoid having to take a max(non-differentiable)
    
   	
   	sel_CN_M=logistic_double(agebins, (selpar_L501_CN-L50_diff), selpar_slope1_CN, selpar_L502_CN, selpar_slope2_CN);
   	sel_CN_F=logistic_double(agebins, selpar_L501_CN, selpar_slope1_CN, selpar_L502_CN, selpar_slope2_CN);  
    sel_MRFSS_keep_M=logistic_double(agebins, (selpar_L501_MRFSS_keep-L50_diff), selpar_slope1_MRFSS_keep, selpar_L502_MRFSS_keep, selpar_slope2_MRFSS_keep);
    sel_MRFSS_keep_F=logistic_double(agebins, selpar_L501_MRFSS_keep, selpar_slope1_MRFSS_keep, selpar_L502_MRFSS_keep, selpar_slope2_MRFSS_keep);
    //sel_MRFSS_keep_M(1)=0.1;
    
     //sel_MRFSS_keep_M(iage)=1./(1.+mfexp(-1.*selpar_slope_MRFSS*(double(agebins(iage))-selpar_L50_MRFSS_keep-L50_diff))); 
     //sel_MRFSS_keep_F(iage)=1./(1.+mfexp(-1.*selpar_slope_MRFSS*(double(agebins(iage))-selpar_L50_MRFSS_keep)));  
   //sel_MRFSS_keep_M=logistic_exponential(agebins, selpar_L501_MRFSS_keep, selpar_slope1_MRFSS_keep, selpar_sigma_MRFSS, selpar_afull_MRFSS);        
   //sel_MRFSS_keep_F=logistic_exponential(agebins, selpar_L501_MRFSS_keep, selpar_slope1_MRFSS_keep, selpar_sigma_MRFSS, selpar_afull_MRFSS);        
   //sel_MRFSS_keep_M=sel_MRFSS_keep_M/max(sel_MRFSS_keep_M);
  //sel_GN_keep_M2=sel_GN_keep_M2/max(sel_GN_keep_M2); //renormalize
  sel_PN_M=sel_PN_M/max(sel_PN_M);
  sel_PN_F=sel_PN_F/max(sel_PN_F);
  //sel_GN_keep_F2=sel_GN_keep_M2;   
 // cout <<"HL male selex"<<sel_HL_keep_M<< endl;
 // cout <<"PN male selex"<<sel_PN_M<< endl;
 // cout <<"GN male selex"<<sel_GN_keep_M<< endl;
 // cout <<"CN male selex"<<sel_CN_M<< endl;
 // cout <<"MRFSS male selex"<<sel_MRFSS_keep_M<< endl;

FUNCTION get_length_at_age_dist
  //compute matrix of length at age, based on the normal distribution
  for (iage=1;iage<=nages;iage++)
  {
    //len_cv(iage)=mfexp(log_len_cv+log_len_cv_dev(iage));
    len_cv(iage)=mfexp(log_len_cv);
    for (ilen=1;ilen<=nlenbins2;ilen++)
    {
      //convert to centimeters for matching
      lenprob_m2(iage,ilen)=(mfexp(-(square(lenbins2(ilen)-0.1*meanlen_m(iage))/
      (2.*square(len_cv(iage)*0.1*meanlen_m(iage)))))/(sqrt2pi*0.1*len_cv(iage)*meanlen_m(iage))); 
      
      lenprob_f2(iage,ilen)=(mfexp(-(square(lenbins2(ilen)-0.1*meanlen_f(iage))/
      (2.*square(len_cv(iage)*0.1*meanlen_f(iage)))))/(sqrt2pi*0.1*len_cv(iage)*meanlen_f(iage)));
    }
    lenprob_m2(iage)/=sum(lenprob_m2(iage)); //standardize to account for truncated normal (i.e., no sizes<0)
    lenprob_f2(iage)/=sum(lenprob_f2(iage)); //standardize to account for truncated normal (i.e., no sizes<0)
    for(ilen=1;ilen<=nlenbins;ilen++){
      lenprob_m(iage,ilen)=lenprob_m2(iage,ilen);
      lenprob_f(iage,ilen)=lenprob_f2(iage,ilen);
    }
    for(ilen=nlenbins+1;ilen<=nlenbins2;ilen++){
      lenprob_m(iage,nlenbins)+=lenprob_m2(iage,ilen);
      lenprob_f(iage,nlenbins)+=lenprob_f2(iage,ilen);
    }
  }

FUNCTION get_spr_F0
  reprod=elem_prod(maturity_f,wgt_f);  
  //at mdyr, apply half this yr's mortality, half next yr's
  N_spr_F0(1)=prop_f_a0*mfexp(-1.0*M(1)/2.0);  //at midpt of year
  N_bpr_F0(1)=1.;
  for (iage=2; iage<=nages; iage++)
  {
    //N_spr_F0(iage)=N_spr_F0(iage-1)*mfexp(-1.0*(M(iage-1));
    N_spr_F0(iage)=N_spr_F0(iage-1)*mfexp(-1.0*(M(iage-1)/2.0 + M(iage)/2.0)); 
    N_bpr_F0(iage)=N_bpr_F0(iage-1)*mfexp(-1.0*(M(iage-1)));    
  }
  N_spr_F0(nages)=N_spr_F0(nages)/(1.0-mfexp(-1.0*M(nages))); //plus group (sum of geometric series)
  N_bpr_F0(nages)=N_bpr_F0(nages)/(1.0-mfexp(-1.0*M(nages)));
  
  spr_F0=sum(elem_prod(N_spr_F0,reprod));
  bpr_F0=prop_f_a0*sum(elem_prod(N_bpr_F0,wgt_f))+(1.-prop_f_a0)*sum(elem_prod(N_bpr_F0,wgt_m)); 
  B0=R0*bpr_F0;   

FUNCTION get_mortality
  fullF=0.0;

  for (iyear=styr; iyear<=endyr; iyear++) {
    if(iyear<styr_HL_L){
      F_HL_out(iyear)=0.0; 
    }
    else{
      F_HL_out(iyear)=mfexp(log_avg_F_HL+log_F_dev_HL(iyear));
    }
    F_HL_F(iyear)=sel_HL_keep_F*F_HL_out(iyear);
    F_HL_M(iyear)=sel_HL_keep_M*F_HL_out(iyear);
    fullF(iyear)+=F_HL_out(iyear);
    
    if(iyear<styr_PN_L){
      F_PN_out(iyear)=0.0;
    }
    else{
      F_PN_out(iyear)=mfexp(log_avg_F_PN+log_F_dev_PN(iyear));
    }
    F_PN_F(iyear)=sel_PN_F*F_PN_out(iyear);
    F_PN_M(iyear)=sel_PN_M*F_PN_out(iyear);
    fullF(iyear)+=F_PN_out(iyear);
    
    if(iyear<styr_GN_L){
      F_GN_out(iyear)=0.0;  
    }
    else {
      F_GN_out(iyear)=mfexp(log_avg_F_GN+log_F_dev_GN(iyear));
    }
       
      F_GN_F(iyear)=sel_GN_keep_F*F_GN_out(iyear);
      F_GN_M(iyear)=sel_GN_keep_M*F_GN_out(iyear);
      fullF(iyear)+=F_GN_out(iyear);
    
    if(iyear<styr_CN_L){
      F_CN_out(iyear)=0.0;
    }
    else{
      F_CN_out(iyear)=mfexp(log_avg_F_CN+log_F_dev_CN(iyear));
    }
    F_CN_F(iyear)=sel_CN_F*F_CN_out(iyear);
    F_CN_M(iyear)=sel_CN_M*F_CN_out(iyear);
    fullF(iyear)+=F_CN_out(iyear);
    
    if(iyear<styr_MRFSS_L){
      F_MRFSS_out(iyear)=0.0;
    }
    else{
      F_MRFSS_out(iyear)=mfexp(log_avg_F_MRFSS+log_F_dev_MRFSS(iyear));
    }
    F_MRFSS_F(iyear)=sel_MRFSS_keep_F*F_MRFSS_out(iyear);
    F_MRFSS_M(iyear)=sel_MRFSS_keep_M*F_MRFSS_out(iyear);
    fullF(iyear)+=F_MRFSS_out(iyear);
    //cout<<"Rec F"<<F_MRFSS_M<<endl;
       
    //discards
    
    if(iyear<styr_MRFSS_D){
      F_MRFSS_D_out(iyear)=0.0;
    }
    else{
      F_MRFSS_D_out(iyear)=mfexp(log_avg_F_MRFSS_D+log_F_dev_MRFSS_D(iyear));
    }
    fullF(iyear)+=F_MRFSS_D_out(iyear);
    F_MRFSS_D_M(iyear)=sel_MRFSS_D_M*F_MRFSS_D_out(iyear);  
    F_MRFSS_D_F(iyear)=sel_MRFSS_D_F*F_MRFSS_D_out(iyear);
  
    if(iyear<styr_shrimp_B){
      F_shrimp_out(iyear)=0.0;
    }
    else{
      F_shrimp_out(iyear)=mfexp(log_avg_F_shrimp+log_F_dev_shrimp(iyear));
    //F_shrimp_out(iyear)=mfexp(log_q_shrimp_B+log_obs_shrimp_B_effort(iyear));
    }
    fullF(iyear)+=F_shrimp_out(iyear);
    F_shrimp(iyear)=sel_shrimp*F_shrimp_out(iyear);  
  }

     
  for (iyear=styr; iyear<=endyr; iyear++) {
    F_F(iyear)=F_HL_F(iyear);
    F_F(iyear)+=F_PN_F(iyear);
    F_F(iyear)+=F_GN_F(iyear);
    F_F(iyear)+=F_CN_F(iyear);
    F_F(iyear)+=F_MRFSS_F(iyear);
    F_F(iyear)+=F_MRFSS_D_F(iyear);
    F_F(iyear)+=F_shrimp(iyear);
    Z_F(iyear)=M+F_F(iyear);
    
    F_M(iyear)=F_HL_M(iyear);
    F_M(iyear)+=F_PN_M(iyear);
    F_M(iyear)+=F_GN_M(iyear);
    F_M(iyear)+=F_CN_M(iyear);
    F_M(iyear)+=F_MRFSS_M(iyear);
    F_M(iyear)+=F_MRFSS_D_M(iyear);
    F_M(iyear)+=F_shrimp(iyear);
    Z_M(iyear)=M+F_M(iyear);    
  }
  F_F_hist=sel_historical_F*F_hist;
  Z_hist=F_F_hist+M;
  
FUNCTION get_bias_corr
  //var_rec_dev=norm2(log_dev_N_rec(styr_rec_dev,(endyr-2))-sum(log_dev_N_rec(styr_rec_dev,(endyr-2)))
  //            /(nyrs_rec-2.0))/(nyrs_rec-3.0); //sample variance from yrs styr_rec_dev-2005 
  //if (set_BiasCor <= 0.0) {BiasCor=mfexp(var_rec_dev/2.0);}   //bias correction            
  //else {BiasCor=set_BiasCor;}
  	
  	 //may exclude last BiasCor_exclude_yrs yrs bc constrained or lack info to estimate
  //var_rec_dev=norm2(log_rec_dev(styr_rec_dev,(endyr-BiasCor_exclude_yrs))-
  //            sum(log_rec_dev(styr_rec_dev,(endyr-BiasCor_exclude_yrs)))
  //            /(nyrs_rec-BiasCor_exclude_yrs))/(nyrs_rec-BiasCor_exclude_yrs-1.0);
  var_rec_dev=norm2(log_rec_dev(styr_rec_dev,(endyr-2))-sum(log_rec_dev(styr_rec_dev,(endyr-2)))
              /(nyrs_rec-2.0))/(nyrs_rec-3.0);  //sample variance from yrs styr_rec_dev - 2009
                           
  //if (set_BiasCor <= 0.0) {BiasCor=mfexp(var_rec_dev/2.0);}   //bias correction
  rec_sigma_sq=square(rec_sigma);
  if (set_BiasCor <= 0.0) {BiasCor=mfexp(rec_sigma_sq/2.0);}   //bias correction               
  else {BiasCor=set_BiasCor;}

FUNCTION get_numbers_at_age
//Initial age
  S0=spr_F0*R0;
  
  //Assume equilibrium age structure for first year dependent on historical F
  N_spr_F_init(1)=prop_f_a0;
  N_spr_F_init_mdyr(1)=N_spr_F_init(1)*mfexp((-1.*(Z_hist(1)))/2.0);
  for (iage=2; iage<=nages; iage++)
  {
    N_spr_F_init(iage)=N_spr_F_init(iage-1)*mfexp(-1.*(Z_hist(iage-1)));
    N_spr_F_init_mdyr(iage)=N_spr_F_init(iage)*mfexp((-1.*(Z_hist(iage)))/2.0);
  }
  N_spr_F_init(nages)=N_spr_F_init(nages)/(1.0-mfexp(-1.*(Z_hist(nages))));
  N_spr_F_init_mdyr(nages)=N_spr_F_init(nages)*mfexp((-1.*(Z_hist(nages)))/2.0); 
  spr_F_init=sum(elem_prod(N_spr_F_init_mdyr,reprod));

  R1=(R0/((5.0*steep-1.0)*spr_F_init))*(BiasCor*4.0*steep*spr_F_init-spr_F0*(1.0-steep));
  //R1=(R0/((5.0*steep-1.0)*spr_F_init))*(4.0*steep*spr_F_init-spr_F0*(1.0-steep));  //take out bias correction in first year because recruitment deterministic at the beginning
  if(R1<0.0)
    {
     R1=1.0;
     } 
  // cout<<"R1"<<R1<<endl;
  N_M(styr,1)=(1.-prop_f_a0)*R1;  
  N_F(styr,1)=prop_f_a0*R1;             
  N_M_mdyr(styr,1)=N_M(styr,1)*mfexp(-1.*Z_M(styr,1)/2.0);
  N_F_mdyr(styr,1)=N_F(styr,1)*mfexp(-1.*Z_F(styr,1)/2.0);
  for (iage=2; iage<=nages; iage++)
  {
    N_M(styr,iage)=N_M(styr,iage-1)*mfexp(-1.*Z_hist(iage-1));
    N_M_mdyr(styr,iage)=N_M(styr,iage)*mfexp(-1.*Z_M(styr,iage)/2.0);
    N_F(styr,iage)=N_F(styr,iage-1)*mfexp(-1.*Z_hist(iage-1));
    N_F_mdyr(styr,iage)=N_F(styr,iage)*mfexp(-1.*Z_F(styr,iage)/2.0);
  }
  //plus group calculation
  N_M(styr,nages)=N_M(styr,nages)/(1.-mfexp(-1.*Z_hist(nages)));
  N_M_mdyr(styr,nages)=N_M(styr,nages)*mfexp(-1.*Z_M(styr,nages)/2.0);  
  N_F(styr,nages)=N_F(styr,nages)/(1.-mfexp(-1.*Z_hist(nages)));
  N_F_mdyr(styr,nages)=N_F(styr,nages)*mfexp(-1.*Z_F(styr,nages)/2.0); 
  
  SSB(styr)=sum(elem_prod(N_F_mdyr(styr),reprod));
  B(styr)=elem_prod(N_M(styr),wgt_m)+elem_prod(N_F(styr),wgt_f);
  totB(styr)=sum(B(styr));
    
//Rest of years 
  for (iyear=styr; iyear<endyr; iyear++)
  {
    
    N_F(iyear+1)(2,nages)=++elem_prod(N_F(iyear)(1,nages-1),(mfexp(-1.*Z_F(iyear)(1,nages-1))));
    N_F(iyear+1,nages)+=N_F(iyear,nages)*mfexp(-1.*Z_F(iyear,nages));//plus group
    N_F_mdyr(iyear+1)(2,nages)=elem_prod(N_F(iyear+1)(2,nages),(mfexp(-1.*(Z_F(iyear+1)(2,nages))/2.0))); //mdyr  
    N_F_mdyr(iyear+1,1)=0.0;       
    N_M(iyear+1)(2,nages)=++elem_prod(N_M(iyear)(1,nages-1),(mfexp(-1.*Z_M(iyear)(1,nages-1))));
    N_M(iyear+1,nages)+=N_M(iyear,nages)*mfexp(-1.*Z_M(iyear,nages));//plus group
    N_M_mdyr(iyear+1)(2,nages)=elem_prod(N_M(iyear+1)(2,nages),(mfexp(-1.*(Z_M(iyear+1)(2,nages))/2.0))); //mdyr 
    N_M_mdyr(iyear+1,1)=0.0;   
    SSB(iyear+1)=sum(elem_prod(N_F_mdyr(iyear+1),reprod));
    
    for(iage=1;iage<=nages;iage++){   
      prop_f(iyear,iage)=N_F(iyear,iage)/(N_F(iyear,iage)+N_M(iyear,iage)+dzero_dum);
    }
    
    if(iyear<(styr_rec_dev-1)) //recruitment follows S-R curve exactly
    {
        //add 0.00001 to avoid log(zero)
        N_F(iyear+1,1)=prop_f_a0*BiasCor*mfexp(log(((0.8*R0*steep*SSB(iyear+1))/(0.2*R0*spr_F0*      
            (1.0-steep)+(steep-0.2)*SSB(iyear+1)))+dzero_dum));
    }
    else{
        N_F(iyear+1,1)=prop_f_a0*mfexp(log(((0.8*R0*steep*SSB(iyear+1))/(0.2*R0*spr_F0*
            (1.0-steep)+(steep-0.2)*SSB(iyear+1)))+dzero_dum)+log_rec_dev(iyear+1)); 
      // cout << "reading out Fem abundance" <<N_F(endyr+1,1) << endl;        
    }
  
   
    N_M(iyear+1,1)=N_F(iyear+1,1)*(1./prop_f_a0-1.);  //this is how it works out algebraically
    N_M_mdyr(iyear+1,1)=N_M(iyear+1,1)*mfexp(-1.*(Z_M(iyear+1,1)/2));
    N_F_mdyr(iyear+1,1)=N_F(iyear+1,1)*mfexp(-1.*(Z_F(iyear+1,1)/2)); 
    B(iyear+1)=elem_prod(N_F(iyear+1),wgt_f)+elem_prod(N_M(iyear+1),wgt_m);
    totB(iyear+1)=sum(B(iyear+1));
  }
    
    //last year (projection) has no recruitment variability
    N_M(endyr+1)(2,nages)=++elem_prod(N_M(endyr)(1,nages-1),(mfexp(-1.*Z_M(endyr)(1,nages-1))));
    N_M(endyr+1,nages)+=N_M(endyr,nages)*mfexp(-1.*Z_M(endyr,nages));//plus group
    N_M_mdyr(endyr+1)(1,nages)=elem_prod(N_M(endyr+1)(1,nages),(mfexp(-1.*(Z_M(endyr)(1,nages))/2.0))); //mdyr 
    N_F(endyr+1)(2,nages)=++elem_prod(N_F(endyr)(1,nages-1),(mfexp(-1.*Z_F(endyr)(1,nages-1))));
    N_F(endyr+1,nages)+=N_F(endyr,nages)*mfexp(-1.*Z_F(endyr,nages));//plus group      
    N_F_mdyr(endyr+1)(1,nages)=elem_prod(N_F(endyr+1)(1,nages),(mfexp(-1.*(Z_F(endyr)(1,nages))/2.0))); //mdyr            
    SSB_extra=sum(elem_prod(N_F_mdyr(endyr+1),reprod));
    N_F(endyr+1,1)=prop_f_a0*mfexp(log(((0.8*R0*steep*SSB_extra)/(0.2*R0*spr_F0*
                 (1.0-steep)+(steep-0.2)*SSB_extra))+dzero_dum));
    N_M(endyr+1,1)=(1.-prop_f_a0)*mfexp(log(((0.8*R0*steep*SSB_extra)/(0.2*R0*spr_F0*
                 (1.0-steep)+(steep-0.2)*SSB_extra))+dzero_dum));
    B(endyr+1)=elem_prod(N_M(endyr+1),wgt_m)+elem_prod(N_F(endyr+1),wgt_f);
    totB(endyr+1)=sum(B(endyr+1));    

//Recruitment time series
  rec=column(N_F,1)+column(N_M,1);

//Benchmark parameters
  S1S0=SSB(styr)/S0;
  popstatus=SSB(endyr)/S0;
  
FUNCTION get_landings_numbers //Baranov catch eqn
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      C_HL_M(iyear,iage)=N_M(iyear,iage)*F_HL_M(iyear,iage)*
        (1.-mfexp(-1.*Z_M(iyear,iage)))/Z_M(iyear,iage);
     //  cout << "Handline catches, Male" <<C_HL_M << endl;
      C_HL_F(iyear,iage)=N_F(iyear,iage)*F_HL_F(iyear,iage)*
        (1.-mfexp(-1.*Z_F(iyear,iage)))/Z_F(iyear,iage);
      C_PN_M(iyear,iage)=N_M(iyear,iage)*F_PN_M(iyear,iage)*
        (1.-mfexp(-1.*Z_M(iyear,iage)))/Z_M(iyear,iage);
      C_PN_F(iyear,iage)=N_F(iyear,iage)*F_PN_F(iyear,iage)*
        (1.-mfexp(-1.*Z_F(iyear,iage)))/Z_F(iyear,iage);
      C_GN_M(iyear,iage)=N_M(iyear,iage)*F_GN_M(iyear,iage)*
        (1.-mfexp(-1.*Z_M(iyear,iage)))/Z_M(iyear,iage);
      C_GN_F(iyear,iage)=N_F(iyear,iage)*F_GN_F(iyear,iage)*
        (1.-mfexp(-1.*Z_F(iyear,iage)))/Z_F(iyear,iage);
      //  cout << "C_GN_F" <<C_GN_F << endl;
      C_CN_M(iyear,iage)=N_M(iyear,iage)*F_CN_M(iyear,iage)*
        (1.-mfexp(-1.*Z_M(iyear,iage)))/Z_M(iyear,iage);
      //  cout << "C_GN_F" <<C_GN_F << endl;
      C_CN_F(iyear,iage)=N_F(iyear,iage)*F_CN_F(iyear,iage)*
        (1.-mfexp(-1.*Z_F(iyear,iage)))/Z_F(iyear,iage);
      C_MRFSS_M(iyear,iage)=N_M(iyear,iage)*F_MRFSS_M(iyear,iage)*
        (1.-mfexp(-1.*Z_M(iyear,iage)))/Z_M(iyear,iage);
      //  cout << "MRFSS catches, Male" <<C_MRFSS_M << endl;
      C_MRFSS_F(iyear,iage)=N_F(iyear,iage)*F_MRFSS_F(iyear,iage)*
        (1.-mfexp(-1.*Z_F(iyear,iage)))/Z_F(iyear,iage);
      //    cout << "MRFSS catches, Female" <<C_MRFSS_F << endl;
      C_shrimp_M(iyear,iage)=N_M(iyear,iage)*F_shrimp(iyear,iage)*
        (1.-mfexp(-1.*Z_M(iyear,iage)))/Z_M(iyear,iage);
      //    cout << "Shrimp bycatch, Male" <<C_shrimp_M << endl;
      C_shrimp_F(iyear,iage)=N_F(iyear,iage)*F_shrimp(iyear,iage)*
        (1.-mfexp(-1.*Z_F(iyear,iage)))/Z_F(iyear,iage);
    }
  }

FUNCTION get_landings_wgt
//---Predicted landings------------------------
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    L_HL(iyear)=elem_prod(C_HL_M(iyear),wgt_land)+elem_prod(C_HL_F(iyear),wgt_land);           //in mt
    L_HL_klb(iyear)=elem_prod(C_HL_F(iyear),wgt_klb_land)+elem_prod(C_HL_M(iyear),wgt_klb_land); //in 1000 lb
 
    L_PN(iyear)=elem_prod(C_PN_M(iyear),wgt_land)+elem_prod(C_PN_F(iyear),wgt_land);   //in mt
    L_PN_klb(iyear)=elem_prod(C_PN_F(iyear),wgt_klb_land)+elem_prod(C_PN_M(iyear),wgt_klb_land); //in 1000 lb
    
    L_GN(iyear)=elem_prod(C_GN_M(iyear),wgt_land)+elem_prod(C_GN_F(iyear),wgt_land); //in mt
    L_GN_klb(iyear)=elem_prod(C_GN_F(iyear),wgt_klb_land)+elem_prod(C_GN_M(iyear),wgt_klb_land); //in 1000 lb
    
    L_CN(iyear)=elem_prod(C_CN_M(iyear),wgt_land)+elem_prod(C_CN_F(iyear),wgt_land); //in mt
    L_CN_klb(iyear)=elem_prod(C_CN_F(iyear),wgt_klb_land)+elem_prod(C_CN_M(iyear),wgt_klb_land); //in 1000 lb
    
    L_MRFSS(iyear)=elem_prod(C_MRFSS_M(iyear),wgt_land)+elem_prod(C_MRFSS_F(iyear),wgt_land); //in mt
    L_MRFSS_klb(iyear)=elem_prod(C_MRFSS_F(iyear),wgt_klb_land)+elem_prod(C_MRFSS_M(iyear),wgt_klb_land); //in 1000 lb
    
    L_shrimp(iyear)=elem_prod(C_shrimp_F(iyear),wgt_land)+elem_prod(C_shrimp_M(iyear),wgt_land); //in mt
  }
 // for (iyear=styr_MRFSS_LL; iyear<+endyr_MRFSS_LL; iyear++)
 // {   L_MRFSS(iyear)=elem_prod(C_MRFSS_M(iyear),wgt_m)+elem_prod(C_MRFSS_F(iyear),wgt_f); //in mt
 //     L_MRFSS_klb(iyear)=elem_prod(C_MRFSS_F(iyear),wgt_klb_f)+elem_prod(C_MRFSS_M(iyear),wgt_klb_m); //in 1000 lb
 // }
 // for (iyear=styr_shrimp_B; iyear<+endyr_shrimp_B; iyear++) 
 //     L_shrimp(iyear)=elem_prod(C_shrimp_F(iyear),wgt_f)+elem_prod(C_shrimp_M(iyear),wgt_m); //in mt
 // }
  for (iyear=styr_HL_L; iyear<=endyr_HL_L; iyear++){
     pred_HL_L(iyear)=sum(L_HL_klb(iyear));
  }
  for (iyear=styr_PN_L; iyear<=endyr_PN_L; iyear++){
     pred_PN_L(iyear)=sum(L_PN_klb(iyear));
  }
  for (iyear=styr_GN_L; iyear<=endyr_GN_L; iyear++){
     pred_GN_L(iyear)=sum(L_GN_klb(iyear));
  }
  for (iyear=styr_CN_L; iyear<=endyr_CN_L; iyear++){
     pred_CN_L(iyear)=sum(L_CN_klb(iyear));
  }
  for(iyear=styr_MRFSS_L;iyear<=endyr_MRFSS_L;iyear++){
     pred_MRFSS_L(iyear)=sum(C_MRFSS_F(iyear)+C_MRFSS_M(iyear))/1000.0; //in 1000's (numbers)
  }  
  for(iyear=styr_shrimp_B;iyear<=endyr_shrimp_B;iyear++){
     pred_shrimp_B(iyear)=sum(C_shrimp_F(iyear)+C_shrimp_M(iyear))/1000.0; //in 1000's (numbers)
  }    
  //log_pred_shrimp_B=log(pred_shrimp_B(styr_shrimp_B_fit,endyr_shrimp_B_fit));
  //pred_shrimp_geomean=mfexp(sum(log_pred_shrimp_B)/(endyr_shrimp_B_fit-styr_shrimp_B_fit+1.0));

FUNCTION get_dead_discards //Baranov catch eqn  
  //dead discards at age (number fish) 

  for (iyear=styr_MRFSS_D; iyear<=endyr_MRFSS_D; iyear++){
    for (iage=1; iage<=nages; iage++){
      C_MRFSS_D_F(iyear,iage)=N_F(iyear,iage)*F_MRFSS_D_F(iyear,iage)*
        (1.-mfexp(-1.*Z_F(iyear,iage)))/Z_F(iyear,iage);
      C_MRFSS_D_M(iyear,iage)=N_M(iyear,iage)*F_MRFSS_D_M(iyear,iage)*
        (1.-mfexp(-1.*Z_M(iyear,iage)))/Z_M(iyear,iage);
    }
    L_MRFSS_D(iyear)=elem_prod(C_MRFSS_D_F(iyear),wgt_land)+elem_prod(C_MRFSS_D_M(iyear),wgt_land);        //discards in 1000lb whole weight     
    pred_MRFSS_D(iyear)=(sum(C_MRFSS_D_M(iyear))+sum(C_MRFSS_D_F(iyear)))/1000.0; //pred annual dead discards in 1000s
  }
       
FUNCTION get_indices
 //FL hand lines cpue
  for (iyear=styr_FL_HL_cpue; iyear<=endyr_FL_HL_cpue; iyear++)
  {   //index in whole wgt (lb) units, wgt_klb in 1000 lb, but the multiplier (1000) is absorbed by q
      N_FL_HL(iyear)=elem_prod(elem_prod(N_F_mdyr(iyear),sel_wgted_tot_F),wgt_klb_f)+
       elem_prod(elem_prod(N_M_mdyr(iyear),sel_wgted_tot_F),wgt_klb_m); 
      pred_FL_HL_cpue(iyear)=mfexp(log_q_FL_HL)*sum(N_FL_HL(iyear));
  }
  //MRFSS cpue  //includes discards
  for (iyear=styr_MRFSS_cpue; iyear<=endyr_MRFSS_cpue; iyear++)
  {   //index in 1000's (numbers)
      N_MRFSS(iyear)=elem_prod(N_F_mdyr(iyear),sel_MRFSS_keep_F)+elem_prod(N_F_mdyr(iyear),sel_MRFSS_D_F)+
       elem_prod(N_M_mdyr(iyear),sel_MRFSS_keep_M)+elem_prod(N_M_mdyr(iyear),sel_MRFSS_D_M); 
      pred_MRFSS_cpue(iyear)=mfexp(log_q_MRFSS)*sum(N_MRFSS(iyear));
  }

   //SEAMAP YOY cpue  
  for (iyear=styr_SMAP_YOY_cpue; iyear<=endyr_SMAP_YOY_cpue; iyear++)
  {   //index in 1000's (numbers)
      pred_SMAP_YOY_cpue(iyear)=mfexp(log_q_SMAP_YOY)*(N_F(iyear,1)+N_M(iyear,1));
  }
  //SEAMAP 1YR cpue
  //for (iyear=styr_SMAP_1YR_cpue; iyear<=endyr_SMAP_1YR_cpue; iyear++)
  //{   //index in 1000's (numbers)
  //    pred_SMAP_1YR_cpue(iyear)=mfexp(log_q_SMAP_1YR)*(N_F(iyear,2)+N_M(iyear,2));
  //}
  
//FUNCTION get_length_comps
//Hand lines
//  for (iyear=1;iyear<=nyr_HL_lenc;iyear++)
// {
//    pred_HL_lenc(iyear)=(C_HL_F(yrs_HL_lenc(iyear))*lenprob_f+C_HL_M(yrs_HL_lenc(iyear))*lenprob_m)/(sum(C_HL_F(yrs_HL_lenc(iyear))+C_HL_M(yrs_HL_lenc(iyear)))+dzero_dum);
//  } 
////Pound Nets
//  for (iyear=1;iyear<=nyr_PN_lenc;iyear++)
//  {
//    iyear2=yrs_PN_lenc(iyear);
//    pred_PN_lenc(iyear)=(C_PN_F(iyear2)*lenprob_f+C_PN_M(iyear2)*lenprob_m)/sum(C_PN_F(iyear2)+C_PN_M(iyear2));
//  }
////Gillnet
//  for (iyear=1;iyear<=nyr_GN_lenc;iyear++)
//  {
//    iyear2=yrs_GN_lenc(iyear);
//    pred_GN_lenc(iyear)=(C_GN_F(iyear2)*lenprob_f+C_GN_M(iyear2)*lenprob_m)/sum(C_GN_F(iyear2)+C_GN_M(iyear2));
//  }
//  
////Castnets
//  for (iyear=1;iyear<=nyr_CN_lenc;iyear++)
//  {
//    iyear2=yrs_CN_lenc(iyear);
//    pred_CN_lenc(iyear)=(C_CN_F(iyear2)*lenprob_f+C_CN_M(iyear2)*lenprob_m)/sum(C_CN_F(iyear2)+C_CN_M(iyear2));
//  }
////MRFSS
//  for (iyear=styr_MRFSS_lenc;iyear<=endyr_MRFSS_lenc;iyear++)
//  {
//    pred_MRFSS_lenc(iyear)=(C_MRFSS_F(iyear)*lenprob_f+C_MRFSS_M(iyear)*lenprob_m)/sum(C_MRFSS_F(iyear)+C_MRFSS_M(iyear));
//  }
//    
FUNCTION get_age_comps
   //Hand lines
  for (iyear=1;iyear<=nyr_HL_agec;iyear++)
  {
    pred_HL_agec(iyear)=(C_HL_F(yrs_HL_agec(iyear))+C_HL_M(yrs_HL_agec(iyear)))/
                            (sum(C_HL_F(yrs_HL_agec(iyear))+C_HL_M(yrs_HL_agec(iyear)))+dzero_dum);
  }
  pred_HL_agec=pred_HL_agec*age_error_matrix;
 //  cout << "HL age comps" <<pred_HL_agec << endl;

  //Pound nets
  for (iyear=1;iyear<=nyr_PN_agec;iyear++)
  {
    pred_PN_agec(iyear)=(C_PN_F(yrs_PN_agec(iyear))+C_PN_M(yrs_PN_agec(iyear)))/
                             sum(C_PN_F(yrs_PN_agec(iyear))+C_PN_M(yrs_PN_agec(iyear)));
  }
  pred_PN_agec=pred_PN_agec*age_error_matrix;
 // cout << "PN age comps" <<pred_PN_agec << endl;
 //Gill nets
  for (iyear=1;iyear<=nyr_GN_agec;iyear++)
  {
    pred_GN_agec(iyear)=(C_GN_F(yrs_GN_agec(iyear))+C_GN_M(yrs_GN_agec(iyear)))/
                             sum(C_GN_F(yrs_GN_agec(iyear))+C_GN_M(yrs_GN_agec(iyear)));
  }  
//  cout << "GN age comps" <<pred_GN_agec << endl;
  //Cast nets
  for (iyear=1;iyear<=nyr_CN_agec;iyear++)
  {
    pred_CN_agec(iyear)=(C_CN_F(yrs_CN_agec(iyear))+C_CN_M(yrs_CN_agec(iyear)))/
                             sum(C_CN_F(yrs_CN_agec(iyear))+C_CN_M(yrs_CN_agec(iyear)));
  }
//  cout << "CN age comps" <<pred_CN_agec << endl;
  //MRFSS
  for (iyear=1;iyear<=nyr_MRFSS_agec;iyear++)
  {
    pred_MRFSS_agec(iyear)=(C_MRFSS_F(yrs_MRFSS_agec(iyear))+C_MRFSS_M(yrs_MRFSS_agec(iyear)))/
                             sum(C_MRFSS_F(yrs_MRFSS_agec(iyear))+C_MRFSS_M(yrs_MRFSS_agec(iyear)));
  }
  //cout << "MRFSS age comps" <<pred_MRFSS_agec << endl;
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_sel_weighted_current //get average of most recent 3 years selectivity for MSY calcs
  F_temp_sum=0.0;
  F_temp_sum+=mfexp((3.0*log_avg_F_HL+sum(log_F_dev_HL(endyr-2,endyr)))/3.0);
  F_temp_sum+=mfexp((3.0*log_avg_F_PN+sum(log_F_dev_PN(endyr-2,endyr)))/3.0);
  F_temp_sum+=mfexp((3.0*log_avg_F_GN+sum(log_F_dev_GN(endyr-2,endyr)))/3.0);
  F_temp_sum+=mfexp((3.0*log_avg_F_CN+sum(log_F_dev_CN(endyr-2,endyr)))/3.0);
  F_temp_sum+=mfexp((3.0*log_avg_F_MRFSS+sum(log_F_dev_MRFSS(endyr-2,endyr)))/3.0);

  F_temp_sum+=mfexp((3.0*log_avg_F_MRFSS_D+sum(log_F_dev_MRFSS_D(endyr-2,endyr)))/3.0);

  F_temp_sum+=mfexp((3.0*log_avg_F_shrimp+sum(log_F_dev_shrimp(endyr-2,endyr)))/3.0);
     
  F_HL_prop=mfexp((3.0*log_avg_F_HL+sum(log_F_dev_HL(endyr-2,endyr)))/3.0)/F_temp_sum;
  F_PN_prop=mfexp((3.0*log_avg_F_PN+sum(log_F_dev_PN(endyr-2,endyr)))/3.0)/F_temp_sum;
  F_GN_prop=mfexp((3.0*log_avg_F_GN+sum(log_F_dev_GN(endyr-2,endyr)))/3.0)/F_temp_sum;
  F_CN_prop=mfexp((3.0*log_avg_F_CN+sum(log_F_dev_CN(endyr-2,endyr)))/3.0)/F_temp_sum;
  F_MRFSS_prop=mfexp((3.0*log_avg_F_MRFSS+sum(log_F_dev_MRFSS(endyr-2,endyr)))/3.0)/F_temp_sum;
  
  F_MRFSS_D_prop=mfexp((3.0*log_avg_F_MRFSS_D+sum(log_F_dev_MRFSS_D(endyr-2,endyr)))/3.0)/F_temp_sum;

  F_shrimp_prop=mfexp((3.0*log_avg_F_shrimp+sum(log_F_dev_shrimp(endyr-2,endyr)))/3.0)/F_temp_sum;
  
  sel_wgted_L_F=F_HL_prop*sel_HL_keep_F+F_PN_prop*sel_PN_F+F_GN_prop*sel_GN_keep_F+F_CN_prop*sel_CN_F+F_MRFSS_prop*sel_MRFSS_keep_F;  
  sel_wgted_D_F=F_MRFSS_D_prop*sel_MRFSS_D_F+F_shrimp_prop*sel_shrimp;
  sel_wgted_tot_F=sel_wgted_L_F+sel_wgted_D_F;
  
  sel_wgted_L_M=F_HL_prop*sel_HL_keep_M+F_PN_prop*sel_PN_M+F_GN_prop*sel_GN_keep_M+F_CN_prop*sel_CN_M+F_MRFSS_prop*sel_MRFSS_keep_M;  
  sel_wgted_D_M=F_MRFSS_D_prop*sel_MRFSS_D_M+F_shrimp_prop*sel_shrimp;
  sel_wgted_tot_M=sel_wgted_L_M+sel_wgted_D_M;
                
  //max_sel_wgted_tot_F=max(sel_wgted_tot_F);              
  //sel_wgted_tot_F/=max_sel_wgted_tot_F;
 //sel_wgted_L_F/=max_sel_wgted_tot_F; //landings sel bumped up by same amount as total sel              
  //sel_wgted_D_F/=max_sel_wgted_tot_F;                 
  
  //max_sel_wgted_tot_M=max(sel_wgted_tot_M);  
  //sel_wgted_tot_M/=max_sel_wgted_tot_M;
  //sel_wgted_L_M/=max_sel_wgted_tot_M; //landings sel bumped up by same amount as total sel              
  //sel_wgted_D_M/=max_sel_wgted_tot_M;    
  
FUNCTION get_msy

  //fill in Fs for per-recruit stuff
  F_msy.fill_seqadd(0,.0001);  //step size should be of inverse dimension to n_iter_msy
    
  //compute values as functions of F
  for(ff=1; ff<=n_iter_msy; ff++){		//commented out 'int' before ff=1
    //int ff=1001;
    //uses fishery-weighted F's
    Z_age_msy_F=0.0;
    F_L_age_msy_F=0.0;
    F_D_age_msy_F=0.0;
    Z_age_msy_M=0.0;
    F_L_age_msy_M=0.0;
    F_D_age_msy_M=0.0;
      
    F_L_age_msy_F=F_msy(ff)*sel_wgted_L_F;
    F_D_age_msy_F=F_msy(ff)*sel_wgted_D_F;
    Z_age_msy_F=M+F_L_age_msy_F+F_D_age_msy_F;   
    F_L_age_msy_M=F_msy(ff)*sel_wgted_L_M;
    F_D_age_msy_M=F_msy(ff)*sel_wgted_D_M;
    Z_age_msy_M=M+F_L_age_msy_M+F_D_age_msy_M;          
    
    N_age_msy_F(1)=prop_f_a0;                            
    for (iage=2; iage<=nages; iage++)
    {
      N_age_msy_F(iage)=N_age_msy_F(iage-1)*mfexp(-1.*Z_age_msy_F(iage-1));
    }
    N_age_msy_F(nages)=N_age_msy_F(nages)/(1.0-mfexp(-1.*Z_age_msy_F(nages)));
    N_age_msy_mdyr_F(1,(nages-1))=elem_prod(N_age_msy_F(1,(nages-1)),
                                   mfexp((-1.*Z_age_msy_F(1,(nages-1)))/2.0));                 
    N_age_msy_mdyr_F(nages)=(N_age_msy_mdyr_F(nages-1)*
                          (mfexp(-1.*(Z_age_msy_F(nages-1)/2 + Z_age_msy_F(nages)/2) )))
                          /(1.0-mfexp(-1.*Z_age_msy_F(nages)));
    
    N_age_msy_M(1)=1.-prop_f_a0;                            
    for (iage=2; iage<=nages; iage++)
    {
      N_age_msy_M(iage)=N_age_msy_M(iage-1)*mfexp(-1.*Z_age_msy_M(iage-1));
    }
    N_age_msy_M(nages)=N_age_msy_M(nages)/(1.0-mfexp(-1.*Z_age_msy_M(nages)));
    N_age_msy_mdyr_M(1,(nages-1))=elem_prod(N_age_msy_M(1,(nages-1)),
                                   mfexp((-1.*Z_age_msy_M(1,(nages-1)))/2.0));                 
    N_age_msy_mdyr_M(nages)=(N_age_msy_mdyr_M(nages-1)*
                          (mfexp(-1.*(Z_age_msy_M(nages-1)/2 + Z_age_msy_M(nages)/2) )))
                          /(1.0-mfexp(-1.*Z_age_msy_M(nages)));                 
    
    for(iage=1;iage<=nages;iage++){   
      prop_f_F0(iage)=N_age_msy_F(iage)/(N_age_msy_F(iage)+N_age_msy_M(iage));  //not really at F=0 here; just using vector
    }
    
    spr_msy(ff)=sum(elem_prod(N_age_msy_mdyr_F,reprod));       

    //Compute equilibrium values of R (including bias correction), SSB and Yield at each F
    R_eq(ff)=(R0/((5.0*steep-1.0)*spr_msy(ff)))*
                 (BiasCor*4.0*steep*spr_msy(ff)-spr_F0*(1.0-steep));
    if (R_eq(ff)<dzero_dum) {R_eq(ff)=dzero_dum;}    
    N_age_msy_F*=R_eq(ff); //proportion female/male already accounted for
    N_age_msy_mdyr_F*=R_eq(ff);
    N_age_msy_M*=R_eq(ff);
    N_age_msy_mdyr_M*=R_eq(ff);
    
    for (iage=1; iage<=nages; iage++){
      C_age_msy_F(iage)=N_age_msy_F(iage)*(F_L_age_msy_F(iage)/Z_age_msy_F(iage))*
                      (1.-mfexp(-1.0*Z_age_msy_F(iage)));
      C_age_msy_M(iage)=N_age_msy_M(iage)*(F_L_age_msy_M(iage)/Z_age_msy_M(iage))*
                      (1.-mfexp(-1.0*Z_age_msy_M(iage)));
      D_age_msy_F(iage)=N_age_msy_F(iage)*(F_D_age_msy_F(iage)/Z_age_msy_F(iage))*
                      (1.-mfexp(-1.0*Z_age_msy_F(iage)));
      D_age_msy_M(iage)=N_age_msy_M(iage)*(F_D_age_msy_M(iage)/Z_age_msy_M(iage))*
                      (1.-mfexp(-1.0*Z_age_msy_M(iage)));
    }
    
    SSB_eq(ff)=sum(elem_prod(N_age_msy_mdyr_F,reprod));
    B_eq(ff)=sum(elem_prod(N_age_msy_F,wgt_f))+sum(elem_prod(N_age_msy_M,wgt_m));
    L_eq(ff)=sum(elem_prod(C_age_msy_F,wgt_land))+sum(elem_prod(C_age_msy_M,wgt_land));  
    E_eq(ff)=(sum(C_age_msy_F(E_age_st,nages))+sum(C_age_msy_M(E_age_st,nages)));
    E_eq(ff)/=(sum(N_age_msy_F(E_age_st,nages))+sum(N_age_msy_M(E_age_st,nages)));
    D_eq(ff)=(sum(D_age_msy_F)+sum(D_age_msy_M));//1000.0;
    
  }
  
  msy_out=max(L_eq);
  
  for(ff=1; ff<=n_iter_msy; ff++)
  {
   if(L_eq(ff) == msy_out) 
      {    
        SSB_msy_out=SSB_eq(ff);
        B_msy_out=B_eq(ff);
        R_msy_out=R_eq(ff);
        D_msy_out=D_eq(ff);
        E_msy_out=E_eq(ff);
        F_msy_out=F_msy(ff);  
        spr_msy_out=spr_msy(ff);      
      }
  }

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_miscellaneous_stuff
  
  sigma_rec_dev=sqrt(var_rec_dev+0.0001); //sample SD of predicted residuals (may not equal rec_sigma)  
  
  //compute total catch-at-age and landings
  C_total=C_HL_F+C_HL_M+C_PN_F+C_PN_M+C_GN_M+C_GN_F+C_CN_M+C_CN_F+C_MRFSS_M+C_MRFSS_F; //catch in number fish
  L_total=L_HL+L_PN+L_GN+L_CN+L_MRFSS;   //landings in mt whole weight
  D_total=L_MRFSS_D; //total discards in mt whole weight
  B_total=L_shrimp;
  
  //compute exploitation rate of age E_age_st +
  for(iyear=styr; iyear<=endyr; iyear++)
  {
    E(iyear)=(sum(C_total(iyear)(E_age_st,nages)))/(sum(N_F(iyear)(E_age_st,nages))+sum(N_M(iyear)(E_age_st,nages))); 
    L_total_yr(iyear)=sum(L_total(iyear));
    B_total_yr(iyear)=sum(B_total(iyear));
    D_total_yr(iyear)=sum(D_total(iyear));
  }
  
  steep_sd=steep;
  fullF_sd=fullF;
  E_sd=E;
  
  if(E_msy_out>0)
    {
      EdE_msy=E/E_msy_out;
      EdE_msy_end=EdE_msy(endyr);
    }
  if(F_msy_out>0)
    {
      FdF_msy=fullF/F_msy_out;
      FdF_msy_end=FdF_msy(endyr);
    }
  if(SSB_msy_out>0)
    {
      SdSSB_msy=SSB/SSB_msy_out;
      SdSSB_msy_end=SdSSB_msy(endyr);
    }  

   //fill in log recruitment deviations for yrs they are nonzero
   for(iyear=styr_rec_dev; iyear<=endyr; iyear++)
   { 
    log_dev_R(iyear)=log_rec_dev(iyear);
   }

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_per_recruit_stuff

  //static per-recruit stuff
 
  for(iyear=styr; iyear<=endyr; iyear++)
  {
    N_age_spr_F(1)=prop_f_a0;
    for(iage=2; iage<=nages; iage++)
    {
      N_age_spr_F(iage)=N_age_spr_F(iage-1)*mfexp(-1.*Z_F(iyear,iage-1));
    }
    N_age_spr_F(nages)=N_age_spr_F(nages)/(1.0-mfexp(-1.*Z_F(iyear,nages)));    
    N_age_spr_mdyr_F(1,(nages-1))=elem_prod(N_age_spr_F(1,(nages-1)),
                                mfexp(-1.*Z_F(iyear)(1,(nages-1))/2.0));
    N_age_spr_mdyr_F(nages)=(N_age_spr_mdyr_F(nages-1)*
                          (mfexp(-1.*(Z_F(iyear)(nages-1)/2.0 + Z_F(iyear)(nages)/2.0) )))
                          /(1.0-mfexp(-1.*Z_F(iyear)(nages)));           
    spr_static(iyear)=sum(elem_prod(N_age_spr_mdyr_F,reprod))/spr_F0;
  }
  
  //fill in Fs for per-recruit stuff
  F_spr.fill_seqadd(0,.001);
  
  //compute SSB/R and YPR as functions of F
  for(int ff=1; ff<=n_iter_spr; ff++)
  {
    //uses fishery-weighted F's, same as in MSY calculations
    Z_age_spr_F=0.0;
    F_L_age_spr_F=0.0;
    Z_age_spr_M=0.0;
    F_L_age_spr_M=0.0;

    F_L_age_spr_F=F_spr(ff)*sel_wgted_L_F;
    F_L_age_spr_M=F_spr(ff)*sel_wgted_L_M;
    
    Z_age_spr_F=M+F_L_age_spr_F+F_spr(ff)*sel_wgted_D_F;
    Z_age_spr_M=M+F_L_age_spr_M+F_spr(ff)*sel_wgted_D_M;

    N_age_spr_F(1)=prop_f_a0;
    N_age_spr_M(1)=1.-prop_f_a0;
    for (iage=2; iage<=nages; iage++){
      N_age_spr_F(iage)=N_age_spr_F(iage-1)*mfexp(-1.*Z_age_spr_F(iage-1));
      N_age_spr_M(iage)=N_age_spr_M(iage-1)*mfexp(-1.*Z_age_spr_M(iage-1));
    }
    N_age_spr_F(nages)=N_age_spr_F(nages)/(1-mfexp(-1.*Z_age_spr_F(nages)));
    N_age_spr_M(nages)=N_age_spr_M(nages)/(1-mfexp(-1.*Z_age_spr_M(nages)));

    N_age_spr_mdyr_F(1,(nages-1))=elem_prod(N_age_spr_F(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr_F(1,(nages-1)))/2.0));                 
    N_age_spr_mdyr_F(nages)=(N_age_spr_mdyr_F(nages-1)*
                          (mfexp(-1.*(Z_age_spr_F(nages-1)/2 + Z_age_spr_F(nages)/2) )))
                          /(1.0-mfexp(-1.*Z_age_spr_F(nages)));
    
    for(iage=1;iage<=nages;iage++){   
      prop_f_F0(iage)=N_age_spr_F(iage)/(N_age_spr_F(iage)+N_age_spr_M(iage));  //not really at F=0 here; just using vector
    }
    spr_spr(ff)=sum(elem_prod(N_age_spr_mdyr_F,reprod));
    L_spr(ff)=0.0;
    for (iage=1; iage<=nages; iage++)
    {
      C_age_spr_F(iage)=N_age_spr_F(iage)*(F_L_age_spr_F(iage)/Z_age_spr_F(iage))*
                      (1.-mfexp(-1.*Z_age_spr_F(iage)));
      C_age_spr_M(iage)=N_age_spr_M(iage)*(F_L_age_spr_M(iage)/Z_age_spr_M(iage))*
                      (1.-mfexp(-1.*Z_age_spr_M(iage)));
      L_spr(ff)+=(C_age_spr_M(iage)*wgt_m(iage)+C_age_spr_F(iage)*wgt_land(iage));
    }
    E_spr(ff)=(sum(C_age_spr_M(E_age_st,nages))+sum(C_age_spr_F(E_age_st,nages)))/
              (sum(N_age_spr_M(E_age_st,nages))+sum(N_age_spr_F(E_age_st,nages)));
    
  }
  
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   

FUNCTION evaluate_objective_function
  fval=0.0;
  fval_unwgt=0.0;

//---likelihoods---------------------------
	//	fval=square(x_dum-3.0);

//---Indices-------------------------------
  f_FL_HL_cpue=0.0;
  for (iyear=styr_FL_HL_cpue; iyear<=endyr_FL_HL_cpue; iyear++)
  {
    f_FL_HL_cpue+=square(log((pred_FL_HL_cpue(iyear)+dzero_dum)/
        (obs_FL_HL_cpue(iyear)+dzero_dum)))/(2.0*square(FL_HL_cpue_cv(iyear)));
  }
  fval+=w_I_FL_HL*f_FL_HL_cpue;
  fval_unwgt+=f_FL_HL_cpue;
   
  f_MRFSS_cpue=0.0;
  for (iyear=styr_MRFSS_cpue; iyear<=endyr_MRFSS_cpue; iyear++)
  {
    f_MRFSS_cpue+=square(log((pred_MRFSS_cpue(iyear)+dzero_dum)/
        (obs_MRFSS_cpue(iyear)+dzero_dum)))/(2.0*square(MRFSS_cpue_cv(iyear)));
  }
  //cout<<w_I_MRFSS<<endl<<f_MRFSS_cpue<<endl<<endl;
  fval+=w_I_MRFSS*f_MRFSS_cpue;
  fval_unwgt+=f_MRFSS_cpue;

  f_SMAP_YOY_cpue=0.0;
  for (iyear=styr_SMAP_YOY_cpue; iyear<=endyr_SMAP_YOY_cpue; iyear++)
  {
    f_SMAP_YOY_cpue+=square(log((pred_SMAP_YOY_cpue(iyear)+dzero_dum)/
        (obs_SMAP_YOY_cpue(iyear)+dzero_dum)))/(2.0*square(SMAP_YOY_cpue_cv(iyear)));
  }
  fval+=w_I_SMAP_YOY*f_SMAP_YOY_cpue;
  fval_unwgt+=f_SMAP_YOY_cpue;

  //f_SMAP_1YR_cpue=0.0;
  //for (iyear=styr_SMAP_1YR_cpue; iyear<=endyr_SMAP_1YR_cpue; iyear++)
  //{
  //  f_SMAP_1YR_cpue+=square(log((pred_SMAP_1YR_cpue(iyear)+dzero_dum)/
  //      (obs_SMAP_1YR_cpue(iyear)+dzero_dum)))/(2.0*square(SMAP_1YR_cpue_cv(iyear)));
  //}
  //fval+=w_I_SMAP_1YR*f_SMAP_1YR_cpue;
  //fval_unwgt+=f_SMAP_1YR_cpue;
  
//---Landings------------------------------- 

  f_HL_L=0.0; //in 1000s total pounds
  for (iyear=styr_HL_L; iyear<=endyr_HL_L; iyear++)
  {
      f_HL_L+=square(log((pred_HL_L(iyear)+dzero_dum)/
        (obs_HL_L(iyear)+dzero_dum)))/(2.0*square(HL_L_cv(iyear)));
  }
  fval+=w_L*f_HL_L;
  fval_unwgt+=f_HL_L;
  
  f_PN_L=0.0; //in 1000s total pounds
  for (iyear=styr_PN_L; iyear<=endyr_PN_L; iyear++)
  {
       f_PN_L+=square(log((pred_PN_L(iyear)+dzero_dum)/
        (obs_PN_L(iyear)+dzero_dum)))/(2.0*square(PN_L_cv(iyear)));
  }
  fval+=w_L*f_PN_L;
  fval_unwgt+=f_PN_L;  

  f_GN_L=0.0; //in 1000s total pounds
  for (iyear=styr_GN_L; iyear<=endyr_GN_L; iyear++)
  {
      f_GN_L+=square(log((pred_GN_L(iyear)+dzero_dum)/
        (obs_GN_L(iyear)+dzero_dum)))/(2.0*square(GN_L_cv(iyear)));
  }
  fval+=w_L*f_GN_L;
  fval_unwgt+=f_GN_L;  
  
  f_CN_L=0.0; //in 1000s total pounds
  for (iyear=styr_CN_L; iyear<=endyr_CN_L; iyear++)
  {
    f_CN_L+=square(log((pred_CN_L(iyear)+dzero_dum)/
        (obs_CN_L(iyear)+dzero_dum)))/(2.0*square(CN_L_cv(iyear)));
  }
  fval+=w_L*f_CN_L;
  fval_unwgt+=f_CN_L;  
  
  f_MRFSS_L=0.0; //in 1000s - numbers
  for (iyear=styr_MRFSS_L; iyear<=endyr_MRFSS_L; iyear++)
  {
      f_MRFSS_L+=square(log((pred_MRFSS_L(iyear)+dzero_dum)/
        (obs_MRFSS_L(iyear)+dzero_dum)))/(2.0*square(MRFSS_L_cv(iyear)));
  }
  fval+=w_L*f_MRFSS_L;
  fval_unwgt+=f_MRFSS_L;  
  
//---Discards & Bycatch-------------------------------
  
  f_MRFSS_D=0.0; //in 1000s (numbers)
  for (iyear=styr_MRFSS_D; iyear<=endyr_MRFSS_D; iyear++)
  {
      f_MRFSS_D+=square(log((pred_MRFSS_D(iyear)+dzero_dum)/
        (obs_MRFSS_D(iyear)+dzero_dum)))/(2.0*square(MRFSS_D_cv(iyear)));
  }
  fval+=w_D*f_MRFSS_D;
  fval_unwgt+=f_MRFSS_D;
  
  f_shrimp_B=0.0; 
  for (iyear=styr_shrimp_B; iyear<=endyr_shrimp_B; iyear++)
  {
      f_shrimp_B+=square(log((pred_shrimp_B(iyear)+dzero_dum)/
        (obs_shrimp_B(iyear)+dzero_dum)))/(2.0*square(shrimp_B_cv(iyear)));
  //      f_shrimp_B+=square(log((pred_shrimp_B(iyear)+dzero_dum)/
  //        (obs_shrimp_geomean+dzero_dum)))/(2.0*square(shrimp_B_cv(iyear)));
  }
  //f_shrimp_B=square(obs_shrimp_geomean-pred_shrimp_geomean);
  
  fval+=w_D*f_shrimp_B;
  //fval+=w_shrimp*f_shrimp_B;
  fval_unwgt+=f_shrimp_B; 

//---Length comps-------------------------------

//  f_HL_lenc=0.0;
//  f_HL_lenc=lk_robust_multinomial(nsamp_HL_lenc, pred_HL_lenc, obs_HL_lenc, nyr_HL_lenc, double(nlenbins), minSS_HL_lenc, w_lc_HL);
//  fval+=f_HL_lenc;
//  fval_unwgt+=f_HL_lenc;
//  for (iyear=1; iyear<=nyr_HL_lenc; iyear++)
//  {
//    f_HL_lenc-=nsamp_HL_lenc(iyear)*
//        sum( elem_prod((obs_HL_lenc(iyear)+dzero_dum),
//            log(elem_div((pred_HL_lenc(iyear)+dzero_dum),
//               (obs_HL_lenc(iyear)+dzero_dum)))));  
//  }
//  fval+=w_lc_HL*f_HL_lenc;
//  fval_unwgt+=f_HL_lenc;
//  
//  f_PN_lenc=0.0;
//  f_PN_lenc=lk_robust_multinomial(nsamp_PN_lenc, pred_PN_lenc, obs_PN_lenc, nyr_PN_lenc, double(nlenbins), minSS_PN_lenc, w_lc_PN);
//  fval+=f_PN_lenc;
//  fval_unwgt+=f_PN_lenc;
//  for (iyear=1; iyear<=nyr_PN_lenc; iyear++)
//  {
//    f_PN_lenc-=nsamp_PN_lenc(iyear)*
//        sum( elem_prod((obs_PN_lenc(iyear)+dzero_dum),
//            log(elem_div((pred_PN_lenc(iyear)+dzero_dum),
//               (obs_PN_lenc(iyear)+dzero_dum)))));  
//  }
//  fval+=w_lc_PN*f_PN_lenc;
//  fval_unwgt+=f_PN_lenc;
//
//  f_GN_lenc=0.0;
//  f_GN_lenc=lk_robust_multinomial(nsamp_GN_lenc, pred_GN_lenc, obs_GN_lenc, nyr_GN_lenc, double(nlenbins), minSS_GN_lenc, w_lc_GN);
//  fval+=f_GN_lenc;
//  fval_unwgt+=f_GN_lenc;
//  for (iyear=1; iyear<=nyr_GN_lenc; iyear++)
//  {
//    f_GN_lenc-=nsamp_GN_lenc(iyear)*
//        sum( elem_prod((obs_GN_lenc(iyear)+dzero_dum),
//            log(elem_div((pred_GN_lenc(iyear)+dzero_dum),
//               (obs_GN_lenc(iyear)+dzero_dum)))));
//  }
//  fval+=w_lc_GN*f_GN_lenc;
//  fval_unwgt+=f_GN_lenc;
//	f_CN_lenc=0.0;
//	f_CN_lenc=lk_robust_multinomial(nsamp_CN_lenc, pred_CN_lenc, obs_CN_lenc, nyr_CN_lenc, double(nlenbins), minSS_CN_lenc, w_lc_CN);
//  fval+=f_CN_lenc;
//  fval_unwgt+=f_CN_lenc;
//  for (iyear=1; iyear<=nyr_CN_lenc; iyear++)
//  {
//    f_CN_lenc-=nsamp_CN_lenc(iyear)*
//         sum(elem_prod((obs_CN_lenc(iyear)+dzero_dum),
//             log(elem_div((pred_CN_lenc(iyear)+dzero_dum), 
//                (obs_CN_lenc(iyear)+dzero_dum)))));
//  }
//  fval+=w_lc_CN*f_CN_lenc;
//  fval_unwgt+=f_CN_lenc;
//
//	f_MRFSS_lenc=0.0;
//	f_MRFSS_lenc=lk_robust_multinomial(nsamp_MRFSS_lenc, pred_MRFSS_lenc, obs_MRFSS_lenc, nyr_MRFSS_lenc, double(nlenbins), minSS_MRFSS_lenc, w_lc_MRFSS);
//  fval+=f_MRFSS_lenc;
//  fval_unwgt+=f_MRFSS_lenc;
//  for (iyear=styr_MRFSS_lenc; iyear<=endyr_MRFSS_lenc; iyear++)
//  {
//    f_MRFSS_lenc-=nsamp_MRFSS_lenc(iyear)*
//         sum(elem_prod((obs_MRFSS_lenc(iyear)+dzero_dum),
//             log(elem_div((pred_MRFSS_lenc(iyear)+dzero_dum), 
//                (obs_MRFSS_lenc(iyear)+dzero_dum)))));
//  }
//  fval+=w_lc_MRFSS*f_MRFSS_lenc;
//  fval_unwgt+=f_MRFSS_lenc;
///---Age comps-------------------------------
  f_HL_agec=0.0;
  f_HL_agec=lk_robust_multinomial(nsamp_HL_agec, pred_HL_agec, obs_HL_agec, nyr_HL_agec, double(nages), minSS_HL_agec, w_ac_HL);
  fval+=f_HL_agec;
  fval_unwgt+=f_HL_agec;
  
//  for (iyear=1; iyear<=nyr_HL_agec; iyear++)
//  {
//    f_HL_agec-=nsamp_HL_agec(iyear)*
//        sum( elem_prod((obs_HL_agec(iyear)+dzero_dum),
//            log(elem_div((pred_HL_agec(iyear)+dzero_dum),
//               (obs_HL_agec(iyear)+dzero_dum)))));  
//  }
//  fval+=w_ac_HL*f_HL_agec;
    
  f_PN_agec=0.0;
 // f_PN_agec=lk_robust_multinomial(nsamp_PN_agec, pred_PN_agec, obs_PN_agec, nyr_PN_agec, double(nages), minSS_PN_agec, w_ac_PN);
 // fval+=f_PN_agec;
 // fval_unwgt+=f_PN_agec;
  
  for (iyear=1; iyear<=nyr_PN_agec; iyear++)
  {
    f_PN_agec-=nsamp_PN_agec(iyear)*
        sum( elem_prod((obs_PN_agec(iyear)+dzero_dum),
            log(elem_div((pred_PN_agec(iyear)+dzero_dum),
               (obs_PN_agec(iyear)+dzero_dum)))));  
  }
  fval+=w_ac_PN*f_PN_agec;
  fval_unwgt+=f_PN_agec; 

  f_GN_agec=0.0;
  f_GN_agec=lk_robust_multinomial(nsamp_GN_agec, pred_GN_agec, obs_GN_agec, nyr_GN_agec, double(nages), minSS_GN_agec, w_ac_GN);
  fval+=f_GN_agec;
  fval_unwgt+=f_GN_agec;
  
//  for (iyear=1; iyear<=nyr_GN_agec; iyear++)
//  {
//    f_GN_agec-=nsamp_GN_agec(iyear)*
//        sum( elem_prod((obs_GN_agec(iyear)+dzero_dum),
//            log(elem_div((pred_GN_agec(iyear)+dzero_dum),
//               (obs_GN_agec(iyear)+dzero_dum)))));  
//  }
//  fval+=w_ac_GN*f_GN_agec;
//  fval_unwgt+=f_GN_agec;
    
  f_CN_agec=0.0;
  f_CN_agec=lk_robust_multinomial(nsamp_CN_agec, pred_CN_agec, obs_CN_agec, nyr_CN_agec, double(nages), minSS_CN_agec, w_ac_CN);
  fval+=f_CN_agec;
  fval_unwgt+=f_CN_agec;
  
//  for (iyear=1; iyear<=nyr_CN_agec; iyear++)
//  {
//    f_CN_agec-=nsamp_CN_agec(iyear)*
//        sum( elem_prod((obs_CN_agec(iyear)+dzero_dum),
//            log(elem_div((pred_CN_agec(iyear)+dzero_dum),
//               (obs_CN_agec(iyear)+dzero_dum)))));  
//  }
//  fval+=w_ac_CN*f_CN_agec;
//  fval_unwgt+=f_CN_agec;
  
  f_MRFSS_agec=0.0;
  f_MRFSS_agec=lk_robust_multinomial(nsamp_MRFSS_agec, pred_MRFSS_agec, obs_MRFSS_agec, nyr_MRFSS_agec, double(nages), minSS_MRFSS_agec, w_ac_MRFSS);
  fval+=f_MRFSS_agec;
  fval_unwgt+=f_MRFSS_agec;
  
//  for (iyear=styr_MRFSS_agec; iyear<=endyr_MRFSS_agec; iyear++)
//  {
//    f_MRFSS_agec-=nsamp_MRFSS_agec(iyear)*
//        sum( elem_prod((obs_MRFSS_agec(iyear)+dzero_dum),
//            log(elem_div((pred_MRFSS_agec(iyear)+dzero_dum),
//               (obs_MRFSS_agec(iyear)+dzero_dum)))));  
//  }
//  fval+=w_ac_MRFSS*f_MRFSS_agec;
//  fval_unwgt+=f_MRFSS_agec;
  
//-----------Constraints and penalties--------------------------------
  f_N_dev=0.0;
    f_N_dev=pow(log_rec_dev(styr_rec_dev),2);
  for(iyear=(styr_rec_dev+1); iyear<=endyr; iyear++)
  {
    f_N_dev+=pow(log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1),2);
  }
  fval+=w_R*f_N_dev;

//  f_N_dev_early=0.0;
//  f_N_dev_early=norm2(log_dev_N_rec(styr_rec_dev,(styr_rec_dev+5)));
//  fval+=w_R_init*f_N_dev_early;

  f_N_dev_end=0.0; //last 3 yrs
    f_N_dev_end=norm2(log_rec_dev(endyr-2,endyr));
  fval+=w_R_end*f_N_dev_end;

//  f_B1dB0_constraint=0.0;
//    f_B1dB0_constraint=square(totB(styr)/B0-B1dB0);
//  fval+=w_B1dB0*f_B1dB0_constraint;

  f_Fend_constraint=0.0; //last 3 yrs
  f_Fend_constraint=norm2(first_difference(fullF(endyr-2,endyr)));
  fval+=w_F*f_Fend_constraint;

  f_fullF_constraint=0.0;
  for (iyear=styr; iyear<=endyr; iyear++)
    {
     if (fullF(iyear)>3.0)
     {
     f_fullF_constraint+=square(fullF(iyear)-3.0);
     }
    }
  fval+=w_fullF*f_fullF_constraint;
  
//  f_cvlen_diff_constraint=0.0;
//    f_cvlen_diff_constraint=norm2(first_difference(log_len_cv_dev));
//  fval+=w_cvlen_diff*f_cvlen_diff_constraint;
//  
//  f_cvlen_dev_constraint=0.0;
//    f_cvlen_dev_constraint=norm2(log_len_cv_dev);  
//  fval+=w_cvlen_dev*f_cvlen_dev_constraint;

//  cout << "fval = " << fval << "  fval_unwgt = " << fval_unwgt << endl;
  //cout << "avg MRFSS " << log_avg_F_MRFSS <<endl;;
  
  //cout<<"pred ac"<<endl<<Pred_HL_agec<<endl<<endl<<"obs ac"<<endl<<endl<<obs_HL_agec;
//---Priors---------------------------------------------------
//neg_log_prior arguments: estimate, prior, variance, pdf type
//Variance input as a negative value is considered to be CV in arithmetic space (CV=-1 implies loose prior) 
//pdf type 1=none, 2=lognormal, 3=normal, 4=beta 
  f_priors=0.0; 
  //f_priors+=neg_log_prior(len_cv_val, set_len_cv, square(set_len_cv_se), 2); 
  //f_priors+=neg_log_prior(steep, set_steep, square(set_steep_se), 3); 
  //f_priors+=neg_log_prior(rec_sigma,set_rec_sigma,square(set_rec_sigma_se),3);   
  //f_priors+=neg_log_prior(R_autocorr,set_R_autocorr, 1.0, 1);
  //f_priors+=neg_log_prior(q_rate, set_q_rate, 0.00001+square(set_q_rate), 2);  
  //f_priors+=neg_log_prior(F_init, set_F_init, -1.0 , 2);
  //f_priors+=neg_log_prior(M_constant, set_M_constant, square(set_M_constant_se), 2);	
  
  //f_priors+=neg_log_prior(selpar_L50_HL, set_selpar_L50_HL, 1.0, 3);  
  f_priors+=neg_log_prior(selpar_slope_HL,set_selpar_slope_HL, -0.25, 3);
  f_priors+=neg_log_prior(selpar_L501_PN, set_selpar_L501_PN, -0.15, 3);  
  f_priors+=neg_log_prior(selpar_slope1_PN,set_selpar_slope1_PN, -0.15, 3);
  //f_priors+=neg_log_prior(selpar_L50_GN, set_selpar_L50_GN, 1.0, 3);  
  f_priors+=neg_log_prior(selpar_slope_GN,set_selpar_slope_GN, -0.15, 3);
  f_priors+=neg_log_prior(selpar_L502_CN, set_selpar_L502_CN, 1.0, 3);  
  f_priors+=neg_log_prior(selpar_slope1_CN,set_selpar_slope1_CN, -0.25, 3); 
  f_priors+=neg_log_prior(selpar_slope2_CN,set_selpar_slope2_CN, -0.25, 3);
  
  f_priors+=neg_log_prior(selpar_L502_MRFSS_keep, set_selpar_L502_MRFSS_keep, -0.25, 3);  
  //f_priors+=neg_log_prior(selpar_slope2_MRFSS_keep,set_selpar_slope2_MRFSS_keep, -0.5, 2);
  //f_priors+=neg_log_prior(selpar_L501_MRFSS_keep, set_selpar_L501_MRFSS_keep, -0.25, 2);    
  f_priors+=neg_log_prior(selpar_slope1_MRFSS_keep,set_selpar_slope1_MRFSS_keep, -0.25, 3);
  f_priors+=neg_log_prior(selpar_slope2_MRFSS_keep,set_selpar_slope2_MRFSS_keep, -0.25, 3);    
  f_priors+=neg_log_prior(selpar_L502_PN, set_selpar_L502_PN, -0.25, 3);
  f_priors+=neg_log_prior(selpar_slope2_PN, set_selpar_slope2_PN, -0.25, 3);
  //f_priors+=neg_log_prior(selpar_min_PN, set_selpar_min_PN, -1.0, 1);
  //f_priors+=neg_log_prior(selpar_afull_cL, set_selpar_afull_cL, -1.0, 1);    
   
  fval+=f_priors;
  
  if(!last_phase())
  {
    for (iyear=styr; iyear<=endyr; iyear++)
    {
      if(fullF(iyear)>1.0)
      {
        fval+=10*(mfexp(fullF(iyear)-1.0)-1.0);
      }
    }
  }
  //cout << "fval = " << fval << "  fval_data = " << fval_data << endl;
  //cout << endl;
  //cout << "f_cL_U = " << f_cL_cpue << " f_mm_U = " << f_mm_cpue << " f_cL_L = " << f_cL_L << " f_cH_L = " << f_cH_L 
  //<< " f_rA_L = " << f_rA_L << endl;
  //cout << "f_cL_lenc = " << f_cL_lenc << " f_cH_lenc = " << f_cH_lenc << " f_rA_lenc = " << f_rA_lenc << " f_mm_lenc = " << f_mm_lenc << endl;
  //cout << "f_cL_agec = " << f_cL_agec << " f_cH_agec = " << f_cH_agec << " f_mm_agec = " << f_mm_agec << endl; 
  //cout << "f_rec_dev = " << f_rec_dev << " f_rec_dev_early = " << f_rec_dev_early << " f_rec_dev_end = " << f_rec_dev_end <<  " f_rec_hist_dev = " << f_rec_historic_dev <<  " f_cL_RW = " << f_cL_RW_cpue << endl;
  //cout << endl;
//----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------                                                                           
FUNCTION dvar_vector logistic_exponential(const dvar_vector& ages, const dvariable& L50, const dvariable& slope, const dvariable& sigma, const dvariable& joint)
  //ages=vector of ages, L50=age at 50% sel (ascending limb), slope=rate of increase, sigma=controls rate of descent (descending)                               
  //joint=age to join curves                                                                                                                                    
  RETURN_ARRAYS_INCREMENT();                                                                                                                                    
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());                                                                                                         
  Sel_Tmp=1.0;                                                                                                                                                  
  for (iage=1; iage<=nages; iage++)                                                                                                                             
  {                                                                                                                                                             
   if (ages(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope*(ages(iage)-L50)));}                                                                             
   if (ages(iage)>joint){Sel_Tmp(iage)=mfexp(-1.*square((ages(iage)-joint)/sigma));}                                                                            
  }                                                                                                                                                             
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);                                                                                                                                 
  RETURN_ARRAYS_DECREMENT();                                                                                                                                    
  return Sel_Tmp;   
  
//---------------------------------------------------------------------------------------  
//Logistic function: 2 parameters                                                                                                                                                
//FUNCTION dvar_vector logistic(const dvar_vector& ages, const dvariable& L50, const dvariable& slope)                                                                           
//  //ages=vector of ages, L50=age at 50% selectivity, slope=rate of increase                                                                                                    
//  RETURN_ARRAYS_INCREMENT();                                                                                                                                                   
//  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());                                                                                                                        
//  Sel_Tmp=1./(1.+mfexp(-1.*slope*(ages-L50))); //logistic;                                                                                                                     
//  RETURN_ARRAYS_DECREMENT();                                                                                                                                                   
//  return Sel_Tmp;                                                                                                                                                              
                                                                                                                                                                                 
//-----------------------------------------------------------------------------------                                                                                            
//Logistic function: 4 parameters                                                                                                                                                
FUNCTION dvar_vector logistic_double(const dvar_vector& ages, const dvariable& L501, const dvariable& slope1, const dvariable& L502, const dvariable& slope2)                  
  //ages=vector of ages, L50=age at 50% selectivity, slope=rate of increase, L502=age at 50% decrease additive to L501, slope2=slope of decrease                               
  RETURN_ARRAYS_INCREMENT();                                                                                                                                                   
                                                                                                                                                                               
                                                                                                                                                                               
                                                                                                                                            
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=elem_prod( (1./(1.+mfexp(-1.*slope1*(ages-L501)))),(1.-(1./(1.+mfexp(-1.*slope2*(ages-(L501+L502)))))) );     
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

//-----------------------------------------------------------------------------------
//Jointed logistic function: 6 parameters (increasing and decreasing logistics joined at peak selectivity)
//FUNCTION dvar_vector logistic_joint(const dvar_vector& ages, const dvariable& L501, const dvariable& slope1, const dvariable& L502, const dvariable& slope2, const dvariable& satval, const dvariable& joint)
//  //ages=vector of ages, L501=age at 50% sel (ascending limb), slope1=rate of increase,L502=age at 50% sel (descending), slope1=rate of increase (ascending), 
//  //satval=saturation value of descending limb, joint=location in age vector to join curves (may equal age or age + 1 if age-0 is included)
//  RETURN_ARRAYS_INCREMENT();
//  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
//  Sel_Tmp=1.0; 
//  for (iage=1; iage<=nages; iage++)
//  {
//   if (double(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope1*(ages(iage)-L501)));}  
//   if (double(iage)>joint){Sel_Tmp(iage)=1.0-(1.0-satval)/(1.+mfexp(-1.*slope2*(ages(iage)-L502)));}  
//  }  
//  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
//  RETURN_ARRAYS_DECREMENT();
//  return Sel_Tmp;

//-----------------------------------------------------------------------------------  
//Double Gaussian function: 6 parameters (as in SS3)
FUNCTION dvar_vector gaussian_double(const dvar_vector& ages, const dvariable& peak, const dvariable& top, const dvariable& ascwid, const dvariable& deswid, const dvariable& init, const dvariable& final)
  //ages=vector of ages, peak=ascending inflection location (as logistic), top=width of plateau, ascwid=ascent width (as log(width))
  //deswid=descent width (as log(width))
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step1(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step2(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step3(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step4(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step5(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step6(ages.indexmin(),ages.indexmax());
  dvar_vector pars_tmp(1,6); dvar_vector sel_tmp_iq(1,2);
  
  pars_tmp(1)=peak;
  pars_tmp(2)=peak+1.0+(0.99*ages(nages)-peak-1.0)/(1.0+mfexp(-top));
  pars_tmp(3)=mfexp(ascwid);
  pars_tmp(4)=mfexp(deswid);
  pars_tmp(5)=1.0/(1.0+mfexp(-init));
  pars_tmp(6)=1.0/(1.0+mfexp(-final));
       
  sel_tmp_iq(1)=mfexp(-(square(ages(1)-pars_tmp(1))/pars_tmp(3)));
  sel_tmp_iq(2)=mfexp(-(square(ages(nages)-pars_tmp(2))/pars_tmp(4)));
  
  sel_step1=mfexp(-(square(ages-pars_tmp(1))/pars_tmp(3)));
  sel_step2=pars_tmp(5)+(1.0-pars_tmp(5))*(sel_step1-sel_tmp_iq(1))/(1.0-sel_tmp_iq(1));  
  sel_step3=mfexp(-(square(ages-pars_tmp(2))/pars_tmp(4)));
  sel_step4=1.0+(pars_tmp(6)-1.0)*(sel_step3-1.0)/(sel_tmp_iq(2)-1.0);
  sel_step5=1.0/ (1.0+mfexp(-(20.0* elem_div((ages-pars_tmp(1)), (1.0+sfabs(ages-pars_tmp(1)))) )));
  sel_step6=1.0/(1.0+mfexp(-(20.0*elem_div((ages-pars_tmp(2)),(1.0+sfabs(ages-pars_tmp(2)))) )));  

  Sel_Tmp=elem_prod(sel_step2,(1.0-sel_step5))+ 
          elem_prod(sel_step5,((1.0-sel_step6)+ elem_prod(sel_step4,sel_step6)) ); 
 
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
    
//-----------------------------------------------------------------------------------    
//Spawner-recruit function (Beverton-Holt)
//FUNCTION dvariable SR_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& SSB)
//  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, SSB=spawning biomass
//  RETURN_ARRAYS_INCREMENT();
//  dvariable Recruits_Tmp;
//  Recruits_Tmp=((0.8*R0*h*SSB)/(0.2*R0*spr_F0*(1.0-h)+(h-0.2)*SSB));
//  RETURN_ARRAYS_DECREMENT();
//  return Recruits_Tmp;

//-----------------------------------------------------------------------------------
//compute multinomial effective sample size for a single yr
//FUNCTION dvariable multinom_eff_N(const dvar_vector& pred_comp, const dvar_vector& obs_comp)
//  //pred_comp=vector of predicted comps, obscomp=vector of observed comps
//  dvariable EffN_Tmp; dvariable numer; dvariable denom;
//  RETURN_ARRAYS_INCREMENT();
//  numer=sum( elem_prod(pred_comp,(1.0-pred_comp)) );
//  denom=sum( square(obs_comp-pred_comp) );
//  if (denom>0.0) {EffN_Tmp=numer/denom;}
//  else {EffN_Tmp=-missing;}                            
//  RETURN_ARRAYS_DECREMENT();
//  return EffN_Tmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: lognormal
//FUNCTION dvariable lk_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
//  //pred=vector of predicted vals, obs=vector of observed vals, cv=vector of CVs in arithmetic space, wgt_dat=constant scaling of CVs
//  
//  RETURN_ARRAYS_INCREMENT();
//  dvariable LkvalTmp;
//  dvar_vector var(cv.indexmin(),cv.indexmax()); //variance in log space
//  var=log(1.0+square(cv/wgt_dat));   // convert cv in arithmetic space to variance in log space
//  LkvalTmp=sum(0.5*elem_div(square(log(elem_div((pred+0.00001),(obs+0.00001)))),var) );
//  RETURN_ARRAYS_DECREMENT();
//  return LkvalTmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: multinomial
//FUNCTION dvariable lk_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const double& minSS, const dvariable& wgt_dat)
//  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, minSS=min N threshold, wgt_dat=scaling of N's
//  RETURN_ARRAYS_INCREMENT();
//  dvariable LkvalTmp;
//  LkvalTmp=0.0;
//  for (int ii=1; ii<=ncomp; ii++)
//  {if (nsamp(ii)>=minSS)
//    {LkvalTmp-=wgt_dat*nsamp(ii)*sum(elem_prod((obs_comp(ii)+0.00001),
//               log(elem_div((pred_comp(ii)+0.00001), (obs_comp(ii)+0.00001)))));
//    }
//  }  
//  RETURN_ARRAYS_DECREMENT();
//  return LkvalTmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: multinomial
FUNCTION dvariable lk_robust_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& wgt_dat)
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  LkvalTmp=0.0;
  dvar_matrix Eprime=elem_prod((1.0-obs_comp), obs_comp)+0.1/mbin; //E' of Francis 2011, p.1131  
  dvar_vector nsamp_wgt=nsamp*wgt_dat;
  //cout<<nsamp_wgt<<endl;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {LkvalTmp+= sum(0.5*log(Eprime(ii))-log(0.00001+mfexp(elem_div((-square(obs_comp(ii)-pred_comp(ii))) , (Eprime(ii)*2.0/nsamp_wgt(ii)) ))) );
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: priors
FUNCTION  dvariable neg_log_prior(dvariable pred, const double& prior, dvariable var, int pdf)
  //prior=prior point estimate, var=variance (if negative, treated as CV in arithmetic space), pred=predicted value, pdf=prior type (1=none, 2=lognormal, 3=normal, 4=beta)
    dvariable LkvalTmp;
    dvariable alpha, beta, ab_iq;
    LkvalTmp=0.0;
    // compute generic pdf's
    switch(pdf) {
        case 1: //option to turn off prior
          LkvalTmp=0.0;
          break;
        case 2: // lognormal
          if(prior<=0.0) cout << "YIKES: Don't use a lognormal distn for a negative prior" << endl;
          else if(pred<=0) LkvalTmp=huge_number;
          else {
            if(var<0.0) var=log(1.0+var*var) ;      // convert cv to variance on log scale
            LkvalTmp= 0.5*( square(log(pred/prior))/var + log(var) );
          }
        break;
        case 3: // normal
          if(var<0.0 && prior!=0.0) var=square(var*prior);       // convert cv to variance on observation scale
          else if(var<0.0 && prior==0.0) var=-var;               // cv not really appropriate if prior value equals zero
          LkvalTmp= 0.5*( square(pred-prior)/var + log(var) );
          break;
        case 4: // beta
          if(var<0.0) var=square(var*prior);          // convert cv to variance on observation scale
          if(prior<=0.0 || prior>=1.0) cout << "YIKES: Don't use a beta distn for a prior outside (0,1)" << endl;
          ab_iq=prior*(1.0-prior)/var - 1.0; alpha=prior*ab_iq; beta=(1.0-prior)*ab_iq;
          if(pred>=0 && pred<=1) LkvalTmp= (1.0-alpha)*log(pred)+(1.0-beta)*log(1.0-pred)-gammln(alpha+beta)+gammln(alpha)+gammln(beta);
          else LkvalTmp=huge_number;
          break;
        default: // no such prior pdf currently available
          cout << "The prior must be either 1(lognormal), 2(normal), or 3(beta)." << endl;
          cout << "Presently it is " << pdf << endl;
          exit(0);
    }
    return LkvalTmp;
 


  
REPORT_SECTION
//  cout<<"start report"<<endl;
  get_sel_weighted_current();
//  cout<<"got sel weighted"<<endl;
  get_msy();
//  cout<<"got msy"<<endl;
  get_miscellaneous_stuff();
//  cout<<"got misc stuff"<<endl;
  get_per_recruit_stuff();
//  cout<<"got per recruit"<<endl;     
  cout << "BC Fmsy=" << F_msy_out<< "   BC SSBmsy=" << SSB_msy_out <<endl;
  cout<<"Pop status="<<SSB(endyr)/SSB_msy_out<<endl;
  cout<<"SSB last year "<<SSB(endyr);
  cout<<"SSB msy "<<SSB_msy_out<<endl;
  cout << "var_rec_resid="<<var_rec_dev<<endl;
////		cout << "x_dum="<<x_dum<<endl; 

  report << "TotalLikelihood " << fval << endl;
  report<<" "<<endl;

  report << "Bias-corrected (BC) MSY stuff" << endl;
  report << "BC Fmsy " << F_msy_out << endl;
  report << "BC Emsy " << E_msy_out << endl;
  report << "BC SSBmsy " << SSB_msy_out << endl;
  report << "BC Rmsy " << R_msy_out << endl;  
  report << "BC Bmsy " << B_msy_out << endl;  
  report << "BC MSY " << msy_out << endl;
  report << "BC F/Fmsy " << fullF/F_msy_out << endl;
  report << "BC E/Emsy " << E/E_msy_out << endl;  
  report << "BC SSB/SSBmsy " << SSB/SSB_msy_out << endl;
  report << "BC B/Bmsy " << totB/B_msy_out << endl;   
  report << "BC Yield/MSY " << L_total_yr/msy_out <<endl; 
  report << "BC F(2010)/Fmsy " << fullF(endyr)/F_msy_out << endl;
  report << "BC E(2010)/Emsy " << E(endyr)/E_msy_out << endl;
  report << "BC SSB(2010)/SSBmsy " << SSB(endyr)/SSB_msy_out << endl; 
  report << "BC Predicted Landings(2010)/MSY " << L_total_yr(endyr)/msy_out <<endl;   
  report  << " "<<endl;

  report << "Mortality and growth" << endl;
  report << "M "<<M<<endl;
  //report << "Linf,males="<<Linf_m << "   K, males=" <<K_m<<"   t0_males="<< t0_m<<endl;
  report << "mean length, females " << meanlen_f << endl;
  //report << "Linf,males="<<Linf_f << "   K, males=" <<K_f<<"   t0_fales="<< t0_f<<endl;
  report << "mean length, males " << meanlen_m << endl;
  report << "cv length " << len_cv << endl;
  report << "wgt, males" << wgt_m << endl;
  report << "wgt, females" << wgt_f << endl;
  report<<" "<<endl;

  report << "Stock-Recruit " << endl;
  report << "R0= " << R0 << endl;
  report << "Steepness= " << steep << endl;
  report << "spr_F0= " << spr_F0 << endl;
  report << "Recruits(R) " << rec << endl;
  report << "VirginSSB " << S0 << endl;
  report << "SSB(styr)/VirginSSB " << S1S0 << endl;
  report << "SSB(2010)/VirginSSB " << popstatus << endl;
  report << "SSB " << SSB << endl;
  report << "Biomass " << totB << endl;  
  report << "log recruit deviations (styr_rec_dev-2011)   " << log_rec_dev(styr_rec_dev,2011) <<endl;
  report << "variance of log rec dev (select yrs)  "<<var_rec_dev<<endl;
  report << "autocorrelation " << R_autocorr <<endl;
  report<<" "<<endl;

  report << "Exploitation rate (1958-2007)" << endl;
 report << E << endl;
 report << "Fully-selected F (1958-2007)" << endl;
 report << fullF << endl;
 report << "Comm hand lines F" << endl;
 report << F_HL_out << endl;
 report << "Poundnet F" << endl;
 report << F_PN_out << endl;
 report << "Gillnet F" << endl;
 report << F_GN_out << endl;
 report << "Castnet F" << endl;
 report << F_CN_out << endl;
 report << "MRFSS F" << endl;
 report << F_MRFSS_out << endl;
 report << "Shrimp Bycatch F"<<endl;
 report << F_shrimp_out << endl;
 
 report << "selpar_L50_HL_keep" <<endl;
 report << selpar_L50_HL_keep <<endl;
 report << "selpar_slope_HL" <<endl;
 report << selpar_slope_HL <<endl;
 report << "selpar_L501_PN" <<endl;
 report << selpar_L501_PN <<endl;
 report << "selpar_slope1_PN" <<endl;
 report << selpar_slope1_PN <<endl;
 report << "selpar_L50_GN_keep" <<endl;
 report << selpar_L50_GN_keep <<endl;
 report << "selpar_slope_GN" <<endl;
 report << selpar_slope_GN <<endl;
 report << "selpar_L501_CN" <<endl;
 report << selpar_L501_CN <<endl;
 report << "selpar_slope1_CN" <<endl;
 report << selpar_slope1_CN <<endl;   
 report << "selpar_L502_CN" <<endl;
 report << selpar_L502_CN <<endl;
 report << "selpar_slope2_CN" <<endl;
 report << selpar_slope2_CN <<endl; 
 report << "selpar_L50_MRFSS_keep" <<endl;
 report << selpar_L501_MRFSS_keep <<endl;
 report << "selpar_slope_MRFSS_keep" <<endl;
 report << selpar_slope1_MRFSS_keep <<endl;
 report << "Hand lines selectivity - females" << endl;
 report << sel_HL_keep_F << endl;
 report << "Hand lines selectivity - males" << endl;
 report << sel_HL_keep_M << endl;
 report << "Poundnet selectivity - females" << endl;
 report << sel_PN_F << endl;
 report << "Poundnet selectivity - males" << endl;
 report << sel_PN_M << endl; 
 report << "Gillnet selectivity - pre-1995- females" << endl;
 report << sel_GN_keep_F << endl;
 report << "Gillnet selectivity - pre-1995- males" << endl;
 report << sel_GN_keep_M << endl;
 report << "Gillnet selectivity - post-1995- females" << endl;
// report << sel_GN_keep_F2 << endl;
 report << "Gillnet selectivity - post-1995- males" << endl;
// report << sel_GN_keep_M2 << endl;
 report << "Castnet selectivity - females" << endl;
 report << sel_CN_F << endl;
 report << "Castnet selectivity - males" << endl;
 report << sel_CN_M << endl;
 report << "MRFSS selectivity - females" << endl;
 report << sel_MRFSS_keep_F << endl;
 report << "MRFSS selectivity - males" << endl;
 report << sel_MRFSS_keep_M << endl;
 
 report << "mean log F - HL" << endl;
 report << log_avg_F_HL << endl;
 report << "log F deviations - HL" << endl;
 report << log_F_dev_HL << endl;
 report << "mean log F - PN" << endl;
 report << log_avg_F_PN << endl;
 report << "log F deviations - PN" << endl;
 report << log_F_dev_PN << endl;
 report << "mean log F - GN" << endl;
 report << log_avg_F_GN << endl;
 report << "log F deviations - GN" << endl;
 report << log_F_dev_GN << endl;  
 report << "mean log F - CN" << endl;
 report << log_avg_F_CN << endl;
 report << "log F deviations - CN" << endl;
 report << log_F_dev_CN << endl; 
 report << "mean log F - MRFSS" << endl;
 report << log_avg_F_MRFSS << endl;
 report << "log F deviations - MRFSS" << endl;
 report << log_F_dev_MRFSS << endl; 
 report << "mean log F - shrimp" << endl;
 report << log_avg_F_shrimp << endl;
 report << "log F deviations - shrimp" << endl;
 report << log_F_dev_shrimp << endl; 
 
 report << "mean log F - MRFSS - Discards" << endl;
 report << log_avg_F_MRFSS_D << endl;
 report << "log F deviations - MRFSS - Discards" << endl;
 report << log_F_dev_MRFSS_D << endl; 
 
 report << "Obs  FL_HL U "<<obs_FL_HL_cpue << endl;
 report << "pred FL_HL U "<<pred_FL_HL_cpue << endl;
 report << "Obs  MRFSS U "<<obs_MRFSS_cpue << endl;
 report << "pred MRFSS U "<<pred_MRFSS_cpue << endl;
 report << "Obs  SMAP_YOY U "<<obs_SMAP_YOY_cpue << endl;
 report << "pred SMAP_YOY U "<<pred_SMAP_YOY_cpue << endl;
 //report << "Obs  SMAP_1YR U "<<obs_SMAP_1YR_cpue << endl;
 //report << "pred SMAP_1YR U "<<pred_SMAP_1YR_cpue << endl;

 report << "Obs  HL landings (1000 lb) "<<obs_HL_L << endl;
 report << "pred HL landings (1000 lb) "<<pred_HL_L << endl;
 report << "Obs  PN landings (1000 lb) "<<obs_PN_L << endl;
 report << "pred PN landings (1000 lb) "<<pred_PN_L << endl; 
 report << "Obs  GN landings (1000 lb) "<<obs_GN_L << endl;
 report << "pred GN landings (1000 lb) "<<pred_GN_L << endl;
 report << "Obs  CN landings (1000 lb) "<<obs_CN_L << endl;
 report << "pred CN landings (1000 lb) "<<pred_CN_L << endl;
 report << "Obs  MRFSS landings (1000's) "<<obs_MRFSS_L << endl;
 report << "pred MRFSS landings (1000's) "<<pred_MRFSS_L << endl;
 report << "Obs shrimp bycatch (1000's) "<<obs_shrimp_B<<endl;
 report << "pred shrimp bycatch (1000's) "<<pred_shrimp_B<<endl;

 report << "Obs MRFSS discards (1000's) "<<obs_MRFSS_D<<endl;
 report << "pred MRFSS discards (1000's) "<<pred_MRFSS_D<<endl;
 #include "sm_make_Robject14.cxx"   // write the S-compatible report

  
