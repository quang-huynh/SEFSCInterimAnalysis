//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##
//##  SEDAR56 Standard Assessment: Black Sea Bass, 2017 - 8
//##
//##  NMFS, Beaufort Lab, Sustainable Fisheries Branch
//##
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

DATA_SECTION

!!cout << "Starting Beaufort Assessment Model" << endl;
!!cout << endl;
!!cout << "                BAM!" << endl;
!!cout << endl;

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: set-up section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

// Starting and ending year of the model (year data starts)
init_int styr;
init_int endyr;
//!!cout << styr << endl;
//!!cout << "endyr" << endyr << endl; 
//Starting year to estimate recruitment deviation from S-R curve
init_int styr_rec_dev;
//possible 3 phases of constraints on recruitment deviations
init_int endyr_rec_phase1;
init_int endyr_rec_phase2;

//Reg blocks -- 3 possible periods of comm size regs: styr-83 no restrictions, 1984-98 8-inch TL, 1999-2010 10-in TL
//           -- 4 possible periods of recr size regs: styr-83 no restrictions, 1984-98 8-inch TL, 1999-2006 10-in TL, 2007-2010 12-in TL 
init_int endyr_period1;
init_int endyr_period2;
init_int endyr_recr_period3;
init_int endyr_period4;
//first yr commercial fisheries were closed due to quotas
init_int styr_comm_closed;  
//size limits
init_number sizelim1;//limit_8in;   //8 inch limit in mm
init_number sizelim2//limit_10in;  //10 inch limit in mm
init_number sizelim3//limit_11in;  //11 inch limit in mm
init_number sizelim4//limit_12in;  //12 inch limit in mm
init_number sizelim5//limit_13;    //13 inch limit in mm
init_number limit_disc; //max size applied to discards in block one, prior to fed regs

//Total number of ages
init_int nages;

// Vector of ages for age bins
init_vector agebins(1,nages);

//number assessment years
number nyrs;
number nyrs_rec;
//this section MUST BE INDENTED!!!
 LOCAL_CALCS
   nyrs=endyr-styr+1.;
   //nyrs_rec=endyr-styr_rec_dev+1.;
   nyrs_rec=endyr_rec_phase2-styr_rec_dev+1.;
 END_CALCS

//Total number of length bins for each matrix and length bins used to compute mass in largest bin (plus group)
init_int nlenbins;       //used to match data
init_int nlenbins_plus;  //used to compute density of largest bin (plus group)
init_number lenbins_width;  //width of length bins (mm)

//Vector of lengths for length bins (mm)(midpoint) and bins used in computation of plus group
init_ivector lenbins(1,nlenbins);
init_ivector lenbins_plus(1,nlenbins_plus);
int nlenbins_all;    //largest size class used to compute average lengths and weights

//this section MUST BE INDENTED!!!
 LOCAL_CALCS
   nlenbins_all=nlenbins+nlenbins_plus;   
 END_CALCS
  

//Max F used in spr and msy calcs
init_number max_F_spr_msy;
//Total number of iterations for spr calcs
init_int n_iter_spr;
//Total number of iterations for msy calcs
init_int n_iter_msy;
 
 LOCAL_CALCS
		n_iter_msy=n_iter_spr; 
 END_CALCS

//Number years at end of time series over which to average sector F's, for weighted selectivities
init_int selpar_n_yrs_wgted;
//bias correction (set to 1.0 for no bias correction or a negative value to compute from rec variance)
init_number set_BiasCor;
//exclude these years from end of time series for computing bias correction
init_number BiasCor_exclude_yrs;

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: observed data section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//#############################################################################
//###MARMAP bft/fst#################
//CPUE
init_int styr_Mbft_cpue;
init_int endyr_Mbft_cpue;
init_vector obs_Mbft_cpue(styr_Mbft_cpue,endyr_Mbft_cpue);   //Observed CPUE
init_vector Mbft_cpue_cv(styr_Mbft_cpue,endyr_Mbft_cpue);    //CV of cpue

// Length Compositions (1 cm bins)
init_int nyr_Mbft_lenc;
init_ivector yrs_Mbft_lenc(1,nyr_Mbft_lenc);
init_vector nsamp_Mbft_lenc(1,nyr_Mbft_lenc);
init_vector nfish_Mbft_lenc(1,nyr_Mbft_lenc);
init_matrix obs_Mbft_lenc(1,nyr_Mbft_lenc,1,nlenbins);

// Age Compositions 
init_int nyr_Mbft_agec;
init_ivector yrs_Mbft_agec(1,nyr_Mbft_agec);
init_vector nsamp_Mbft_agec(1,nyr_Mbft_agec);
init_vector nfish_Mbft_agec(1,nyr_Mbft_agec);
init_matrix obs_Mbft_agec(1,nyr_Mbft_agec,1,nages);

//#############################################################################
//###SERFS Video#################
//CPUE
//init_int styr_Vid_cpue;
//init_int endyr_Vid_cpue;
//init_vector obs_Vid_cpue(styr_Vid_cpue,endyr_Vid_cpue);   //Observed CPUE
//init_vector Vid_cpue_cv(styr_Vid_cpue,endyr_Vid_cpue);    //CV of cpue

//#############################################################################
//###MARMAP cvt/CVID#################
//CPUE
init_int styr_Mcvt_cpue;
init_int endyr_Mcvt_cpue;
init_vector obs_Mcvt_cpue(styr_Mcvt_cpue,endyr_Mcvt_cpue);   //Observed CPUE
init_vector Mcvt_cpue_cv(styr_Mcvt_cpue,endyr_Mcvt_cpue);    //CV of cpue
// !!cout << "CVT CPUE" << obs_Mcvt_cpue << endl;   
// Age Compositions 
init_int nyr_Mcvt_agec;
init_ivector yrs_Mcvt_agec(1,nyr_Mcvt_agec);
init_vector nsamp_Mcvt_agec(1,nyr_Mcvt_agec);
init_vector nfish_Mcvt_agec(1,nyr_Mcvt_agec);
init_matrix obs_Mcvt_agec(1,nyr_Mcvt_agec,1,nages);
//  !!cout << "CVT age comps" << obs_Mcvt_agec << endl; 
//#############################################################################
//###################Commercial Hook and Line fishery #########################
//CPUE
init_int styr_cL_cpue;                                             
init_int endyr_cL_cpue;                                            
init_vector obs_cL_cpue(styr_cL_cpue,endyr_cL_cpue);//Observed CPUE
init_vector cL_cpue_cv(styr_cL_cpue,endyr_cL_cpue); //CV of cpue
 
// Landings  (1000 lb whole weight)
init_int styr_cL_L;
init_int endyr_cL_L;
init_vector obs_cL_L(styr_cL_L,endyr_cL_L);   //vector of observed landings by year 
init_vector cL_L_cv(styr_cL_L,endyr_cL_L);    //vector of CV of landings by year

// Discards (1000 fish)
init_int styr_cL_D;
init_int endyr_cL_D;
init_vector obs_cL_released(styr_cL_D,endyr_cL_D); //vector of observed releases by year, gets multiplied by discard mortality for fitting
init_vector cL_D_cv(styr_cL_D,endyr_cL_D);         //vector of CV of discards by year
// Discards (1000 fish) during closed season
init_int styr_cL_closed_D;
init_int endyr_cL_closed_D;
init_vector obs_cL_closed_released(styr_cL_closed_D,endyr_cL_closed_D); //vector of observed releases by year, gets multiplied by discard mortality for fitting

// Length Compositions (1 cm bins)
init_int nyr_cL_lenc;
//!!cout << "nyr cL lenc" << nyr_cL_lenc << endl;
init_ivector yrs_cL_lenc(1,nyr_cL_lenc);
//!!cout << "yrs_cL_lenc" << yrs_cL_lenc << endl;
init_vector nsamp_cL_lenc(1,nyr_cL_lenc);
//!!cout << "nsamp_cL_lenc" << nsamp_cL_lenc << endl;
init_vector nfish_cL_lenc(1,nyr_cL_lenc);
//!!cout << "nfish_cL_lenc" << nfish_cL_lenc << endl;
init_matrix obs_cL_lenc(1,nyr_cL_lenc,1,nlenbins);
//!!cout << "obs_cL_lenc" << obs_cL_lenc << endl;
// Age Compositions 
init_int nyr_cL_agec;
//!!cout << "nyr cL agec" << nyr_cL_agec << endl;
init_ivector yrs_cL_agec(1,nyr_cL_agec);
init_vector nsamp_cL_agec(1,nyr_cL_agec);
init_vector nfish_cL_agec(1,nyr_cL_agec);
init_matrix obs_cL_agec(1,nyr_cL_agec,1,nages);


//#############################################################################
//##Commercial pot (+ other) fleet 
// Landings (1000 lb whole weight)
init_int styr_cP_L;
init_int endyr_cP_L;
init_vector obs_cP_L(styr_cP_L,endyr_cP_L);
init_vector cP_L_cv(styr_cP_L,endyr_cP_L);    //vector of CV of landings by year

// Discards (1000 fish)
init_int styr_cP_D;
init_int endyr_cP_D;
init_vector obs_cP_released(styr_cP_D,endyr_cP_D); //vector of observed releases by year, gets multiplied by discard mortality for fitting
init_vector cP_D_cv(styr_cP_D,endyr_cP_D);         //vector of CV of discards by year
// Discards (1000 fish) during closed season
init_int styr_cP_closed_D;
init_int endyr_cP_closed_D;
init_vector obs_cP_closed_released(styr_cP_closed_D,endyr_cP_closed_D); //vector of observed releases by year, gets multiplied by discard mortality for fitting

// Length Compositions (1 cm bins)
init_int nyr_cP_lenc;
init_ivector yrs_cP_lenc(1,nyr_cP_lenc);
init_vector nsamp_cP_lenc(1,nyr_cP_lenc);
init_vector nfish_cP_lenc(1,nyr_cP_lenc);
init_matrix obs_cP_lenc(1,nyr_cP_lenc,1,nlenbins);
//!!cout << "nyr cP lenc" << nyr_cP_lenc << endl;
init_int nyr_cP_lenc_pool;     //years and weights to pool predicted cP length comps to match pooled observations
init_ivector yrs_cP_lenc_pool(1,nyr_cP_lenc_pool);
init_vector nsamp_cP_lenc_pool(1,nyr_cP_lenc_pool);

// Age Compositions
init_int nyr_cP_agec;
init_ivector yrs_cP_agec(1,nyr_cP_agec);
init_vector nsamp_cP_agec(1,nyr_cP_agec);
init_vector nfish_cP_agec(1,nyr_cP_agec);
init_matrix obs_cP_agec(1,nyr_cP_agec,1,nages);

//#############################################################################
//#############################################################################
//##Commercial Trawl fleet 
// Landings (1000 lb whole weight)
init_int styr_cT_L;
init_int endyr_cT_L;
init_vector obs_cT_L(styr_cT_L,endyr_cT_L);
init_vector cT_L_cv(styr_cT_L,endyr_cT_L);    //vector of CV of landings by year


//#############################################################################
//################################Headboat fleet ########################################
//CPUE
init_int styr_HB_cpue;
init_int endyr_HB_cpue;
init_vector obs_HB_cpue(styr_HB_cpue,endyr_HB_cpue);//Observed CPUE
init_vector HB_cpue_cv(styr_HB_cpue,endyr_HB_cpue); //CV of cpue
//###HBD index (headboat discards from at sea observer program#################
//init_int styr_HBD_cpue;
//init_int endyr_HBD_cpue;
//init_vector obs_HBD_cpue(styr_HBD_cpue,endyr_HBD_cpue);   //Observed CPUE
//init_vector HBD_cpue_cv(styr_HBD_cpue,endyr_HBD_cpue);    //CV of cpue
// Landings (1000 lb)
init_int styr_HB_L;
init_int endyr_HB_L;
init_vector obs_HB_L(styr_HB_L,endyr_HB_L);
init_vector HB_L_cv(styr_HB_L,endyr_HB_L);
// Discards (1000s)
init_int styr_HB_D;
init_int endyr_HB_D;
init_vector obs_HB_released(styr_HB_D,endyr_HB_D);  //vector of observed releases by year, multiplied by discard mortality for fitting  
init_vector HB_D_cv(styr_HB_D,endyr_HB_D);          //vector of CV of discards by year
// Length Compositions (1 cm bins) of landings
init_int nyr_HB_lenc;
init_ivector yrs_HB_lenc(1,nyr_HB_lenc);
init_vector nsamp_HB_lenc(1,nyr_HB_lenc);
init_vector nfish_HB_lenc(1,nyr_HB_lenc);
init_matrix obs_HB_lenc(1,nyr_HB_lenc,1,nlenbins);
//!!cout << "nyr HB lenc" << nyr_HB_lenc << endl;
// Age compositions of landings
init_int nyr_HB_agec;
init_ivector yrs_HB_agec(1,nyr_HB_agec);
init_vector nsamp_HB_agec(1,nyr_HB_agec);
init_vector nfish_HB_agec(1,nyr_HB_agec);
init_matrix obs_HB_agec(1,nyr_HB_agec,1,nages);
//Length Compositions (1 cm bins) of HB discards
init_int nyr_HB_D_lenc;
init_ivector yrs_HB_D_lenc(1,nyr_HB_D_lenc);
init_vector nsamp_HB_D_lenc(1,nyr_HB_D_lenc);
init_vector nfish_HB_D_lenc(1,nyr_HB_D_lenc);
init_matrix obs_HB_D_lenc(1,nyr_HB_D_lenc,1,nlenbins);
//!!cout << "nyr HB_D lenc" << nyr_HB_D_lenc << endl;

//#############################################################################
//############################mrip recreational fleet #################################
// Landings (1000 lb)
init_int styr_mrip_L;
init_int endyr_mrip_L;
init_vector obs_mrip_L(styr_mrip_L,endyr_mrip_L);
init_vector mrip_L_cv(styr_mrip_L,endyr_mrip_L);
// Discards (1000s)
init_int styr_mrip_D;
init_int endyr_mrip_D;
init_vector obs_mrip_released(styr_mrip_D,endyr_mrip_D); //vector of observed releases by year, multiplied by discard mortality for fitting 
init_vector mrip_D_cv(styr_mrip_D,endyr_mrip_D);         //vector of CV of discards by year
// Length Compositions (1 cm bins)
init_int nyr_mrip_lenc;
init_ivector yrs_mrip_lenc(1,nyr_mrip_lenc);
init_vector nsamp_mrip_lenc(1,nyr_mrip_lenc);
init_vector nfish_mrip_lenc(1,nyr_mrip_lenc);
init_matrix obs_mrip_lenc(1,nyr_mrip_lenc,1,nlenbins);
//init_int nyr_mrip_lenc_pool;     //years and weights to pool predicted mrip length comps to match pooled observations
//init_ivector yrs_mrip_lenc_pool(1,nyr_mrip_lenc_pool);
//init_vector nsamp_mrip_lenc_pool(1,nyr_mrip_lenc_pool);
// Age Compositions 
//init_int nyr_mrip_agec;
//init_ivector yrs_mrip_agec(1,nyr_mrip_agec);
//init_vector nsamp_mrip_agec(1,nyr_mrip_agec);
//init_vector nfish_mrip_agec(1,nyr_mrip_agec);
//init_matrix obs_mrip_agec(1,nyr_mrip_agec,1,nages);

//Discard mortality constants
init_number set_Dmort_HL;   //handline (commercial)
init_number set_Dmort_HB_HL; //headboat-specific hook and line
init_number set_Dmort_GR_HL; //charterboat and private hook and line
init_number set_Dmort_cP1;  //pots 1.5 inch panel
init_number set_Dmort_cP2;  //pots 2.0 inch panel
//!!cout <<"discard mortality for pots with 2 inch panel"<<set_Dmort_cP2<< endl;
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: parameter section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##################Single Parameter values and initial guesses #################################

// Von Bert parameters in TL mm 
init_vector set_Linf(1,7);
init_vector set_K(1,7);
init_vector set_t0(1,7);
init_vector set_len_cv(1,7);
//init_number set_len_cv_se(1,7);
//Scalar used only for computing MSST.
init_vector set_M_constant(1,7);     //age-independent: used only for MSST and to scale age dependent M, prior if M is estimated

//Standard errors of von bert params
//init_number set_Linf_se;
//init_number set_K_se;
//init_number set_t0_se;
//init_number set_len_cv_se;

//Spawner-recruit parameters (Initial guesses or fixed values)
init_vector set_steep(1,7);         //recruitment steepness
//init_number set_steep_se;      //SE of recruitment steepness
//init_int steep_prior_pdf;      //(1=none, 2=lognormal, 3=normal, 4=beta)
init_vector set_log_R0(1,7);        //recruitment R0
init_vector set_R_autocorr(1,7);    //recruitment autocorrelation
init_vector set_rec_sigma(1,7);     //recruitment standard deviation in log space
//init_number set_rec_sigma_se;  //SE of recruitment standard deviation in log space
//init_int rec_sigma_prior_pdf;      //(1=none, 2=lognormal, 3=normal, 4=beta)

init_vector set_log_dm_Mbft_lc(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_Mcvt_lc(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_cL_lc(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_cP_lc(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_HB_lc(1,7);  //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_HB_D_lc(1,7);   //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_mrip_lc(1,7);   //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_Mbft_ac(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_Mcvt_ac(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_cL_ac(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_cP_ac(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_HB_ac(1,7);   //Dirichlet-multinomial overdispersion parameter
//init_vector set_log_dm_mrip_ac(1,7);   //Dirichlet-multinomial overdispersion parameter

//Initial guesses or fixed values of estimated selectivity parameters

init_vector set_selpar_A50_Mbft(1,7);
init_vector set_selpar_slope_Mbft(1,7);

init_vector set_selpar_A50_Mcvt(1,7);
init_vector set_selpar_slope_Mcvt(1,7);

init_vector set_selpar_A50_Vid(1,7);
init_vector set_selpar_slope_Vid(1,7);

init_vector set_selpar_A50_cL2(1,7);
init_vector set_selpar_slope_cL2(1,7);
init_vector set_selpar_A50_cL3(1,7);
init_vector set_selpar_slope_cL3(1,7);
init_vector set_selpar_A50_cL4(1,7);
init_vector set_selpar_slope_cL4(1,7);

init_vector set_selpar_A50_cP2(1,7);
init_vector set_selpar_slope_cP2(1,7);
init_vector set_selpar_A50_cP3(1,7);
init_vector set_selpar_slope_cP3(1,7);
init_vector set_selpar_A50_cP4(1,7);
init_vector set_selpar_slope_cP4(1,7);

init_vector set_selpar_A50_HB1(1,7);
init_vector set_selpar_slope_HB1(1,7);
init_vector set_selpar_A50_HB2(1,7);
init_vector set_selpar_slope_HB2(1,7);
init_vector set_selpar_A50_HB3(1,7);
init_vector set_selpar_slope_HB3(1,7);
init_vector set_selpar_A50_HB4(1,7);
init_vector set_selpar_slope_HB4(1,7);
init_vector set_selpar_A50_HB5(1,7);
init_vector set_selpar_slope_HB5(1,7);

init_vector set_selpar_A50_mrip1(1,7);
init_vector set_selpar_slope_mrip1(1,7);
init_vector set_selpar_A50_mrip2(1,7);
init_vector set_selpar_slope_mrip2(1,7);
init_vector set_selpar_A50_mrip3(1,7);
init_vector set_selpar_slope_mrip3(1,7);
init_vector set_selpar_A50_mrip4(1,7);
init_vector set_selpar_slope_mrip4(1,7);
init_vector set_selpar_A50_mrip5(1,7);
init_vector set_selpar_slope_mrip5(1,7);

init_vector set_selpar_Age0_HB_D_logit(1,7);
init_vector set_selpar_Age1_HB_D_logit(1,7);
init_vector set_selpar_Age2_HB_D_logit(1,7);

init_vector set_selpar_A50_HBD4(1,7);
init_vector set_selpar_slope_HBD4(1,7);
init_vector set_selpar_A502_HBD4(1,7);
init_vector set_selpar_slope2_HBD4(1,7);
init_vector set_selpar_A50_HBD5(1,7);
init_vector set_selpar_slope_HBD5(1,7);
init_vector set_selpar_A502_HBD5(1,7);
init_vector set_selpar_slope2_HBD5(1,7);
	!!cout <<set_selpar_slope2_HBD5<<endl;
//--index catchability------------------------------------------------------------------------------------------------------------
init_vector set_logq_Mbft(1,7);    //catchability coefficient (log) for Blackfish trap
init_vector set_logq_Mcvt(1,7);    //catchability coefficient (log) for MARMAP chevron trap
//init_vector set_logq_Vid(1,7);    //catchability coefficient (log) for SERFS video
init_vector set_logq_cL(1,7);      //catchability coefficient (log) for commercial logbook index
init_vector set_logq_HB(1,7);      //catchability coefficient (log) for the headboat index
//init_vector set_logq_HBD(1,7);     //catchability coefficient (log) for HBD

////--mean F's in log space--------------------------------
init_vector set_log_avg_F_cL(1,7);
init_vector set_log_avg_F_cP(1,7);
init_vector set_log_avg_F_cT(1,7);
init_vector set_log_avg_F_HB(1,7);
init_vector set_log_avg_F_mrip(1,7);
////--discard F's-----------------------
init_vector set_log_avg_F_comm_D(1,7);
init_vector set_log_avg_F_HB_D(1,7);
init_vector set_log_avg_F_mrip_D(1,7);

////Dev's-------------------------------------
init_vector set_log_F_dev_cL(1,3);
//!!cout<<"cL F devs"<<set_log_F_dev_cL<<endl;
init_vector set_log_F_dev_cP(1,3);
//!!cout<<"cP F devs"<<set_log_F_dev_cP<<endl;
init_vector set_log_F_dev_cT(1,3);
//!!cout<<"cT F devs"<<set_log_F_dev_cT<<endl;
init_vector set_log_F_dev_HB(1,3);
//!!cout<<"HB F devs"<<set_log_F_dev_HB<<endl;
init_vector set_log_F_dev_mrip(1,3);
//!!cout<<"mrip F devs"<<set_log_F_dev_mrip<<endl;

init_vector set_log_F_dev_comm_D(1,3);
init_vector set_log_F_dev_HB_D(1,3);
init_vector set_log_F_dev_mrip_D(1,3);

init_vector set_log_RWq_dev(1,3);
init_vector set_log_rec_dev(1,3);
init_vector set_log_Nage_dev(1,3);

init_vector set_log_F_dev_cL_vals(styr_cL_L,endyr_cL_L);
//!!cout<<"cL F devs by year"<<set_log_F_dev_cL_vals<<endl;
init_vector set_log_F_dev_cP_vals(styr_cP_L,endyr_cP_L);
//!!cout<<"cP F devs by year"<<set_log_F_dev_cP_vals<<endl;
//init_vector set_log_F_dev_cT_vals(styr_cT_L,endyr_cT_L);
//!!cout<<"cT F devs by year"<<set_log_F_dev_cT_vals<<endl;
init_vector set_log_F_dev_HB_vals(styr_HB_L,endyr_HB_L);
//!!cout<<"HB F devs by year"<<set_log_F_dev_HB_vals<<endl;
init_vector set_log_F_dev_mrip_vals(styr_mrip_L,endyr_mrip_L);
//!!cout<<"mrip F devs by year"<<set_log_F_dev_mrip_vals<<endl;
init_vector set_log_F_dev_comm_D_vals(styr_cL_D,endyr_cL_D);
//!!cout<<"cL_D F devs by year"<<set_log_F_dev_comm_D_vals<<endl;
init_vector set_log_F_dev_HB_D_vals(styr_HB_D,endyr_HB_D);
//!!cout<<"HBD F devs by year"<<set_log_F_dev_HB_D_vals<<endl;
init_vector set_log_F_dev_mrip_D_vals(styr_mrip_D,endyr_mrip_D);
//!!cout<<"mrip_D F devs by year"<<set_log_F_dev_mrip_D_vals<<endl;
init_vector set_log_rec_dev_vals(styr_rec_dev,endyr_rec_phase2);
//!!cout<<"Rec devs by year"<<set_log_rec_dev_vals<<endl;
init_vector set_log_Nage_dev_vals(2,nages);   
//!!cout<<"Nage devs by year"<<set_log_Nage_dev_vals<<endl;

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: likelihood weights section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//--weights for likelihood components-------------------------------------------------------------------------------
init_number set_w_L;
init_number set_w_D;

init_number set_w_lc_Mbft; //Won't need with the dirichlet
init_number set_w_lc_cL;
init_number set_w_lc_cP;
init_number set_w_lc_HB;
init_number set_w_lc_HB_D;
init_number set_w_lc_mrip;

init_number set_w_ac_Mbft;  //Won't need with the dirichlet
init_number set_w_ac_Mcvt;
init_number set_w_ac_cL;
init_number set_w_ac_cP;
init_number set_w_ac_HB;
//init_number set_w_ac_mrip;

init_number set_w_I_Mbft;
init_number set_w_I_Mcvt;
//init_number set_w_I_Vid;
//!!cout << "Video index weight" << set_w_I_Vid<< endl;
init_number set_w_I_cL;
init_number set_w_I_HB;
//init_number set_w_I_HBD;

init_number set_w_rec;             //for fitting S-R curve
init_number set_w_rec_early;       //additional constraint on early years recruitment
init_number set_w_rec_end;         //additional constraint on ending years recruitment 
init_number set_w_fullF;           //penalty for any Fapex>3(removed in final phase of optimization)
init_number set_w_Ftune;           //weight applied to tuning F (removed in final phase of optimization)
//init_number set_w_cvlen_dev;         //penalty on cv deviations at age
//init_number set_w_cvlen_diff;       //penalty on first difference of cv deviations at age

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: miscellaneous stuff section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//TL(mm)-weight(whole weight in kg) relationship: W=aL^b
init_number wgtpar_a;
init_number wgtpar_b;

//weight(whole weight)- fecundity (units=eggs/batch) relationship: log(y)=a+bW
init_number fecpar_a;
init_number fecpar_b;
init_number fecpar_batches;  //number of annual batches may be important if age-dependent, otherwise just a scalar
init_number fecpar_scale;    //used for scaling annual egg production (10^X eggs)

//Female maturity and proportion female at age
init_vector maturity_f_obs(1,nages);            //proportion females mature at age
init_vector prop_f_obs(1,nages);                //proportion female at age
init_number spawn_time_frac; //time of year of peak spawning, as a fraction of the year

init_vector set_M(1,nages);     //age-dependent: used in model
init_number max_obs_age;        //max observed age, used to scale M

//value or initial guess for recreational HB and mrip historic landings multiplicative bias, last yr the bias applies  
//this feature is not currently implemented in the likelihood fcn
init_number set_L_hb_bias;
init_number set_L_mrip_bias;
init_number set_L_comm_bias;
init_number endyr_L_HB_bias;
init_number endyr_L_mrip_bias;
init_number endyr_L_comm_bias;

//rate of increase on q
init_int set_q_rate_phase;  //value sets estimation phase of rate increase, negative value turns it off
init_number set_q_rate;
//density dependence on fishery q's 
init_int set_q_DD_phase;      //value sets estimation phase of random walk, negative value turns it off
init_number set_q_DD_beta;    //value of 0.0 is density indepenent
init_number set_q_DD_beta_se;
init_int set_q_DD_stage;      //age to begin counting biomass, should be near full exploitation

//random walk on fishery q's 
init_int set_q_RW_phase;         //value sets estimation phase of random walk, negative value turns it off
init_number set_q_RW_cL_var;     //assumed variance of RW q
init_number set_q_RW_HB_var;     //assumed variance of RW q
//init_number set_q_RW_HBD_var;    //assumed variance of RW q

//init_vector set_F_init_ratio;  //defines initialization F as a ratio of that from first several yrs of assessment
init_number set_F_init_ratio

//Tune Fapex (tuning removed in final year of optimization)
init_number set_Ftune;
init_int set_Ftune_yr;

//threshold sample sizes for including length comps, age comps, respectively 
init_number minSS_lenc;
init_number minSS_agec;

//maximum allowable annual sample sizes for length comps, age comps, respectively
init_number maxSS_lenc;
init_number maxSS_agec;

//ageing error matrix (columns are true ages, rows are ages as read for age comps: columns should sum to one)
init_matrix age_error(1,nages,1,nages);

//proportion of length comp mass below size limit considered when matching length comp
//note: these need length comp and age comp data to be estimable
init_number set_p_lenc_cL2; 
init_number set_p_lenc_cL3; 
init_number set_p_lenc_cP2;
init_number set_p_lenc_cP3;
init_number set_p_lenc_cT2;
init_number set_p_lenc_cT3;
init_number set_p_lenc_HB2;
init_number set_p_lenc_HB3;
init_number set_p_lenc_HB4;
init_number set_p_lenc_HB5; 
init_number set_p_lenc_mrip2;  
init_number set_p_lenc_mrip3;
init_number set_p_lenc_mrip4;
init_number set_p_lenc_mrip5;

init_number set_p_lenc_comm_D2;
init_number set_p_lenc_comm_D3;
init_number set_p_lenc_comm_D4;
init_number set_p_lenc_HB_D2;
init_number set_p_lenc_HB_D3;
init_number set_p_lenc_HB_D4;
init_number set_p_lenc_HB_D5;
init_number set_p_lenc_mrip_D1;
init_number set_p_lenc_mrip_D2;
init_number set_p_lenc_mrip_D3;
init_number set_p_lenc_mrip_D4;
init_number set_p_lenc_mrip_D5;

// #######Indexing integers for year(iyear), age(iage),length(ilen) ###############
int iyear;
int iage;
int ilen;
int ff;

number sqrt2pi;
number g2mt;                    //conversion of grams to metric tons 
number g2kg;                    //conversion of grams to kg   
number g2klb;                   //conversion of grams to 1000 lb   
number mt2klb;                  //conversion of metric tons to 1000 lb
number mt2lb;                   //conversion of metric tons to lb
number dzero;                   //small additive constant to prevent division by zero
number huge_number;             //huge number, to avoid irregular parameter space
number onehalf;                 //0.5

init_number end_of_data_file;
//this section MUST BE INDENTED!!!
 LOCAL_CALCS
   if(end_of_data_file!=999)
   {
     for(iyear=1; iyear<=1000; iyear++)
     {
       cout << "*** WARNING: Data File NOT READ CORRECTLY ****" << endl;
       cout << "" <<endl;
     }
   }
   else
   {
    cout << "Data File read correctly" << endl;
   } 
 END_CALCS   


//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PARAMETER_SECTION
LOCAL_CALCS
  const double Linf_LO=set_Linf(2); const double Linf_HI=set_Linf(3); const double Linf_PH=set_Linf(4);
  const double K_LO=set_K(2); const double K_HI=set_K(3); const double K_PH=set_K(4);
  const double t0_LO=set_t0(2); const double t0_HI=set_t0(3); const double t0_PH=set_t0(4);  
  const double len_cv_LO=set_len_cv(2); const double len_cv_HI=set_len_cv(3); const double len_cv_PH=set_len_cv(4); 
     
  const double M_constant_LO=set_M_constant(2); const double M_constant_HI=set_M_constant(3); const double M_constant_PH=set_M_constant(4);        
  const double steep_LO=set_steep(2); const double steep_HI=set_steep(3); const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); const double log_R0_HI=set_log_R0(3); const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); const double R_autocorr_HI=set_R_autocorr(3); const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); const double rec_sigma_HI=set_rec_sigma(3); const double rec_sigma_PH=set_rec_sigma(4);
  
  const double log_dm_Mbft_lc_LO=set_log_dm_Mbft_lc(2); const double log_dm_Mbft_lc_HI=set_log_dm_Mbft_lc(3); const double log_dm_Mbft_lc_PH=set_log_dm_Mbft_lc(4);
  const double log_dm_Mcvt_lc_LO=set_log_dm_Mcvt_lc(2); const double log_dm_Mcvt_lc_HI=set_log_dm_Mcvt_lc(3); const double log_dm_Mcvt_lc_PH=set_log_dm_Mcvt_lc(4);
  const double log_dm_HB_lc_LO=set_log_dm_HB_lc(2); const double log_dm_HB_lc_HI=set_log_dm_HB_lc(3); const double log_dm_HB_lc_PH=set_log_dm_HB_lc(4);
  const double log_dm_mrip_lc_LO=set_log_dm_mrip_lc(2); const double log_dm_mrip_lc_HI=set_log_dm_mrip_lc(3); const double log_dm_mrip_lc_PH=set_log_dm_mrip_lc(4);
  const double log_dm_HB_D_lc_LO=set_log_dm_HB_D_lc(2); const double log_dm_HB_D_lc_HI=set_log_dm_HB_D_lc(3); const double log_dm_HB_D_lc_PH=set_log_dm_HB_D_lc(4);
  const double log_dm_cL_lc_LO=set_log_dm_cL_lc(2); const double log_dm_cL_lc_HI=set_log_dm_cL_lc(3); const double log_dm_cL_lc_PH=set_log_dm_cL_lc(4);
  const double log_dm_cP_lc_LO=set_log_dm_cP_lc(2); const double log_dm_cP_lc_HI=set_log_dm_cP_lc(3); const double log_dm_cP_lc_PH=set_log_dm_cP_lc(4);
  const double log_dm_Mbft_ac_LO=set_log_dm_Mbft_ac(2); const double log_dm_Mbft_ac_HI=set_log_dm_Mbft_ac(3); const double log_dm_Mbft_ac_PH=set_log_dm_Mbft_ac(4);
  const double log_dm_Mcvt_ac_LO=set_log_dm_Mcvt_ac(2); const double log_dm_Mcvt_ac_HI=set_log_dm_Mcvt_ac(3); const double log_dm_Mcvt_ac_PH=set_log_dm_Mcvt_ac(4);
  const double log_dm_HB_ac_LO=set_log_dm_HB_ac(2); const double log_dm_HB_ac_HI=set_log_dm_HB_ac(3); const double log_dm_HB_ac_PH=set_log_dm_HB_ac(4);
  //const double log_dm_mrip_ac_LO=set_log_dm_mrip_ac(2); const double log_dm_mrip_ac_HI=set_log_dm_mrip_ac(3); const double log_dm_mrip_ac_PH=set_log_dm_mrip_ac(4);
  const double log_dm_cL_ac_LO=set_log_dm_cL_ac(2); const double log_dm_cL_ac_HI=set_log_dm_cL_ac(3); const double log_dm_cL_ac_PH=set_log_dm_cL_ac(4);
  const double log_dm_cP_ac_LO=set_log_dm_cP_ac(2); const double log_dm_cP_ac_HI=set_log_dm_cP_ac(3); const double log_dm_cP_ac_PH=set_log_dm_cP_ac(4);
  
  const double selpar_A50_Mbft_LO=set_selpar_A50_Mbft(2); const double selpar_A50_Mbft_HI=set_selpar_A50_Mbft(3); const double selpar_A50_Mbft_PH=set_selpar_A50_Mbft(4);
  const double selpar_slope_Mbft_LO=set_selpar_slope_Mbft(2); const double selpar_slope_Mbft_HI=set_selpar_slope_Mbft(3); const double selpar_slope_Mbft_PH=set_selpar_slope_Mbft(4);
  
  const double selpar_A50_Mcvt_LO=set_selpar_A50_Mcvt(2); const double selpar_A50_Mcvt_HI=set_selpar_A50_Mcvt(3); const double selpar_A50_Mcvt_PH=set_selpar_A50_Mcvt(4);
  const double selpar_slope_Mcvt_LO=set_selpar_slope_Mcvt(2); const double selpar_slope_Mcvt_HI=set_selpar_slope_Mcvt(3); const double selpar_slope_Mcvt_PH=set_selpar_slope_Mcvt(4);

  const double selpar_A50_cL2_LO=set_selpar_A50_cL2(2); const double selpar_A50_cL2_HI=set_selpar_A50_cL2(3); const double selpar_A50_cL2_PH=set_selpar_A50_cL2(4);
  const double selpar_slope_cL2_LO=set_selpar_slope_cL2(2); const double selpar_slope_cL2_HI=set_selpar_slope_cL2(3); const double selpar_slope_cL2_PH=set_selpar_slope_cL2(4);
  const double selpar_A50_cL3_LO=set_selpar_A50_cL3(2); const double selpar_A50_cL3_HI=set_selpar_A50_cL3(3); const double selpar_A50_cL3_PH=set_selpar_A50_cL3(4);
  const double selpar_slope_cL3_LO=set_selpar_slope_cL3(2); const double selpar_slope_cL3_HI=set_selpar_slope_cL3(3); const double selpar_slope_cL3_PH=set_selpar_slope_cL3(4);
  const double selpar_A50_cL4_LO=set_selpar_A50_cL4(2); const double selpar_A50_cL4_HI=set_selpar_A50_cL4(3); const double selpar_A50_cL4_PH=set_selpar_A50_cL4(4);
  const double selpar_slope_cL4_LO=set_selpar_slope_cL4(2); const double selpar_slope_cL4_HI=set_selpar_slope_cL4(3); const double selpar_slope_cL4_PH=set_selpar_slope_cL4(4);
  
  const double selpar_A50_cP2_LO=set_selpar_A50_cP2(2); const double selpar_A50_cP2_HI=set_selpar_A50_cP2(3); const double selpar_A50_cP2_PH=set_selpar_A50_cP2(4);
  const double selpar_slope_cP2_LO=set_selpar_slope_cP2(2); const double selpar_slope_cP2_HI=set_selpar_slope_cP2(3); const double selpar_slope_cP2_PH=set_selpar_slope_cP2(4);
  const double selpar_A50_cP3_LO=set_selpar_A50_cP3(2); const double selpar_A50_cP3_HI=set_selpar_A50_cP3(3); const double selpar_A50_cP3_PH=set_selpar_A50_cP3(4);
  const double selpar_slope_cP3_LO=set_selpar_slope_cP3(2); const double selpar_slope_cP3_HI=set_selpar_slope_cP3(3); const double selpar_slope_cP3_PH=set_selpar_slope_cP3(4);
  const double selpar_A50_cP4_LO=set_selpar_A50_cP4(2); const double selpar_A50_cP4_HI=set_selpar_A50_cP4(3); const double selpar_A50_cP4_PH=set_selpar_A50_cP4(4);
  const double selpar_slope_cP4_LO=set_selpar_slope_cP4(2); const double selpar_slope_cP4_HI=set_selpar_slope_cP4(3); const double selpar_slope_cP4_PH=set_selpar_slope_cP4(4);
  
  const double selpar_A50_HB1_LO=set_selpar_A50_HB1(2); const double selpar_A50_HB1_HI=set_selpar_A50_HB1(3); const double selpar_A50_HB1_PH=set_selpar_A50_HB1(4);
  const double selpar_slope_HB1_LO=set_selpar_slope_HB1(2); const double selpar_slope_HB1_HI=set_selpar_slope_HB1(3); const double selpar_slope_HB1_PH=set_selpar_slope_HB1(4);
  const double selpar_A50_HB2_LO=set_selpar_A50_HB2(2); const double selpar_A50_HB2_HI=set_selpar_A50_HB2(3); const double selpar_A50_HB2_PH=set_selpar_A50_HB2(4);
  const double selpar_slope_HB2_LO=set_selpar_slope_HB2(2); const double selpar_slope_HB2_HI=set_selpar_slope_HB2(3); const double selpar_slope_HB2_PH=set_selpar_slope_HB2(4);
  const double selpar_A50_HB3_LO=set_selpar_A50_HB3(2); const double selpar_A50_HB3_HI=set_selpar_A50_HB3(3); const double selpar_A50_HB3_PH=set_selpar_A50_HB3(4);
  const double selpar_slope_HB3_LO=set_selpar_slope_HB3(2); const double selpar_slope_HB3_HI=set_selpar_slope_HB3(3); const double selpar_slope_HB3_PH=set_selpar_slope_HB3(4);
  const double selpar_A50_HB4_LO=set_selpar_A50_HB4(2); const double selpar_A50_HB4_HI=set_selpar_A50_HB4(3); const double selpar_A50_HB4_PH=set_selpar_A50_HB4(4);
  const double selpar_slope_HB4_LO=set_selpar_slope_HB4(2); const double selpar_slope_HB4_HI=set_selpar_slope_HB4(3); const double selpar_slope_HB4_PH=set_selpar_slope_HB4(4);
  const double selpar_A50_HB5_LO=set_selpar_A50_HB5(2); const double selpar_A50_HB5_HI=set_selpar_A50_HB5(3); const double selpar_A50_HB5_PH=set_selpar_A50_HB5(4);
  const double selpar_slope_HB5_LO=set_selpar_slope_HB5(2); const double selpar_slope_HB5_HI=set_selpar_slope_HB5(3); const double selpar_slope_HB5_PH=set_selpar_slope_HB5(4);
  
  const double selpar_A50_mrip1_LO=set_selpar_A50_mrip1(2); const double selpar_A50_mrip1_HI=set_selpar_A50_mrip1(3); const double selpar_A50_mrip1_PH=set_selpar_A50_mrip1(4);
  const double selpar_slope_mrip1_LO=set_selpar_slope_mrip1(2); const double selpar_slope_mrip1_HI=set_selpar_slope_mrip1(3); const double selpar_slope_mrip1_PH=set_selpar_slope_mrip1(4);
  const double selpar_A50_mrip2_LO=set_selpar_A50_mrip2(2); const double selpar_A50_mrip2_HI=set_selpar_A50_mrip2(3); const double selpar_A50_mrip2_PH=set_selpar_A50_mrip2(4);
  const double selpar_slope_mrip2_LO=set_selpar_slope_mrip2(2); const double selpar_slope_mrip2_HI=set_selpar_slope_mrip2(3); const double selpar_slope_mrip2_PH=set_selpar_slope_mrip2(4);
  const double selpar_A50_mrip3_LO=set_selpar_A50_mrip3(2); const double selpar_A50_mrip3_HI=set_selpar_A50_mrip3(3); const double selpar_A50_mrip3_PH=set_selpar_A50_mrip3(4);
  const double selpar_slope_mrip3_LO=set_selpar_slope_mrip3(2); const double selpar_slope_mrip3_HI=set_selpar_slope_mrip3(3); const double selpar_slope_mrip3_PH=set_selpar_slope_mrip3(4);
  const double selpar_A50_mrip4_LO=set_selpar_A50_mrip4(2); const double selpar_A50_mrip4_HI=set_selpar_A50_mrip4(3); const double selpar_A50_mrip4_PH=set_selpar_A50_mrip4(4);
  const double selpar_slope_mrip4_LO=set_selpar_slope_mrip4(2); const double selpar_slope_mrip4_HI=set_selpar_slope_mrip4(3); const double selpar_slope_mrip4_PH=set_selpar_slope_mrip4(4);
  const double selpar_A50_mrip5_LO=set_selpar_A50_mrip5(2); const double selpar_A50_mrip5_HI=set_selpar_A50_mrip5(3); const double selpar_A50_mrip5_PH=set_selpar_A50_mrip5(4);
  const double selpar_slope_mrip5_LO=set_selpar_slope_mrip5(2); const double selpar_slope_mrip5_HI=set_selpar_slope_mrip5(3); const double selpar_slope_mrip5_PH=set_selpar_slope_mrip5(4);
    
  const double selpar_Age0_HB_D_logit_LO=set_selpar_Age0_HB_D_logit(2); const double selpar_Age0_HB_D_logit_HI=set_selpar_Age0_HB_D_logit(3); const double selpar_Age0_HB_D_logit_PH=set_selpar_Age0_HB_D_logit(4);
  const double selpar_Age1_HB_D_logit_LO=set_selpar_Age1_HB_D_logit(2); const double selpar_Age1_HB_D_logit_HI=set_selpar_Age1_HB_D_logit(3); const double selpar_Age1_HB_D_logit_PH=set_selpar_Age1_HB_D_logit(4);
  const double selpar_Age2_HB_D_logit_LO=set_selpar_Age2_HB_D_logit(2); const double selpar_Age2_HB_D_logit_HI=set_selpar_Age2_HB_D_logit(3); const double selpar_Age2_HB_D_logit_PH=set_selpar_Age2_HB_D_logit(4);
  
  const double selpar_A50_HBD4_LO=set_selpar_A50_HBD4(2); 
  const double selpar_A50_HBD4_HI=set_selpar_A50_HBD4(3); 
  const double selpar_A50_HBD4_PH=set_selpar_A50_HBD4(4);
  const double selpar_slope_HBD4_LO=set_selpar_slope_HBD4(2); 
  const double selpar_slope_HBD4_HI=set_selpar_slope_HBD4(3); 
  const double selpar_slope_HBD4_PH=set_selpar_slope_HBD4(4);
  const double selpar_A502_HBD4_LO=set_selpar_A502_HBD4(2); 
  const double selpar_A502_HBD4_HI=set_selpar_A502_HBD4(3); 
  const double selpar_A502_HBD4_PH=set_selpar_A502_HBD4(4);
  const double selpar_slope2_HBD4_LO=set_selpar_slope2_HBD4(2); 
  const double selpar_slope2_HBD4_HI=set_selpar_slope2_HBD4(3); 
  const double selpar_slope2_HBD4_PH=set_selpar_slope2_HBD4(4);
  const double selpar_A50_HBD5_LO=set_selpar_A50_HBD5(2); 
  const double selpar_A50_HBD5_HI=set_selpar_A50_HBD5(3); 
  const double selpar_A50_HBD5_PH=set_selpar_A50_HBD5(4);
  const double selpar_slope_HBD5_LO=set_selpar_slope_HBD5(2); 
  const double selpar_slope_HBD5_HI=set_selpar_slope_HBD5(3); 
  const double selpar_slope_HBD5_PH=set_selpar_slope_HBD5(4);
  const double selpar_A502_HBD5_LO=set_selpar_A502_HBD5(2); 
  const double selpar_A502_HBD5_HI=set_selpar_A502_HBD5(3); 
  const double selpar_A502_HBD5_PH=set_selpar_A502_HBD5(4);
  const double selpar_slope2_HBD5_LO=set_selpar_slope2_HBD5(2); 
  const double selpar_slope2_HBD5_HI=set_selpar_slope2_HBD5(3); 
  const double selpar_slope2_HBD5_PH=set_selpar_slope2_HBD5(4);
  
  const double log_q_Mbft_LO=set_logq_Mbft(2); const double log_q_Mbft_HI=set_logq_Mbft(3); const double log_q_Mbft_PH=set_logq_Mbft(4);
  const double log_q_HB_LO=set_logq_HB(2); const double log_q_HB_HI=set_logq_HB(3); const double log_q_HB_PH=set_logq_HB(4);
  const double log_q_Mcvt_LO=set_logq_Mcvt(2); const double log_q_Mcvt_HI=set_logq_Mcvt(3); const double log_q_Mcvt_PH=set_logq_Mcvt(4);
  //const double log_q_Vid_LO=set_logq_Vid(2); const double log_q_Vid_HI=set_logq_Vid(3); const double log_q_Vid_PH=set_logq_Vid(4);
  const double log_q_cL_LO=set_logq_cL(2); const double log_q_cL_HI=set_logq_cL(3); const double log_q_cL_PH=set_logq_cL(4);
  //const double log_q_HBD_LO=set_logq_HBD(2); const double log_q_HBD_HI=set_logq_HBD(3); const double log_q_HBD_PH=set_logq_HBD(4);
  
  const double log_avg_F_cL_LO=set_log_avg_F_cL(2); const double log_avg_F_cL_HI=set_log_avg_F_cL(3); const double log_avg_F_cL_PH=set_log_avg_F_cL(4);
  const double log_avg_F_cP_LO=set_log_avg_F_cP(2); const double log_avg_F_cP_HI=set_log_avg_F_cP(3); const double log_avg_F_cP_PH=set_log_avg_F_cP(4);
  const double log_avg_F_cT_LO=set_log_avg_F_cT(2); const double log_avg_F_cT_HI=set_log_avg_F_cT(3); const double log_avg_F_cT_PH=set_log_avg_F_cT(4);
  const double log_avg_F_HB_LO=set_log_avg_F_HB(2); const double log_avg_F_HB_HI=set_log_avg_F_HB(3); const double log_avg_F_HB_PH=set_log_avg_F_HB(4); 
  const double log_avg_F_mrip_LO=set_log_avg_F_mrip(2); const double log_avg_F_mrip_HI=set_log_avg_F_mrip(3); const double log_avg_F_mrip_PH=set_log_avg_F_mrip(4); 
  
  const double log_avg_F_comm_D_LO=set_log_avg_F_comm_D(2); const double log_avg_F_comm_D_HI=set_log_avg_F_comm_D(3); const double log_avg_F_comm_D_PH=set_log_avg_F_comm_D(4);
  const double log_avg_F_HB_D_LO=set_log_avg_F_HB_D(2); const double log_avg_F_HB_D_HI=set_log_avg_F_HB_D(3); const double log_avg_F_HB_D_PH=set_log_avg_F_HB_D(4); 
  const double log_avg_F_mrip_D_LO=set_log_avg_F_mrip_D(2); const double log_avg_F_mrip_D_HI=set_log_avg_F_mrip_D(3); const double log_avg_F_mrip_D_PH=set_log_avg_F_mrip_D(4); 
  
  //-dev vectors-----------------------------------------------------------------------------------------------------------  
  const double log_F_dev_cL_LO=set_log_F_dev_cL(1); const double log_F_dev_cL_HI=set_log_F_dev_cL(2); const double log_F_dev_cL_PH=set_log_F_dev_cL(3);  
  const double log_F_dev_cP_LO=set_log_F_dev_cP(1); const double log_F_dev_cP_HI=set_log_F_dev_cP(2); const double log_F_dev_cP_PH=set_log_F_dev_cP(3);   
  const double log_F_dev_cT_LO=set_log_F_dev_cT(1); const double log_F_dev_cT_HI=set_log_F_dev_cT(2); const double log_F_dev_cT_PH=set_log_F_dev_cT(3);   
  const double log_F_dev_HB_LO=set_log_F_dev_HB(1); const double log_F_dev_HB_HI=set_log_F_dev_HB(2); const double log_F_dev_HB_PH=set_log_F_dev_HB(3);   
  const double log_F_dev_mrip_LO=set_log_F_dev_mrip(1); const double log_F_dev_mrip_HI=set_log_F_dev_mrip(2); const double log_F_dev_mrip_PH=set_log_F_dev_mrip(3);   
  
  const double log_F_dev_comm_D_LO=set_log_F_dev_comm_D(1); const double log_F_dev_comm_D_HI=set_log_F_dev_comm_D(2); const double log_F_dev_comm_D_PH=set_log_F_dev_comm_D(3);   
  const double log_F_dev_HB_D_LO=set_log_F_dev_HB_D(1); const double log_F_dev_HB_D_HI=set_log_F_dev_HB_D(2); const double log_F_dev_HB_D_PH=set_log_F_dev_HB_D(3);   
  const double log_F_dev_mrip_D_LO=set_log_F_dev_mrip_D(1); const double log_F_dev_mrip_D_HI=set_log_F_dev_mrip_D(2); const double log_F_dev_mrip_D_PH=set_log_F_dev_mrip_D(3);   
  //const double F_init_ratio_LO=set_F_init_ratio(1); const double F_init_ratio_HI=set_F_init_ratio(2); const double F_init_ratio_PH=set_F_init_ratio(3); 

  const double log_RWq_LO=set_log_RWq_dev(1); const double log_RWq_HI=set_log_RWq_dev(2); const double log_RWq_PH=set_log_RWq_dev(3);  
  
  const double log_rec_dev_LO=set_log_rec_dev(1); const double log_rec_dev_HI=set_log_rec_dev(2); const double log_rec_dev_PH=set_log_rec_dev(3);          
  const double log_Nage_dev_LO=set_log_Nage_dev(1); const double log_Nage_dev_HI=set_log_Nage_dev(2); const double log_Nage_dev_PH=set_log_Nage_dev(3);          
  
 END_CALCS
 
////--------------Growth---------------------------------------------------------------------------
 
  //Population growth parms and conversions
  init_bounded_number Linf(Linf_LO,Linf_HI,Linf_PH);
  init_bounded_number K(K_LO,K_HI,K_PH);
  init_bounded_number t0(t0_LO,t0_HI,t0_PH);
  init_bounded_number len_cv_val(len_cv_LO,len_cv_HI,len_cv_PH);  
  vector Linf_out(1,8);
  vector K_out(1,8);
  vector t0_out(1,8);
  vector len_cv_val_out(1,8);
  
  ////init_bounded_number Linf(300,800,3);
  ////init_bounded_number K(0.05,0.5,3);
  ////init_bounded_number t0(-1.5,-0.01,3);
  //number Linf;
  //number K;
  //number t0;  
  //init_bounded_number len_cv_val(0.07,0.3,4);
  ////  init_bounded_dev_vector log_len_cv_dev(1,nages,-2,2,-4)
  ////number len_cv_val;
  //vector len_sd(1,nages);
  //vector len_cv(1,nages);  
  
  vector meanlen_TL(1,nages);   //mean Total length (mm) at age
  
  vector wgt_g(1,nages);        //whole wgt in g
  vector wgt_kg(1,nages);       //whole wgt in kg
  vector wgt_mt(1,nages);       //whole wgt in mt
  vector wgt_klb(1,nages);      //whole wgt in 1000 lb
  vector wgt_lb(1,nages);       //whole wgt in lb
  vector fecundity(1,nages);    //fecundity at age, perhaps scaled  

  matrix len_cL_mm(styr,endyr,1,nages);       //mean length at age of cL landings in mm (may differ from popn mean)    
  matrix wgt_cL_klb(styr,endyr,1,nages);      //whole wgt of cL landings in 1000 lb  
  matrix len_cP_mm(styr,endyr,1,nages);       //mean length at age of cP landings in mm (may differ from popn mean)    
  matrix wgt_cP_klb(styr,endyr,1,nages);      //whole wgt of cP landings in 1000 lb  
  matrix len_cT_mm(styr,endyr,1,nages);       //mean length at age of cT landings in mm (may differ from popn mean)    
  matrix wgt_cT_klb(styr,endyr,1,nages);      //whole wgt of cT landings in 1000 lb  
  matrix len_HB_mm(styr,endyr,1,nages);       //mean length at age of HB landings in mm (may differ from popn mean)    
  matrix wgt_HB_klb(styr,endyr,1,nages);      //whole wgt of HB landings in 1000 lb  
  matrix len_mrip_mm(styr,endyr,1,nages);     //mean length at age of mrip landings in mm (may differ from popn mean)    
  matrix wgt_mrip_klb(styr,endyr,1,nages);    //whole wgt of mrip landings in 1000 lb  

  matrix len_comm_D_mm(styr,endyr,1,nages);    //mean length at age of cL discards in mm (may differ from popn mean)    
  matrix wgt_comm_D_klb(styr,endyr,1,nages);   //whole wgt of cL discards in 1000 lb  
  matrix len_HB_D_mm(styr,endyr,1,nages);      //mean length at age of cL discards in mm (may differ from popn mean)    
  matrix wgt_HB_D_klb(styr,endyr,1,nages);     //whole wgt of cL discards in 1000 lb  
  matrix len_mrip_D_mm(styr,endyr,1,nages);    //mean length at age of cL discards in mm (may differ from popn mean)    
  matrix wgt_mrip_D_klb(styr,endyr,1,nages);   //whole wgt of cL discards in 1000 lb  
 
  matrix lenprob(1,nages,1,nlenbins);           //distn of size at age (age-length key, 1 cm bins) in population
    
  matrix lenprob_plus(1,nages,1,nlenbins_plus); //used to compute mass in last length bin (a plus group)   
  matrix lenprob_all(1,nages,1,nlenbins_all);   //extended lenprob
  vector lenbins_all(1,nlenbins_all);
  
  number zscore_len;                            //standardized normal values used for computing lenprob
  vector cprob_lenvec(1,nlenbins);              //cumulative probabilities used for computing lenprob
  number zscore_lzero;                          //standardized normal values for length = 0
  number cprob_lzero;                           //length probability mass below zero, used for computing lenprob
  
  //matrices below are used to match length comps
  matrix lenprob_Mbft(1,nages,1,nlenbins);    //distn of size at age of Mbft
  matrix lenprob_cL1(1,nages,1,nlenbins);     //distn of size at age in cL block 1
  matrix lenprob_cL2(1,nages,1,nlenbins);     //distn of size at age in cL block 2
  matrix lenprob_cL3(1,nages,1,nlenbins);     //distn of size at age in cL block 3
  matrix lenprob_cP1(1,nages,1,nlenbins);     //distn of size at age in cP block 1
  matrix lenprob_cP2(1,nages,1,nlenbins);     //distn of size at age in cP block 2
  matrix lenprob_cP3(1,nages,1,nlenbins);     //distn of size at age in cP block 3
  matrix lenprob_cT1(1,nages,1,nlenbins);     //distn of size at age in cT block 1
  matrix lenprob_cT2(1,nages,1,nlenbins);     //distn of size at age in cT block 2
  matrix lenprob_HB1(1,nages,1,nlenbins);     //distn of size at age in HB block 1
  matrix lenprob_HB2(1,nages,1,nlenbins);     //distn of size at age in HB block 2
  matrix lenprob_HB3(1,nages,1,nlenbins);     //distn of size at age in HB block 3
  matrix lenprob_HB4(1,nages,1,nlenbins);     //distn of size at age in HB block 4 
  matrix lenprob_HB5(1,nages,1,nlenbins);     //distn of size at age in HB block 5    
  matrix lenprob_mrip1(1,nages,1,nlenbins);   //distn of size at age in mrip block 2
  matrix lenprob_mrip2(1,nages,1,nlenbins);   //distn of size at age in mrip block 2
  matrix lenprob_mrip3(1,nages,1,nlenbins);   //distn of size at age in mrip block 3
  matrix lenprob_mrip4(1,nages,1,nlenbins);   //distn of size at age in mrip block 4
  matrix lenprob_mrip5(1,nages,1,nlenbins);   //distn of size at age in mrip block 5

  matrix lenprob_comm_D2(1,nages,1,nlenbins);    //distn of size at age in comm discards comm block 2   
  matrix lenprob_comm_D3(1,nages,1,nlenbins);    //distn of size at age in comm discards comm block 3
  matrix lenprob_comm_D4(1,nages,1,nlenbins);    //distn of size at age in comm discards comm block 4
  matrix lenprob_HB_D2(1,nages,1,nlenbins);      //distn of size at age in HB discards rec block 2 
  matrix lenprob_HB_D3(1,nages,1,nlenbins);      //distn of size at age in HB discards rec block 3
  matrix lenprob_HB_D4(1,nages,1,nlenbins);      //distn of size at age in HB discards rec block 4
  matrix lenprob_HB_D5(1,nages,1,nlenbins);      //distn of size at age in HB discards rec block 5
  matrix lenprob_mrip_D1(1,nages,1,nlenbins);    //distn of size at age in mrip discards rec block 1
  matrix lenprob_mrip_D2(1,nages,1,nlenbins);    //distn of size at age in mrip discards rec block 2
  matrix lenprob_mrip_D3(1,nages,1,nlenbins);    //distn of size at age in mrip discards rec block 3
  matrix lenprob_mrip_D4(1,nages,1,nlenbins);    //distn of size at age in mrip discards rec block 4
  matrix lenprob_mrip_D5(1,nages,1,nlenbins);    //distn of size at age in mrip discards rec block 5

  //matrices below used to compute mean weights
  matrix lenprob_cL1_all(1,nages,1,nlenbins_all);     //distn of size at age in cL block 1
  matrix lenprob_cL2_all(1,nages,1,nlenbins_all);     //distn of size at age in cL block 2
  matrix lenprob_cL3_all(1,nages,1,nlenbins_all);     //distn of size at age in cL block 3
  matrix lenprob_cP1_all(1,nages,1,nlenbins_all);     //distn of size at age in cP block 1  
  matrix lenprob_cP2_all(1,nages,1,nlenbins_all);     //distn of size at age in cP block 2
  matrix lenprob_cP3_all(1,nages,1,nlenbins_all);     //distn of size at age in cP block 3
  matrix lenprob_cT1_all(1,nages,1,nlenbins_all);     //distn of size at age in cT block 1  
  matrix lenprob_cT2_all(1,nages,1,nlenbins_all);     //distn of size at age in cT block 2
  matrix lenprob_HB1_all(1,nages,1,nlenbins_all);     //distn of size at age in HB block 1
  matrix lenprob_HB2_all(1,nages,1,nlenbins_all);     //distn of size at age in HB block 2
  matrix lenprob_HB3_all(1,nages,1,nlenbins_all);     //distn of size at age in HB block 3
  matrix lenprob_HB4_all(1,nages,1,nlenbins_all);     //distn of size at age in HB block 4
  matrix lenprob_HB5_all(1,nages,1,nlenbins_all);     //distn of size at age in HB block 5
  matrix lenprob_mrip1_all(1,nages,1,nlenbins_all);   //distn of size at age in mrip block 1
  matrix lenprob_mrip2_all(1,nages,1,nlenbins_all);   //distn of size at age in mrip block 2
  matrix lenprob_mrip3_all(1,nages,1,nlenbins_all);   //distn of size at age in mrip block 3
  matrix lenprob_mrip4_all(1,nages,1,nlenbins_all);   //distn of size at age in mrip block 4
  matrix lenprob_mrip5_all(1,nages,1,nlenbins_all);   //distn of size at age in mrip block 5

  matrix lenprob_comm_D2_all(1,nages,1,nlenbins_all);   //distn of size at age in cL discards comm block 2   
  matrix lenprob_comm_D3_all(1,nages,1,nlenbins_all);   //distn of size at age in cL discards comm block 3
  matrix lenprob_comm_D4_all(1,nages,1,nlenbins_all);   //distn of size at age in cL discards comm block 4
  matrix lenprob_HB_D2_all(1,nages,1,nlenbins_all);     //distn of size at age in HB discards rec block 2 
  matrix lenprob_HB_D3_all(1,nages,1,nlenbins_all);     //distn of size at age in HB discards rec block 3
  matrix lenprob_HB_D4_all(1,nages,1,nlenbins_all);     //distn of size at age in HB discards rec block 4
  matrix lenprob_HB_D5_all(1,nages,1,nlenbins_all);     //distn of size at age in HB discards rec block 5
  matrix lenprob_mrip_D1_all(1,nages,1,nlenbins_all);   //distn of size at age in mrip discards rec block 1
  matrix lenprob_mrip_D2_all(1,nages,1,nlenbins_all);   //distn of size at age in mrip discards rec block 2
  matrix lenprob_mrip_D3_all(1,nages,1,nlenbins_all);   //distn of size at age in mrip discards rec block 3
  matrix lenprob_mrip_D4_all(1,nages,1,nlenbins_all);   //distn of size at age in mrip discards rec block 4
  matrix lenprob_mrip_D5_all(1,nages,1,nlenbins_all);   //distn of size at age in mrip discards rec block 5

  vector len_sd(1,nages);
  vector len_cv(1,nages); //for fishgraph 
  
//----Predicted length and age compositions
  matrix pred_Mbft_lenc(1,nyr_Mbft_lenc,1,nlenbins);
  matrix pred_cL_lenc(1,nyr_cL_lenc,1,nlenbins);
  matrix pred_cP_lenc(1,nyr_cP_lenc,1,nlenbins);
  matrix pred_HB_lenc(1,nyr_HB_lenc,1,nlenbins);
  matrix pred_HB_D_lenc(1,nyr_HB_D_lenc,1,nlenbins);
  matrix pred_mrip_lenc(1,nyr_mrip_lenc,1,nlenbins);  
  
  matrix L_cP_num_pool(1,nyr_cP_lenc,1,nages);          //landings (numbers) at age pooled for length comps
  matrix L_cP_num_pool_yr(1,nyr_cP_lenc_pool,1,nages);  //scaled and weighted landings (numbers) for pooling length comps
  //matrix L_mrip_num_pool(1,nyr_mrip_lenc,1,nages);          //landings (numbers) at age pooled for length comps
  //matrix L_mrip_num_pool_yr(1,nyr_mrip_lenc_pool,1,nages);  //scaled and weighted landings (numbers) for pooling length comps
  
//  //##p_lenc_fishery pars require age comp and length comp data for estimation
  number p_lenc_cL2;
  number p_lenc_cL3;  
  number p_lenc_cP2;
  number p_lenc_cP3;  
  number p_lenc_cT2;
  number p_lenc_cT3;  
  number p_lenc_HB2;
  number p_lenc_HB3;
  number p_lenc_HB4;
  number p_lenc_HB5;
  number p_lenc_mrip2;
  number p_lenc_mrip3;
  number p_lenc_mrip4;
  number p_lenc_mrip5;
    
  number p_lenc_comm_D2;
  number p_lenc_comm_D3;
  number p_lenc_comm_D4;
  number p_lenc_HB_D2; 
  number p_lenc_HB_D3;
  number p_lenc_HB_D4;
  number p_lenc_HB_D5;  
  number p_lenc_mrip_D1; 
  number p_lenc_mrip_D2; 
  number p_lenc_mrip_D3;
  number p_lenc_mrip_D4;
  number p_lenc_mrip_D5;
  
  matrix pred_Mbft_agec(1,nyr_Mbft_agec,1,nages);
  matrix ErrorFree_Mbft_agec(1,nyr_Mbft_agec,1,nages);
  matrix pred_Mcvt_agec(1,nyr_Mcvt_agec,1,nages);
  matrix ErrorFree_Mcvt_agec(1,nyr_Mcvt_agec,1,nages);  
  matrix pred_cL_agec(1,nyr_cL_agec,1,nages);
  matrix ErrorFree_cL_agec(1,nyr_cL_agec,1,nages);
  matrix pred_cP_agec(1,nyr_cP_agec,1,nages);
  matrix ErrorFree_cP_agec(1,nyr_cP_agec,1,nages);  
  matrix pred_HB_agec(1,nyr_HB_agec,1,nages);
  matrix ErrorFree_HB_agec(1,nyr_HB_agec,1,nages);
  //matrix pred_mrip_agec(1,nyr_mrip_agec,1,nages);
  //matrix ErrorFree_mrip_agec(1,nyr_mrip_agec,1,nages);
 
//effective sample size applied in multinomial distributions
  vector nsamp_Mbft_lenc_allyr(styr,endyr);
  vector nsamp_cL_lenc_allyr(styr,endyr);
  vector nsamp_cP_lenc_allyr(styr,endyr);
  vector nsamp_HB_lenc_allyr(styr,endyr);
  vector nsamp_HB_D_lenc_allyr(styr,endyr);
  vector nsamp_mrip_lenc_allyr(styr,endyr);
  vector nsamp_Mbft_agec_allyr(styr,endyr);
  vector nsamp_Mcvt_agec_allyr(styr,endyr);
  vector nsamp_cL_agec_allyr(styr,endyr);
  vector nsamp_cP_agec_allyr(styr,endyr);  
  vector nsamp_HB_agec_allyr(styr,endyr);
  //vector nsamp_mrip_agec_allyr(styr,endyr);

//Nfish used in MCB analysis (not used in fitting)
  vector nfish_Mbft_lenc_allyr(styr,endyr);
  vector nfish_cL_lenc_allyr(styr,endyr);
  vector nfish_cP_lenc_allyr(styr,endyr);
  vector nfish_HB_lenc_allyr(styr,endyr);
  vector nfish_HB_D_lenc_allyr(styr,endyr);
  vector nfish_mrip_lenc_allyr(styr,endyr);
  vector nfish_Mbft_agec_allyr(styr,endyr);
  vector nfish_Mcvt_agec_allyr(styr,endyr);
  vector nfish_cL_agec_allyr(styr,endyr);
  vector nfish_cP_agec_allyr(styr,endyr);  
  vector nfish_HB_agec_allyr(styr,endyr);
  //vector nfish_mrip_agec_allyr(styr,endyr);
  
//Computed effective sample size for output (not used in fitting)
  vector neff_Mbft_lenc_allyr_out(styr,endyr);
  vector neff_cL_lenc_allyr_out(styr,endyr);
  vector neff_cP_lenc_allyr_out(styr,endyr);
  vector neff_HB_lenc_allyr_out(styr,endyr);
  vector neff_HB_D_lenc_allyr_out(styr,endyr);
  vector neff_mrip_lenc_allyr_out(styr,endyr);
  vector neff_Mbft_agec_allyr_out(styr,endyr);
  vector neff_Mcvt_agec_allyr_out(styr,endyr);
  vector neff_cL_agec_allyr_out(styr,endyr);
  vector neff_cP_agec_allyr_out(styr,endyr);  
  vector neff_HB_agec_allyr_out(styr,endyr);
  //vector neff_mrip_agec_allyr_out(styr,endyr);


//-----Population-----------------------------------------------------------------------------------
  matrix N(styr,endyr,1,nages);             //Population numbers by year and age at start of yr
  matrix N_mdyr(styr,endyr,1,nages);        //Population numbers by year and age at mdpt of yr: used for comps and cpue
  matrix N_spawn(styr,endyr,1,nages);       //Population numbers by year and age at peaking spawning: used for SSB  
  //init_bounded_vector log_Nage_dev(2,nages,-5,3,1); //log deviations on initial abundance at age
  init_bounded_vector log_Nage_dev(2,nages,log_Nage_dev_LO,log_Nage_dev_HI,log_Nage_dev_PH);
  //vector log_Nage_dev(2,nages);
  vector log_Nage_dev_output(1,nages);      //used in output. equals zero for first age
  matrix B(styr,endyr,1,nages);             //Population biomass by year and age at start of yr
  vector totB(styr,endyr);                  //Total biomass by year
  vector totN(styr,endyr);                  //Total abundance by year
  vector SSB(styr,endyr);                   //Total spawning biomass by year (scaled popn fecundity)
  vector MatFemB(styr,endyr);               //Total spawning biomass by year (total mature female biomass)  
  vector rec(styr,endyr);                   //Recruits by year
  vector prop_f(1,nages);                   //Proportion female by age
  vector maturity_f(1,nages);               //Proportion of female mature at age
  vector reprod(1,nages);                   //vector used to compute spawning biomass (scaled popn fecundity)
  vector reprod2(1,nages);                  //vector used to compute mature female biomass 

////---Stock-Recruit Function (Beverton-Holt, steepness parameterization)----------
  //init_bounded_number log_R0(13,20,1);        //log(virgin Recruitment)
  //number log_R0;
  init_bounded_number log_R0(log_R0_LO,log_R0_HI,log_R0_PH);        //log(virgin Recruitment)
  vector log_R0_out(1,8);
  number R0;                                  //virgin recruitment
  //init_bounded_number steep(0.21,0.991,3);    //steepness
  init_bounded_number steep(steep_LO,steep_HI,steep_PH); //steepness
  vector steep_out(1,8);
  //number steep;  //uncomment to fix steepness, comment line directly above
  //init_bounded_number rec_sigma(0.1,1.5,4);  //sd recruitment residuals
  init_bounded_number rec_sigma(rec_sigma_LO,rec_sigma_HI,rec_sigma_PH);  //sd recruitment residuals  
  vector rec_sigma_out(1,8);
  init_bounded_number R_autocorr(R_autocorr_LO,R_autocorr_HI,R_autocorr_PH);  //autocorrelation in SR  
  vector R_autocorr_out(1,8);
  
  number rec_sigma_sq;                        //square of rec_sigma      
  number rec_sigma_sqd2;                      //square of rec_sigma divided by two
  number rec_logL_add;                        //additive term in -logL term   
    
  //init_bounded_dev_vector log_rec_dev(styr_rec_dev,endyr,-3,3,2);  //log recruitment deviations
  //vector log_rec_dev(styr_rec_dev,endyr);
  init_bounded_dev_vector log_rec_dev(styr_rec_dev,endyr_rec_phase2,log_rec_dev_LO,log_rec_dev_HI,log_rec_dev_PH);
  vector log_rec_dev_output(styr_rec_dev,endyr);             //used in output. equals zero except for yrs in log_rec_dev
  vector log_rec_dev_out(styr_rec_dev,endyr_rec_phase2);  //used in output for bound checking

  number var_rec_dev;                                //variance of log recruitment deviations, from yrs with unconstrainted S-R(XXXX-XXXX)
  number sigma_rec_dev;                              //sample SD of log residuals (may not equal rec_sigma 
  number BiasCor;                               //Bias correction in equilibrium recruits
  //number R_autocorr;
  number S0;                                    //equal to spr_F0*R0 = virgin SSB
  number B0;                                    //equal to bpr_F0*R0 = virgin B  
  number R1;                                    //Recruits in styr
  number R_virgin;                              //unfished recruitment with bias correction
  vector SdS0(styr,endyr);                      //SSB / virgin SSB

  init_bounded_number log_dm_Mbft_lc(log_dm_Mbft_lc_LO,log_dm_Mbft_lc_HI,log_dm_Mbft_lc_PH);
  init_bounded_number log_dm_Mcvt_lc(log_dm_Mcvt_lc_LO,log_dm_Mcvt_lc_HI,log_dm_Mcvt_lc_PH);
  init_bounded_number log_dm_cL_lc(log_dm_cL_lc_LO,log_dm_cL_lc_HI,log_dm_cL_lc_PH);
  init_bounded_number log_dm_cP_lc(log_dm_cP_lc_LO,log_dm_cP_lc_HI,log_dm_cP_lc_PH);
  init_bounded_number log_dm_HB_lc(log_dm_HB_lc_LO,log_dm_HB_lc_HI,log_dm_HB_lc_PH);
  init_bounded_number log_dm_HB_D_lc(log_dm_HB_D_lc_LO,log_dm_HB_D_lc_HI,log_dm_HB_D_lc_PH);
  init_bounded_number log_dm_mrip_lc(log_dm_mrip_lc_LO,log_dm_mrip_lc_HI,log_dm_mrip_lc_PH);
  init_bounded_number log_dm_Mbft_ac(log_dm_Mbft_ac_LO,log_dm_Mbft_ac_HI,log_dm_Mbft_ac_PH);
  init_bounded_number log_dm_Mcvt_ac(log_dm_Mcvt_ac_LO,log_dm_Mcvt_ac_HI,log_dm_Mcvt_ac_PH);
  init_bounded_number log_dm_cL_ac(log_dm_cL_ac_LO,log_dm_cL_ac_HI,log_dm_cL_ac_PH);
  init_bounded_number log_dm_cP_ac(log_dm_cP_ac_LO,log_dm_cP_ac_HI,log_dm_cP_ac_PH);
  init_bounded_number log_dm_HB_ac(log_dm_HB_ac_LO,log_dm_HB_ac_HI,log_dm_HB_ac_PH);
  //init_bounded_number log_dm_mrip_ac(log_dm_mrip_ac_LO,log_dm_mrip_ac_HI,log_dm_mrip_ac_PH);
  vector log_dm_Mbft_lc_out(1,8);
  vector log_dm_Mcvt_lc_out(1,8);
  vector log_dm_cL_lc_out(1,8);
  vector log_dm_cP_lc_out(1,8);
  vector log_dm_HB_lc_out(1,8);
  vector log_dm_mrip_lc_out(1,8);
  vector log_dm_HB_D_lc_out(1,8);
  vector log_dm_Mbft_ac_out(1,8);
  vector log_dm_Mcvt_ac_out(1,8);
  vector log_dm_cL_ac_out(1,8);
  vector log_dm_cP_ac_out(1,8);
  vector log_dm_HB_ac_out(1,8);
  //vector log_dm_mrip_ac_out(1,8);
  
//-----------------------------------------------------------------------------------------------------------------------------------------------
//---Selectivity-------------------------------------------------------------------------

//MARMAP Mbft selectivity -------------------------------------------------------------------------
  matrix sel_Mbft(styr,endyr,1,nages);
  vector sel_Mbft_vec(1,nages);
  
  init_bounded_number selpar_A50_Mbft(selpar_A50_Mbft_LO,selpar_A50_Mbft_HI,selpar_A50_Mbft_PH);
  init_bounded_number selpar_slope_Mbft(selpar_slope_Mbft_LO,selpar_slope_Mbft_HI,selpar_slope_Mbft_PH);
   
  vector selpar_A50_Mbft_out(1,8);
  vector selpar_slope_Mbft_out(1,8);
  
//MARMAP Mcvt selectivity -------------------------------------------------------------------------
  matrix sel_Mcvt(styr,endyr,1,nages);
  vector sel_Mcvt_vec(1,nages);
  
  init_bounded_number selpar_A50_Mcvt(selpar_A50_Mcvt_LO,selpar_A50_Mcvt_HI,selpar_A50_Mcvt_PH);
  init_bounded_number selpar_slope_Mcvt(selpar_slope_Mcvt_LO,selpar_slope_Mcvt_HI,selpar_slope_Mcvt_PH);
   
  vector selpar_A50_Mcvt_out(1,8);
  vector selpar_slope_Mcvt_out(1,8);
  
//Commercial handline selectivity-------------------------------------------------
  matrix sel_cL(styr,endyr,1,nages);  
  //vector sel_cL_1(1,nages); //sel in period 1 assumed equal to period 2
  vector sel_cL_2(1,nages); //sel in period 2
  vector sel_cL_3(1,nages); //sel in period 3 
  vector sel_cL_4(1,nages); //sel in period 4 
  
  init_bounded_number selpar_A50_cL2(selpar_A50_cL2_LO,selpar_A50_cL2_HI,selpar_A50_cL2_PH);
  init_bounded_number selpar_slope_cL2(selpar_slope_cL2_LO,selpar_slope_cL2_HI,selpar_slope_cL2_PH);
  init_bounded_number selpar_A50_cL3(selpar_A50_cL3_LO,selpar_A50_cL3_HI,selpar_A50_cL3_PH);
  init_bounded_number selpar_slope_cL3(selpar_slope_cL3_LO,selpar_slope_cL3_HI,selpar_slope_cL3_PH);
  init_bounded_number selpar_A50_cL4(selpar_A50_cL4_LO,selpar_A50_cL4_HI,selpar_A50_cL4_PH);
  init_bounded_number selpar_slope_cL4(selpar_slope_cL4_LO,selpar_slope_cL4_HI,selpar_slope_cL4_PH);
 
  vector selpar_A50_cL2_out(1,8);
  vector selpar_slope_cL2_out(1,8);
  vector selpar_A50_cL3_out(1,8);
  vector selpar_slope_cL3_out(1,8);
  vector selpar_A50_cL4_out(1,8);
  vector selpar_slope_cL4_out(1,8);
     
//commercial discards (handline + pots)
  matrix sel_comm_D(styr,endyr,1,nages); 
  vector sel_comm_D_2(1,nages);         //sel in period 2
  vector sel_comm_D_3(1,nages);         //sel in period 3 
  vector sel_comm_D_4(1,nages);         //sel in period 4  
  vector sel_comm_D_quota3(1,nages);    //sel in period 3 when quotas were in place (2009,2010) Also in 2011-2012   
 
//values used for weighting selex and avg weights of discards during yrs with quotas
  number Dopen_cL; number Dclosed_cL; number Lopen_cL;  
  number Dopen_cP; number Dclosed_cP; number Lopen_cP;
  number D_sum_cLcP; 
  number Dprop_comm_sel_D; number Dprop_comm_sel_cL; number Dprop_comm_sel_cP; 
         
//Commercial pots selectivity -------------------------------------------------            
  matrix sel_cP(styr,endyr,1,nages); 
  //  vector sel_cP_1(1,nages); //sel vector in period 1 assumed equal to period 2 
  vector sel_cP_2(1,nages);   //sel vector in period 2  
  vector sel_cP_3(1,nages);   //sel vector in period 3  
  vector sel_cP_4(1,nages);   //sel vector in period 3    
  
  init_bounded_number selpar_A50_cP2(selpar_A50_cP2_LO,selpar_A50_cP2_HI,selpar_A50_cP2_PH);
  init_bounded_number selpar_slope_cP2(selpar_slope_cP2_LO,selpar_slope_cP2_HI,selpar_slope_cP2_PH);
  init_bounded_number selpar_A50_cP3(selpar_A50_cP3_LO,selpar_A50_cP3_HI,selpar_A50_cP3_PH);
  init_bounded_number selpar_slope_cP3(selpar_slope_cP3_LO,selpar_slope_cP3_HI,selpar_slope_cP3_PH);
  init_bounded_number selpar_A50_cP4(selpar_A50_cP4_LO,selpar_A50_cP4_HI,selpar_A50_cP4_PH);
  init_bounded_number selpar_slope_cP4(selpar_slope_cP4_LO,selpar_slope_cP4_HI,selpar_slope_cP4_PH);
 
  vector selpar_A50_cP2_out(1,8);
  vector selpar_slope_cP2_out(1,8);
  vector selpar_A50_cP3_out(1,8);
  vector selpar_slope_cP3_out(1,8);
  vector selpar_A50_cP4_out(1,8);
  vector selpar_slope_cP4_out(1,8);

//Commercial trawl selectivity -------------------------------------------------            
  matrix sel_cT(styr,endyr,1,nages);  //mirrors comm pot sel

//Headboat selectivity -------------------------------------------------
  matrix sel_HB(styr,endyr,1,nages);  
  vector sel_HB_1(1,nages); //sel in period 1 
  vector sel_HB_2(1,nages); //sel in period 2
  vector sel_HB_3(1,nages); //sel in period 3 
  vector sel_HB_4(1,nages); //sel in period 4  
  vector sel_HB_5(1,nages); //sel in period 5   
  
  init_bounded_number selpar_A50_HB1(selpar_A50_HB1_LO,selpar_A50_HB1_HI,selpar_A50_HB1_PH);
  init_bounded_number selpar_slope_HB1(selpar_slope_HB1_LO,selpar_slope_HB1_HI,selpar_slope_HB1_PH);
  init_bounded_number selpar_A50_HB2(selpar_A50_HB2_LO,selpar_A50_HB2_HI,selpar_A50_HB2_PH);
  init_bounded_number selpar_slope_HB2(selpar_slope_HB2_LO,selpar_slope_HB2_HI,selpar_slope_HB2_PH);    
  init_bounded_number selpar_A50_HB3(selpar_A50_HB3_LO,selpar_A50_HB3_HI,selpar_A50_HB3_PH);
  init_bounded_number selpar_slope_HB3(selpar_slope_HB3_LO,selpar_slope_HB3_HI,selpar_slope_HB3_PH);
  init_bounded_number selpar_A50_HB4(selpar_A50_HB4_LO,selpar_A50_HB4_HI,selpar_A50_HB4_PH);
  init_bounded_number selpar_slope_HB4(selpar_slope_HB4_LO,selpar_slope_HB4_HI,selpar_slope_HB4_PH);
  init_bounded_number selpar_A50_HB5(selpar_A50_HB5_LO,selpar_A50_HB5_HI,selpar_A50_HB5_PH);
  init_bounded_number selpar_slope_HB5(selpar_slope_HB5_LO,selpar_slope_HB5_HI,selpar_slope_HB5_PH);
  
  vector selpar_A50_HB1_out(1,8);
  vector selpar_slope_HB1_out(1,8);
  vector selpar_A50_HB2_out(1,8);
  vector selpar_slope_HB2_out(1,8);
  vector selpar_A50_HB3_out(1,8);
  vector selpar_slope_HB3_out(1,8);
  vector selpar_A50_HB4_out(1,8);
  vector selpar_slope_HB4_out(1,8);
  vector selpar_A50_HB5_out(1,8);
  vector selpar_slope_HB5_out(1,8);
  
  number selpar_Age0_HB_D;
  number selpar_Age1_HB_D;
  number selpar_Age2_HB_D;
  
  //---headboat discards--------------------------------- 
  matrix sel_HB_D(styr,endyr,1,nages); 
  vector sel_HB_D_1(1,nages); //sel in period 1, assumed equal to period 2
  vector sel_HB_D_2(1,nages); //sel in period 2
  vector sel_HB_D_3(1,nages); //sel in period 3 
  vector sel_HB_D_4(1,nages); //sel in period 4  
  vector sel_HB_D_5(1,nages); //sel in period 5 
  
  init_bounded_number selpar_A50_HBD4(selpar_A50_HBD4_LO,selpar_A50_HBD4_HI,selpar_A50_HBD4_PH);
  init_bounded_number selpar_slope_HBD4(selpar_slope_HBD4_LO,selpar_slope_HBD4_HI,selpar_slope_HBD4_PH);
  init_bounded_number selpar_A502_HBD4(selpar_A502_HBD4_LO,selpar_A502_HBD4_HI,selpar_A502_HBD4_PH);
  init_bounded_number selpar_slope2_HBD4(selpar_slope2_HBD4_LO,selpar_slope2_HBD4_HI,selpar_slope2_HBD4_PH);
  
  init_bounded_number selpar_A50_HBD5(selpar_A50_HBD5_LO,selpar_A50_HBD5_HI,selpar_A50_HBD5_PH);
  init_bounded_number selpar_slope_HBD5(selpar_slope_HBD5_LO,selpar_slope_HBD5_HI,selpar_slope_HBD5_PH);
  init_bounded_number selpar_A502_HBD5(selpar_A502_HBD5_LO,selpar_A502_HBD5_HI,selpar_A502_HBD5_PH);
  init_bounded_number selpar_slope2_HBD5(selpar_slope2_HBD5_LO,selpar_slope2_HBD5_HI,selpar_slope2_HBD5_PH);
  
  vector selpar_A50_HBD4_out(1,8);
  vector selpar_slope_HBD4_out(1,8);
  vector selpar_A502_HBD4_out(1,8);
  vector selpar_slope2_HBD4_out(1,8);
  vector selpar_A50_HBD5_out(1,8);
  vector selpar_slope_HBD5_out(1,8);
  vector selpar_A502_HBD5_out(1,8);
  vector selpar_slope2_HBD5_out(1,8);
  
  vector vecprob_HB_D2(4,nages);     //prob of less than size limit
  vector vecprob_HB_D3(4,nages);     //prob of less than size limit
  vector vecprob_HB_D4(4,nages);     //prob of less than size limit
  vector vecprob_HB_D5(4,nages);     //prob of less than size limit

  init_bounded_number selpar_Age0_HB_D_logit(selpar_Age0_HB_D_logit_LO,selpar_Age0_HB_D_logit_HI,selpar_Age0_HB_D_logit_PH);
  init_bounded_number selpar_Age1_HB_D_logit(selpar_Age1_HB_D_logit_LO,selpar_Age1_HB_D_logit_HI,selpar_Age1_HB_D_logit_PH);
  init_bounded_number selpar_Age2_HB_D_logit(selpar_Age2_HB_D_logit_LO,selpar_Age2_HB_D_logit_HI,selpar_Age2_HB_D_logit_PH);
  
  vector selpar_Age0_HB_D_logit_out(1,8);                                  
  vector selpar_Age1_HB_D_logit_out(1,8);                                  
  vector selpar_Age2_HB_D_logit_out(1,8); 
  
  vector prob_belowsizelim_block1(1,nages);
  vector prob_belowsizelim_block2(1,nages);
  vector prob_belowsizelim_block3(1,nages);
  vector prob_belowsizelim_block4(1,nages);
  vector prob_belowsizelim_block5(1,nages);
  
  number zscore_lsizelim1;
  number zscore_lsizelim2;
  number zscore_lsizelim3;
  number zscore_lsizelim4;
  number zscore_lsizelim5;
  
  number cprob_lsizelim1;
  number cprob_lsizelim2;
  number cprob_lsizelim3;
  number cprob_lsizelim4;
  number cprob_lsizelim5;
  
    //mrip selectivity  -------------------------------------------------
  matrix sel_mrip(styr,endyr,1,nages);  
  matrix sel_mrip_D(styr,endyr,1,nages); 
  vector sel_mrip1(1,nages); //sel in period 1 
  vector sel_mrip2(1,nages); //sel in period 2
  vector sel_mrip3(1,nages); //sel in period 3 
  vector sel_mrip4(1,nages); //sel in period 4  
  vector sel_mrip5(1,nages); //sel in period 5   
     
  init_bounded_number selpar_A50_mrip1(selpar_A50_mrip1_LO,selpar_A50_mrip1_HI,selpar_A50_mrip1_PH);
  init_bounded_number selpar_slope_mrip1(selpar_slope_mrip1_LO,selpar_slope_mrip1_HI,selpar_slope_mrip1_PH);
  init_bounded_number selpar_A50_mrip2(selpar_A50_mrip2_LO,selpar_A50_mrip2_HI,selpar_A50_mrip2_PH);
  init_bounded_number selpar_slope_mrip2(selpar_slope_mrip2_LO,selpar_slope_mrip2_HI,selpar_slope_mrip2_PH);    
  init_bounded_number selpar_A50_mrip3(selpar_A50_mrip3_LO,selpar_A50_mrip3_HI,selpar_A50_mrip3_PH);
  init_bounded_number selpar_slope_mrip3(selpar_slope_mrip3_LO,selpar_slope_mrip3_HI,selpar_slope_mrip3_PH);
  init_bounded_number selpar_A50_mrip4(selpar_A50_mrip4_LO,selpar_A50_mrip4_HI,selpar_A50_mrip4_PH);
  init_bounded_number selpar_slope_mrip4(selpar_slope_mrip4_LO,selpar_slope_mrip4_HI,selpar_slope_mrip4_PH);
  init_bounded_number selpar_A50_mrip5(selpar_A50_mrip5_LO,selpar_A50_mrip5_HI,selpar_A50_mrip5_PH);
  init_bounded_number selpar_slope_mrip5(selpar_slope_mrip5_LO,selpar_slope_mrip5_HI,selpar_slope_mrip5_PH);
  
  vector selpar_A50_mrip1_out(1,8);
  vector selpar_slope_mrip1_out(1,8);
  vector selpar_A50_mrip2_out(1,8);
  vector selpar_slope_mrip2_out(1,8);
  vector selpar_A50_mrip3_out(1,8);
  vector selpar_slope_mrip3_out(1,8);
  vector selpar_A50_mrip4_out(1,8);
  vector selpar_slope_mrip4_out(1,8);
  vector selpar_A50_mrip5_out(1,8);
  vector selpar_slope_mrip5_out(1,8);
  
//effort-weighted, recent selectivities
  vector sel_wgted_L(1,nages);  //toward landings
  vector sel_wgted_D(1,nages);  //toward discards  
  vector sel_wgted_tot(1,nages);//toward Z, landings plus deads discards

//-----------------------------------------------------------------------------------------------------------------------------------------------
//-------CPUE Predictions--------------------------------
  vector pred_Mbft_cpue(styr_Mbft_cpue,endyr_Mbft_cpue);       //predicted Mbft U (fish/trap-hour)
  matrix N_Mbft(styr_Mbft_cpue,endyr_Mbft_cpue,1,nages);       //used to compute Mbft index
  vector pred_Mcvt_cpue(styr_Mcvt_cpue,endyr_Mcvt_cpue);       //predicted Mcvt U (fish/trap-hour)
  matrix N_Mcvt(styr_Mcvt_cpue,endyr_Mcvt_cpue,1,nages);       //used to compute Mcvt index
  //vector pred_Vid_cpue(styr_Vid_cpue,endyr_Vid_cpue);       //predicted Mcvt U (fish/trap-hour)
  //matrix N_Vid(styr_Vid_cpue,endyr_Vid_cpue,1,nages);       //used to compute Mcvt index
  vector pred_cL_cpue(styr_cL_cpue,endyr_cL_cpue);             //predicted cL U (pounds/hook-hour)
  matrix N_cL(styr_cL_cpue,endyr_cL_cpue,1,nages);             //used to compute cL index
  vector pred_HB_cpue(styr_HB_cpue,endyr_HB_cpue);             //predicted HB U (pounds/hour)
  matrix N_HB(styr_HB_cpue,endyr_HB_cpue,1,nages);             //used to compute HB index
  //vector pred_HBD_cpue(styr_HBD_cpue,endyr_HBD_cpue);          //predicted HBD U (fish/angler-hour)
  //matrix N_HBD(styr_HBD_cpue,endyr_HBD_cpue,1,nages);          //used to compute HBD index

//---Catchability (CPUE q's)----------------------------------------------------------
  //init_bounded_number log_q_Mbft(-20,-10,1);
  //init_bounded_number log_q_Mcvt(-20,-10,1);
  //init_bounded_number log_q_Vid(-20,-10,1);  
  //init_bounded_number log_q_cL(-20,-5,1);
  //init_bounded_number log_q_HB(-20,-5,1);
  //init_bounded_number log_q_HBD(-20,-10,1);
  init_bounded_number q_rate(0.001,0.1,set_q_rate_phase);
  
  init_bounded_number log_q_Mbft(log_q_Mbft_LO,log_q_Mbft_HI,log_q_Mbft_PH);
  init_bounded_number log_q_Mcvt(log_q_Mcvt_LO,log_q_Mcvt_HI,log_q_Mcvt_PH);
  init_bounded_number log_q_cL(log_q_cL_LO,log_q_cL_HI,log_q_cL_PH);
  init_bounded_number log_q_HB(log_q_HB_LO,log_q_HB_HI,log_q_HB_PH);
  //init_bounded_number log_q_HBD(log_q_HBD_LO,log_q_HBD_HI,log_q_HBD_PH);
  //init_bounded_number log_q_Vid(log_q_Vid_LO,log_q_Vid_HI,log_q_Vid_PH);
  
  vector log_q_Mbft_out(1,8);
  vector log_q_Mcvt_out(1,8);
  vector log_q_cL_out(1,8);
  vector log_q_HB_out(1,8);
  //vector log_q_HBD_out(1,8);
  //vector log_q_Vid_out(1,8);
  
  //number q_rate;
  vector q_rate_fcn_cL(styr_cL_cpue,endyr_cL_cpue);       //increase due to technology creep (saturates in 2003)
  vector q_rate_fcn_HB(styr_HB_cpue,endyr_HB_cpue);         //increase due to technology creep (saturates in 2003)
  //vector q_rate_fcn_HBD(styr_HBD_cpue,endyr_HBD_cpue);      //increase due to technology creep (saturates in 2003)
  
  //init_bounded_number q_DD_beta(0.1,0.9,set_q_DD_phase);  
  number q_DD_beta;
  vector q_DD_fcn(styr,endyr);    //density dependent function as a multiple of q (scaled a la Katsukawa and Matsuda. 2003)
  number B0_q_DD;                 //B0 of ages q_DD_age plus
  vector B_q_DD(styr,endyr);      //annual biomass of ages q_DD_age plus

  //init_bounded_vector q_RW_log_dev_cL(styr_cL_cpue,endyr_cL_cpue-1,-3.0,3.0,set_q_RW_phase);
  //init_bounded_vector q_RW_log_dev_HB(styr_HB_cpue,endyr_HB_cpue-1,-3.0,3.0,set_q_RW_phase);
  //init_bounded_vector q_RW_log_dev_HBD(styr_HBD_cpue,endyr_HBD_cpue-1,-3.0,3.0,set_q_RW_phase);
 
  //init_bounded_vector q_RW_log_dev_cL(styr_cL_cpue,endyr_cL_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  //init_bounded_vector q_RW_log_dev_HB(styr_HB_cpue,endyr_HB_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  //init_bounded_vector q_RW_log_dev_HBD(styr_HBD_cpue,endyr_HBD_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  init_bounded_vector q_RW_log_dev_cL(styr_cL_cpue,endyr_cL_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  init_bounded_vector q_RW_log_dev_HB(styr_HB_cpue,endyr_HB_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  //init_bounded_vector q_RW_log_dev_HBD(styr_HBD_cpue,endyr_HBD_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
 
  vector q_cL(styr_cL_cpue,endyr_cL_cpue);
  vector q_HB(styr_HB_cpue,endyr_HB_cpue);
  //vector q_HBD(styr_HBD_cpue,endyr_HBD_cpue);

//----------------------------------------------------------------------------------------------------------------------------------------------- 
//---Landings Bias for recreational landings------------------------------------------------------------------
  //init_bounded_number L_mrip_bias(0.1,10.0,3);  
  number L_hb_bias;
  number L_mrip_bias;
  number L_comm_bias;

//---Landings in numbers (total or 1000 fish) and in wgt (klb)--------------------------------------------------
  matrix L_cL_num(styr,endyr,1,nages);  //landings (numbers) at age
  matrix L_cL_klb(styr,endyr,1,nages);  //landings (1000 lb whole weight) at age
  vector pred_cL_L_knum(styr,endyr);    //yearly landings in 1000 fish summed over ages
  vector pred_cL_L_klb(styr,endyr);     //yearly landings in 1000 lb summed over ages

  matrix L_cP_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cP_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  vector pred_cP_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages   
  vector pred_cP_L_klb(styr,endyr);      //yearly landings in 1000 lb summed over ages 

  matrix L_cT_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cT_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  vector pred_cT_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages   
  vector pred_cT_L_klb(styr,endyr);      //yearly landings in 1000 lb summed over ages 

  matrix L_HB_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_HB_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  vector pred_HB_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages  
  vector pred_HB_L_klb(styr,endyr);      //yearly landings in 1000 lb summed over ages
  vector obs_HB_L_wbias(styr,endyr);     //yearly landings observed, perhaps adjusted for multiplicative bias

  matrix L_mrip_num(styr,endyr,1,nages);  //landings (numbers) at age
  matrix L_mrip_klb(styr,endyr,1,nages);  //landings (1000 lb whole weight) at age    
  vector pred_mrip_L_knum(styr,endyr);    //yearly landings in 1000 fish summed over ages  
  vector pred_mrip_L_klb(styr,endyr);     //yearly landings in 1000 lb summed over ages
  vector obs_mrip_L_wbias(styr,endyr);    //yearly landings observed, perhaps adjusted for multiplicative bias

  matrix L_total_num(styr,endyr,1,nages);//total landings in number at age
  matrix L_total_klb(styr,endyr,1,nages);//landings in klb at age 
  vector L_total_knum_yr(styr,endyr);    //total landings in 1000 fish by yr summed over ages  
  vector L_total_klb_yr(styr,endyr);     //total landings (klb) by yr summed over ages

//---Dead discards in numbers (total or 1000 fish) and in wgt (klb) --------------------------------------------------
  matrix D_comm_num(styr,endyr,1,nages);   //discards (numbers) at age
  matrix D_comm_klb(styr,endyr,1,nages);   //discards (1000 lb) at age
  vector pred_comm_D_knum(styr,endyr);     //yearly discards summed over ages
  vector pred_comm_D_klb(styr,endyr);      //yearly discards in klb summed over ages
  vector obs_comm_D(styr_cL_D,endyr_cL_D); //observed releases multiplied by discard mortality
  vector obs_cL_D(styr_cL_D,endyr_cL_D);   //observed releases multiplied by discard mortality
  vector obs_cP_D(styr_cL_D,endyr_cL_D);   //observed releases multiplied by discard mortality
  vector comm_D_cv(styr_cL_D,endyr_cL_D);  //CVs for fitting combined discards
  
  matrix D_HB_num(styr,endyr,1,nages);     //discards (numbers) at age
  matrix D_HB_klb(styr,endyr,1,nages);     //discards (1000 lb) at age
  vector pred_HB_D_knum(styr,endyr);       //yearly discards summed over ages
  vector pred_HB_D_klb(styr,endyr);        //yearly discards in klb summed over ages
  vector obs_HB_D(styr_HB_D,endyr_HB_D);   //observed releases multiplied by discard mortality

  matrix D_mrip_num(styr,endyr,1,nages);   //discards (numbers) at age
  matrix D_mrip_klb(styr,endyr,1,nages);   //discards (1000 lb) at age
  vector pred_mrip_D_knum(styr,endyr);     //yearly discards summed over ages
  vector pred_mrip_D_klb(styr,endyr);      //yearly discards in klb summed over ages
  vector obs_mrip_D(styr_mrip_D,endyr_mrip_D); //observed releases multiplied by discard mortality
  
  matrix D_total_num(styr,endyr,1,nages);  //total discards in number at age
  matrix D_total_klb(styr,endyr,1,nages);  //discards in klb at age 
  vector D_total_knum_yr(styr,endyr);      //total discards in 1000 fish by yr summed over ages  
  vector D_total_klb_yr(styr,endyr);       //total discards (klb) by yr summed over ages

////---MSY calcs----------------------------------------------------------------------------
  number F_cL_prop;       //proportion of F_sum attributable to hal, last X=selpar_n_yrs_wgted yrs, used for avg body weights
  number F_cP_prop;       //proportion of F_sum attributable to pots, last X yrs
  number F_HB_prop;       //proportion of F_sum attributable to headboat, last X yrs
  number F_mrip_prop;     //proportion of F_sum attributable to mrip, last X yrs
  number F_comm_D_prop;   //proportion of F_sum attributable to comm discards, last X yrs
  number F_HB_D_prop;     //proportion of F_sum attributable to headboat discards, last X yrs
  number F_mrip_D_prop;   //proportion of F_sum attributable to mrip discards, last X yrs
  number F_temp_sum;      //sum of geom mean Fsum's in last X yrs, used to compute F_fishery_prop

  vector F_end(1,nages);
  vector F_end_L(1,nages);
  vector F_end_D(1,nages);    
  number F_end_apex;
  
  number SSB_msy_out;           //SSB (popn fecudity) at msy
  number F_msy_out;             //F at msy
  number msy_klb_out;           //max sustainable yield (1000 lb)
  number msy_knum_out;          //max sustainable yield (1000 fish)  
  number B_msy_out;             //total biomass at MSY 
  number R_msy_out;             //equilibrium recruitment at F=Fmsy
  number D_msy_knum_out;        //equilibrium dead discards (1000 fish) at F=Fmsy
  number D_msy_klb_out;         //equilibrium dead discards (1000 lb) at F=Fmsy  
  number spr_msy_out;           //spr at F=Fmsy

  vector N_age_msy(1,nages);         //numbers at age for MSY calculations: beginning of yr
  vector N_age_msy_spawn(1,nages);   //numbers at age for MSY calculations: time of peak spawning
  vector L_age_msy(1,nages);         //catch at age for MSY calculations
  vector Z_age_msy(1,nages);         //total mortality at age for MSY calculations
  vector D_age_msy(1,nages);         //discard mortality (dead discards) at age for MSY calculations
  vector F_L_age_msy(1,nages);       //fishing mortality landings (not discards) at age for MSY calculations
  vector F_D_age_msy(1,nages);       //fishing mortality of discards at age for MSY calculations
  vector F_msy(1,n_iter_msy);        //values of full F to be used in equilibrium calculations
  vector spr_msy(1,n_iter_msy);      //reproductive capacity-per-recruit values corresponding to F values in F_msy
  vector R_eq(1,n_iter_msy);         //equilibrium recruitment values corresponding to F values in F_msy
  vector L_eq_klb(1,n_iter_msy);     //equilibrium landings(klb) values corresponding to F values in F_msy
  vector L_eq_knum(1,n_iter_msy);    //equilibrium landings(1000 fish) values corresponding to F values in F_msy
  vector SSB_eq(1,n_iter_msy);       //equilibrium reproductive capacity values corresponding to F values in F_msy
  vector B_eq(1,n_iter_msy);         //equilibrium biomass values corresponding to F values in F_msy
  vector D_eq_klb(1,n_iter_msy);     //equilibrium discards (klb) corresponding to F values in F_msy
  vector D_eq_knum(1,n_iter_msy);    //equilibrium discards (1000s) corresponding to F values in F_msy

  vector FdF_msy(styr,endyr);
  vector SdSSB_msy(styr,endyr);
  number SdSSB_msy_end;
  number FdF_msy_end;
  number FdF_msy_end_mean;          //geometric mean of last 3 yrs  

  vector wgt_wgted_L_klb(1,nages);  //fishery-weighted average weight at age of landings
  vector wgt_wgted_D_klb(1,nages);  //fishery-weighted average weight at age of discards  
  number wgt_wgted_L_denom;         //used in intermediate calculations
  number wgt_wgted_D_denom;         //used in intermediate calculations
  
  number iter_inc_msy;               //increments used to compute msy, equals 1/(n_iter_msy-1)
  
////--------Mortality------------------------------------------------------------------

// Stuff immediately below used only if M is estimated
//  //init_bounded_number M_constant(0.1,0.2,1);  //age-indpendent: used only for MSST
//  vector Mscale_ages(1,max_obs_age);
//  vector Mscale_len(1,max_obs_age);
//  vector Mscale_wgt_g(1,max_obs_age); 
//  vector M_lorenzen(1,max_obs_age);   
//  number cum_surv_1plus;  

  vector M(1,nages);                         //age-dependent natural mortality
  //number M_constant;                         //age-indpendent: used only for MSST
  init_bounded_number M_constant(M_constant_LO,M_constant_HI,M_constant_PH);   //age-indpendent: used only for MSST
  vector M_constant_out(1,8);
 
  number smsy2msstM;                         //scales Smsy to get msst using (1-M). Used only in output.
  number smsy2msst75;                        //scales Smsy to get msst using 75%. Used only in output.  
  
  matrix F(styr,endyr,1,nages);
  vector Fsum(styr,endyr);                   //Full fishing mortality rate by year
  vector Fapex(styr,endyr);                  //Max across ages, fishing mortality rate by year (may differ from Fsum bc of dome-shaped sel 
//  sdreport_vector fullF_sd(styr,endyr);
  matrix Z(styr,endyr,1,nages);

  //init_bounded_number log_avg_F_cL(-10,0.0,1); 
  //init_bounded_dev_vector log_F_dev_cL(styr_cL_L,endyr_cL_L,-10.0,5.0,2);  
  
  init_bounded_number log_avg_F_cL(log_avg_F_cL_LO,log_avg_F_cL_HI,log_avg_F_cL_PH);
  vector log_avg_F_cL_out(1,8);  
  init_bounded_dev_vector log_F_dev_cL(styr_cL_L,endyr_cL_L,log_F_dev_cL_LO,log_F_dev_cL_HI,log_F_dev_cL_PH);
  vector log_F_dev_cL_out(styr_cL_L,endyr_cL_L);
  matrix F_cL(styr,endyr,1,nages);
  vector F_cL_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cL;
  number log_F_dev_end_cL;  

  //init_bounded_number log_avg_F_cP(-10,0.0,1);
  //init_bounded_dev_vector log_F_dev_cP(styr_cP_L,endyr_cP_L,-10.0,5.0,2);
  init_bounded_number log_avg_F_cP(log_avg_F_cP_LO,log_avg_F_cP_HI,log_avg_F_cP_PH);
  vector log_avg_F_cP_out(1,8);  
  init_bounded_dev_vector log_F_dev_cP(styr_cP_L,endyr_cP_L,log_F_dev_cP_LO,log_F_dev_cP_HI,log_F_dev_cP_PH);
  vector log_F_dev_cP_out(styr_cP_L,endyr_cP_L);
  matrix F_cP(styr,endyr,1,nages);
  vector F_cP_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cP;
  number log_F_dev_end_cP;  

  //init_bounded_number log_avg_F_cT(-10,0.0,1);
  //init_bounded_dev_vector log_F_dev_cT(styr_cT_L,endyr_cT_L,-10.0,5.0,2);
  init_bounded_number log_avg_F_cT(log_avg_F_cT_LO,log_avg_F_cT_HI,log_avg_F_cT_PH);
  vector log_avg_F_cT_out(1,8);  
  init_bounded_dev_vector log_F_dev_cT(styr_cT_L,endyr_cT_L,log_F_dev_cT_LO,log_F_dev_cT_HI,log_F_dev_cT_PH);
  vector log_F_dev_cT_out(styr_cT_L,endyr_cT_L);
  matrix F_cT(styr,endyr,1,nages);
  vector F_cT_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cT;
  number log_F_dev_end_cT;  
 
  //init_bounded_number log_avg_F_HB(-10.0,0.0,1);  
  //init_bounded_dev_vector log_F_dev_HB(styr_HB_L,endyr_HB_L,-10.0,5.0,2);
  init_bounded_number log_avg_F_HB(log_avg_F_HB_LO,log_avg_F_HB_HI,log_avg_F_HB_PH);
  vector log_avg_F_HB_out(1,8); 
  init_bounded_dev_vector log_F_dev_HB(styr_HB_L,endyr_HB_L,log_F_dev_HB_LO,log_F_dev_HB_HI,log_F_dev_HB_PH);    
  vector log_F_dev_HB_out(styr_HB_L,endyr_HB_L);
  matrix F_HB(styr,endyr,1,nages);
  vector F_HB_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_init_HB;
  number log_F_dev_end_HB;    

  //init_bounded_number log_avg_F_mrip(-10,0.0,1);
  //init_bounded_dev_vector log_F_dev_mrip(styr_mrip_L,endyr_mrip_L,-10.0,5.0,2);
  init_bounded_number log_avg_F_mrip(log_avg_F_mrip_LO,log_avg_F_mrip_HI,log_avg_F_mrip_PH);
  vector log_avg_F_mrip_out(1,8); 
  init_bounded_dev_vector log_F_dev_mrip(styr_mrip_L,endyr_mrip_L,log_F_dev_mrip_LO,log_F_dev_mrip_HI,log_F_dev_mrip_PH);    
  vector log_F_dev_mrip_out(styr_mrip_L,endyr_mrip_L);
  matrix F_mrip(styr,endyr,1,nages);
  vector F_mrip_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_mrip;
  number log_F_dev_end_mrip;  
  
  init_bounded_number F_init_ratio(0.1,1.5,1); //scales initial F, which is geometric mean first three yrs
  //init_bounded_number F_init_ratio(F_init_ratio_LO,F_init_ratio_HI,F_init_ratio_PH); //number F_init_ratio;
  //vector F_init_ratio_out(1,8); 

  //--Discard mortality stuff------------------------------------------------------------------------------
  //init_bounded_number log_avg_F_comm_D(-10.0,0.0,1);
  //init_bounded_dev_vector log_F_dev_comm_D(styr_cL_D,endyr_cL_D,-10.0,5.0,2);
  init_bounded_number log_avg_F_comm_D(log_avg_F_comm_D_LO,log_avg_F_comm_D_HI,log_avg_F_comm_D_PH);
  vector log_avg_F_comm_D_out(1,8);  
  init_bounded_dev_vector log_F_dev_comm_D(styr_cL_D,endyr_cL_D,log_F_dev_comm_D_LO,log_F_dev_comm_D_HI,log_F_dev_comm_D_PH);
  vector log_F_dev_comm_D_out(styr_cL_D,endyr_cL_D);
  matrix F_comm_D(styr,endyr,1,nages);
  vector F_comm_D_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_comm_D2;   //avg log deviations in reg period 2 (for estimation 1984-1992, prior to data) 
  number log_F_dev_end_comm_D;  
  
  //init_bounded_number log_avg_F_HB_D(-10.0,0.0,1);
  //init_bounded_dev_vector log_F_dev_HB_D(styr_HB_D,endyr_HB_D,-10.0,5.0,2);
  init_bounded_number log_avg_F_HB_D(log_avg_F_HB_D_LO,log_avg_F_HB_D_HI,log_avg_F_HB_D_PH);
  vector log_avg_F_HB_D_out(1,8);  
  init_bounded_dev_vector log_F_dev_HB_D(styr_HB_D,endyr_HB_D,log_F_dev_HB_D_LO,log_F_dev_HB_D_HI,log_F_dev_HB_D_PH);
  vector log_F_dev_HB_D_out(styr_HB_D,endyr_HB_D);
  matrix F_HB_D(styr,endyr,1,nages);
  vector F_HB_D_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_end_HB_D;  
       
  //init_bounded_number log_avg_F_mrip_D(-10.0,0.0,1);
  //init_bounded_dev_vector log_F_dev_mrip_D(styr_mrip_D,endyr_mrip_D,-10.0,5.0,2);
  init_bounded_number log_avg_F_mrip_D(log_avg_F_mrip_D_LO,log_avg_F_mrip_D_HI,log_avg_F_mrip_D_PH);
  vector log_avg_F_mrip_D_out(1,8);  
  init_bounded_dev_vector log_F_dev_mrip_D(styr_mrip_D,endyr_mrip_D,log_F_dev_mrip_D_LO,log_F_dev_mrip_D_HI,log_F_dev_mrip_D_PH);
  vector log_F_dev_mrip_D_out(styr_mrip_D,endyr_mrip_D);
  matrix F_mrip_D(styr,endyr,1,nages);
  vector F_mrip_D_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_mrip_D;
  number log_F_dev_end_mrip_D;  

  number Dmort_HL;
  number Dmort_HB_HL;
  number Dmort_GR_HL;
  number Dmort_cP1;
  number Dmort_cP2;  

//---Per-recruit stuff----------------------------------------------------------------------------------
  vector N_age_spr(1,nages);         //numbers at age for SPR calculations: beginning of year
  vector N_age_spr_spawn(1,nages);   //numbers at age for SPR calculations: at time of peak spawning  
  vector L_age_spr(1,nages);         //catch at age for SPR calculations
  vector Z_age_spr(1,nages);         //total mortality at age for SPR calculations
  vector spr_static(styr,endyr);     //vector of static SPR values by year
  vector F_L_age_spr(1,nages);       //fishing mortality of landings (not discards) at age for SPR calculations
  vector F_spr(1,n_iter_spr);        //values of full F to be used in per-recruit calculations
  vector spr_spr(1,n_iter_spr);      //reproductive capacity-per-recruit values corresponding to F values in F_spr
  vector L_spr(1,n_iter_spr);        //landings(lb)-per-recruit (ypr) values corresponding to F values in F_spr

  vector N_spr_F0(1,nages);          //Used to compute spr at F=0: at time of peak spawning
  vector N_bpr_F0(1,nages);          //Used to compute bpr at F=0: at start of year  
  vector N_spr_initial(1,nages);     //Initial spawners per recruit at age given initial F
  vector N_initial_eq(1,nages);      //Initial equilibrium abundance at age
  vector F_initial(1,nages);         //initial F at age
  vector Z_initial(1,nages);         //initial Z at age
  number spr_initial;                //initial spawners per recruit
  number spr_F0;                     //Spawning biomass per recruit at F=0
  number bpr_F0;                     //Biomass per recruit at F=0

  number iter_inc_spr;               //increments used to compute msy, equals max_F_spr_msy/(n_iter_spr-1)

////-------SDNR output-----------------------------------------------------------------------------
 
  number sdnr_lc_Mbft;
  number sdnr_lc_Mcvt;
  number sdnr_lc_cL;
  number sdnr_lc_cP;
  number sdnr_lc_HB;  
  number sdnr_lc_HB_D;  
  number sdnr_lc_mrip;  
 
  number sdnr_ac_Mbft;
  number sdnr_ac_Mcvt;
  number sdnr_ac_cL; 
  number sdnr_ac_cP;  
  number sdnr_ac_HB;
  //number sdnr_ac_mrip;
  
  number sdnr_I_Mbft;
  number sdnr_I_Mcvt;
  number sdnr_I_cL;
  number sdnr_I_HB;  
  //number sdnr_I_Vid;    
  
//-------Objective function components-----------------------------------------------------------------------------
  number w_L;
  number w_D;

  number w_lc_Mbft;
  number w_lc_cL;
  number w_lc_cP;
  number w_lc_HB;
  number w_lc_HB_D;
  number w_lc_mrip;

  number w_ac_Mbft; 
  number w_ac_Mcvt; 
  number w_ac_cL; 
  number w_ac_cP;  
  number w_ac_HB;
  //number w_ac_mrip;
  
  number w_I_Mbft;  
  number w_I_Mcvt;
  //number w_I_Vid;
  number w_I_cL;
  number w_I_HB;
  //number w_I_HBD;
  
  number w_rec;
  number w_rec_early;
  number w_rec_end;
  number w_fullF;  
  number w_Ftune;
//  number w_cvlen_dev;
//  number w_cvlen_diff;
  
  number f_Mbft_cpue;
  number f_Mcvt_cpue;
  number f_Vid_cpue;
  number f_cL_cpue;
  number f_HB_cpue;
  //number f_HBD_cpue;
 
  number f_cL_L;  
  number f_cP_L;
  number f_cT_L;
  number f_HB_L;
  number f_mrip_L;

  number f_comm_D; 
  number f_HB_D;
  number f_mrip_D;

  number f_Mbft_lenc;   
  number f_cL_lenc;
  number f_cP_lenc;
  number f_HB_lenc;
  number f_HB_D_lenc;
  number f_mrip_lenc;

  number f_Mbft_agec;
  number f_Mcvt_agec;  
  number f_cL_agec;  
  number f_cP_agec;    
  number f_HB_agec; 
  number f_mrip_agec;
  
  number f_cL_RW_cpue; //random walk component of indices
  number f_HB_RW_cpue;
  number f_HBD_RW_cpue;   
  
//Penalties and constraints. Not all are used.
  number f_rec_dev;                //weight on recruitment deviations to fit S-R curve
  number f_rec_dev_early;          //extra weight on deviations in first recruitment stanza
  number f_rec_dev_end;            //extra weight on deviations in first recruitment stanza
  number f_Ftune;                  //penalty for tuning F in Ftune yr.  Not applied in final optimization phase.
  number f_fullF_constraint;       //penalty for Fapex>X
  number f_priors;                  //prior information on parameters
  //number f_cvlen_dev_constraint; //deviation penalty on cv's of length at age
  //number f_cvlen_diff_constraint;//first diff penalty on cv's of length at age
  
  objective_function_value fval;
  number fval_data;
  number grad_max;
  
//--Dummy variables ----
  number denom;                   //denominator used in some calculations
  number numer;                   //numerator used in some calculations
     
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
INITIALIZATION_SECTION


//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
GLOBALS_SECTION
  #include "admodel.h"          // Include AD class definitions
  #include "admb2r.cpp"    // Include S-compatible output functions (needs preceding)
	time_t start,finish;
	long hour,minute,second;	
	double elapsed_time;
	
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
RUNTIME_SECTION
 maximum_function_evaluations 10000, 20000, 30000, 50000, 100000, 100000, 100000;
 convergence_criteria 1e-2, 1e-2,1e-2, 1e-4;
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PRELIMINARY_CALCS_SECTION

//// Set values of fixed parameters or set initial guess of estimated parameters
 sqrt2pi=sqrt(2.*3.14159265);
 g2mt=0.000001;         //conversion of grams to metric tons
 g2kg=0.001;            //conversion of grams to kg 
 mt2klb=2.20462;        //conversion of metric tons to 1000 lb 
 mt2lb=mt2klb*1000.0;   //conversion of metric tons to lb
 g2klb=g2mt*mt2klb;     //conversion of grams to 1000 lb 
 dzero=0.00001;         
 huge_number=1.0e+10;   
 onehalf=0.5;

  Dmort_HL=set_Dmort_HL;
  Dmort_HB_HL=set_Dmort_HB_HL;
  Dmort_GR_HL=set_Dmort_GR_HL;
  Dmort_cP1=set_Dmort_cP1;
  Dmort_cP2=set_Dmort_cP2;

  //values used for weighting selex and avg weight of comm discards in yrs with quotas
  //geometric mean of last three yrs
  //avg weights of landings were near 1 lb, so those values are left in weight
  Dopen_cL=Dmort_HL*pow((obs_cL_released(2014)*obs_cL_released(2015)*obs_cL_released(endyr_cL_D)),(1.0/3.0));
  Dclosed_cL=Dmort_HL*((obs_cL_closed_released(2013)+obs_cL_closed_released(2014)+obs_cL_closed_released(2015)+obs_cL_closed_released(endyr_cL_D))/4.0);
  //Dclosed_cL=Dmort_HL*pow((obs_cL_closed_released(2014)*obs_cL_closed_released(2015)*obs_cL_closed_released(endyr_cL_D)),(1.0/3.0));
  Lopen_cL=Dmort_HL*pow((obs_cL_L(2014)*obs_cL_L(2015)*obs_cL_L(endyr_cL_L)),(1.0/3.0));
  
  Dopen_cP=Dmort_cP2*pow((obs_cP_released(2014)*obs_cP_released(2015)*obs_cP_released(endyr_cP_D)),(1.0/3.0));
  Dclosed_cP=Dmort_cP2*(obs_cP_closed_released(2013)+(obs_cP_closed_released(2014)+obs_cP_closed_released(2015)+obs_cP_closed_released(endyr_cP_D))/4.0);
  //Dclosed_cP=Dmort_cP2*pow((obs_cP_closed_released(2014)*obs_cP_closed_released(2015)*obs_cP_closed_released(endyr_cP_D)),(1.0/3.0));
  Lopen_cP=Dmort_cP2*pow((obs_cP_L(2014)*obs_cP_L(2015)*obs_cP_L(endyr_cP_L)),(1.0/3.0));
  cout << "Dclosed_cP" <<Dclosed_cP<< endl;
  cout << "Dclosed_cL" <<Dclosed_cL<< endl;
  cout << "Dopen_cP" <<Dopen_cP<< endl;
  cout << "Dopen_cL" <<Dopen_cL<< endl;
  
  D_sum_cLcP=Dopen_cL+Dclosed_cL+Dopen_cP+Dclosed_cP;

  Dprop_comm_sel_D=(Dopen_cL + Dopen_cP + Dclosed_cL*(Dopen_cL/(Dopen_cL+Lopen_cL)) +
                    Dclosed_cP*(Dopen_cP/(Dopen_cP+Lopen_cP)))/D_sum_cLcP; 
  Dprop_comm_sel_cL=Dclosed_cL*(Lopen_cL/(Dopen_cL+Lopen_cL))/D_sum_cLcP; 
  Dprop_comm_sel_cP=Dclosed_cP*(Lopen_cP/(Dopen_cP+Lopen_cP))/D_sum_cLcP; 
  cout << Dprop_comm_sel_cP <<Dprop_comm_sel_cP<< endl;
  //cout<<"prop cP"<<Dprop_comm_sel_cP;
  //cout<<"prop cL"<<Dprop_comm_sel_cL;
  //cout<<"overall prop"<<Dprop_comm_sel_D;
  //discards values for fitting, include discard mortality
  
  obs_cL_D=Dmort_HL*obs_cL_released;
  obs_cL_D(styr_cL_closed_D,endyr_cL_closed_D)+=Dmort_HL*obs_cL_closed_released;
  
  obs_cP_D(styr_cP_D,2006)=Dmort_cP1*obs_cP_released(styr_cP_D,2006);
  obs_cP_D(2007,endyr_cP_D)=Dmort_cP2*obs_cP_released(2007,endyr_cP_D);
  obs_cP_D(styr_cP_closed_D,endyr_cP_closed_D)+=Dmort_cP2*obs_cP_closed_released;
  
  obs_comm_D=obs_cL_D+obs_cP_D;
  
  obs_HB_D=Dmort_HB_HL*obs_HB_released;
  obs_mrip_D=Dmort_GR_HL*obs_mrip_released;
  
  comm_D_cv=cL_D_cv;
 
  Linf=set_Linf(1);
  K=set_K(1);
  t0=set_t0(1);

//  age_limit_8in=t0-log(1.0-limit_8in/Linf)/K; //age at size limit: 8" limit;
//  age_limit_10in=t0-log(1.0-limit_10in/Linf)/K; //age at size limit: 10" limit;
//  age_limit_12in=t0-log(1.0-limit_12in/Linf)/K; //age at size limit: 12" limit;

  M=set_M; 
  M_constant=set_M_constant(1);
//  for (iage=1;iage<=max_obs_age;iage++){Mscale_ages(iage)=iage;}
  
  
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);
  
  log_dm_Mbft_lc=set_log_dm_Mbft_lc(1);
  log_dm_Mcvt_lc=set_log_dm_Mcvt_lc(1);
  log_dm_cL_lc=set_log_dm_cL_lc(1);
  log_dm_cP_lc=set_log_dm_cP_lc(1);
  log_dm_HB_lc=set_log_dm_HB_lc(1);
  log_dm_HB_D_lc=set_log_dm_HB_D_lc(1);
  log_dm_mrip_lc=set_log_dm_mrip_lc(1);
  
  log_dm_Mbft_ac=set_log_dm_Mbft_ac(1);
  log_dm_Mcvt_ac=set_log_dm_Mcvt_ac(1);
  log_dm_cL_ac=set_log_dm_cL_ac(1);
  log_dm_cP_ac=set_log_dm_cP_ac(1);
  log_dm_HB_ac=set_log_dm_HB_ac(1);
  //log_dm_mrip_ac=set_log_dm_mrip_ac(1);
    
  log_q_Mbft=set_logq_Mbft(1);  
  log_q_Mcvt=set_logq_Mcvt(1);
  log_q_cL=set_logq_cL(1);
  log_q_HB=set_logq_HB(1);
  //log_q_HBD=set_logq_HBD(1);
  q_rate=set_q_rate;
  q_rate_fcn_cL=1.0;
  q_rate_fcn_HB=1.0;
  //q_rate_fcn_HBD=1.0;
  q_DD_beta=set_q_DD_beta;
  q_DD_fcn=1.0;
  q_RW_log_dev_cL.initialize();
  q_RW_log_dev_HB.initialize();
  //q_RW_log_dev_HBD.initialize();
 	//cout <<"q_rate"<< set_q_rate_phase << endl;
  if (set_q_rate_phase<0 & q_rate!=0.0)
  {
      for (iyear=styr_cL_cpue; iyear<=endyr_cL_cpue; iyear++)
      {   if (iyear>styr_cL_cpue & iyear <=2003) 
          {//q_rate_fcn_cL(iyear)=(1.0+q_rate)*q_rate_fcn_cL(iyear-1); //compound
             q_rate_fcn_cL(iyear)=(1.0+(iyear-styr_cL_cpue)*q_rate)*q_rate_fcn_cL(styr_cL_cpue);  //linear
          }
          if (iyear>2003) {q_rate_fcn_cL(iyear)=q_rate_fcn_cL(iyear-1);} 
      }   
      for (iyear=styr_HB_cpue; iyear<=endyr_HB_cpue; iyear++)
      {   if (iyear>styr_HB_cpue & iyear <=2003) 
          {//q_rate_fcn_HB(iyear)=(1.0+q_rate)*q_rate_fcn_HB(iyear-1); //compound
             q_rate_fcn_HB(iyear)=(1.0+(iyear-styr_HB_cpue)*q_rate)*q_rate_fcn_HB(styr_HB_cpue);  //linear
          }
          if (iyear>2003) {q_rate_fcn_HB(iyear)=q_rate_fcn_HB(iyear-1);} 
      }  
      //for (iyear=styr_HBD_cpue; iyear<=endyr_HBD_cpue; iyear++)
      //{   if (iyear>styr_HBD_cpue & iyear <=2003) 
       //   {//q_rate_fcn_HBD(iyear)=(1.0+q_rate)*q_rate_fcn_HBD(iyear-1); //compound
       //      q_rate_fcn_HBD(iyear)=(1.0+(iyear-styr_HBD_cpue)*q_rate)*q_rate_fcn_HBD(styr_HBD_cpue);  //linear
       //   }
       //   if (iyear>2003) {q_rate_fcn_HBD(iyear)=q_rate_fcn_HBD(iyear-1);} 
      //}
  } //end q_rate conditional      


  L_hb_bias=set_L_hb_bias;
  L_mrip_bias=set_L_mrip_bias;
  L_comm_bias=set_L_comm_bias;

  w_L=set_w_L;
  w_D=set_w_D;
  
  w_lc_Mbft=set_w_lc_Mbft;
  w_lc_cL=set_w_lc_cL;
  w_lc_cP=set_w_lc_cP;
  w_lc_HB=set_w_lc_HB;
  w_lc_HB_D=set_w_lc_HB_D;
  w_lc_mrip=set_w_lc_mrip;
  
  w_ac_Mbft=set_w_ac_Mbft;
  w_ac_Mcvt=set_w_ac_Mcvt;
  w_ac_cL=set_w_ac_cL;
  w_ac_cP=set_w_ac_cP;
  w_ac_HB=set_w_ac_HB;
  //w_ac_mrip=set_w_ac_mrip;  

  w_I_Mcvt=set_w_I_Mcvt;
  w_I_Mbft=set_w_I_Mbft;
  //w_I_Vid=set_w_I_Vid;
  w_I_cL=set_w_I_cL;
  w_I_HB=set_w_I_HB;
  //w_I_HBD=set_w_I_HBD;

  w_rec=set_w_rec;
  w_fullF=set_w_fullF;
  w_rec_early=set_w_rec_early;
  w_rec_end=set_w_rec_end;
  w_Ftune=set_w_Ftune;
  //w_cvlen_dev=set_w_cvlen_dev;
  //w_cvlen_diff=set_w_cvlen_diff;

  log_avg_F_cL=set_log_avg_F_cL(1);  
  log_avg_F_cP=set_log_avg_F_cP(1);
  log_avg_F_cT=set_log_avg_F_cT(1);  
  log_avg_F_HB=set_log_avg_F_HB(1);
  log_avg_F_mrip=set_log_avg_F_mrip(1);  
  F_init_ratio=set_F_init_ratio;

  
  log_avg_F_comm_D=set_log_avg_F_comm_D(1); 
  log_avg_F_HB_D=set_log_avg_F_HB_D(1);
  log_avg_F_mrip_D=set_log_avg_F_mrip_D(1);  
 
  len_cv_val=set_len_cv(1);

  log_R0=set_log_R0(1);

  selpar_A50_Mbft= set_selpar_A50_Mbft(1);
  selpar_slope_Mbft=set_selpar_slope_Mbft(1); 

  selpar_A50_Mcvt=set_selpar_A50_Mcvt(1);
  selpar_slope_Mcvt=set_selpar_slope_Mcvt(1); 

  selpar_A50_cL2=set_selpar_A50_cL2(1);
  selpar_slope_cL2=set_selpar_slope_cL2(1); 
  selpar_A50_cL3=set_selpar_A50_cL3(1);
  selpar_slope_cL3=set_selpar_slope_cL3(1); 
  selpar_A50_cL4=set_selpar_A50_cL4(1);
  selpar_slope_cL4=set_selpar_slope_cL4(1); 

  selpar_A50_cP2=set_selpar_A50_cP2(1);
  selpar_slope_cP2=set_selpar_slope_cP2(1); 
  selpar_A50_cP3=set_selpar_A50_cP3(1);
  selpar_slope_cP3=set_selpar_slope_cP3(1); 
  selpar_A50_cP4=set_selpar_A50_cP4(1);
  selpar_slope_cP4=set_selpar_slope_cP4(1);

  selpar_A50_HB1=set_selpar_A50_HB1(1);
  selpar_slope_HB1=set_selpar_slope_HB1(1); 
  selpar_A50_HB2=set_selpar_A50_HB2(1);
  selpar_slope_HB2=set_selpar_slope_HB2(1); 
  selpar_A50_HB3=set_selpar_A50_HB3(1);
  selpar_slope_HB3=set_selpar_slope_HB3(1); 
  selpar_A50_HB4=set_selpar_A50_HB4(1);
  selpar_slope_HB4=set_selpar_slope_HB4(1); 
  selpar_A50_HB5=set_selpar_A50_HB5(1);
  selpar_slope_HB5=set_selpar_slope_HB5(1); 

  selpar_A50_mrip1=set_selpar_A50_mrip1(1);
  selpar_slope_mrip1=set_selpar_slope_mrip1(1); 
  selpar_A50_mrip2=set_selpar_A50_mrip2(1);
  selpar_slope_mrip2=set_selpar_slope_mrip2(1); 
  selpar_A50_mrip3=set_selpar_A50_mrip3(1);
  selpar_slope_mrip3=set_selpar_slope_mrip3(1); 
  selpar_A50_mrip4=set_selpar_A50_mrip4(1);
  selpar_slope_mrip4=set_selpar_slope_mrip4(1); 
  selpar_A50_mrip5=set_selpar_A50_mrip5(1);
  selpar_slope_mrip5=set_selpar_slope_mrip5(1); 


  selpar_Age0_HB_D_logit=set_selpar_Age0_HB_D_logit(1);  
  selpar_Age1_HB_D_logit=set_selpar_Age1_HB_D_logit(1); 
  selpar_Age2_HB_D_logit=set_selpar_Age2_HB_D_logit(1);

  selpar_A50_HBD4=set_selpar_A50_HBD4(1);
  selpar_slope_HBD4=set_selpar_slope_HBD4(1); 
  selpar_A502_HBD4=set_selpar_A502_HBD4(1);
  selpar_slope2_HBD4=set_selpar_slope2_HBD4(1); 
  selpar_A50_HBD5=set_selpar_A50_HBD5(1);
  selpar_slope_HBD5=set_selpar_slope_HBD5(1);
  selpar_A502_HBD5=set_selpar_A502_HBD5(1);
  selpar_slope2_HBD5=set_selpar_slope2_HBD5(1);
   
 SSB_msy_out=0.0;

 iter_inc_msy=max_F_spr_msy/(n_iter_msy-1);
 iter_inc_spr=max_F_spr_msy/(n_iter_spr-1); 

 maturity_f=maturity_f_obs;
 prop_f=prop_f_obs;

 p_lenc_cL2=set_p_lenc_cL2;
 p_lenc_cL3=set_p_lenc_cL3; 
 p_lenc_cP2=set_p_lenc_cP2;
 p_lenc_cP3=set_p_lenc_cP3; 
 p_lenc_cT2=set_p_lenc_cT2;
 p_lenc_cT3=set_p_lenc_cT3;  
 p_lenc_HB2=set_p_lenc_HB2;
 p_lenc_HB3=set_p_lenc_HB3;
 p_lenc_HB4=set_p_lenc_HB4; 
 p_lenc_HB5=set_p_lenc_HB5; 
 p_lenc_mrip2=set_p_lenc_mrip2;
 p_lenc_mrip3=set_p_lenc_mrip3;
 p_lenc_mrip4=set_p_lenc_mrip4;
 p_lenc_mrip5=set_p_lenc_mrip5;
  
 p_lenc_comm_D2=set_p_lenc_comm_D2; 
 p_lenc_comm_D3=set_p_lenc_comm_D3;
 p_lenc_comm_D4=set_p_lenc_comm_D4;
 p_lenc_HB_D2=set_p_lenc_HB_D2;
 p_lenc_HB_D3=set_p_lenc_HB_D3; 
 p_lenc_HB_D4=set_p_lenc_HB_D4;
 p_lenc_HB_D5=set_p_lenc_HB_D5; 
 p_lenc_mrip_D1=set_p_lenc_mrip_D1; 
 p_lenc_mrip_D2=set_p_lenc_mrip_D2;
 p_lenc_mrip_D3=set_p_lenc_mrip_D3; 
 p_lenc_mrip_D4=set_p_lenc_mrip_D4; 
 p_lenc_mrip_D5=set_p_lenc_mrip_D5;
 
 lenbins_all(1,nlenbins)=lenbins(1,nlenbins);  
 for (iyear=1;iyear<=nlenbins_plus; iyear++) {lenbins_all(nlenbins+iyear)=lenbins_plus(iyear);} 

 //multiplicative bias for early rec data
 obs_HB_L_wbias(styr_HB_L,endyr_L_HB_bias)=L_hb_bias*obs_HB_L(styr_HB_L,endyr_L_HB_bias); 
 obs_HB_L_wbias((endyr_L_HB_bias+1),endyr_HB_L)=obs_HB_L((endyr_L_HB_bias+1),endyr_HB_L); 

 obs_mrip_L_wbias(styr_mrip_L,endyr_L_mrip_bias)=L_hb_bias*obs_mrip_L(styr_mrip_L,endyr_L_mrip_bias); 
 obs_mrip_L_wbias((endyr_L_mrip_bias+1),endyr_mrip_L)=obs_mrip_L((endyr_L_mrip_bias+1),endyr_mrip_L); 

 
//Fill in sample sizes of comps, possibly sampled in nonconsec yrs 
//Used primarily for output in R object   
  nsamp_Mbft_lenc_allyr=missing;//"missing" defined in admb2r.cpp
  nsamp_cL_lenc_allyr=missing;  
  nsamp_cP_lenc_allyr=missing;
  nsamp_HB_lenc_allyr=missing;
  nsamp_HB_D_lenc_allyr=missing;
  nsamp_mrip_lenc_allyr=missing;
  nsamp_Mbft_agec_allyr=missing;
  nsamp_Mcvt_agec_allyr=missing;  
  nsamp_cL_agec_allyr=missing;
  nsamp_cP_agec_allyr=missing;      
  nsamp_HB_agec_allyr=missing;
  //nsamp_mrip_agec_allyr=missing;

  nfish_Mbft_lenc_allyr=missing;//"missing" defined in admb2r.cpp
  nfish_cL_lenc_allyr=missing;  
  nfish_cP_lenc_allyr=missing;
  nfish_HB_lenc_allyr=missing;
  nfish_HB_D_lenc_allyr=missing;
  nfish_mrip_lenc_allyr=missing;
  nfish_Mbft_agec_allyr=missing;
  nfish_Mcvt_agec_allyr=missing;  
  nfish_cL_agec_allyr=missing;
  nfish_cP_agec_allyr=missing;      
  nfish_HB_agec_allyr=missing;
  //nfish_mrip_agec_allyr=missing;
          

      for (iyear=1; iyear<=nyr_Mbft_lenc; iyear++)
         {  if (nsamp_Mbft_lenc(iyear)>maxSS_lenc)
             {nsamp_Mbft_lenc(iyear)=maxSS_lenc;}
            if (nsamp_Mbft_lenc(iyear)>=minSS_lenc)
             {nsamp_Mbft_lenc_allyr(yrs_Mbft_lenc(iyear))=nsamp_Mbft_lenc(iyear);
              nfish_Mbft_lenc_allyr(yrs_Mbft_lenc(iyear))=nfish_Mbft_lenc(iyear);}}
            
      for (iyear=1; iyear<=nyr_cL_lenc; iyear++)
         {  if (nsamp_cL_lenc(iyear)>maxSS_lenc)
             {nsamp_cL_lenc(iyear)=maxSS_lenc;}
            if (nsamp_cL_lenc(iyear)>=minSS_lenc)
             {nsamp_cL_lenc_allyr(yrs_cL_lenc(iyear))=nsamp_cL_lenc(iyear);
              nfish_cL_lenc_allyr(yrs_cL_lenc(iyear))=nfish_cL_lenc(iyear);}}
         
      for (iyear=1; iyear<=nyr_cP_lenc; iyear++)
         {  if (nsamp_cP_lenc(iyear)>maxSS_lenc)
             {nsamp_cP_lenc(iyear)=maxSS_lenc;}
            if (nsamp_cP_lenc(iyear)>=minSS_lenc)
             {nsamp_cP_lenc_allyr(yrs_cP_lenc(iyear))=nsamp_cP_lenc(iyear);
              nfish_cP_lenc_allyr(yrs_cP_lenc(iyear))=nfish_cP_lenc(iyear);}}
      
      for (iyear=1; iyear<=nyr_HB_lenc; iyear++)
         {  if (nsamp_HB_lenc(iyear)>maxSS_lenc)
             {nsamp_HB_lenc(iyear)=maxSS_lenc;}
            if (nsamp_HB_lenc(iyear)>=minSS_lenc)
             {nsamp_HB_lenc_allyr(yrs_HB_lenc(iyear))=nsamp_HB_lenc(iyear);
              nfish_HB_lenc_allyr(yrs_HB_lenc(iyear))=nfish_HB_lenc(iyear);}}

      for (iyear=1; iyear<=nyr_HB_D_lenc; iyear++)
         {  if (nsamp_HB_D_lenc(iyear)>maxSS_lenc)
             {nsamp_HB_D_lenc(iyear)=maxSS_lenc;}
            if (nsamp_HB_D_lenc(iyear)>=minSS_lenc)
             {nsamp_HB_D_lenc_allyr(yrs_HB_D_lenc(iyear))=nsamp_HB_D_lenc(iyear);
              nfish_HB_D_lenc_allyr(yrs_HB_D_lenc(iyear))=nfish_HB_D_lenc(iyear);}}
         
      for (iyear=1; iyear<=nyr_mrip_lenc; iyear++)
         {  if (nsamp_mrip_lenc(iyear)>maxSS_lenc)
             {nsamp_mrip_lenc(iyear)=maxSS_lenc;}
            if (nsamp_mrip_lenc(iyear)>=minSS_lenc)
             {nsamp_mrip_lenc_allyr(yrs_mrip_lenc(iyear))=nsamp_mrip_lenc(iyear);
              nfish_mrip_lenc_allyr(yrs_mrip_lenc(iyear))=nfish_mrip_lenc(iyear);}}

                  
      for (iyear=1; iyear<=nyr_Mbft_agec; iyear++)
         {  if (nsamp_Mbft_agec(iyear)>maxSS_agec)
             {nsamp_Mbft_agec(iyear)=maxSS_agec;}
            if (nsamp_Mbft_agec(iyear)>=minSS_agec)
             {nsamp_Mbft_agec_allyr(yrs_Mbft_agec(iyear))=nsamp_Mbft_agec(iyear);
              nfish_Mbft_agec_allyr(yrs_Mbft_agec(iyear))=nfish_Mbft_agec(iyear);}} 

      for (iyear=1; iyear<=nyr_Mcvt_agec; iyear++)
         {  if (nsamp_Mcvt_agec(iyear)>maxSS_agec)
             {nsamp_Mcvt_agec(iyear)=maxSS_agec;}
            if (nsamp_Mcvt_agec(iyear)>=minSS_agec)
             {nsamp_Mcvt_agec_allyr(yrs_Mcvt_agec(iyear))=nsamp_Mcvt_agec(iyear);
              nfish_Mcvt_agec_allyr(yrs_Mcvt_agec(iyear))=nfish_Mcvt_agec(iyear);}} 

      for (iyear=1; iyear<=nyr_cL_agec; iyear++)
         {  if (nsamp_cL_agec(iyear)>maxSS_agec)
             {nsamp_cL_agec(iyear)=maxSS_agec;}
            if (nsamp_cL_agec(iyear)>=minSS_agec)
             {nsamp_cL_agec_allyr(yrs_cL_agec(iyear))=nsamp_cL_agec(iyear);
              nfish_cL_agec_allyr(yrs_cL_agec(iyear))=nfish_cL_agec(iyear);}} 
            
      for (iyear=1; iyear<=nyr_cP_agec; iyear++)
         {  if (nsamp_cP_agec(iyear)>maxSS_agec)
             {nsamp_cP_agec(iyear)=maxSS_agec;}
            if (nsamp_cP_agec(iyear)>=minSS_agec)
             {nsamp_cP_agec_allyr(yrs_cP_agec(iyear))=nsamp_cP_agec(iyear);
              nfish_cP_agec_allyr(yrs_cP_agec(iyear))=nfish_cP_agec(iyear);}}
             
      for (iyear=1; iyear<=nyr_HB_agec; iyear++)
         {  if (nsamp_HB_agec(iyear)>maxSS_agec)
             {nsamp_HB_agec(iyear)=maxSS_agec;}
            if (nsamp_HB_agec(iyear)>=minSS_agec)
             {nsamp_HB_agec_allyr(yrs_HB_agec(iyear))=nsamp_HB_agec(iyear);
              nfish_HB_agec_allyr(yrs_HB_agec(iyear))=nfish_HB_agec(iyear);}}
             
      //for (iyear=1; iyear<=nyr_mrip_agec; iyear++)
      //   {  if (nsamp_mrip_agec(iyear)>maxSS_agec)
      //       {nsamp_mrip_agec(iyear)=maxSS_agec;}
      //      if (nsamp_mrip_agec(iyear)>=minSS_agec)
      //       {nsamp_mrip_agec_allyr(yrs_mrip_agec(iyear))=nsamp_mrip_agec(iyear);
      //        nfish_mrip_agec_allyr(yrs_mrip_agec(iyear))=nfish_mrip_agec(iyear);}}  

//fill in Fs for msy and per-recruit analyses
  F_msy(1)=0.0;  
    for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
    for (ff=2;ff<=n_iter_spr;ff++){F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}


//fill in F's, Catch matrices, and log rec dev with zero's
  F_cL.initialize(); 
  F_cP.initialize(); 
  F_cT.initialize(); 
  F_HB.initialize(); 
  F_mrip.initialize();
  
  F_comm_D.initialize(); 
  F_HB_D.initialize(); 
  F_mrip_D.initialize();  
      
  L_cL_num.initialize(); 
  L_cP_num.initialize(); 
  L_cT_num.initialize(); 
  L_HB_num.initialize(); 
  L_mrip_num.initialize();
  
  D_comm_num.initialize(); 
  D_HB_num.initialize(); 
  D_mrip_num.initialize();
      
  F_cL_out.initialize(); 
  F_cP_out.initialize(); 
  F_cT_out.initialize(); 
  F_HB_out.initialize(); 
  F_mrip_out.initialize();

  F_comm_D_out.initialize(); 
  F_HB_D_out.initialize(); 
  F_mrip_D_out.initialize();
  
  pred_cL_L_klb.initialize();
  pred_cP_L_klb.initialize();
  pred_cT_L_klb.initialize(); 
  pred_HB_L_klb.initialize(); 
  pred_mrip_L_klb.initialize();
  pred_cL_L_knum.initialize();
  pred_cP_L_knum.initialize();
  pred_cT_L_knum.initialize(); 
  pred_HB_L_knum.initialize(); 
  pred_mrip_L_knum.initialize();
  
  pred_comm_D_klb.initialize(); 
  pred_HB_D_klb.initialize();
  pred_mrip_D_klb.initialize();
  pred_comm_D_knum.initialize(); 
  pred_HB_D_knum.initialize(); 
  pred_mrip_D_knum.initialize();
 
  sel_Mcvt.initialize(); 
  sel_Mbft.initialize();  
  sel_cL.initialize(); 
  sel_cP.initialize(); 
  sel_cT.initialize();  
  sel_HB.initialize(); 
  sel_mrip.initialize();  
  sel_comm_D.initialize(); 
  sel_HB_D.initialize(); 
  sel_mrip_D.initialize();  
  
  prob_belowsizelim_block1.initialize();
  prob_belowsizelim_block2.initialize();
  prob_belowsizelim_block3.initialize();
  prob_belowsizelim_block4.initialize();
  prob_belowsizelim_block5.initialize();
  
  log_rec_dev_output.initialize();
  log_Nage_dev_output.initialize();
  log_rec_dev.initialize();
  log_Nage_dev.initialize();
  
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
TOP_OF_MAIN_SECTION
  time(&start);
  arrmblsize=20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(1600);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);
  

//>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PROCEDURE_SECTION

 R0=mfexp(log_R0);
 
 //cout<<"start"<<endl;

 //get_M_at_age(); //Needed only if M is estimated
 
 get_length_weight_at_age(); 
 //cout << "got length, weight, fecundity transitions" <<endl;
 get_reprod();
// cout << "got reprod" << endl;
 get_length_at_age_dist(); 
// cout<< "got predicted length at age distribution"<<endl;
 get_weight_at_age_landings();
// cout<< "got weight at age of landings"<<endl; 
 get_spr_F0();
// cout << "got F0 spr" << endl;
 get_selectivity(); 
// cout << "got selectivity" << endl;
 get_mortality(); 
// cout << "got mortalities" << endl;
 get_bias_corr(); 
// cout<< "got recruitment bias correction" << endl;
 get_numbers_at_age(); 
// cout << "got numbers at age" << endl;
 get_landings_numbers();
// cout << "got catch at age" << endl;
 get_landings_wgt();
// cout << "got landings" << endl;
 get_dead_discards();
 //cout << "got discards" << endl;
 get_catchability_fcns();
 //cout << "got catchability_fcns" << endl;
 get_indices();
 //cout << "got indices" << endl;
 get_length_comps();
 //cout<< "got length comps"<< endl;
 get_age_comps();
// cout<< "got age comps"<< endl;
 evaluate_objective_function();
 //cout << "objective function calculations complete" << endl;


//FUNCTION get_M_at_age
//    Mscale_len=Linf*(1.0-mfexp(-K*(Mscale_ages-t0+0.5))); 
//    Mscale_wgt_g=wgtpar_a*pow(Mscale_len,wgtpar_b); 
//    M_lorenzen=3.69*pow(Mscale_wgt_g,-0.305);
//    cum_surv_1plus=mfexp(-max_obs_age*M_constant);  
//    M=M_lorenzen(1,nages)*(-log(cum_surv_1plus)/sum(M_lorenzen(1,max_obs_age)));    


FUNCTION get_length_weight_at_age
  //compute mean length (mm) and weight (whole) at age
    meanlen_TL=Linf*(1.0-mfexp(-K*(agebins-t0+0.5)));    //total length in mm
    wgt_g=wgtpar_a*pow(meanlen_TL,wgtpar_b);             //wgt in grams 
    wgt_kg=g2kg*wgt_g;                                   //wgt in kilograms 
    wgt_mt=g2mt*wgt_g;                                   //mt of whole wgt: g2mt converts g to mt
    wgt_klb=mt2klb*wgt_mt;                               //1000 lb of whole wgt
    wgt_lb=mt2lb*wgt_mt;                                 //1000 lb of whole wgt
    fecundity=fecpar_batches*mfexp(fecpar_a+wgt_g*fecpar_b)/fecpar_scale;    //fecundity at age, scaled

FUNCTION get_reprod 
   //reprod is product of stuff going into reproductive capacity calcs
   reprod=elem_prod(elem_prod(prop_f,maturity_f),fecundity);  
   reprod2=elem_prod(elem_prod(prop_f,maturity_f),wgt_mt);  
 
FUNCTION get_length_at_age_dist
  
//compute matrix of length at age, based on the normal distribution
  //len_cv=len_cv_val; 
  //len_sd=elem_prod(len_cv, meanlen_TL);
  for (iage=1;iage<=nages;iage++)
   {len_cv(iage)=len_cv_val;
    len_sd(iage)=meanlen_TL(iage)*len_cv(iage);
    zscore_lzero=(0.0-meanlen_TL(iage))/len_sd(iage); 
	cprob_lzero=cumd_norm(zscore_lzero);	 
    //lenprob_all(iage)/=sum(lenprob_all(iage)); 
	//cout<<"len_cv"<<len_cv<<endl;
//first length bin
	//population
    zscore_len=((lenbins(1)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage);
    cprob_lenvec(1)=cumd_norm(zscore_len);          //includes any probability mass below zero
    lenprob(iage,1)=cprob_lenvec(1)-cprob_lzero;    //removes any probability mass below zero	
    //cout<<"zscore_len"<<zscore_len<<endl;
	//First size limit 8" Period 1 for both commercial and recreational discards (1984 through 1998)
	zscore_lsizelim1=(sizelim1-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim1=cumd_norm(zscore_lsizelim1);   //includes any probability mass below zero
	prob_belowsizelim_block1(iage)=	cprob_lsizelim1-cprob_lzero; //removes any probability mass below zero
	//cout<<"zscore_lsizelim1"<<zscore_lsizelim1<<endl;
	//cout<<"prob_belowsizelim_block1"<<prob_belowsizelim_block1<<endl;
	//Second size limit 10" Period 2 for both commercial and recreational discards and Period 3 for recreational landings (1999 through 2012 for commercial) (1999 through 2006 for recreational)
	zscore_lsizelim2=(sizelim2-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim2=cumd_norm(zscore_lsizelim2);   //includes any probability mass below zero
	prob_belowsizelim_block2(iage)=	cprob_lsizelim2-cprob_lzero; //removes any probability mass below zero
	//cout<<"zscore_lsizelim2"<<zscore_lsizelim2<<endl;
	//cout<<"prob_belowsizelim_block2"<<prob_belowsizelim_block2<<endl;
	//Third size limit 11" Period 3 for commercial (2013 through terminal year)
	zscore_lsizelim3=(sizelim3-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim3=cumd_norm(zscore_lsizelim3);                 //includes any probability mass below zero
	prob_belowsizelim_block3(iage)=	cprob_lsizelim3-cprob_lzero; //removes any probability mass below zero
	//cout<<"zscore_lsizelim3"<<zscore_lsizelim3<<endl;
	//cout<<"prob_belowsizelim_block3"<<prob_belowsizelim_block3<<endl;
	//Fourth size limit 12" Period 3 for recreational discards and Period 4 for recreational landings (2007 through 2012)
	zscore_lsizelim4=(sizelim4-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim4=cumd_norm(zscore_lsizelim4);                 //includes any probability mass below zero
	prob_belowsizelim_block4(iage)=	cprob_lsizelim4-cprob_lzero; //removes any probability mass below zero
	//cout<<"zscore_lsizelim4"<<zscore_lsizelim4<<endl;
	//cout<<"prob_belowsizelim_block4"<<prob_belowsizelim_block4<<endl;
	//Fifth size limit 13" Period 4 for recreational discards and Period 5 for recreational landings (2013 through the terminal year)
	zscore_lsizelim5=(sizelim5-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim5=cumd_norm(zscore_lsizelim5);                 //includes any probability mass below zero
	prob_belowsizelim_block5(iage)=	cprob_lsizelim5-cprob_lzero; //removes any probability mass below zero
	//cout<<"zscore_lsizelim5"<<zscore_lsizelim5<<endl;
	//cout<<"prob_belowsizelim_block5"<<prob_belowsizelim_block5<<endl;

	//most other length bins  
    //population
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len=((lenbins(ilen)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage); 
		cprob_lenvec(ilen)=cumd_norm(zscore_len);
        lenprob(iage,ilen)=cprob_lenvec(ilen)-cprob_lenvec(ilen-1);
      }
	//last length bin is a plus group
	//population
    zscore_len=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage); 
	lenprob(iage,nlenbins)=1.0-cumd_norm(zscore_len);
      lenprob(iage)=lenprob(iage)/(1.0-cprob_lzero);  //renormalize to account for any prob mass below size=0
 }
   //cout<<"lenprob"<<lenprob<<endl;	  
  // for (ilen=1;ilen<=nlenbins_all;ilen++)
  //  { lenprob_all(iage,ilen)=(mfexp(-(square(lenbins_all(ilen)-meanlen_TL(iage))/
  //    (2.*square(len_sd(iage)))))/(sqrt2pi*len_sd(iage)));
  //  } 
	 
//andardize to approximate integration and to account for truncated normal (i.e., no sizes<smallest)
//for (ilen=1;ilen<=nlenbins;ilen++) {lenprob(iage,ilen)=lenprob_all(iage,ilen);}
//for (ilen=nlenbins+1;ilen<=nlenbins_all;ilen++){lenprob(iage)(nlenbins)=lenprob(iage)(nlenbins)+lenprob_all(iage)(ilen);} //plus group
 
  
  

  //fishery/fleet specific length probs, assumed normal prior to size limits
  lenprob_Mbft=lenprob;
  
  lenprob_cL1=lenprob; 
  lenprob_cL2=lenprob;//_all;    //values may be adjusted based on size limit
  lenprob_cL3=lenprob;//_all;    //values may be adjusted based on size limit 
                         ;//
  lenprob_cP1=lenprob;   ;//
  lenprob_cP2=lenprob;//_all;    //values may be adjusted based on size limit
  lenprob_cP3=lenprob;//_all;    //values may be adjusted based on size limit
                         ;//
  lenprob_cT1=lenprob;   ;//
  lenprob_cT2=lenprob;//_all;    //values may be adjusted based on size limit
                         ;//
  lenprob_HB1=lenprob;   ;//
  lenprob_HB2=lenprob;//_all;    //values may be adjusted based on size limit 
  lenprob_HB3=lenprob;//_all;    //values may be adjusted based on size limit 
  lenprob_HB4=lenprob;//_all;    //values may be adjusted based on size limit
  lenprob_HB5=lenprob;//_all;    //values may be adjusted based on size limit  
    
  lenprob_mrip1=lenprob; 
  lenprob_mrip2=lenprob;//_all;   //values may be adjusted based on size limit 
  lenprob_mrip3=lenprob;//_all;   //values may be adjusted based on size limit 
  lenprob_mrip4=lenprob;//_all;   //values may be adjusted based on size limit 
  lenprob_mrip5=lenprob;//_all;   //values may be adjusted based on size limit 
  
  lenprob_comm_D2=lenprob;//_all; //values may be adjusted based on size limit 
  lenprob_comm_D3=lenprob;//_all; //values may be adjusted based on size limit 
  lenprob_comm_D4=lenprob;//_all; //values may be adjusted based on size limit 
  lenprob_HB_D2=lenprob;//_all;   //values may be adjusted based on size limit 
  lenprob_HB_D3=lenprob;//_all;   //values may be adjusted based on size limit 
  lenprob_HB_D4=lenprob;//_all;   //values may be adjusted based on size limit
  lenprob_HB_D5=lenprob;//_all;   //values may be adjusted based on size limit   
  lenprob_mrip_D1=lenprob;//_all; //values may be adjusted based on size limit 
  lenprob_mrip_D2=lenprob;//_all; //values may be adjusted based on size limit 
  lenprob_mrip_D3=lenprob;//_all; //values may be adjusted based on size limit 
  lenprob_mrip_D4=lenprob;//_all; //values may be adjusted based on size limit
  lenprob_mrip_D5=lenprob;//_all; //values may be adjusted based on size limit  
  
  //for (iage=1;iage<=nages;iage++)
  //{
  //  for (ilen=1;ilen<=nlenbins_all;ilen++)
  //  {
  //      if (lenbins_all(ilen) < limit_8in)  //Landings block two
  //      {
  //        lenprob_cL2_all(iage,ilen)=p_lenc_cL2*lenprob_all(iage,ilen);
  //        lenprob_cP2_all(iage,ilen)=p_lenc_cP2*lenprob_all(iage,ilen); 
  //        lenprob_cT2_all(iage,ilen)=p_lenc_cT2*lenprob_all(iage,ilen);          
  //        lenprob_HB2_all(iage,ilen)=p_lenc_HB2*lenprob_all(iage,ilen);
  //        lenprob_mrip2_all(iage,ilen)=p_lenc_mrip2*lenprob_all(iage,ilen);
  //       
  //      }
  //      if (lenbins_all(ilen) < limit_10in) //Landings block three
  //      {
  //        lenprob_cL3_all(iage,ilen)=p_lenc_cL3*lenprob_all(iage,ilen);
  //        lenprob_cP3_all(iage,ilen)=p_lenc_cP3*lenprob_all(iage,ilen);         
  //        lenprob_HB3_all(iage,ilen)=p_lenc_HB3*lenprob_all(iage,ilen);
  //        lenprob_mrip3_all(iage,ilen)=p_lenc_mrip3*lenprob_all(iage,ilen);
  //       }
  //      if (lenbins_all(ilen) < limit_12in) //Landings block four
  //      {
  //        lenprob_HB4_all(iage,ilen)=p_lenc_HB4*lenprob_all(iage,ilen);
  //        lenprob_mrip4_all(iage,ilen)=p_lenc_mrip4*lenprob_all(iage,ilen);
  //       }
	//	if (lenbins_all(ilen) < limit_13in) //Landings block five
  //      {
  //        lenprob_HB5_all(iage,ilen)=p_lenc_HB5*lenprob_all(iage,ilen);
  //        lenprob_mrip5_all(iage,ilen)=p_lenc_mrip5*lenprob_all(iage,ilen);
  //       }
  //
  //
  //      if (lenbins_all(ilen) > limit_disc) //Discards block 1
  //      { lenprob_mrip_D1_all(iage,ilen)=p_lenc_mrip_D1*lenprob_all(iage,ilen);}
  //
  //      if (lenbins_all(ilen) > limit_8in) //Discards block two
  //      { 
  //        lenprob_comm_D2_all(iage,ilen)=p_lenc_comm_D2*lenprob_all(iage,ilen);
  //        lenprob_HB_D2_all(iage,ilen)=p_lenc_HB_D2*lenprob_all(iage,ilen); 
  //        lenprob_mrip_D2_all(iage,ilen)=p_lenc_mrip_D2*lenprob_all(iage,ilen);
  //      }
  //
  //      if (lenbins_all(ilen) > limit_10in) //Discards block three
  //      { 
  //        lenprob_comm_D3_all(iage,ilen)=p_lenc_comm_D3*lenprob_all(iage,ilen);
  //        lenprob_HB_D3_all(iage,ilen)=p_lenc_HB_D3*lenprob_all(iage,ilen);
  //        lenprob_mrip_D3_all(iage,ilen)=p_lenc_mrip_D3*lenprob_all(iage,ilen);          
  //      }
  //
  //      if (lenbins_all(ilen) > limit_12in) //Discards block four
  //      { 
  //        lenprob_HB_D4_all(iage,ilen)=p_lenc_HB_D4*lenprob_all(iage,ilen);
  //        lenprob_mrip_D4_all(iage,ilen)=p_lenc_mrip_D4*lenprob_all(iage,ilen);          
  //      }
  //
  //
  //  }  //end ilen loop
  //  
  //  
  //  if (iage>=4)  //compute prior to standardizing
  //      {vecprob_HB_D2(iage)=sum(lenprob_HB_D2_all(iage));
  //       vecprob_HB_D3(iage)=sum(lenprob_HB_D3_all(iage));
  //       vecprob_HB_D4(iage)=sum(lenprob_HB_D4_all(iage));
	//	 vecprob_HB_D5(iage)=sum(lenprob_HB_D5_all(iage));
  //      }
//   
//   lenprob_cL2_all(iage)/=sum(lenprob_cL2_all(iage));       //standardize
//   lenprob_cL3_all(iage)/=sum(lenprob_cL3_all(iage));       //standardize
//   lenprob_cP2_all(iage)/=sum(lenprob_cP2_all(iage));       //standardize
//   lenprob_cP3_all(iage)/=sum(lenprob_cP3_all(iage));       //standardize
//   lenprob_cT2_all(iage)/=sum(lenprob_cT2_all(iage));       //standardize
//   lenprob_HB2_all(iage)/=sum(lenprob_HB2_all(iage));       //standardize
//   lenprob_HB3_all(iage)/=sum(lenprob_HB3_all(iage));       //standardize
//   lenprob_HB4_all(iage)/=sum(lenprob_HB4_all(iage));       //standardize
//	lenprob_HB5_all(iage)/=sum(lenprob_HB5_all(iage));       //standardize
//   lenprob_mrip2_all(iage)/=sum(lenprob_mrip2_all(iage));   //standardize
//   lenprob_mrip3_all(iage)/=sum(lenprob_mrip3_all(iage));   //standardize
//   lenprob_mrip4_all(iage)/=sum(lenprob_mrip4_all(iage));   //standardize
//	lenprob_mrip5_all(iage)/=sum(lenprob_mrip5_all(iage));   //standardize
//
//   lenprob_comm_D2_all(iage)/=sum(lenprob_comm_D2_all(iage));   //standardize
//   lenprob_comm_D3_all(iage)/=sum(lenprob_comm_D3_all(iage));   //standardize
//   lenprob_HB_D2_all(iage)/=sum(lenprob_HB_D2_all(iage));       //standardize
//   lenprob_HB_D3_all(iage)/=sum(lenprob_HB_D3_all(iage));       //standardize
//   lenprob_HB_D4_all(iage)/=sum(lenprob_HB_D4_all(iage));       //standardize
//   lenprob_mrip_D2_all(iage)/=sum(lenprob_mrip_D2_all(iage));   //standardize
//   lenprob_mrip_D3_all(iage)/=sum(lenprob_mrip_D3_all(iage));   //standardize
//   lenprob_mrip_D4_all(iage)/=sum(lenprob_mrip_D4_all(iage));   //standardize
//
//   for (ilen=1;ilen<=nlenbins;ilen++) 
//     {lenprob_cL2(iage,ilen)=lenprob_cL2_all(iage,ilen);
//      lenprob_cL3(iage,ilen)=lenprob_cL3_all(iage,ilen);
//      lenprob_cP2(iage,ilen)=lenprob_cP2_all(iage,ilen);
//      lenprob_cP3(iage,ilen)=lenprob_cP3_all(iage,ilen);
//      lenprob_cT2(iage,ilen)=lenprob_cT2_all(iage,ilen);
//      lenprob_HB2(iage,ilen)=lenprob_HB2_all(iage,ilen);
//      lenprob_HB3(iage,ilen)=lenprob_HB3_all(iage,ilen);
//      lenprob_HB4(iage,ilen)=lenprob_HB4_all(iage,ilen); 
//	   lenprob_HB5(iage,ilen)=lenprob_HB5_all(iage,ilen);         
//      lenprob_mrip2(iage,ilen)=lenprob_mrip2_all(iage,ilen);
//      lenprob_mrip3(iage,ilen)=lenprob_mrip3_all(iage,ilen);
//      lenprob_mrip4(iage,ilen)=lenprob_mrip4_all(iage,ilen);
//	   lenprob_mrip5(iage,ilen)=lenprob_mrip5_all(iage,ilen);
//      lenprob_comm_D2(iage,ilen)=lenprob_comm_D2_all(iage,ilen);
//      lenprob_comm_D3(iage,ilen)=lenprob_comm_D3_all(iage,ilen);
//      lenprob_HB_D2(iage,ilen)=lenprob_HB_D2_all(iage,ilen);
//      lenprob_HB_D3(iage,ilen)=lenprob_HB_D3_all(iage,ilen);
//      lenprob_HB_D4(iage,ilen)=lenprob_HB_D4_all(iage,ilen);       
//      lenprob_mrip_D2(iage,ilen)=lenprob_mrip_D2_all(iage,ilen);
//      lenprob_mrip_D3(iage,ilen)=lenprob_mrip_D3_all(iage,ilen); 
//      lenprob_mrip_D4(iage,ilen)=lenprob_mrip_D4_all(iage,ilen);                           
//     }
//  
//   for (ilen=nlenbins+1;ilen<=nlenbins_all;ilen++) //plus group
//     {
//     lenprob_cL2(iage)(nlenbins)=lenprob_cL2(iage)(nlenbins)+lenprob_cL2_all(iage)(ilen);
//     lenprob_cL3(iage)(nlenbins)=lenprob_cL3(iage)(nlenbins)+lenprob_cL3_all(iage)(ilen);
//     lenprob_cP2(iage)(nlenbins)=lenprob_cP2(iage)(nlenbins)+lenprob_cP2_all(iage)(ilen);
//     lenprob_cP3(iage)(nlenbins)=lenprob_cP3(iage)(nlenbins)+lenprob_cP3_all(iage)(ilen);
//     lenprob_cT2(iage)(nlenbins)=lenprob_cT2(iage)(nlenbins)+lenprob_cT2_all(iage)(ilen);
//     lenprob_HB2(iage)(nlenbins)=lenprob_HB2(iage)(nlenbins)+lenprob_HB2_all(iage)(ilen);
//     lenprob_HB3(iage)(nlenbins)=lenprob_HB3(iage)(nlenbins)+lenprob_HB3_all(iage)(ilen);
//     lenprob_HB4(iage)(nlenbins)=lenprob_HB4(iage)(nlenbins)+lenprob_HB4_all(iage)(ilen);
//	  lenprob_HB5(iage)(nlenbins)=lenprob_HB5(iage)(nlenbins)+lenprob_HB5_all(iage)(ilen);
//     lenprob_mrip2(iage)(nlenbins)=lenprob_mrip2(iage)(nlenbins)+lenprob_mrip2_all(iage)(ilen);
//     lenprob_mrip3(iage)(nlenbins)=lenprob_mrip3(iage)(nlenbins)+lenprob_mrip3_all(iage)(ilen);
//     lenprob_mrip4(iage)(nlenbins)=lenprob_mrip4(iage)(nlenbins)+lenprob_mrip4_all(iage)(ilen);
//	  lenprob_mrip5(iage)(nlenbins)=lenprob_mrip5(iage)(nlenbins)+lenprob_mrip5_all(iage)(ilen);
//     lenprob_comm_D2(iage)(nlenbins)=lenprob_comm_D2(iage)(nlenbins)+lenprob_comm_D2_all(iage)(ilen);
//     lenprob_comm_D3(iage)(nlenbins)=lenprob_comm_D3(iage)(nlenbins)+lenprob_comm_D3_all(iage)(ilen);
//     lenprob_HB_D2(iage)(nlenbins)=lenprob_HB_D2(iage)(nlenbins)+lenprob_HB_D2_all(iage)(ilen);
//     lenprob_HB_D3(iage)(nlenbins)=lenprob_HB_D3(iage)(nlenbins)+lenprob_HB_D3_all(iage)(ilen);
//     lenprob_HB_D4(iage)(nlenbins)=lenprob_HB_D4(iage)(nlenbins)+lenprob_HB_D4_all(iage)(ilen);
//     lenprob_mrip_D2(iage)(nlenbins)=lenprob_mrip_D2(iage)(nlenbins)+lenprob_mrip_D2_all(iage)(ilen);
//     lenprob_mrip_D3(iage)(nlenbins)=lenprob_mrip_D3(iage)(nlenbins)+lenprob_mrip_D3_all(iage)(ilen);
//     lenprob_mrip_D4(iage)(nlenbins)=lenprob_mrip_D4(iage)(nlenbins)+lenprob_mrip_D4_all(iage)(ilen);
//     } 
//       
// } //end iage loop
 
FUNCTION get_weight_at_age_landings
  //fleets under identical size limits are set equal at end of fcn
  for (iyear=styr; iyear<=endyr_period1; iyear++)
  {
    len_cL_mm(iyear)=meanlen_TL;
    wgt_cL_klb(iyear)=wgt_klb;      
    //len_cP_mm(iyear)=meanlen_TL;  
    //wgt_cP_klb(iyear)=wgt_klb;  
    //len_cT_mm(iyear)=meanlen_TL;  
    //wgt_cT_klb(iyear)=wgt_klb;  
    
    //len_HB_mm(iyear)=meanlen_TL;
    //wgt_HB_klb(iyear)=wgt_klb;  
    len_mrip_mm(iyear)=meanlen_TL;
    wgt_mrip_klb(iyear)=wgt_klb;  

    for (iage=1;iage<=nages; iage++)
    {
    len_mrip_D_mm(iyear,iage)=sum(elem_prod(lenprob_mrip_D2(iage),lenbins)); //assumes same size distn in period 1 as in period 2       
    }
    wgt_mrip_D_klb(iyear)=g2klb*wgtpar_a*pow(len_mrip_D_mm(iyear),wgtpar_b);
  } // end iyear loop
  
  
  for (iyear=(endyr_period1+1); iyear<=endyr_period2; iyear++)
  {
    for (iage=1;iage<=nages; iage++)
    {
    len_cL_mm(iyear,iage)=sum(elem_prod(lenprob_cL2(iage),lenbins)); 
    //len_cP_mm(iyear,iage)=sum(elem_prod(lenprob_cP2_all(iage),lenbins_all)); 
    //len_cT_mm(iyear,iage)=sum(elem_prod(lenprob_cT2_all(iage),lenbins_all));     
    //len_HB_mm(iyear,iage)=sum(elem_prod(lenprob_HB2_all(iage),lenbins_all));     
    len_mrip_mm(iyear,iage)=sum(elem_prod(lenprob_mrip2(iage),lenbins)); 
    len_comm_D_mm(iyear,iage)=sum(elem_prod(lenprob_comm_D2(iage),lenbins));
    //len_HB_D_mm(iyear,iage)=sum(elem_prod(lenprob_HB_D2_all(iage),lenbins_all));        
    len_mrip_D_mm(iyear,iage)=sum(elem_prod(lenprob_mrip_D2(iage),lenbins));             
    }
    wgt_cL_klb(iyear)=g2klb*wgtpar_a*pow(len_cL_mm(iyear),wgtpar_b);
    //wgt_cP_klb(iyear)=g2klb*wgtpar_a*pow(len_cP_mm(iyear),wgtpar_b);
    //wgt_cT_klb(iyear)=g2klb*wgtpar_a*pow(len_cT_mm(iyear),wgtpar_b);    
    //wgt_HB_klb(iyear)=g2klb*wgtpar_a*pow(len_HB_mm(iyear),wgtpar_b);
    wgt_mrip_klb(iyear)=g2klb*wgtpar_a*pow(len_mrip_mm(iyear),wgtpar_b);
    wgt_comm_D_klb(iyear)=g2klb*wgtpar_a*pow(len_comm_D_mm(iyear),wgtpar_b);  
    //wgt_HB_D_klb(iyear)=g2klb*wgtpar_a*pow(len_HB_D_mm(iyear),wgtpar_b);
    wgt_mrip_D_klb(iyear)=g2klb*wgtpar_a*pow(len_mrip_D_mm(iyear),wgtpar_b);               
  }

  for (iyear=(endyr_period2+1); iyear<=endyr; iyear++) //comm only
  {
    for (iage=1;iage<=nages; iage++)
    {
    len_cL_mm(iyear,iage)=sum(elem_prod(lenprob_cL3(iage),lenbins)); 
    //len_cP_mm(iyear,iage)=sum(elem_prod(lenprob_cP3_all(iage),lenbins_all));      
    len_comm_D_mm(iyear,iage)=sum(elem_prod(lenprob_comm_D3(iage),lenbins));
    }
    wgt_cL_klb(iyear)=g2klb*wgtpar_a*pow(len_cL_mm(iyear),wgtpar_b);
    //wgt_cP_klb(iyear)=g2klb*wgtpar_a*pow(len_cP_mm(iyear),wgtpar_b);   
    wgt_comm_D_klb(iyear)=g2klb*wgtpar_a*pow(len_comm_D_mm(iyear),wgtpar_b);  
  }    
  
  for (iyear=(endyr_period2+1); iyear<=endyr_recr_period3; iyear++) //rec only
  {
    for (iage=1;iage<=nages; iage++)
    {
    //len_HB_mm(iyear,iage)=sum(elem_prod(lenprob_HB3_all(iage),lenbins_all));     
    len_mrip_mm(iyear,iage)=sum(elem_prod(lenprob_mrip3(iage),lenbins)); 
    //len_HB_D_mm(iyear,iage)=sum(elem_prod(lenprob_HB_D3_all(iage),lenbins_all));        
    len_mrip_D_mm(iyear,iage)=sum(elem_prod(lenprob_mrip_D3(iage),lenbins));             
    }
    //wgt_HB_klb(iyear)=g2klb*wgtpar_a*pow(len_HB_mm(iyear),wgtpar_b);
    wgt_mrip_klb(iyear)=g2klb*wgtpar_a*pow(len_mrip_mm(iyear),wgtpar_b);
    //wgt_HB_D_klb(iyear)=g2klb*wgtpar_a*pow(len_HB_D_mm(iyear),wgtpar_b);
    wgt_mrip_D_klb(iyear)=g2klb*wgtpar_a*pow(len_mrip_D_mm(iyear),wgtpar_b);               
  }
  
  for (iyear=(endyr_recr_period3+1); iyear<=endyr; iyear++) //rec only
  {
    for (iage=1;iage<=nages; iage++)
    {
    //len_HB_mm(iyear,iage)=sum(elem_prod(lenprob_HB4_all(iage),lenbins_all));     
    len_mrip_mm(iyear,iage)=sum(elem_prod(lenprob_mrip4(iage),lenbins)); 
    //len_HB_D_mm(iyear,iage)=sum(elem_prod(lenprob_HB_D4_all(iage),lenbins_all));        
    len_mrip_D_mm(iyear,iage)=sum(elem_prod(lenprob_mrip_D4(iage),lenbins));             
    }
    //wgt_HB_klb(iyear)=g2klb*wgtpar_a*pow(len_HB_mm(iyear),wgtpar_b);
    wgt_mrip_klb(iyear)=g2klb*wgtpar_a*pow(len_mrip_mm(iyear),wgtpar_b);
    //wgt_HB_D_klb(iyear)=g2klb*wgtpar_a*pow(len_HB_D_mm(iyear),wgtpar_b);
    wgt_mrip_D_klb(iyear)=g2klb*wgtpar_a*pow(len_mrip_D_mm(iyear),wgtpar_b);               
  }  
   
   for (iyear=(endyr_period4+1); iyear<=endyr; iyear++) //rec only
  {
    for (iage=1;iage<=nages; iage++)
    {
    //len_HB_mm(iyear,iage)=sum(elem_prod(lenprob_HB4_all(iage),lenbins_all));     
    len_mrip_mm(iyear,iage)=sum(elem_prod(lenprob_mrip5(iage),lenbins)); 
    //len_HB_D_mm(iyear,iage)=sum(elem_prod(lenprob_HB_D4_all(iage),lenbins_all));        
    len_mrip_D_mm(iyear,iage)=sum(elem_prod(lenprob_mrip_D4(iage),lenbins));             
    }
    //wgt_HB_klb(iyear)=g2klb*wgtpar_a*pow(len_HB_mm(iyear),wgtpar_b);
    wgt_mrip_klb(iyear)=g2klb*wgtpar_a*pow(len_mrip_mm(iyear),wgtpar_b);
    //wgt_HB_D_klb(iyear)=g2klb*wgtpar_a*pow(len_HB_D_mm(iyear),wgtpar_b);
    wgt_mrip_D_klb(iyear)=g2klb*wgtpar_a*pow(len_mrip_D_mm(iyear),wgtpar_b);               
  }  
  
    //identical fleets set equal here (for speed)      
    len_cP_mm=len_cL_mm;  wgt_cP_klb=wgt_cL_klb;  
    len_cT_mm=len_cL_mm;  wgt_cT_klb=wgt_cL_klb;  
    len_HB_mm=len_mrip_mm;  wgt_HB_klb=wgt_mrip_klb; 
    len_HB_D_mm=len_mrip_D_mm;  wgt_HB_D_klb=wgt_mrip_D_klb; 

  for (iyear=styr_comm_closed; iyear<=endyr; iyear++) //overwrite last two yrs comm discards, accnt for quotas
  { len_comm_D_mm(iyear)=Dprop_comm_sel_D*len_comm_D_mm(iyear)+Dprop_comm_sel_cL*len_cL_mm(iyear)+
                         Dprop_comm_sel_cP*len_cP_mm(iyear);
    wgt_comm_D_klb(iyear)=g2klb*wgtpar_a*pow(len_comm_D_mm(iyear),wgtpar_b);                      
  } 
    
FUNCTION get_spr_F0
  //at mdyr, apply half this yr's mortality, half next yr's
  N_spr_F0(1)=1.0*mfexp(-1.0*M(1)*spawn_time_frac); //at peak spawning time
  N_bpr_F0(1)=1.0;      //at start of year
  for (iage=2; iage<=nages; iage++)
  {
    //N_spr_F0(iage)=N_spr_F0(iage-1)*mfexp(-1.0*(M(iage-1));
    N_spr_F0(iage)=N_spr_F0(iage-1)*
                   mfexp(-1.0*(M(iage-1)*(1.0-spawn_time_frac) + M(iage)*spawn_time_frac)); 
    N_bpr_F0(iage)=N_bpr_F0(iage-1)*mfexp(-1.0*(M(iage-1)));    
  }
  N_spr_F0(nages)=N_spr_F0(nages)/(1.0-mfexp(-1.0*M(nages))); //plus group (sum of geometric series)
  N_bpr_F0(nages)=N_bpr_F0(nages)/(1.0-mfexp(-1.0*M(nages)));
   
  spr_F0=sum(elem_prod(N_spr_F0,reprod));
  bpr_F0=sum(elem_prod(N_bpr_F0,wgt_mt));    
 

FUNCTION get_selectivity

// ------- compute landings selectivities by period

  //---flat-topped sels---------------------------
  sel_Mbft_vec=logistic(agebins, selpar_A50_Mbft, selpar_slope_Mbft);
  sel_Mcvt_vec=logistic(agebins, selpar_A50_Mcvt, selpar_slope_Mcvt);

  sel_cL_2=logistic(agebins, selpar_A50_cL2, selpar_slope_cL2);
  sel_cL_3=logistic(agebins, selpar_A50_cL3, selpar_slope_cL3);
  sel_cL_4=logistic(agebins, selpar_A50_cL4, selpar_slope_cL4);

  sel_cP_2=logistic(agebins, selpar_A50_cP2, selpar_slope_cP2);
  sel_cP_3=logistic(agebins, selpar_A50_cP3, selpar_slope_cP3);
  sel_cP_4=logistic(agebins, selpar_A50_cP4, selpar_slope_cP4);
  
  sel_HB_1=logistic(agebins, selpar_A50_HB1, selpar_slope_HB1);
  sel_HB_2=logistic(agebins, selpar_A50_HB2, selpar_slope_HB2);
  sel_HB_3=logistic(agebins, selpar_A50_HB3, selpar_slope_HB3);
  sel_HB_4=logistic(agebins, selpar_A50_HB4, selpar_slope_HB4);
  sel_HB_5=logistic(agebins, selpar_A50_HB5, selpar_slope_HB5);

 // selpar_A50_mrip1=selpar_A50_HB1;
 // selpar_slope_mrip1=selpar_slope_HB1;
 // selpar_A50_mrip2=selpar_A50_HB2;
 // selpar_slope_mrip2=selpar_slope_HB2;
 // selpar_A50_mrip3=selpar_A50_HB3;
 // selpar_slope_mrip3=selpar_slope_HB3;    
 // selpar_A50_mrip4=selpar_A50_HB4;
//  selpar_slope_mrip4=selpar_slope_HB4;
  
  sel_mrip1=logistic(agebins, selpar_A50_mrip1, selpar_slope_mrip1);
  sel_mrip2=logistic(agebins, selpar_A50_mrip2, selpar_slope_mrip2);
  sel_mrip3=logistic(agebins, selpar_A50_mrip3, selpar_slope_mrip3);
  sel_mrip4=logistic(agebins, selpar_A50_mrip4, selpar_slope_mrip4);
  sel_mrip5=logistic(agebins, selpar_A50_mrip5, selpar_slope_mrip5);


//-----------fill in years--------------------------------------------
  
  //Period 1:   
  for (iyear=styr; iyear<=endyr_period1; iyear++)
  {
     sel_Mbft(iyear)=sel_Mbft_vec;
     sel_Mcvt(iyear)=sel_Mcvt_vec;
     sel_cL(iyear)=sel_cL_2; //commercial handline sel mirrors period 2
     sel_cP(iyear)=sel_cP_2; //commercial handline sel mirrors period 2
     sel_HB(iyear)=sel_HB_1; 
     sel_mrip(iyear)=sel_mrip1;      
  }

  //Period 2: 
  for (iyear=endyr_period1+1; iyear<=endyr_period2; iyear++)
  {     
     sel_Mbft(iyear)=sel_Mbft_vec;
     sel_Mcvt(iyear)=sel_Mcvt_vec;
     sel_cL(iyear)=sel_cL_2;
     sel_cP(iyear)=sel_cP_2; 
     sel_HB(iyear)=sel_HB_2;
     sel_mrip(iyear)=sel_mrip2;     
  }

  //Period 3 
  for (iyear=endyr_period2+1; iyear<=endyr; iyear++)
  {
     sel_Mbft(iyear)=sel_Mbft_vec;
     sel_Mcvt(iyear)=sel_Mcvt_vec;
     sel_cL(iyear)=sel_cL_3; 
     sel_cP(iyear)=sel_cP_3;     
     sel_HB(iyear)=sel_HB_3;
     sel_mrip(iyear)=sel_mrip3;     
  }   
  //Period 4: rec only, overwrites yrs calculated for period 3
  for (iyear=endyr_recr_period3; iyear<endyr_period4; iyear++)
  {
     sel_HB(iyear)=sel_HB_4;
     sel_mrip(iyear)=sel_mrip4;     
  }   
   //Period 5: Comm and rec, overwrites yrs previously calculated
   for (iyear=endyr_period4; iyear<=endyr; iyear++)
  {
     sel_cL(iyear)=sel_cL_4; 
     sel_cP(iyear)=sel_cP_4;
	 sel_HB(iyear)=sel_HB_5;
     sel_mrip(iyear)=sel_mrip5;     
  }   
  //set selectivities that mirror others
  sel_cT=sel_cP;  
  //sel_mrip=sel_HB;
  
  
//---Discard selectivities---------------------------------    
//---------------------------------------------------------   
  //mrip and comm mirror headboat discard selectivity

  selpar_Age0_HB_D=1.0/(1.0+mfexp(-selpar_Age0_HB_D_logit));
  selpar_Age1_HB_D=1.0/(1.0+mfexp(-selpar_Age1_HB_D_logit)); 
  selpar_Age2_HB_D=1.0/(1.0+mfexp(-selpar_Age2_HB_D_logit));
  ///cout<<"selpar_Age0_HB_D"<<selpar_Age0_HB_D<<endl;
  ///cout<<"selpar_Age1_HB_D"<<selpar_Age1_HB_D<<endl;
  ///cout<<"selpar_Age2_HB_D"<<selpar_Age2_HB_D<<endl;

  //cout << "selpar_slope2_HBD5" << selpar_slope2_HBD5 << endl;
  sel_HB_D_4=logistic_exponential(agebins, selpar_A50_HBD4, selpar_slope_HBD4, selpar_A502_HBD4, selpar_slope2_HBD4);//logistic_double(agebins, selpar_A50_HBD4, selpar_slope_HBD4, selpar_A502_HBD4, selpar_slope2_HBD4); 
  sel_HB_D_5=logistic_exponential(agebins, selpar_A50_HBD5, selpar_slope_HBD5, selpar_A502_HBD5, selpar_slope2_HBD5);//logistic_double(agebins, selpar_A50_HBD5, selpar_slope_HBD5, selpar_A502_HBD5, selpar_slope2_HBD5); 
   
 //Assume same sel of age 0's across periods 
  sel_HB_D_2(1)=selpar_Age0_HB_D; 
  sel_HB_D_3(1)=selpar_Age0_HB_D; 
  //sel_HB_D_4(1)=selpar_Age0_HB_D;
  //sel_HB_D_5(1)=selpar_Age0_HB_D;
  sel_comm_D_3(1)=selpar_Age0_HB_D;
 //Assume same sel of age 1's across periods 
  sel_HB_D_2(2)=selpar_Age1_HB_D; 
  sel_HB_D_3(2)=selpar_Age1_HB_D; 
  //sel_HB_D_4(2)=selpar_Age1_HB_D;
  //sel_HB_D_5(2)=selpar_Age1_HB_D;
  sel_comm_D_3(2)=selpar_Age1_HB_D;
 //Assume same sel of age 2's across periods 
  sel_HB_D_2(3)=selpar_Age2_HB_D; 
  sel_HB_D_3(3)=selpar_Age2_HB_D; 
  //sel_HB_D_4(3)=selpar_Age2_HB_D;
  //sel_HB_D_5(3)=selpar_Age2_HB_D;  
  sel_comm_D_3(3)=selpar_Age2_HB_D;
 //Assume full sel at age 3 across periods 
  sel_HB_D_2(4)=1.0; 
  sel_HB_D_3(4)=1.0; 
  //sel_HB_D_4(4)=1.0;
  //sel_HB_D_5(4)=1.0; 
  sel_comm_D_3(4)=1.0;  
   
  for (iage=5; iage<=nages; iage++)
      {sel_HB_D_2(iage)=prob_belowsizelim_block1(iage);
       sel_HB_D_3(iage)=prob_belowsizelim_block2(iage);
	   sel_comm_D_3(iage)=prob_belowsizelim_block3(iage);
       //sel_HB_D_4(iage)=prob_belowsizelim_block4(iage);
	   //sel_HB_D_4(iage)=prob_belowsizelim_block5(iage);
       }
	//cout<<"sel_HB_D_2"<<sel_HB_D_2<<endl;
	//cout<<"sel_HB_D_3"<<sel_HB_D_3<<endl;
	//cout<<"sel_HB_D_4"<<sel_HB_D_4<<endl;
	//cout<<"sel_HB_D_5"<<sel_HB_D_5<<endl;
	//cout<<"selpar_slope2_HBD4"<<selpar_slope2_HBD4<<endl;
	//cout<<"selpar_slope2_HBD5"<<selpar_slope2_HBD5<<endl;
  //Period 1: assumed same as in period 1, no commercial discards
  for (iyear=styr; iyear<=endyr_period1; iyear++)
      {sel_HB_D(iyear)=sel_HB_D_2;}

  //Period 2: 
  for (iyear=endyr_period1+1; iyear<=endyr_period2; iyear++)
      {sel_HB_D(iyear)=sel_HB_D_2;
       sel_comm_D(iyear)=sel_HB_D_2;
      }
	//cout<<"sel_HB_D after period 2 loop"<<sel_HB_D<<endl;
  //Period 3: Starts in 1999
  for (iyear=endyr_period2+1; iyear<=endyr; iyear++)
  {sel_HB_D(iyear)=sel_HB_D_3;
   sel_comm_D(iyear)=sel_HB_D_3;}  
	//cout<<"sel_HB_D after period 3 loop"<<sel_HB_D<<endl;
  

  //Period 4: Starts in 2007, hb and mrip only, overwrites last few yrs calculated for period 3
  for (iyear=endyr_recr_period3+1; iyear<=endyr_period4; iyear++)
  {sel_HB_D(iyear)=sel_HB_D_4;
   sel_comm_D(iyear)=sel_HB_D_3;
   
  }
  //Discard quota: wghted average,, overwrites last few yrs calculated for period 3 styr comm closed=2009
  for (iyear=styr_comm_closed; iyear<=endyr_period4; iyear++)
  {sel_comm_D(iyear)=Dprop_comm_sel_D*sel_HB_D_3 + Dprop_comm_sel_cL*sel_cL_2 +
                     Dprop_comm_sel_cP*sel_cP_2;
   sel_comm_D(iyear)=sel_comm_D(iyear)/max(sel_comm_D(iyear));} 
   
  //cout<<"sel_HB_D after period 4 loop"<<sel_HB_D<<endl;
  //Period 5: Starts in 2013 
  for (iyear=endyr_period4+1; iyear<=endyr; iyear++)
  {sel_HB_D(iyear)=sel_HB_D_5;
   sel_comm_D(iyear)=sel_comm_D_3;
   //Discard quota: wghted average,, overwrites last few yrs calculated for period 3 styr comm closed=2009
   sel_comm_D(iyear)=Dprop_comm_sel_D*sel_comm_D_3 + Dprop_comm_sel_cL*sel_cL_2 +
                     Dprop_comm_sel_cP*sel_cP_2;
   sel_comm_D(iyear)=sel_comm_D(iyear)/max(sel_comm_D(iyear));} 
  
   //cout<<"sel_HB_D after period 5 loop"<<sel_HB_D<<endl;
   //cout<<"sel_HB_D"<<sel_HB_D_4;
   //cout<<"sel_HB_4"<<sel_HB_4;
   //Applied the number of days open to the landings selectivity and the number of days
   //Rec fishery closure in 2011 (10/11) time block 4
   //sel_HB_D(2011)(4,nages)=(182.0/365.0)*sel_HB_4(4,nages) + (183.0/365.0)*sel_HB_D_4(4,nages);  //Rec fishery closure in 2012 (11/12) time block 4
   //sel_HB_D(2012)(4,nages)=(96.0/365.0)*sel_HB_4(4,nages) + (269.0/365.0)*sel_HB_D_4(4,nages);   //Rec fishery closure in 2013 (12/13) time block 4 
   //sel_HB_D(2013)(4,nages)=(214.0/365.0)*sel_HB_5(4,nages) + (151.0/365.0)*sel_HB_D_5(4,nages);  //Rec fishery closure in 2014 (13/14) time block 5 
   
   //sel_HB_D(2011)=sel_HB_D(2011)/max(sel_HB_D(2011));
   //sel_HB_D(2012)=sel_HB_D(2012)/max(sel_HB_D(2012));
   //sel_HB_D(2013)=sel_HB_D(2013)/max(sel_HB_D(2013));
    //sel_HB_D(2011)(4,nages)=sel_HB_D(2011)(4,nages)/max(sel_HB_D(2011)(4,nages));
    //sel_HB_D(2012)(4,nages)=sel_HB_D(2012)(4,nages)/max(sel_HB_D(2012)(4,nages));
    //sel_HB_D(2013)(4,nages)=sel_HB_D(2013)(4,nages)/max(sel_HB_D(2013)(4,nages));
   //cout<<"sel_HB_D after open/closed adjustment"<<sel_HB_D<<endl;   
   //cout<<"sel_HB_D"<<sel_HB_D_4;
   //cout<<"sel_HB_term"<<sel_HB_D;
   //mrip discard selectivity same as headboat;
   sel_mrip_D=sel_HB_D;

FUNCTION get_mortality
  Fsum.initialize();
  Fapex.initialize();
  F.initialize();
  //initialization F is avg from first 3 yrs of observed landings 
  log_F_dev_init_cL=sum(log_F_dev_cL(styr_cL_L,(styr_cL_L+2)))/3.0;
  log_F_dev_init_cP=sum(log_F_dev_cL(styr_cP_L,(styr_cP_L+2)))/3.0;
  log_F_dev_init_cT=sum(log_F_dev_cL(styr_cT_L,(styr_cT_L+2)))/3.0;
  log_F_init_HB=sum(log_F_dev_HB(styr_HB_L,(styr_HB_L+2)))/3.0;
  log_F_dev_init_mrip=sum(log_F_dev_mrip(styr_mrip_L,(styr_mrip_L+2)))/3.0;
  log_F_dev_init_mrip_D=sum(log_F_dev_mrip_D(styr_mrip_D,(styr_mrip_D+2)))/3.0;
  
  log_F_dev_comm_D2=sum(log_F_dev_comm_D(styr_cL_D,(styr_cL_D+5)))/6.0; //for comm D 1984-1992
  
  //cout<<styr<<endl;  
  for (iyear=styr; iyear<=endyr; iyear++) 
  {
    //-------------
    if(iyear>=styr_cL_L & iyear<=endyr_cL_L)
    {  F_cL_out(iyear)=mfexp(log_avg_F_cL+log_F_dev_cL(iyear)); //}    
       //if (iyear<styr_cL_L){F_cL_out(iyear)=mfexp(log_avg_F_cL+log_F_dev_init_cL);}        
       F_cL(iyear)=sel_cL(iyear)*F_cL_out(iyear);
       Fsum(iyear)+=F_cL_out(iyear);
    }
    
    //-------------
    if(iyear>=styr_cP_L & iyear<=endyr_cP_L)
    {  F_cP_out(iyear)=mfexp(log_avg_F_cP+log_F_dev_cP(iyear)); //}
       //if (iyear<styr_cP_L) {F_cP_out(iyear)=0.0;}
       F_cP(iyear)=sel_cP(iyear)*F_cP_out(iyear);
       Fsum(iyear)+=F_cP_out(iyear);
    }

    //-------------
    if(iyear>=styr_cT_L & iyear<=endyr_cT_L)
    {  F_cT_out(iyear)=mfexp(log_avg_F_cT+log_F_dev_cT(iyear)); //}
      // if (iyear<styr_cT_L) {F_cT_out(iyear)=0.0;}
       F_cT(iyear)=sel_cT(iyear)*F_cT_out(iyear);
       Fsum(iyear)+=F_cT_out(iyear);
    } 
        
    //-------------
    if(iyear>=styr_HB_L & iyear<=endyr_HB_L)
    {  F_HB_out(iyear)=mfexp(log_avg_F_HB+log_F_dev_HB(iyear));//}
    // if (iyear<styr_HB_L){F_HB_out(iyear)=mfexp(log_avg_F_HB+log_F_init_HB);}
       F_HB(iyear)=sel_HB(iyear)*F_HB_out(iyear);    
       Fsum(iyear)+=F_HB_out(iyear);
    }
    
    //-------------
    if(iyear>=styr_mrip_L & iyear<=endyr_mrip_L)
       {F_mrip_out(iyear)=mfexp(log_avg_F_mrip+log_F_dev_mrip(iyear));}
    if (iyear<styr_mrip_L){F_mrip_out(iyear)=mfexp(log_avg_F_mrip+log_F_dev_init_mrip);}
    F_mrip(iyear)=sel_mrip(iyear)*F_mrip_out(iyear);    
    Fsum(iyear)+=F_mrip_out(iyear);
    
    //discards-------------

    if(iyear>=styr_cL_D & iyear<=endyr_cL_D)
      {F_comm_D_out(iyear)=mfexp(log_avg_F_comm_D+log_F_dev_comm_D(iyear));}
    if(iyear > endyr_period1 & iyear < styr_cL_D)
      {F_comm_D_out(iyear)=mfexp(log_avg_F_comm_D+log_F_dev_comm_D2);}
    F_comm_D(iyear)=sel_comm_D(iyear)*F_comm_D_out(iyear);  
    Fsum(iyear)+=F_comm_D_out(iyear);    


    if(iyear>=styr_HB_D & iyear<=endyr_HB_D)
    {  F_HB_D_out(iyear)=mfexp(log_avg_F_HB_D+log_F_dev_HB_D(iyear)); 
       F_HB_D(iyear)=sel_HB_D(iyear)*F_HB_D_out(iyear); 
       Fsum(iyear)+=F_HB_D_out(iyear); 
    }
    if(iyear<styr_mrip_D)
       {F_mrip_D_out(iyear)=mfexp(log_avg_F_mrip_D+log_F_dev_init_mrip_D);}
    if(iyear>=styr_mrip_D& iyear<=endyr_mrip_D)
       { F_mrip_D_out(iyear)=mfexp(log_avg_F_mrip_D+log_F_dev_mrip_D(iyear));}
    F_mrip_D(iyear)=sel_mrip_D(iyear)*F_mrip_D_out(iyear);  
    Fsum(iyear)+=F_mrip_D_out(iyear);

    //Total F at age
    F(iyear)=F_cL(iyear); //first in additive series (NO +=)
    F(iyear)+=F_cP(iyear);
    F(iyear)+=F_cT(iyear);    
    F(iyear)+=F_HB(iyear);
    F(iyear)+=F_mrip(iyear);
   
    F(iyear)+=F_comm_D(iyear);
    F(iyear)+=F_HB_D(iyear);
    F(iyear)+=F_mrip_D(iyear);
    
    Fapex(iyear)=max(F(iyear));
   // cout <<"Fapex" << Fapex <<endl;
    Z(iyear)=M+F(iyear);
  }  //end iyear
 
 
FUNCTION get_bias_corr
  //may exclude last BiasCor_exclude_yrs yrs bc constrained or lack info to estimate
  var_rec_dev=norm2(log_rec_dev(styr_rec_dev,(endyr-BiasCor_exclude_yrs))-
              sum(log_rec_dev(styr_rec_dev,(endyr-BiasCor_exclude_yrs)))
              /(nyrs_rec-BiasCor_exclude_yrs))/(nyrs_rec-BiasCor_exclude_yrs-1.0); 
			  //cout<<"var rec dev"<<var_rec_dev<<endl;
  //if (set_BiasCor <= 0.0) {BiasCor=mfexp(var_rec_dev/2.0);}   //bias correction
  rec_sigma_sq=square(rec_sigma);
  rec_sigma_sqd2=rec_sigma_sq/2.0;
  if (set_BiasCor <= 0.0) {BiasCor=mfexp(rec_sigma_sqd2);}   //bias correction               
  else {BiasCor=set_BiasCor;}

//  cout<< "get bias corr" << var_rec_dev <<endl;

FUNCTION get_numbers_at_age
//Initialization
  S0=spr_F0*R0;
  R_virgin=(R0/((5.0*steep-1.0)*spr_F0))*(BiasCor*4.0*steep*spr_F0-spr_F0*(1.0-steep));
  B0=bpr_F0*R_virgin;   
  B0_q_DD=R_virgin*sum(elem_prod(N_bpr_F0(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages))); 

  F_initial=sel_cL(styr)*mfexp(log_avg_F_cL+log_F_dev_init_cL)+
            sel_cP(styr)*mfexp(log_avg_F_cP+log_F_dev_init_cP)+
            sel_cT(styr)*mfexp(log_avg_F_cT+log_F_dev_init_cT)+
            sel_HB(styr)*mfexp(log_avg_F_HB+log_F_init_HB)+
            sel_mrip(styr)*mfexp(log_avg_F_mrip+log_F_dev_init_mrip)+
            sel_mrip_D(styr)*mfexp(log_avg_F_mrip_D+log_F_dev_init_mrip_D);
  Z_initial=M+F_init_ratio*F_initial;


//Initial equilibrium age structure
  N_spr_initial(1)=1.0*mfexp(-1.0*Z_initial(1)*spawn_time_frac); //at peak spawning time;
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   mfexp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
//    N_spr_F_init_mdyr(1,(nages-1))=elem_prod(N_spr_initial(1,(nages-1)),
//                                   mfexp((-1.*(M(nages-1)+ F_initial))/2.0));   

  spr_initial=sum(elem_prod(N_spr_initial,reprod));

  if (styr==styr_rec_dev) {R1=(R0/((5.0*steep-1.0)*spr_initial))*
                 (4.0*steep*spr_initial-spr_F0*(1.0-steep));} //without bias correction (deviation added later)
  else {R1=(R0/((5.0*steep-1.0)*spr_initial))*
                 (BiasCor*4.0*steep*spr_initial-spr_F0*(1.0-steep));} //with bias correction                 

  if(R1<10.0) {R1=10.0;} //Avoid negative (or unreasonably low) popn sizes during search algorithm

  
   
//Compute equilibrium age structure for first year
  N_initial_eq(1)=R1;
  for (iage=2; iage<=nages; iage++)
  {
    N_initial_eq(iage)=N_initial_eq(iage-1)*
        mfexp(-1.0*(Z_initial(iage-1)));    
  }
  //plus group calculation
  N_initial_eq(nages)=N_initial_eq(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  
//Add deviations to initial equilibrium N
  N(styr)(2,nages)=elem_prod(N_initial_eq(2,nages),mfexp(log_Nage_dev));
   
  if (styr==styr_rec_dev) {N(styr,1)=N_initial_eq(1)*mfexp(log_rec_dev(styr_rec_dev));}
  else {N(styr,1)=N_initial_eq(1);}
  
  N_mdyr(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*0.5))); //mid year 
  N_spawn(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*spawn_time_frac))); //peak spawning time 

  SSB(styr)=sum(elem_prod(N_spawn(styr),reprod));
  MatFemB(styr)=sum(elem_prod(N_spawn(styr),reprod2));
  B_q_DD(styr)=sum(elem_prod(N(styr)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
    
//Rest of years 
  for (iyear=styr; iyear<endyr; iyear++)
  {
    if(iyear<(styr_rec_dev-1)) //recruitment follows S-R curve exactly
    {
        N(iyear+1,1)=0.0;   //no age 0's mature in SSB calculations, value replaced below for abundance calcs
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages));//plus group
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod)); 
        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod2));   
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));   
        N(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(iyear+1));
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year                    
    }
    
    else   //recruitment follows S-R curve with lognormal deviation
    {       
        N(iyear+1,1)=0.0;   //no age 0's mature in SSB calculations, value replaced below for abundance calcs
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages));//plus group
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));
        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod2));
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
        N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear+1))*mfexp(log_rec_dev(iyear+1));
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        
    }
  } //end iyear
  
    ////last year (projection) cannot compute recuitment of age 0's past terminal yr bc spawning occurs at spawn_time_frac
    //N(endyr+1,1)=SR_func(R0, steep, spr_F0, SSB(endyr));  
    //N(endyr+1)(2,nages)=++elem_prod(N(endyr)(1,nages-1),(mfexp(-1.*Z(endyr)(1,nages-1))));
    //N(endyr+1,nages)+=N(endyr,nages)*mfexp(-1.*Z(endyr,nages));//plus group
    //SSB(endyr+1)=sum(elem_prod(N(endyr+1),reprod));


//Time series of interest
  rec=column(N,1);
  SdS0=SSB/S0;
  
FUNCTION get_landings_numbers //Baranov catch eqn
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      L_cL_num(iyear,iage)=N(iyear,iage)*F_cL(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_cP_num(iyear,iage)=N(iyear,iage)*F_cP(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_cT_num(iyear,iage)=N(iyear,iage)*F_cT(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_HB_num(iyear,iage)=N(iyear,iage)*F_HB(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_mrip_num(iyear,iage)=N(iyear,iage)*F_mrip(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }
                
    pred_cL_L_knum(iyear)=sum(L_cL_num(iyear))/1000.0;
    pred_cP_L_knum(iyear)=sum(L_cP_num(iyear))/1000.0;
    pred_cT_L_knum(iyear)=sum(L_cT_num(iyear))/1000.0;    
    pred_HB_L_knum(iyear)=sum(L_HB_num(iyear))/1000.0;
    pred_mrip_L_knum(iyear)=sum(L_mrip_num(iyear))/1000.0;
  }
 
 
FUNCTION get_landings_wgt

////---Predicted landings------------------------
  for (iyear=styr; iyear<=endyr; iyear++)
  {    
    L_cL_klb(iyear)=elem_prod(L_cL_num(iyear),wgt_cL_klb(iyear));     //in 1000 lb
    L_cP_klb(iyear)=elem_prod(L_cP_num(iyear),wgt_cP_klb(iyear));     //in 1000 lb
    L_cT_klb(iyear)=elem_prod(L_cT_num(iyear),wgt_cT_klb(iyear));     //in 1000 lb    
    L_HB_klb(iyear)=elem_prod(L_HB_num(iyear),wgt_HB_klb(iyear));     //in 1000 lb
    L_mrip_klb(iyear)=elem_prod(L_mrip_num(iyear),wgt_mrip_klb(iyear));  //in 1000 lb    

    pred_cL_L_klb(iyear)=sum(L_cL_klb(iyear));
    pred_cP_L_klb(iyear)=sum(L_cP_klb(iyear));
    pred_cT_L_klb(iyear)=sum(L_cT_klb(iyear));
    pred_HB_L_klb(iyear)=sum(L_HB_klb(iyear));
    pred_mrip_L_klb(iyear)=sum(L_mrip_klb(iyear));
  }
 
   
    
FUNCTION get_dead_discards //Baranov catch eqn
  //dead discards at age (number fish) 

  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_comm_num(iyear,iage)=N(iyear,iage)*F_comm_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      D_HB_num(iyear,iage)=N(iyear,iage)*F_HB_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      D_mrip_num(iyear,iage)=N(iyear,iage)*F_mrip_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);        
    }
    pred_comm_D_knum(iyear)=sum(D_comm_num(iyear))/1000.0;         //pred annual dead discards in 1000s (for matching data)
    pred_comm_D_klb(iyear)=sum(elem_prod(D_comm_num(iyear),wgt_comm_D_klb(iyear)));  //annual dead discards in 1000 lb (for output only)

    pred_HB_D_knum(iyear)=sum(D_HB_num(iyear))/1000.0;             //pred annual dead discards in 1000s (for matching data)
    pred_HB_D_klb(iyear)=sum(elem_prod(D_HB_num(iyear),wgt_HB_D_klb(iyear)));        //annual dead discards in 1000 lb (for output only)

    pred_mrip_D_knum(iyear)=sum(D_mrip_num(iyear))/1000.0;         //pred annual dead discards in 1000s (for matching data)
    pred_mrip_D_klb(iyear)=sum(elem_prod(D_mrip_num(iyear),wgt_mrip_D_klb(iyear)));  //annual dead discards in 1000 lb (for output only)
  }

FUNCTION get_catchability_fcns    
 //Get rate increase if estimated, otherwise fixed above
  if (set_q_rate_phase>0.0)
  {
      for (iyear=styr_cL_cpue; iyear<=endyr_cL_cpue; iyear++)
      {   if (iyear>styr_cL_cpue & iyear <=2003) 
          {//q_rate_fcn_cL(iyear)=(1.0+q_rate)*q_rate_fcn_cL(iyear-1); //compound
             q_rate_fcn_cL(iyear)=(1.0+(iyear-styr_cL_cpue)*q_rate)*q_rate_fcn_cL(styr_cL_cpue);  //linear
          }
          if (iyear>2003) {q_rate_fcn_cL(iyear)=q_rate_fcn_cL(iyear-1);} 
      }   
      for (iyear=styr_HB_cpue; iyear<=endyr_HB_cpue; iyear++)
      {   if (iyear>styr_HB_cpue & iyear <=2003) 
          {//q_rate_fcn_HB(iyear)=(1.0+q_rate)*q_rate_fcn_HB(iyear-1); //compound
             q_rate_fcn_HB(iyear)=(1.0+(iyear-styr_HB_cpue)*q_rate)*q_rate_fcn_HB(styr_HB_cpue);  //linear
          }
          if (iyear>2003) {q_rate_fcn_HB(iyear)=q_rate_fcn_HB(iyear-1);} 
      }  
      //for (iyear=styr_HBD_cpue; iyear<=endyr_HBD_cpue; iyear++)
      //{   if (iyear>styr_HBD_cpue & iyear <=2003) 
      //    {//q_rate_fcn_HBD(iyear)=(1.0+q_rate)*q_rate_fcn_HBD(iyear-1); //compound
      //       q_rate_fcn_HBD(iyear)=(1.0+(iyear-styr_HBD_cpue)*q_rate)*q_rate_fcn_HBD(styr_HBD_cpue);  //linear
      //    }
      //    if (iyear>2003) {q_rate_fcn_HBD(iyear)=q_rate_fcn_HBD(iyear-1);} 
      //}
  } //end q_rate conditional      

 //Get density dependence scalar (=1.0 if density independent model is used)   
  if (q_DD_beta>0.0) 
  {
    B_q_DD+=dzero;
    for (iyear=styr;iyear<=endyr;iyear++)
        {q_DD_fcn(iyear)=pow(B0_q_DD,q_DD_beta)*pow(B_q_DD(iyear),-q_DD_beta);}
          //{q_DD_fcn(iyear)=1.0+4.0/(1.0+mfexp(0.75*(B_q_DD(iyear)-0.1*B0_q_DD))); }
  }  
  
      
FUNCTION get_indices
//---Predicted CPUEs------------------------
 //Survey 1: Mbft
  for (iyear=styr_Mbft_cpue; iyear<=endyr_Mbft_cpue; iyear++)
  {   //index in number units
      N_Mbft(iyear)=elem_prod(N_mdyr(iyear),sel_Mbft(iyear)); 
      pred_Mbft_cpue(iyear)=mfexp(log_q_Mbft)*sum(N_Mbft(iyear));
  }

 //Survey 2: CVID
  for (iyear=styr_Mcvt_cpue; iyear<=endyr_Mcvt_cpue; iyear++)
  {   //index in number units
      N_Mcvt(iyear)=elem_prod(N_mdyr(iyear),sel_Mcvt(iyear)); 
      pred_Mcvt_cpue(iyear)=mfexp(log_q_Mcvt)*sum(N_Mcvt(iyear));
  }

  //Survey 3: Vid
  //for (iyear=styr_Vid_cpue; iyear<=endyr_Vid_cpue; iyear++)
  //{   //index in number units
  //    N_Vid(iyear)=elem_prod(N_mdyr(iyear),sel_Mcvt(iyear)); 
  //    pred_Vid_cpue(iyear)=mfexp(log_q_Vid)*sum(N_Vid(iyear));
  //}
  
 //Commercial handline cpue
  q_cL(styr_cL_cpue)=mfexp(log_q_cL); 
  for (iyear=styr_cL_cpue; iyear<=endyr_cL_cpue; iyear++)
  {   //index in weight units. original index in lb and re-scaled. predicted in klb, but difference is absorbed by q
      N_cL(iyear)=elem_prod(elem_prod(N_mdyr(iyear),sel_cL(iyear)),wgt_cL_klb(iyear)); 
      pred_cL_cpue(iyear)=q_cL(iyear)*q_rate_fcn_cL(iyear)*q_DD_fcn(iyear)*sum(N_cL(iyear));
      if (iyear<endyr_cL_cpue){q_cL(iyear+1)=q_cL(iyear)*mfexp(q_RW_log_dev_cL(iyear));}
  }

 //Headboat cpue
  q_HB(styr_HB_cpue)=mfexp(log_q_HB);
  for (iyear=styr_HB_cpue; iyear<=endyr_HB_cpue; iyear++)
  {   //index in weight units. original index in lb and re-scaled. predicted in klb, but difference is absorbed by q
      N_HB(iyear)=elem_prod(elem_prod(N_mdyr(iyear),sel_HB(iyear)),wgt_HB_klb(iyear)); 
      pred_HB_cpue(iyear)=q_HB(iyear)*q_rate_fcn_HB(iyear)*q_DD_fcn(iyear)*sum(N_HB(iyear));
      if (iyear<endyr_HB_cpue){q_HB(iyear+1)=q_HB(iyear)*mfexp(q_RW_log_dev_HB(iyear));}
  }
  //HBD cpue
  //q_HBD(styr_HBD_cpue)=mfexp(log_q_HBD);
  //for (iyear=styr_HBD_cpue; iyear<=endyr_HBD_cpue; iyear++)
  //{   //index in number units
  //    N_HBD(iyear)=elem_prod(N_mdyr(iyear),sel_HB_D(iyear)); 
  //    pred_HBD_cpue(iyear)=q_HBD(iyear)*q_rate_fcn_HBD(iyear)*q_DD_fcn(iyear)*sum(N_HBD(iyear));
  //    if (iyear<endyr_HBD_cpue){q_HBD(iyear+1)=q_HBD(iyear)*mfexp(q_RW_log_dev_HBD(iyear));}
  //}
  
FUNCTION get_length_comps
  
  //Mbft
  for (iyear=1;iyear<=nyr_Mbft_lenc;iyear++)
      {pred_Mbft_lenc(iyear)=(N_Mbft(yrs_Mbft_lenc(iyear))*lenprob)/sum(N_Mbft(yrs_Mbft_lenc(iyear)));} 

  //Commercial lines
  for (iyear=1;iyear<=nyr_cL_lenc;iyear++) //all yrs within periods 2,3
  {  if (yrs_cL_lenc(iyear)<=endyr_period2)
     {pred_cL_lenc(iyear)=(L_cL_num(yrs_cL_lenc(iyear))*lenprob_cL2)
                          /sum(L_cL_num(yrs_cL_lenc(iyear)));       
     } 
     if (yrs_cL_lenc(iyear)>endyr_period2)
     {pred_cL_lenc(iyear)=(L_cL_num(yrs_cL_lenc(iyear))*lenprob_cL3)
                          /sum(L_cL_num(yrs_cL_lenc(iyear)));       
     }  
  }

  //Commercial pots: pooled all from period 2
  L_cP_num_pool.initialize();
  for (iyear=1;iyear<=nyr_cP_lenc_pool;iyear++)
  {  L_cP_num_pool_yr(iyear)=nsamp_cP_lenc_pool(iyear)*L_cP_num(yrs_cP_lenc_pool(iyear))
                            /sum(L_cP_num(yrs_cP_lenc_pool(iyear)));                             
     if (yrs_cP_lenc_pool(iyear)<=endyr_period2) {L_cP_num_pool(1)+=L_cP_num_pool_yr(iyear);}                             
  } 
  for (iyear=1;iyear<=nyr_cP_lenc;iyear++) //all yrs within periods 2
  {  if (yrs_cP_lenc(iyear)<=endyr_period2)
       {pred_cP_lenc(iyear)=(L_cP_num_pool(iyear)*lenprob_cP2)/sum(L_cP_num_pool(iyear)); } 
	   pred_cP_lenc(iyear)=(L_cP_num(yrs_cP_lenc(iyear))*lenprob_cP3)  //added to calculate comps for last period
						  /sum(L_cP_num(yrs_cP_lenc(iyear)));
  }  
   

 
 //Headboat 
  for (iyear=1;iyear<=nyr_HB_lenc;iyear++)  //all in periods 1,2,3
  {  if (yrs_HB_lenc(iyear)<=endyr_period1)
     {pred_HB_lenc(iyear)=(L_HB_num(yrs_HB_lenc(iyear))*lenprob_HB1)
                          /sum(L_HB_num(yrs_HB_lenc(iyear)));       
     } 
     if (yrs_HB_lenc(iyear)>endyr_period1 & yrs_HB_lenc(iyear)<=endyr_period2)
     {pred_HB_lenc(iyear)=(L_HB_num(yrs_HB_lenc(iyear))*lenprob_HB2)
                          /sum(L_HB_num(yrs_HB_lenc(iyear)));       
     }  
     if (yrs_HB_lenc(iyear)>endyr_period2)
     {pred_HB_lenc(iyear)=(L_HB_num(yrs_HB_lenc(iyear))*lenprob_HB3)
                          /sum(L_HB_num(yrs_HB_lenc(iyear)));       
     }  
  }
 //HB discards 
  for (iyear=1;iyear<=nyr_HB_D_lenc;iyear++) //all yrs within period 3,4
  {  if (yrs_HB_D_lenc(iyear)<=endyr_recr_period3)
     {pred_HB_D_lenc(iyear)=(D_HB_num(yrs_HB_D_lenc(iyear))*lenprob_HB_D3)
                        /sum(D_HB_num(yrs_HB_D_lenc(iyear)));}      
    if (yrs_HB_D_lenc(iyear)>endyr_recr_period3)
     {pred_HB_D_lenc(iyear)=(D_HB_num(yrs_HB_D_lenc(iyear))*lenprob_HB_D4)
                        /sum(D_HB_num(yrs_HB_D_lenc(iyear)));}                      
  }

 
 //MRIP
  for (iyear=1;iyear<=nyr_mrip_lenc;iyear++)  //all in periods 1,2,3
  {  if (yrs_mrip_lenc(iyear)<=endyr_period1)
     {pred_mrip_lenc(iyear)=(L_mrip_num(yrs_mrip_lenc(iyear))*lenprob_mrip1)
                          /sum(L_mrip_num(yrs_mrip_lenc(iyear)));       
     } 
     if (yrs_mrip_lenc(iyear)>endyr_period1 & yrs_mrip_lenc(iyear)<=endyr_period2)
     {pred_mrip_lenc(iyear)=(L_mrip_num(yrs_mrip_lenc(iyear))*lenprob_mrip2)
                          /sum(L_mrip_num(yrs_mrip_lenc(iyear)));       
     }  
     if (yrs_mrip_lenc(iyear)>endyr_period2 & yrs_mrip_lenc(iyear)<=endyr_recr_period3)
     {pred_mrip_lenc(iyear)=(L_mrip_num(yrs_mrip_lenc(iyear))*lenprob_mrip3)
                          /sum(L_mrip_num(yrs_mrip_lenc(iyear)));       
     }  
     if (yrs_mrip_lenc(iyear)>endyr_recr_period3 &
	 yrs_mrip_lenc(iyear)<=endyr_period4)
     {pred_mrip_lenc(iyear)=(L_mrip_num(yrs_mrip_lenc(iyear))*lenprob_mrip4)
                          /sum(L_mrip_num(yrs_mrip_lenc(iyear)));       
     }  
	 if (yrs_mrip_lenc(iyear)>endyr_period4)
     {pred_mrip_lenc(iyear)=(L_mrip_num(yrs_mrip_lenc(iyear))*lenprob_mrip5)
                          /sum(L_mrip_num(yrs_mrip_lenc(iyear)));       
     }  
  } 
 
  
FUNCTION get_age_comps

 //MARMAP bft
  for (iyear=1;iyear<=nyr_Mbft_agec;iyear++)
  {
    ErrorFree_Mbft_agec(iyear)=N_Mbft(yrs_Mbft_agec(iyear))/sum(N_Mbft(yrs_Mbft_agec(iyear)));
    pred_Mbft_agec(iyear)=age_error*ErrorFree_Mbft_agec(iyear);                      
  }

 //MARMAP cvt
  for (iyear=1;iyear<=nyr_Mcvt_agec;iyear++)
  {
    ErrorFree_Mcvt_agec(iyear)=N_Mcvt(yrs_Mcvt_agec(iyear))/sum(N_Mcvt(yrs_Mcvt_agec(iyear)));
    pred_Mcvt_agec(iyear)=age_error*ErrorFree_Mcvt_agec(iyear);                      
  }
    
 //Commercial lines
  for (iyear=1;iyear<=nyr_cL_agec;iyear++)
  {
    ErrorFree_cL_agec(iyear)=L_cL_num(yrs_cL_agec(iyear))/sum(L_cL_num(yrs_cL_agec(iyear)));
    pred_cL_agec(iyear)=age_error*ErrorFree_cL_agec(iyear);                      
  }

 //Commercial pots
  for (iyear=1;iyear<=nyr_cP_agec;iyear++)
  {
    ErrorFree_cP_agec(iyear)=L_cP_num(yrs_cP_agec(iyear))/sum(L_cP_num(yrs_cP_agec(iyear)));
    pred_cP_agec(iyear)=age_error*ErrorFree_cP_agec(iyear);                      
  }
  
 //Headboat
  for (iyear=1;iyear<=nyr_HB_agec;iyear++)
  {
    ErrorFree_HB_agec(iyear)=L_HB_num(yrs_HB_agec(iyear))/sum(L_HB_num(yrs_HB_agec(iyear)));
    pred_HB_agec(iyear)=age_error*ErrorFree_HB_agec(iyear);                    
  }
  
 ////mrip
 // for (iyear=1;iyear<=nyr_mrip_agec;iyear++)
 // {
 //   ErrorFree_mrip_agec(iyear)=L_mrip_num(yrs_mrip_agec(iyear))/sum(L_mrip_num(yrs_mrip_agec(iyear)));
 //   pred_mrip_agec(iyear)=age_error*ErrorFree_mrip_agec(iyear);                    
 // }
  
    
////--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_weighted_current 
  F_temp_sum=0.0;
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_cL+
        sum(log_F_dev_cL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_cP+
        sum(log_F_dev_cP((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_HB+
        sum(log_F_dev_HB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_mrip+
        sum(log_F_dev_mrip((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_comm_D+
        sum(log_F_dev_comm_D((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_HB_D+
        sum(log_F_dev_HB_D((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_mrip_D+
        sum(log_F_dev_mrip_D((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
      
  F_cL_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_cL+
        sum(log_F_dev_cL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_cP_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_cP+
        sum(log_F_dev_cP((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_HB_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_HB+
        sum(log_F_dev_HB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_mrip_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_mrip+
        sum(log_F_dev_mrip((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_comm_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_comm_D+
        sum(log_F_dev_comm_D((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_HB_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_HB_D+
        sum(log_F_dev_HB_D((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_mrip_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_mrip_D+
        sum(log_F_dev_mrip_D((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;

  log_F_dev_end_cL=sum(log_F_dev_cL((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_cP=sum(log_F_dev_cP((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  
  log_F_dev_end_HB=sum(log_F_dev_HB((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_mrip=sum(log_F_dev_mrip((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;

  log_F_dev_end_comm_D=sum(log_F_dev_comm_D((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_HB_D=sum(log_F_dev_HB_D((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_mrip_D=sum(log_F_dev_mrip_D((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;

  F_end_L=sel_cL(endyr)*mfexp(log_avg_F_cL+log_F_dev_end_cL)+
          sel_cP(endyr)*mfexp(log_avg_F_cP+log_F_dev_end_cP)+
          sel_HB(endyr)*mfexp(log_avg_F_HB+log_F_dev_end_HB)+
          sel_mrip(endyr)*mfexp(log_avg_F_mrip+log_F_dev_end_mrip);
                 
  F_end_D=sel_comm_D(endyr)*mfexp(log_avg_F_comm_D+log_F_dev_end_comm_D)+
          sel_HB_D(endyr)*mfexp(log_avg_F_HB_D+log_F_dev_end_HB_D)+
          sel_mrip_D(endyr)*mfexp(log_avg_F_mrip_D+log_F_dev_end_mrip_D);    

  F_end=F_end_L+F_end_D;
  F_end_apex=max(F_end);
  
  sel_wgted_tot=F_end/F_end_apex;
  sel_wgted_L=elem_prod(sel_wgted_tot, elem_div(F_end_L,F_end));
  sel_wgted_D=elem_prod(sel_wgted_tot, elem_div(F_end_D,F_end));
  
  wgt_wgted_L_denom=F_cL_prop+F_cP_prop+F_HB_prop+F_mrip_prop;
  wgt_wgted_L_klb=F_cL_prop/wgt_wgted_L_denom*wgt_cL_klb(endyr)+
              F_cP_prop/wgt_wgted_L_denom*wgt_cP_klb(endyr)+
              F_HB_prop/wgt_wgted_L_denom*wgt_HB_klb(endyr)+
              F_mrip_prop/wgt_wgted_L_denom*wgt_mrip_klb(endyr);                

  wgt_wgted_D_denom=F_comm_D_prop+F_HB_D_prop+F_mrip_D_prop;
  wgt_wgted_D_klb=F_comm_D_prop/wgt_wgted_D_denom*wgt_comm_D_klb(endyr)+
              F_HB_D_prop/wgt_wgted_D_denom*wgt_HB_D_klb(endyr)+
              F_mrip_D_prop/wgt_wgted_D_denom*wgt_mrip_D_klb(endyr);                
  
FUNCTION get_msy
  
  //compute values as functions of F
  for(ff=1; ff<=n_iter_msy; ff++)
  {
    //uses fishery-weighted F's
    Z_age_msy=0.0;
    F_L_age_msy=0.0;
    F_D_age_msy=0.0;
      
    F_L_age_msy=F_msy(ff)*sel_wgted_L;
    F_D_age_msy=F_msy(ff)*sel_wgted_D;
    Z_age_msy=M+F_L_age_msy+F_D_age_msy;         
    
    N_age_msy(1)=1.0;
    for (iage=2; iage<=nages; iage++)
    {
      N_age_msy(iage)=N_age_msy(iage-1)*mfexp(-1.*Z_age_msy(iage-1));
    }
    N_age_msy(nages)=N_age_msy(nages)/(1.0-mfexp(-1.*Z_age_msy(nages)));
    N_age_msy_spawn(1,(nages-1))=elem_prod(N_age_msy(1,(nages-1)),
                                   mfexp((-1.*Z_age_msy(1,(nages-1)))*spawn_time_frac));                 
    N_age_msy_spawn(nages)=(N_age_msy_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_msy(nages-1)*(1.0-spawn_time_frac) + 
                                 Z_age_msy(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_msy(nages)));
                     
    spr_msy(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
        

    //Compute equilibrium values of R (including bias correction), SSB and Yield at each F
    R_eq(ff)=(R0/((5.0*steep-1.0)*spr_msy(ff)))*
                 (BiasCor*4.0*steep*spr_msy(ff)-spr_F0*(1.0-steep));
    if (R_eq(ff)<dzero) {R_eq(ff)=dzero;}    
    N_age_msy*=R_eq(ff);
    N_age_msy_spawn*=R_eq(ff);
    
    for (iage=1; iage<=nages; iage++)
    {
      L_age_msy(iage)=N_age_msy(iage)*(F_L_age_msy(iage)/Z_age_msy(iage))*
                      (1.-mfexp(-1.*Z_age_msy(iage)));
      D_age_msy(iage)=N_age_msy(iage)*(F_D_age_msy(iage)/Z_age_msy(iage))*
                      (1.-mfexp(-1.0*Z_age_msy(iage)));
    }
    
    
    SSB_eq(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
    B_eq(ff)=sum(elem_prod(N_age_msy,wgt_mt));
    L_eq_klb(ff)=sum(elem_prod(L_age_msy,wgt_wgted_L_klb));
    L_eq_knum(ff)=sum(L_age_msy)/1000.0;  
    D_eq_klb(ff)=sum(elem_prod(D_age_msy,wgt_wgted_D_klb));    
    D_eq_knum(ff)=sum(D_age_msy)/1000.0;
  }  //end ff loop
  
  msy_klb_out=max(L_eq_klb);
  
  for(ff=1; ff<=n_iter_msy; ff++)
  {
   if(L_eq_klb(ff) == msy_klb_out) 
      {    
        SSB_msy_out=SSB_eq(ff);
        B_msy_out=B_eq(ff);
        R_msy_out=R_eq(ff);
        msy_knum_out=L_eq_knum(ff);
        D_msy_knum_out=D_eq_knum(ff);
        D_msy_klb_out=D_eq_klb(ff);
        F_msy_out=F_msy(ff);  
        spr_msy_out=spr_msy(ff);      
      }
  }


//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_per_recruit_stuff

  //static per-recruit stuff
 
  for(iyear=styr; iyear<=endyr; iyear++)
  {
    N_age_spr(1)=1.0;
    for(iage=2; iage<=nages; iage++)
    {
      N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z(iyear,iage-1));
    }
    N_age_spr(nages)=N_age_spr(nages)/(1.0-mfexp(-1.*Z(iyear,nages)));    
    N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                mfexp(-1.*Z(iyear)(1,(nages-1))*spawn_time_frac));
    N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z(iyear)(nages-1)*(1.0-spawn_time_frac) + Z(iyear)(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z(iyear)(nages)));           
    spr_static(iyear)=sum(elem_prod(N_age_spr_spawn,reprod))/spr_F0;
  }
  

  //compute SSB/R and YPR as functions of F
  for(ff=1; ff<=n_iter_spr; ff++)
  {
    //uses fishery-weighted F's, same as in MSY calculations
    Z_age_spr=0.0;
    F_L_age_spr=0.0;

    F_L_age_spr=F_spr(ff)*sel_wgted_L;
    
    Z_age_spr=M+F_L_age_spr+F_spr(ff)*sel_wgted_D;

    N_age_spr(1)=1.0;
    for (iage=2; iage<=nages; iage++)
    {
      N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));
    }
    N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
    N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
    N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));
    
    spr_spr(ff)=sum(elem_prod(N_age_spr_spawn,reprod));
    L_spr(ff)=0.0;
    for (iage=1; iage<=nages; iage++)
    {
      L_age_spr(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      L_spr(ff)+=L_age_spr(iage)*wgt_wgted_L_klb(iage)*1000.0; //in lb
    }   
  }

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_miscellaneous_stuff
  
  sigma_rec_dev=sqrt(var_rec_dev); //sample SD of predicted residuals (may not equal rec_sigma)  
  //len_cv=elem_div(len_sd,meanlen_TL);
  
  //compute total landings- and discards-at-age in 1000 fish and klb
  L_total_num.initialize();
  L_total_klb.initialize();
  D_total_num.initialize();
  D_total_klb.initialize();
  L_total_knum_yr.initialize();
  L_total_klb_yr.initialize();  
  D_total_knum_yr.initialize();
  D_total_klb_yr.initialize();
  
  for(iyear=styr; iyear<=endyr; iyear++)
  {
        L_total_klb_yr(iyear)= pred_cL_L_klb(iyear)+pred_cP_L_klb(iyear)+pred_cT_L_klb(iyear)+
                               pred_HB_L_klb(iyear)+pred_mrip_L_klb(iyear);
        L_total_knum_yr(iyear)=pred_cL_L_knum(iyear)+pred_cP_L_knum(iyear)+pred_cT_L_knum(iyear)+
                               pred_HB_L_knum(iyear)+pred_mrip_L_knum(iyear);
                
        D_total_knum_yr(iyear)+=pred_comm_D_knum(iyear);
        D_total_klb_yr(iyear)+=pred_comm_D_klb(iyear);

        D_total_knum_yr(iyear)+=pred_HB_D_knum(iyear);
        D_total_klb_yr(iyear)+=pred_HB_D_klb(iyear);

        D_total_knum_yr(iyear)+=pred_mrip_D_knum(iyear);
        D_total_klb_yr(iyear)+=pred_mrip_D_klb(iyear);

        D_comm_klb(iyear)=elem_prod(D_comm_num(iyear),wgt_comm_D_klb(iyear));     //in 1000 lb
        D_HB_klb(iyear)=elem_prod(D_HB_num(iyear),wgt_HB_D_klb(iyear));         //in 1000 lb
        D_mrip_klb(iyear)=elem_prod(D_mrip_num(iyear),wgt_mrip_D_klb(iyear));   //in 1000 lb 
        
        B(iyear)=elem_prod(N(iyear),wgt_mt);
        totN(iyear)=sum(N(iyear));
        totB(iyear)=sum(B(iyear));              
  }
  
  L_total_num=(L_cL_num+L_cP_num+L_cT_num+L_HB_num+L_mrip_num); //landings at age in number fish
  L_total_klb=L_cL_klb+L_cP_klb+L_cT_klb+L_HB_klb+L_mrip_klb;   //landings at age in klb whole weight

  D_total_num=(D_comm_num+D_HB_num+D_mrip_num);          //discards at age in number fish
  D_total_klb=D_comm_klb+D_HB_klb+D_mrip_klb;            //discards at age in klb whole weight
  
  //B(endyr+1)=elem_prod(N(endyr+1),wgt_mt);
  //totN(endyr+1)=sum(N(endyr+1));
  //totB(endyr+1)=sum(B(endyr+1));  
  
//  steep_sd=steep;
//  fullF_sd=Fsum;
  
  if(F_msy_out>0)
    {
      FdF_msy=Fapex/F_msy_out;
      FdF_msy_end=FdF_msy(endyr_cL_L);
      //cout<<"Fstatus2012"<< FdF_msy_end << endl;
      FdF_msy_end_mean=pow((FdF_msy(endyr_cL_L)*FdF_msy(endyr_cL_L-1)),(1.0/2.0));
      //FdF_msy_end_mean=pow((FdF_msy(endyr)*FdF_msy(endyr-1)*FdF_msy(endyr-2)),(1.0/3.0));
    }
  if(SSB_msy_out>0)
    {
      SdSSB_msy=SSB/SSB_msy_out;
      SdSSB_msy_end=SdSSB_msy(endyr);
    }  

   //fill in log recruitment deviations for yrs they are nonzero
   for(iyear=styr_rec_dev; iyear<=endyr; iyear++)
     {log_rec_dev_output(iyear)=log_rec_dev(iyear);}
   //fill in log Nage deviations for ages they are nonzero (ages2+)
   for(iage=2; iage<=nages; iage++)
   { log_Nage_dev_output(iage)=log_Nage_dev(iage);}


FUNCTION get_effective_sample_sizes
  neff_Mbft_lenc_allyr_out=missing;//"missing" defined in admb2r.cpp
  neff_cL_lenc_allyr_out=missing; 
  neff_cP_lenc_allyr_out=missing;
  neff_HB_lenc_allyr_out=missing;
  neff_HB_D_lenc_allyr_out=missing;
  neff_mrip_lenc_allyr_out=missing;
  neff_Mbft_agec_allyr_out=missing;
  neff_Mcvt_agec_allyr_out=missing;
  neff_cL_agec_allyr_out=missing;
  neff_cP_agec_allyr_out=missing;      
  neff_HB_agec_allyr_out=missing;
  //neff_mrip_agec_allyr_out=missing;
  
      for (iyear=1; iyear<=nyr_Mbft_lenc; iyear++)
         {if (nsamp_Mbft_lenc(iyear)>=minSS_lenc)
            {neff_Mbft_lenc_allyr_out(yrs_Mbft_lenc(iyear))=multinom_eff_N(pred_Mbft_lenc(iyear),obs_Mbft_lenc(iyear));} 
          else {neff_Mbft_lenc_allyr_out(yrs_Mbft_lenc(iyear))=-99;}
         }                  

      for (iyear=1; iyear<=nyr_cL_lenc; iyear++)
         {if (nsamp_cL_lenc(iyear)>=minSS_lenc)
            {neff_cL_lenc_allyr_out(yrs_cL_lenc(iyear))=multinom_eff_N(pred_cL_lenc(iyear),obs_cL_lenc(iyear));} 
          else {neff_cL_lenc_allyr_out(yrs_cL_lenc(iyear))=-99;}
         }                  

      for (iyear=1; iyear<=nyr_cP_lenc; iyear++)
         {if (nsamp_cP_lenc(iyear)>=minSS_lenc)
             {neff_cP_lenc_allyr_out(yrs_cP_lenc(iyear))=multinom_eff_N(pred_cP_lenc(iyear),obs_cP_lenc(iyear));}                            
          else {neff_cP_lenc_allyr_out(yrs_cP_lenc(iyear))=-99;}
         }

      for (iyear=1; iyear<=nyr_HB_lenc; iyear++)
         {if (nsamp_HB_lenc(iyear)>=minSS_lenc)
             {neff_HB_lenc_allyr_out(yrs_HB_lenc(iyear))=multinom_eff_N(pred_HB_lenc(iyear),obs_HB_lenc(iyear));}                            
          else {neff_HB_lenc_allyr_out(yrs_HB_lenc(iyear))=-99;}
         }
                  
      for (iyear=1; iyear<=nyr_HB_D_lenc; iyear++)
         {if (nsamp_HB_D_lenc(iyear)>=minSS_lenc)
            {neff_HB_D_lenc_allyr_out(yrs_HB_D_lenc(iyear))=multinom_eff_N(pred_HB_D_lenc(iyear),obs_HB_D_lenc(iyear));}                            
           else {neff_HB_D_lenc_allyr_out(yrs_HB_D_lenc(iyear))=-99;}
         }

      for (iyear=1; iyear<=nyr_mrip_lenc; iyear++)
         {if (nsamp_mrip_lenc(iyear)>=minSS_lenc)
            {neff_mrip_lenc_allyr_out(yrs_mrip_lenc(iyear))=multinom_eff_N(pred_mrip_lenc(iyear),obs_mrip_lenc(iyear)); }                           
           else {neff_mrip_lenc_allyr_out(yrs_mrip_lenc(iyear))=-99;}
         }
         

      for (iyear=1; iyear<=nyr_Mbft_agec; iyear++)
         {if (nsamp_Mbft_agec(iyear)>=minSS_agec)
            {neff_Mbft_agec_allyr_out(yrs_Mbft_agec(iyear))=multinom_eff_N(pred_Mbft_agec(iyear),obs_Mbft_agec(iyear));}                            
          else {neff_Mbft_agec_allyr_out(yrs_Mbft_agec(iyear))=-99;}
         }

      for (iyear=1; iyear<=nyr_Mcvt_agec; iyear++)
         {if (nsamp_Mcvt_agec(iyear)>=minSS_agec)
            {neff_Mcvt_agec_allyr_out(yrs_Mcvt_agec(iyear))=multinom_eff_N(pred_Mcvt_agec(iyear),obs_Mcvt_agec(iyear));}                            
          else {neff_Mcvt_agec_allyr_out(yrs_Mcvt_agec(iyear))=-99;}
         }

      for (iyear=1; iyear<=nyr_cL_agec; iyear++)
         {if (nsamp_cL_agec(iyear)>=minSS_agec)
            {neff_cL_agec_allyr_out(yrs_cL_agec(iyear))=multinom_eff_N(pred_cL_agec(iyear),obs_cL_agec(iyear));}                            
          else {neff_cL_agec_allyr_out(yrs_cL_agec(iyear))=-99;}
         }

      for (iyear=1; iyear<=nyr_cP_agec; iyear++)
         {if (nsamp_cP_agec(iyear)>=minSS_agec)
            {neff_cP_agec_allyr_out(yrs_cP_agec(iyear))=multinom_eff_N(pred_cP_agec(iyear),obs_cP_agec(iyear));}                            
          else {neff_cP_agec_allyr_out(yrs_cP_agec(iyear))=-99;}
         }
           
      for (iyear=1; iyear<=nyr_HB_agec; iyear++)
         {if (nsamp_HB_agec(iyear)>=minSS_agec)
            {neff_HB_agec_allyr_out(yrs_HB_agec(iyear))=multinom_eff_N(pred_HB_agec(iyear),obs_HB_agec(iyear));}                            
          else {neff_HB_agec_allyr_out(yrs_HB_agec(iyear))=-99;}
         }

      //for (iyear=1; iyear<=nyr_mrip_agec; iyear++)
      //   {if (nsamp_mrip_agec(iyear)>=minSS_agec)
      //      {neff_mrip_agec_allyr_out(yrs_mrip_agec(iyear))=multinom_eff_N(pred_mrip_agec(iyear),obs_mrip_agec(iyear));}                            
      //    else {neff_mrip_agec_allyr_out(yrs_mrip_agec(iyear))=-99;}
      //   }
//FUNCTION get_projection
//  
//    switch(Fproj_switch){
//       case 1: //F=Fcurrent
//          F_reg_proj=Fend_mean;
//          break;
//       case 2: //F=Fmsy
//          F_reg_proj=F_msy_out;
//          break;
//       case 3: //F=F30
//          F_reg_proj=F30_out;
//          break;     
//       case 4: //F=F40
//          F_reg_proj=F40_out;
//          break;          		  
//       default: // no such switch available
//          cout << "Error in input: Projection switch Fproj_switch must be set to 1, 2, 3, or 4." << endl;
//          cout << "Presently it is set to " << Fproj_switch <<"."<< endl;
//          exit(0);          
//   }
//
//  N_proj(styr_proj)=N(endyr+1); //initial conditions computed previously
// 
//  for (iyear=styr_proj; iyear<=endyr_proj; iyear++) //recruitment follows S-R curve (with bias correction) exactly
//  {     
//        if (iyear<styr_regs) {F_proj(iyear)=Fend_mean;}
//		else {F_proj(iyear)=Fproj_mult*F_reg_proj;}
//		
//		FL_age_proj=sel_wgted_L*F_proj(iyear);
//		FD_age_proj=sel_wgted_D*F_proj(iyear);
//		
//        Z_proj(iyear)=M+FL_age_proj+FD_age_proj;
//        N_spawn_proj(iyear)(1,nages)=elem_prod(N_proj(iyear)(1,nages),(mfexp(-1.*(Z_proj(iyear)(1,nages))*spawn_time_frac))); //peak spawning time
//		SSB_proj(iyear)= sum(elem_prod(N_spawn_proj(iyear),reprod));
//        B_proj(iyear)=sum(elem_prod(N_proj(iyear),wgt_mt)); //uses spawning weight
//	     
//		for (iage=1; iage<=nages; iage++)
//			{L_age_proj(iyear,iage)=N_proj(iyear,iage)*FL_age_proj(iage)*(1.-mfexp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
//		     D_age_proj(iyear,iage)=N_proj(iyear,iage)*FD_age_proj(iage)*(1.-mfexp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
//			}          
//        L_knum_proj(iyear)=sum(L_age_proj(iyear))/1000.0;
//	    D_knum_proj(iyear)=sum(D_age_proj(iyear))/1000.0;
//	    L_klb_proj(iyear)=sum(elem_prod(L_age_proj(iyear),wgt_wgted_L_klb));     //in 1000 lb
//        D_klb_proj(iyear)=sum(elem_prod(D_age_proj(iyear),wgt_wgted_D_klb));     //in 1000 lb
//		
//		if (iyear<endyr_proj) {
//			N_proj(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB_proj(iyear),SR_switch);
//			N_proj(iyear+1)(2,nages)=++elem_prod(N_proj(iyear)(1,nages-1),(mfexp(-1.*Z_proj(iyear)(1,nages-1))));
//			N_proj(iyear+1,nages)+=N_proj(iyear,nages)*mfexp(-1.*Z_proj(iyear,nages)); //plus group		
//		}
//  }
//   R_proj=column(N_proj,1);                          
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   


FUNCTION evaluate_objective_function
  fval=0.0;
  fval_data=0.0;

//---likelihoods---------------------------

//---Indices-------------------------------

  f_Mbft_cpue=0.0;
  f_Mbft_cpue=lk_lognormal(pred_Mbft_cpue, obs_Mbft_cpue, Mbft_cpue_cv, w_I_Mbft);
  fval+=f_Mbft_cpue;
  fval_data+=f_Mbft_cpue;

  f_Mcvt_cpue=0.0;
  f_Mcvt_cpue=lk_lognormal(pred_Mcvt_cpue, obs_Mcvt_cpue, Mcvt_cpue_cv, w_I_Mcvt);
  fval+=f_Mcvt_cpue;
  fval_data+=f_Mcvt_cpue;
  
  //f_Vid_cpue=0.0;
  //f_Vid_cpue=lk_lognormal(pred_Vid_cpue, obs_Vid_cpue, Vid_cpue_cv, w_I_Vid);
  //fval+=f_Vid_cpue;
  //fval_data+=f_Vid_cpue;

  f_cL_cpue=0.0;
  f_cL_cpue=lk_lognormal(pred_cL_cpue, obs_cL_cpue, cL_cpue_cv, w_I_cL);
  fval+=f_cL_cpue;
  fval_data+=f_cL_cpue;
  
  f_HB_cpue=0.0;
  f_HB_cpue=lk_lognormal(pred_HB_cpue, obs_HB_cpue, HB_cpue_cv, w_I_HB);
  fval+=f_HB_cpue;
  fval_data+=f_HB_cpue;
  
  //f_HBD_cpue=0.0;
  //f_HBD_cpue=lk_lognormal(pred_HBD_cpue, obs_HBD_cpue, HBD_cpue_cv, w_I_HBD);
  //fval+=f_HBD_cpue;
  //fval_data+=f_HBD_cpue;

////---Landings-------------------------------
 
  //f_cL_L in 1000 lb ww
  f_cL_L=lk_lognormal(pred_cL_L_klb(styr_cL_L,endyr_cL_L), obs_cL_L(styr_cL_L,endyr_cL_L), 
                      cL_L_cv(styr_cL_L,endyr_cL_L), w_L);
  fval+=f_cL_L;
  fval_data+=f_cL_L;  
  //f_cP_L in 1000 lb ww
  f_cP_L=lk_lognormal(pred_cP_L_klb(styr_cP_L,endyr_cP_L), obs_cP_L(styr_cP_L,endyr_cP_L), 
                      cP_L_cv(styr_cP_L,endyr_cP_L), w_L);
  fval+=f_cP_L;
  fval_data+=f_cP_L; 
  //f_cT_L in 1000 lb ww
  f_cT_L=lk_lognormal(pred_cT_L_klb(styr_cT_L,endyr_cT_L), obs_cT_L(styr_cT_L,endyr_cT_L), 
                      cT_L_cv(styr_cT_L,endyr_cT_L), w_L);
  fval+=f_cT_L;
  fval_data+=f_cT_L;
  //f_HB_L in 1000 lb
  f_HB_L=lk_lognormal(pred_HB_L_klb(styr_HB_L,endyr_HB_L), obs_HB_L(styr_HB_L,endyr_HB_L),
                      HB_L_cv(styr_HB_L,endyr_HB_L), w_L);
  fval+=f_HB_L;
  fval_data+=f_HB_L;

  //f_mrip_L in 1000 lb
  f_mrip_L=lk_lognormal(pred_mrip_L_klb(styr_mrip_L,endyr_mrip_L), obs_mrip_L(styr_mrip_L,endyr_mrip_L), 
                       mrip_L_cv(styr_mrip_L,endyr_mrip_L), w_L);
  fval+=f_mrip_L;
  fval_data+=f_mrip_L;


//---Discards-------------------------------

  //f_comm_D in 1000 fish
  f_comm_D=lk_lognormal(pred_comm_D_knum(styr_cL_D,endyr_cL_D), obs_comm_D(styr_cL_D,endyr_cL_D), 
                      comm_D_cv(styr_cL_D,endyr_cL_D), w_D);
  fval+=f_comm_D;
  fval_data+=f_comm_D;
  
  
  //f_HB_D in 1000 fish
  f_HB_D=lk_lognormal(pred_HB_D_knum(styr_HB_D,endyr_HB_D), obs_HB_D(styr_HB_D,endyr_HB_D), 
                         HB_D_cv(styr_HB_D,endyr_HB_D), w_D);
  fval+=f_HB_D;
  fval_data+=f_HB_D;

  //f_mrip_D in 1000 fish
  f_mrip_D=lk_lognormal(pred_mrip_D_knum(styr_mrip_D,endyr_mrip_D), obs_mrip_D(styr_mrip_D,endyr_mrip_D), 
                       mrip_D_cv(styr_mrip_D,endyr_mrip_D), w_D);
  fval+=f_mrip_D;
  fval_data+=f_mrip_D;


//---Length comps-------------------------------

  //f_Mbft_lenc
  f_Mbft_lenc=lk_dirichlet_multinomial(nsamp_Mbft_lenc, pred_Mbft_lenc, obs_Mbft_lenc, nyr_Mbft_lenc, double(nlenbins), minSS_lenc, log_dm_Mbft_lc); //w_lc_Mbft);
  fval+=f_Mbft_lenc;
  fval_data+=f_Mbft_lenc;

  //f_cL_lenc
  f_cL_lenc=lk_dirichlet_multinomial(nsamp_cL_lenc, pred_cL_lenc, obs_cL_lenc, nyr_cL_lenc, double(nlenbins), minSS_lenc, log_dm_cL_lc); //w_lc_cL);
  fval+=f_cL_lenc;
  fval_data+=f_cL_lenc;
  
  //f_cP_lenc
  f_cP_lenc=lk_dirichlet_multinomial(nsamp_cP_lenc, pred_cP_lenc, obs_cP_lenc, nyr_cP_lenc, double(nlenbins), minSS_lenc, log_dm_cP_lc); //w_lc_cP);
  fval+=f_cP_lenc;
  fval_data+=f_cP_lenc;
  
  //f_HB_lenc
  f_HB_lenc=lk_dirichlet_multinomial(nsamp_HB_lenc, pred_HB_lenc, obs_HB_lenc, nyr_HB_lenc, double(nlenbins), minSS_lenc, log_dm_HB_lc); //w_lc_HB);
  fval+=f_HB_lenc;
  fval_data+=f_HB_lenc;  
  
  //f_HB_D_lenc
  f_HB_D_lenc=lk_dirichlet_multinomial(nsamp_HB_D_lenc, pred_HB_D_lenc, obs_HB_D_lenc, nyr_HB_D_lenc, double(nlenbins), minSS_lenc, log_dm_HB_D_lc); //w_lc_HB_D);
  fval+=f_HB_D_lenc;
  fval_data+=f_HB_D_lenc;
    
  //f_mrip_lenc
  f_mrip_lenc=lk_dirichlet_multinomial(nsamp_mrip_lenc, pred_mrip_lenc, obs_mrip_lenc, nyr_mrip_lenc, double(nlenbins), minSS_lenc, log_dm_mrip_lc); //w_lc_mrip);
  fval+=f_mrip_lenc;
  fval_data+=f_mrip_lenc;


//---Age comps-------------------------------

  //f_Mbft_agec
  f_Mbft_agec=lk_dirichlet_multinomial(nsamp_Mbft_agec, pred_Mbft_agec, obs_Mbft_agec, nyr_Mbft_agec, double(nages), minSS_agec, log_dm_Mbft_ac); //w_ac_Mbft);
  fval+=f_Mbft_agec;
  fval_data+=f_Mbft_agec;

  //f_Mcvt_agec
  f_Mcvt_agec=lk_dirichlet_multinomial(nsamp_Mcvt_agec, pred_Mcvt_agec, obs_Mcvt_agec, nyr_Mcvt_agec, double(nages), minSS_agec, log_dm_Mcvt_ac); //w_ac_Mcvt);
  fval+=f_Mcvt_agec;
  fval_data+=f_Mcvt_agec;

  //f_cL_agec
  f_cL_agec=lk_dirichlet_multinomial(nsamp_cL_agec, pred_cL_agec, obs_cL_agec, nyr_cL_agec, double(nages), minSS_agec, log_dm_cL_ac); //w_ac_cL);
  fval+=f_cL_agec;
  fval_data+=f_cL_agec;
  
  //f_cP_agec
  f_cP_agec=lk_dirichlet_multinomial(nsamp_cP_agec, pred_cP_agec, obs_cP_agec, nyr_cP_agec, double(nages), minSS_agec, log_dm_cP_ac); //w_ac_cP);
  fval+=f_cP_agec;
  fval_data+=f_cP_agec;
    
  //f_HB_agec
  f_HB_agec=lk_dirichlet_multinomial(nsamp_HB_agec, pred_HB_agec, obs_HB_agec, nyr_HB_agec, double(nages), minSS_agec, log_dm_HB_ac); //w_ac_HB);
  fval+=f_HB_agec;
  fval_data+=f_HB_agec;  
  
  ////f_mrip_agec
  //f_mrip_agec=lk_dirichlet_multinomial(nsamp_mrip_agec, pred_mrip_agec, obs_mrip_agec, nyr_mrip_agec, double(nages), minSS_agec, log_dm_mrip_ac); //w_ac_mrip);
  //fval+=f_mrip_agec;
  //fval_data+=f_mrip_agec;    
  ////f_mrip_agec=lk_multinomial(nsamp_mrip_agec, pred_mrip_agec, obs_mrip_agec, nyr_mrip_agec, minSS_agec, w_ac_mrip);
  ////fval+=f_mrip_agec;
  ////fval_data+=f_mrip_agec;  
//-----------Constraints and penalties--------------------------------

  f_rec_dev=0.0;
  rec_logL_add=nyrs_rec*log(rec_sigma);
  f_rec_dev=(square(log_rec_dev(styr_rec_dev) + rec_sigma_sqd2)/(2.0*rec_sigma_sq));
  for(iyear=(styr_rec_dev+1); iyear<=endyr; iyear++)
  {f_rec_dev+=(square(log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1) + rec_sigma_sqd2)/
               (2.0*rec_sigma_sq));}
  f_rec_dev+=rec_logL_add;            
  fval+=w_rec*f_rec_dev;

  f_rec_dev_early=0.0; //possible extra constraint on early rec deviations
  if (w_rec_early>0.0)
    { if (styr_rec_dev<endyr_rec_phase1)
        {  
          f_rec_dev_early=(square(log_rec_dev(styr_rec_dev) + rec_sigma_sq/2.0)/(2.0*rec_sigma_sq)) + rec_logL_add;
          for(iyear=(styr_rec_dev+1); iyear<=endyr_rec_phase1; iyear++)
          {f_rec_dev_early+=(square(log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1) + rec_sigma_sqd2)/
                            (2.0*rec_sigma_sq)) + rec_logL_add;}
        }
  fval+=w_rec_early*f_rec_dev_early;
  }
  
  f_rec_dev_end=0.0; //possible extra constraint on ending rec deviations
  if (w_rec_end>0.0)
  { if (endyr_rec_phase2<endyr)
        {  
          for(iyear=(endyr_rec_phase2+1); iyear<=endyr; iyear++)
          {f_rec_dev_end+=(square(log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1) + rec_sigma_sqd2)/
                            (2.0*rec_sigma_sq)) + rec_logL_add;}
        }
      fval+=w_rec_end*f_rec_dev_end;
   }
//	cout<< "log rec devs" << log_rec_dev << endl;

//  f_rec_dev_early=0.0; //possible extra constraint on early rec deviations
//  if (styr_rec_dev<endyr_rec_phase1)
//    {  
//      f_rec_dev_early=pow(log_rec_dev(styr_rec_dev),2);
//      for(iyear=(styr_rec_dev+1); iyear<=endyr_rec_phase1; iyear++)
//      {f_rec_dev_early+=pow((log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1)),2);}
//    }
//  fval+=w_rec_early*f_rec_dev_early;

//  f_rec_dev_end=0.0; //possible extra constraint on ending rec deviations
//  if (endyr_rec_phase2<endyr)
//    {  
//      for(iyear=(endyr_rec_phase2+1); iyear<=endyr; iyear++)
//      {f_rec_dev_end+=pow((log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1)),2);}
//    }
//  fval+=w_rec_end*f_rec_dev_end;

//  f_Ftune=0.0; 
//  if (set_Ftune>0.0 && !last_phase()) {f_Ftune=square(Fapex(set_Ftune_yr)-set_Ftune);}
//  fval+=w_Ftune*f_Ftune;
  
//  //code below contingent on four phases
//  f_fullF_constraint=0.0;
//  if (!last_phase())
//  {for (iyear=styr; iyear<=endyr; iyear++)
//       {if (Fapex(iyear)>3.0){f_fullF_constraint+=mfexp(Fapex(iyear)-3.0);}}
//   if (current_phase()==1) {w_fullF=set_w_fullF;}
//   if (current_phase()==2) {w_fullF=set_w_fullF/10.0;}  
//   if (current_phase()==3) {w_fullF=set_w_fullF/100.0;}
//  }

//  fval+=w_fullF*f_fullF_constraint;
//
//  f_fullF_constraint=0.0;
//  for (iyear=styr; iyear<=endyr; iyear++)
//       {if (Fapex(iyear)>3.0){f_fullF_constraint+=mfexp(Fapex(iyear)-3.0);}}      
//  fval+=w_fullF*f_fullF_constraint;  
    
//  f_cvlen_diff_constraint=0.0;
//    f_cvlen_diff_constraint=norm2(first_difference(log_len_cv_dev));
//  fval+=w_cvlen_diff*f_cvlen_diff_constraint;
//  
//  f_cvlen_dev_constraint=0.0;
//    f_cvlen_dev_constraint=norm2(log_len_cv_dev);  
//  fval+=w_cvlen_dev*f_cvlen_dev_constraint;
  
  
  //applies if initial age structure is estimated
  fval+=norm2(log_Nage_dev); 

  //Random walk components of fishery dependent indices: these components equal zero if RW turned off
  f_cL_RW_cpue=0.0;
  for (iyear=styr_cL_cpue; iyear<endyr_cL_cpue; iyear++)
      {f_cL_RW_cpue+=square(q_RW_log_dev_cL(iyear))/(2.0*set_q_RW_cL_var);}
  fval+=f_cL_RW_cpue;  
  
  f_HB_RW_cpue=0.0;
  for (iyear=styr_HB_cpue; iyear<endyr_HB_cpue; iyear++)
      {f_HB_RW_cpue+=square(q_RW_log_dev_HB(iyear))/(2.0*set_q_RW_HB_var);}      
  fval+=f_HB_RW_cpue;   
  
  //f_HBD_RW_cpue=0.0;
  //for (iyear=styr_HBD_cpue; iyear<endyr_HBD_cpue; iyear++)
  //    {f_HBD_RW_cpue+=square(q_RW_log_dev_HBD(iyear))/(2.0*set_q_RW_HBD_var);}      
  //fval+=f_HBD_RW_cpue;  
  
  
//---Priors---------------------------------------------------
//neg_log_prior arguments: estimate, prior, variance, pdf type
//Variance input as a negative value is considered to be CV in arithmetic space (CV=-1 implies loose prior) 
//pdf type(1=none, 2=lognormal, 3=normal, 4=beta)

  f_priors=0.0; 
  f_priors+=neg_log_prior(steep, set_steep(5), set_steep(6), set_steep(7)); 
//	cout<<f_priors<<endl;
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));
  f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  //f_priors+=neg_log_prior(q_DD_beta, set_q_DD_beta, square(set_q_DD_beta_se), 3);	
  //f_priors+=neg_log_prior(q_rate, set_q_rate, dzero+square(set_q_rate), 3);  
  f_priors+=neg_log_prior(F_init_ratio, set_F_init_ratio, -1.0 , 3);
  //f_priors+=neg_log_prior(M_constant, set_M_constant, square(set_M_constant_se), 3);	

  //f_priors+=neg_log_prior(Linf,set_Linf,square(set_Linf_se),3);
  //f_priors+=neg_log_prior(K,set_K,square(set_K_se),3);
  //f_priors+=neg_log_prior(t0,set_K,square(set_t0_se),3);
  f_priors+=neg_log_prior(len_cv_val,set_len_cv(5),set_len_cv(6),set_len_cv(7));

  f_priors+=neg_log_prior(selpar_A50_Mbft, set_selpar_A50_Mbft(5),set_selpar_A50_Mbft(6),set_selpar_A50_Mbft(7));  
  f_priors+=neg_log_prior(selpar_slope_Mbft, set_selpar_slope_Mbft(5),set_selpar_slope_Mbft(6),set_selpar_slope_Mbft(7));                  

  f_priors+=neg_log_prior(selpar_A50_Mcvt, set_selpar_A50_Mcvt(5),set_selpar_A50_Mcvt(6),set_selpar_A50_Mcvt(7));  
  f_priors+=neg_log_prior(selpar_slope_Mcvt, set_selpar_slope_Mcvt(5),set_selpar_slope_Mcvt(6),set_selpar_slope_Mcvt(7));                  

  f_priors+=neg_log_prior(selpar_A50_cL2, set_selpar_A50_cL2(5),set_selpar_A50_cL2(6),set_selpar_A50_cL2(7));  
  f_priors+=neg_log_prior(selpar_slope_cL2, set_selpar_slope_cL2(5),set_selpar_slope_cL2(6),set_selpar_slope_cL2(7));                  
  f_priors+=neg_log_prior(selpar_A50_cL3, set_selpar_A50_cL3(5),set_selpar_A50_cL3(6),set_selpar_A50_cL3(7));  
  f_priors+=neg_log_prior(selpar_slope_cL3, set_selpar_slope_cL3(5),set_selpar_slope_cL3(6),set_selpar_slope_cL3(7)); 
  f_priors+=neg_log_prior(selpar_A50_cL4, set_selpar_A50_cL4(5),set_selpar_A50_cL4(6),set_selpar_A50_cL4(7));  
  f_priors+=neg_log_prior(selpar_slope_cL4, set_selpar_slope_cL4(5),set_selpar_slope_cL4(6),set_selpar_slope_cL4(7));   

  f_priors+=neg_log_prior(selpar_A50_cP2, set_selpar_A50_cP2(5),set_selpar_A50_cP2(6),set_selpar_A50_cP2(7));  
  f_priors+=neg_log_prior(selpar_slope_cP2, set_selpar_slope_cP2(5),set_selpar_slope_cP2(6),set_selpar_slope_cP2(7));                  
  f_priors+=neg_log_prior(selpar_A50_cP3, set_selpar_A50_cP3(5),set_selpar_A50_cP3(6),set_selpar_A50_cP3(7));  
  f_priors+=neg_log_prior(selpar_slope_cP3, set_selpar_slope_cP3(5),set_selpar_slope_cP3(6),set_selpar_slope_cP3(7)); 
  f_priors+=neg_log_prior(selpar_A50_cP4, set_selpar_A50_cP4(5),set_selpar_A50_cP4(6),set_selpar_A50_cP4(7));  
  f_priors+=neg_log_prior(selpar_slope_cP4, set_selpar_slope_cP4(5),set_selpar_slope_cP4(6),set_selpar_slope_cP4(7));   

  f_priors+=neg_log_prior(selpar_A50_HB1, set_selpar_A50_HB1(5),set_selpar_A50_HB1(6),set_selpar_A50_HB1(7));  
  f_priors+=neg_log_prior(selpar_slope_HB1, set_selpar_slope_HB1(5),set_selpar_slope_HB1(6),set_selpar_slope_HB1(7));                  
  f_priors+=neg_log_prior(selpar_A50_HB2, set_selpar_A50_HB2(5),set_selpar_A50_HB2(6),set_selpar_A50_HB2(7));  
  f_priors+=neg_log_prior(selpar_slope_HB2, set_selpar_slope_HB2(5),set_selpar_slope_HB2(6),set_selpar_slope_HB2(7));                  
  f_priors+=neg_log_prior(selpar_A50_HB3, set_selpar_A50_HB3(5),set_selpar_A50_HB3(6),set_selpar_A50_HB3(7));  
  f_priors+=neg_log_prior(selpar_slope_HB3, set_selpar_slope_HB3(5),set_selpar_slope_HB3(6),set_selpar_slope_HB3(7));      
  f_priors+=neg_log_prior(selpar_A50_HB4, set_selpar_A50_HB4(5),set_selpar_A50_HB4(6),set_selpar_A50_HB4(7)); 
  f_priors+=neg_log_prior(selpar_slope_HB4, set_selpar_slope_HB4(5),set_selpar_slope_HB4(6),set_selpar_slope_HB4(7)); 
  f_priors+=neg_log_prior(selpar_A50_HB5, set_selpar_A50_HB5(5),set_selpar_A50_HB5(6),set_selpar_A50_HB5(7));   
  f_priors+=neg_log_prior(selpar_slope_HB5, set_selpar_slope_HB5(5),set_selpar_slope_HB5(6),set_selpar_slope_HB5(7));   
  
  f_priors+=neg_log_prior(selpar_A50_mrip1, set_selpar_A50_mrip1(5),set_selpar_A50_mrip1(6),set_selpar_A50_mrip1(7));  
  f_priors+=neg_log_prior(selpar_slope_mrip1, set_selpar_slope_mrip1(5),set_selpar_slope_mrip1(6),set_selpar_slope_mrip1(7));                  
  f_priors+=neg_log_prior(selpar_A50_mrip2, set_selpar_A50_mrip2(5),set_selpar_A50_mrip2(6),set_selpar_A50_mrip2(7));  
  f_priors+=neg_log_prior(selpar_slope_mrip2, set_selpar_slope_mrip2(5),set_selpar_slope_mrip2(6),set_selpar_slope_mrip2(7));                  
  f_priors+=neg_log_prior(selpar_A50_mrip3, set_selpar_A50_mrip3(5),set_selpar_A50_mrip3(6),set_selpar_A50_mrip3(7));  
  f_priors+=neg_log_prior(selpar_slope_mrip3, set_selpar_slope_mrip3(5),set_selpar_slope_mrip3(6),set_selpar_slope_mrip3(7));      
  f_priors+=neg_log_prior(selpar_A50_mrip4, set_selpar_A50_mrip4(5),set_selpar_A50_mrip4(6),set_selpar_A50_mrip4(7)); 
  f_priors+=neg_log_prior(selpar_slope_mrip4, set_selpar_slope_mrip4(5),set_selpar_slope_mrip4(6),set_selpar_slope_mrip4(7)); 
  f_priors+=neg_log_prior(selpar_A50_mrip5, set_selpar_A50_mrip5(5),set_selpar_A50_mrip5(6),set_selpar_A50_mrip5(7));   
  f_priors+=neg_log_prior(selpar_slope_mrip5, set_selpar_slope_mrip5(5),set_selpar_slope_mrip5(6),set_selpar_slope_mrip5(7));  

  f_priors+=neg_log_prior(selpar_Age0_HB_D_logit, set_selpar_Age0_HB_D_logit(5),set_selpar_Age0_HB_D_logit(6),set_selpar_Age0_HB_D_logit(7));
  f_priors+=neg_log_prior(selpar_Age1_HB_D_logit, set_selpar_Age1_HB_D_logit(5),set_selpar_Age1_HB_D_logit(6),set_selpar_Age1_HB_D_logit(7));
  f_priors+=neg_log_prior(selpar_Age2_HB_D_logit, set_selpar_Age2_HB_D_logit(5),set_selpar_Age2_HB_D_logit(6), set_selpar_Age2_HB_D_logit(7));
    
  f_priors+=neg_log_prior(selpar_A50_HBD4, set_selpar_A50_HBD4(5),set_selpar_A50_HBD4(6),set_selpar_A50_HBD4(7)); 
  f_priors+=neg_log_prior(selpar_slope_HBD4, set_selpar_slope_HBD4(5),set_selpar_slope_HBD4(6),set_selpar_slope_HBD4(7)); 
  f_priors+=neg_log_prior(selpar_A502_HBD4, set_selpar_A502_HBD4(5),set_selpar_A502_HBD4(6),set_selpar_A502_HBD4(7)); f_priors+=neg_log_prior(selpar_A50_HBD5, set_selpar_A50_HBD5(5),set_selpar_A50_HBD5(6),set_selpar_A50_HBD5(7));   
  f_priors+=neg_log_prior(selpar_slope_HBD5, set_selpar_slope_HBD5(5),set_selpar_slope_HBD5(6),set_selpar_slope_HBD5(7)); 
  f_priors+=neg_log_prior(selpar_A502_HBD5, set_selpar_A502_HBD5(5),set_selpar_A502_HBD5(6),set_selpar_A502_HBD5(7));
  
  fval+=f_priors;

  //cout << "fval = " << fval << "  fval_data = " << fval_data << endl;


//----------------------------------------------------------------------------------
//Logistic function: 2 parameters
FUNCTION dvar_vector logistic(const dvar_vector& ages, const dvariable& A50, const dvariable& slope)
  //ages=vector of ages, A50=age at 50% selectivity, slope=rate of increase
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1./(1.+mfexp(-1.*slope*(ages-A50))); //logistic;  
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

//-----------------------------------------------------------------------------------
//Logistic function: 4 parameters
FUNCTION dvar_vector logistic_double(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2)
  //ages=vector of ages, A50=age at 50% selectivity, slope=rate of increase, A502=age at 50% decrease additive to A501, slope2=slope of decrease
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=elem_prod( (1./(1.+mfexp(-1.*slope1*(ages-A501)))),(1.-(1./(1.+mfexp(-1.*slope2*(ages-(A501+A502)))))) );     
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
//-----------------------------------------------------------------------------------
  //Logistic-exponential: 4 parameters (but 1 is fixed)
FUNCTION dvar_vector logistic_exponential(const dvar_vector& ages, const dvariable& A50, const dvariable& slope, const dvariable& sigma, const dvariable& joint)
  //ages=vector of ages, A50=age at 50% sel (ascending limb), slope=rate of increase, sigma=controls rate of descent (descending)                               
  //joint=age to join curves                                                                                                                                    
  RETURN_ARRAYS_INCREMENT();                                                                                                                                  
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());                                                                                                         
  Sel_Tmp=1.0;                                                                                                                                                  
  for (iage=1; iage<=nages; iage++)                                                                                                                             
  {                                                                                                                                                             
   if (ages(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope*(ages(iage)-A50)));}                                                                             
   if (ages(iage)>joint){Sel_Tmp(iage)=mfexp(-1.*square((ages(iage)-joint)/sigma));}                                                                            
  }                                                                                                                                                             
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);                                                                                                                                 
  RETURN_ARRAYS_DECREMENT();                                                                                                                                    
  return Sel_Tmp;   
  

//-----------------------------------------------------------------------------------
//Jointed logistic function: 6 parameters (increasing and decreasing logistics joined at peak selectivity)
FUNCTION dvar_vector logistic_joint(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2, const dvariable& satval, const dvariable& joint)
  //ages=vector of ages, A501=age at 50% sel (ascending limb), slope1=rate of increase,A502=age at 50% sel (descending), slope1=rate of increase (ascending), 
  //satval=saturation value of descending limb, joint=location in age vector to join curves (may equal age or age + 1 if age-0 is included)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1.0; 
  for (iage=1; iage<=nages; iage++)
  {
   if (double(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope1*(ages(iage)-A501)));}  
   if (double(iage)>joint){Sel_Tmp(iage)=1.0-(1.0-satval)/(1.+mfexp(-1.*slope2*(ages(iage)-A502)));}  
  }  
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

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
FUNCTION dvariable SR_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& SSB)
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, SSB=spawning biomass
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  Recruits_Tmp=((0.8*R0*h*SSB)/(0.2*R0*spr_F0*(1.0-h)+(h-0.2)*SSB));
  RETURN_ARRAYS_DECREMENT();
  return Recruits_Tmp;

//-----------------------------------------------------------------------------------
//compute multinomial effective sample size for a single yr
FUNCTION dvariable multinom_eff_N(const dvar_vector& pred_comp, const dvar_vector& obs_comp)
  //pred_comp=vector of predicted comps, obscomp=vector of observed comps
  dvariable EffN_Tmp; dvariable numer; dvariable denom;
  RETURN_ARRAYS_INCREMENT();
  numer=sum( elem_prod(pred_comp,(1.0-pred_comp)) );
  denom=sum( square(obs_comp-pred_comp) );
  if (denom>0.0) {EffN_Tmp=numer/denom;}
  else {EffN_Tmp=-missing;}                            
  RETURN_ARRAYS_DECREMENT();
  return EffN_Tmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: lognormal
FUNCTION dvariable lk_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
  //pred=vector of predicted vals, obs=vector of observed vals, cv=vector of CVs in arithmetic space, wgt_dat=constant scaling of CVs
  //dzero is small value to avoid log(0) during search
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvar_vector var(cv.indexmin(),cv.indexmax()); //variance in log space
  var=log(1.0+square(cv/wgt_dat));   // convert cv in arithmetic space to variance in log space
  LkvalTmp=sum(0.5*elem_div(square(log(elem_div((pred+dzero),(obs+dzero)))),var) );
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: multinomial
FUNCTION dvariable lk_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const double& minSS, const dvariable& wgt_dat)
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  LkvalTmp=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {LkvalTmp-=wgt_dat*nsamp(ii)*sum(elem_prod((obs_comp(ii)+dzero),
               log(elem_div((pred_comp(ii)+dzero), (obs_comp(ii)+dzero)))));
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

  //-----------------------------------------------------------------------------------
//Likelihood contribution: Dirichlet-multinomial
FUNCTION dvariable lk_dirichlet_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& log_dir_par)
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.00001;
  LkvalTmp=0.0; 
  dvar_vector nsamp_adjust=nsamp*mfexp(log_dir_par);
   //dvar_vector nsamp_adjust=mfexp(log_dir_par);
  for (int ii=1; ii<=ncomp; ii++)
  {
	if (nsamp(ii)>=minSS)
    {
		LkvalTmp-=gammln(nsamp_adjust(ii))-gammln(nsamp(ii)+nsamp_adjust(ii));
		LkvalTmp-=sum(gammln(nsamp(ii)*obs_comp(ii)+nsamp_adjust(ii)*pred_comp(ii)+small_number));
		LkvalTmp+=sum(gammln(nsamp_adjust(ii)*pred_comp(ii)+small_number));
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

// //Likelihood contribution: Dirichlet-multinomial
// FUNCTION dvariable lk_dirichlet_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& log_dir_par)
  // //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold, wgt_dat=scaling of N's
  // RETURN_ARRAYS_INCREMENT();
  // dvariable LkvalTmp;
  // LkvalTmp=0.0; 
  // dvar_vector nsamp_adjust=nsamp*mfexp(log_dir_par);
  // //dvar_vector nsamp_adjust=mfexp(log_dir_par);
  // for (int ii=1; ii<=ncomp; ii++)
  // {
	// if (nsamp(ii)>=minSS)
    // {
		// LkvalTmp-=gammln(nsamp_adjust(ii))-gammln(nsamp(ii)+nsamp_adjust(ii));
		// LkvalTmp-=sum(gammln(nsamp(ii)*obs_comp(ii)+nsamp_adjust(ii)*pred_comp(ii)));
        // LkvalTmp+=sum(gammln(nsamp_adjust(ii)*pred_comp(ii)));		
    // }
  // }  
  // RETURN_ARRAYS_DECREMENT();
  // return LkvalTmp;
  

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
 
//-----------------------------------------------------------------------------------
//SDNR: age comp likelihood (assumes fits are done with the robust multinomial function)
FUNCTION dvariable sdnr_multinomial(const double& ncomp, const dvar_vector& ages, const dvar_vector& nsamp, 
                                    const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const dvariable& wgt_dat)
  //ncomp=number of years of data, ages=vector of ages, nsamp=vector of N's, 
  //pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, wgt_dat=likelihood weight for data source
  RETURN_ARRAYS_INCREMENT();
  dvariable SdnrTmp;
  dvar_vector o(1,ncomp);  
  dvar_vector p(1,ncomp);  
  dvar_vector ose(1,ncomp);  
  dvar_vector res(1,ncomp);
  SdnrTmp=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {
    o(ii)=sum(elem_prod(ages,obs_comp(ii)));
    p(ii)=sum(elem_prod(ages,pred_comp(ii)));
    ose(ii)=sqrt((sum(elem_prod(square(ages),pred_comp(ii)))-square(p(ii)))/(nsamp(ii)*wgt_dat));
  }
  res=elem_div((o-p),ose); 
  SdnrTmp=sqrt(sum(square(res-(sum(res)/ncomp))/(ncomp-1.0))); 
  RETURN_ARRAYS_DECREMENT();
  return SdnrTmp;

//-----------------------------------------------------------------------------------
//SDNR: lognormal likelihood
FUNCTION dvariable sdnr_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
  //nyr=number of years of data, pred=vector of predicted data, obs=vector of observed data, cv=vector of cv's, wgt_dat=likelihood weight for data source
  RETURN_ARRAYS_INCREMENT();
  dvariable SdnrTmp;
  dvariable small_number=0.00001;
  dvariable n;
  dvar_vector res(cv.indexmin(),cv.indexmax());
  SdnrTmp=0.0;
  res=elem_div(log(elem_div(obs+small_number,pred+small_number)),sqrt(log(1+square(cv/wgt_dat))));
  n=cv.indexmax()-cv.indexmin()+1;
  SdnrTmp=sqrt(sum(square(res-(sum(res)/n))/(n-1.0))); 
  RETURN_ARRAYS_DECREMENT();
  return SdnrTmp; 

//-----------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------
REPORT_SECTION

  if (last_phase())  
  {

      cout<<"start report"<<endl;
      get_weighted_current();
      cout<<"got weighted"<<endl;
      get_msy();
      cout<<"got msy"<<endl;
      get_miscellaneous_stuff();
      cout<<"got misc stuff"<<endl;
      get_per_recruit_stuff();
      cout<<"got per recruit"<<endl;  
      get_effective_sample_sizes();
	  cout<<"got effective sample sizes"<<endl; 
      //get_projection();
	  //cout<<"got projection"<<endl;
	  
      grad_max=objective_function_value::pobjfun->gmax;
      time(&finish);
	  elapsed_time=difftime(finish,start);
	  hour=long(elapsed_time)/3600;
	  minute=long(elapsed_time)%3600/60;
	  second=(long(elapsed_time)%3600)%60;
	  cout<<endl<<endl<<"*******************************************"<<endl;
	  cout<<"--Start time: "<<ctime(&start)<<endl;
	  cout<<"--Finish time: "<<ctime(&finish)<<endl;
	  cout<<"--Runtime: ";
	  cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	  cout << "--TotalLikelihood: " << fval << endl;
      cout<<"--Final gradient: "<<objective_function_value::pobjfun->gmax << endl;
	  cout<<"*******************************************"<<endl;
      cout <<endl;     
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;
      cout << "BC Fmsy=" << F_msy_out<< "   BC SSBmsy=" << SSB_msy_out <<endl;
      cout <<"F status="<<FdF_msy_end<<endl;
      cout <<"Pop status="<<SdSSB_msy_end<<endl;
      cout << "h="<<steep<<"   R0="<<R0<<endl;
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;  
       
      report << "TotalLikelihood " << fval << endl;
      report << "N" << endl;
      report << N<<endl;      
      report << "SSB" << endl;       
      report << SSB << endl;
	    

	  //sdnr_lc_Mbft=sdnr_multinomial(nyr_Mbft_lenc, lenbins, nsamp_Mbft_lenc, pred_Mbft_lenc, obs_Mbft_lenc, w_lc_Mbft); 
      //sdnr_lc_Mcvt=sdnr_multinomial(nyr_Mcvt_lenc, lenbins, nsamp_Mcvt_lenc, pred_Mcvt_lenc, obs_Mcvt_lenc, w_lc_Mcvt); 
      //sdnr_lc_cL=sdnr_multinomial(nyr_cL_lenc, lenbins, nsamp_cL_lenc, pred_cL_lenc, obs_cL_lenc, w_lc_cL);  
      //sdnr_lc_cP=sdnr_multinomial(nyr_cP_lenc, lenbins, nsamp_cP_lenc, pred_cP_lenc, obs_cP_lenc, w_lc_cP);  	  
      //sdnr_lc_HB=sdnr_multinomial(nyr_HB_lenc, lenbins, nsamp_HB_lenc, pred_HB_lenc, obs_HB_lenc, w_lc_HB); 
	  //sdnr_lc_HB_D=sdnr_multinomial(nyr_HB_D_lenc, lenbins, nsamp_HB_D_lenc, pred_HB_D_lenc, obs_HB_D_lenc, w_lc_HB_D); 
	  //sdnr_lc_mrip=sdnr_multinomial(nyr_mrip_lenc, lenbins, nsamp_mrip_lenc, pred_mrip_lenc, obs_mrip_lenc, w_lc_mrip); 
       
      //sdnr_ac_Mbft=sdnr_multinomial(nyr_Mbft_agec, agebins_agec, nsamp_Mbft_agec, pred_Mbft_agec, obs_Mbft_agec, w_ac_Mbft);  
      //sdnr_ac_Mcvt=sdnr_multinomial(nyr_Mcvt_agec, agebins_agec, nsamp_Mcvt_agec, pred_Mcvt_agec, obs_Mcvt_agec, w_ac_Mcvt);  
      //sdnr_ac_cL=sdnr_multinomial(nyr_cL_agec, agebins_agec, nsamp_cL_agec, pred_cL_agec, obs_cL_agec, w_ac_cL);  
      //sdnr_ac_cP=sdnr_multinomial(nyr_cP_agec, agebins_agec, nsamp_cP_agec, pred_cP_agec, obs_cP_agec, w_ac_cP);  
      //sdnr_ac_HB=sdnr_multinomial(nyr_HB_agec, agebins_agec, nsamp_HB_agec, pred_HB_agec, obs_HB_agec, w_ac_HB);  
	  //sdnr_ac_mrip=sdnr_multinomial(nyr_mrip_agec, agebins_agec, nsamp_mrip_agec, pred_mrip_agec, obs_mrip_agec, w_ac_mrip);  
      
      sdnr_I_Mbft=sdnr_lognormal(pred_Mbft_cpue, obs_Mbft_cpue, Mbft_cpue_cv, w_I_Mbft);  
      sdnr_I_Mcvt=sdnr_lognormal(pred_Mcvt_cpue, obs_Mcvt_cpue, Mcvt_cpue_cv, w_I_Mcvt);  
      sdnr_I_cL=sdnr_lognormal(pred_cL_cpue, obs_cL_cpue, cL_cpue_cv, w_I_cL);
      sdnr_I_HB=sdnr_lognormal(pred_HB_cpue, obs_HB_cpue, HB_cpue_cv, w_I_HB);
      //sdnr_I_Vid=sdnr_lognormal(pred_Vid_cpue, obs_Vid_cpue, Vid_cpue_cv, w_I_Vid);
            
  
	   Linf_out(8)=Linf; Linf_out(1,7)=set_Linf; 
       K_out(8)=K; K_out(1,7)=set_K;
       t0_out(8)=t0; t0_out(1,7)=set_t0;
       len_cv_val_out(8)=len_cv_val; len_cv_val_out(1,7)=set_len_cv;
	   	   
       log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
       M_constant_out(8)=M_constant; M_constant_out(1,7)=set_M_constant;
       steep_out(8)=steep; steep_out(1,7)=set_steep;
       rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
       R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;
	   
	   log_dm_Mbft_lc_out(8)=log_dm_Mbft_lc; log_dm_Mbft_lc_out(1,7)=set_log_dm_Mbft_lc;
	   log_dm_Mcvt_lc_out(8)=log_dm_Mcvt_lc; log_dm_Mcvt_lc_out(1,7)=set_log_dm_Mcvt_lc;
	   log_dm_cL_lc_out(8)=log_dm_cL_lc; log_dm_cL_lc_out(1,7)=set_log_dm_cL_lc;
	   log_dm_cP_lc_out(8)=log_dm_cP_lc; log_dm_cP_lc_out(1,7)=set_log_dm_cP_lc;
	   log_dm_HB_lc_out(8)=log_dm_HB_lc; log_dm_HB_lc_out(1,7)=set_log_dm_HB_lc;
	   log_dm_HB_D_lc_out(8)=log_dm_HB_D_lc; log_dm_HB_D_lc_out(1,7)=set_log_dm_HB_D_lc;
	   log_dm_mrip_lc_out(8)=log_dm_mrip_lc; log_dm_mrip_lc_out(1,7)=set_log_dm_mrip_lc;
	   
	   log_dm_Mbft_ac_out(8)=log_dm_Mbft_ac; log_dm_Mbft_ac_out(1,7)=set_log_dm_Mbft_ac;
       log_dm_Mcvt_ac_out(8)=log_dm_Mcvt_ac; log_dm_Mcvt_ac_out(1,7)=set_log_dm_Mcvt_ac;
       log_dm_cL_ac_out(8)=log_dm_cL_ac; log_dm_cL_ac_out(1,7)=set_log_dm_cL_ac;
	   log_dm_cP_ac_out(8)=log_dm_cP_ac; log_dm_cP_ac_out(1,7)=set_log_dm_cP_ac;
	   log_dm_HB_ac_out(8)=log_dm_HB_ac; log_dm_HB_ac_out(1,7)=set_log_dm_HB_ac;
	   //log_dm_mrip_ac_out(8)=log_dm_mrip_ac; log_dm_mrip_ac_out(1,7)=set_log_dm_mrip_ac;
	   
	   
       selpar_A50_Mbft_out(8)=selpar_A50_Mbft; selpar_A50_Mbft_out(1,7)=set_selpar_A50_Mbft;
       selpar_slope_Mbft_out(8)=selpar_slope_Mbft; selpar_slope_Mbft_out(1,7)=set_selpar_slope_Mbft;
       
	   selpar_A50_Mcvt_out(8)=selpar_A50_Mcvt; selpar_A50_Mcvt_out(1,7)=set_selpar_A50_Mcvt;
       selpar_slope_Mcvt_out(8)=selpar_slope_Mcvt; selpar_slope_Mcvt_out(1,7)=set_selpar_slope_Mcvt;
       
	   //selpar_A50_cL1_out(8)=selpar_A50_cL1; selpar_A50_cL1_out(1,7)=set_selpar_A50_cL1;
       //selpar_slope_cL1_out(8)=selpar_slope_cL1; selpar_slope_cL1_out(1,7)=set_selpar_slope_cL1;
       selpar_A50_cL2_out(8)=selpar_A50_cL2; selpar_A50_cL2_out(1,7)=set_selpar_A50_cL2;
       selpar_slope_cL2_out(8)=selpar_slope_cL2; selpar_slope_cL2_out(1,7)=set_selpar_slope_cL2;
       selpar_A50_cL3_out(8)=selpar_A50_cL3; selpar_A50_cL3_out(1,7)=set_selpar_A50_cL3;
       selpar_slope_cL3_out(8)=selpar_slope_cL3; selpar_slope_cL3_out(1,7)=set_selpar_slope_cL3;
	   selpar_A50_cL4_out(8)=selpar_A50_cL4; selpar_A50_cL4_out(1,7)=set_selpar_A50_cL4;
       selpar_slope_cL4_out(8)=selpar_slope_cL4; selpar_slope_cL4_out(1,7)=set_selpar_slope_cL4;
	   
	   selpar_A50_cP2_out(8)=selpar_A50_cP2; selpar_A50_cP2_out(1,7)=set_selpar_A50_cP2;
       selpar_slope_cP2_out(8)=selpar_slope_cP2; selpar_slope_cP2_out(1,7)=set_selpar_slope_cP2;
       selpar_A50_cP3_out(8)=selpar_A50_cP3; selpar_A50_cP3_out(1,7)=set_selpar_A50_cP3;
       selpar_slope_cP3_out(8)=selpar_slope_cP3; selpar_slope_cP3_out(1,7)=set_selpar_slope_cP3;
	   selpar_A50_cP4_out(8)=selpar_A50_cP4; selpar_A50_cP4_out(1,7)=set_selpar_A50_cP4;
       selpar_slope_cP4_out(8)=selpar_slope_cP4; selpar_slope_cP4_out(1,7)=set_selpar_slope_cP4;
          
       selpar_A50_HB1_out(8)=selpar_A50_HB1; selpar_A50_HB1_out(1,7)=set_selpar_A50_HB1;
       selpar_slope_HB1_out(8)=selpar_slope_HB1; selpar_slope_HB1_out(1,7)=set_selpar_slope_HB1;
       selpar_A50_HB2_out(8)=selpar_A50_HB2; selpar_A50_HB2_out(1,7)=set_selpar_A50_HB2;
       selpar_slope_HB2_out(8)=selpar_slope_HB2; selpar_slope_HB2_out(1,7)=set_selpar_slope_HB2;
       selpar_A50_HB3_out(8)=selpar_A50_HB3; selpar_A50_HB3_out(1,7)=set_selpar_A50_HB3;
       selpar_slope_HB3_out(8)=selpar_slope_HB3; selpar_slope_HB3_out(1,7)=set_selpar_slope_HB3;
	   selpar_A50_HB4_out(8)=selpar_A50_HB4; selpar_A50_HB4_out(1,7)=set_selpar_A50_HB4;
       selpar_slope_HB4_out(8)=selpar_slope_HB4; selpar_slope_HB4_out(1,7)=set_selpar_slope_HB4;
	   selpar_A50_HB5_out(8)=selpar_A50_HB5; selpar_A50_HB5_out(1,7)=set_selpar_A50_HB5;
       selpar_slope_HB5_out(8)=selpar_slope_HB5; selpar_slope_HB5_out(1,7)=set_selpar_slope_HB5;
	  
	   selpar_A50_mrip1_out(8)=selpar_A50_mrip1; selpar_A50_mrip1_out(1,7)=set_selpar_A50_mrip1;
       selpar_slope_mrip1_out(8)=selpar_slope_mrip1; selpar_slope_mrip1_out(1,7)=set_selpar_slope_mrip1;
       selpar_A50_mrip2_out(8)=selpar_A50_mrip2; selpar_A50_mrip2_out(1,7)=set_selpar_A50_mrip2;
       selpar_slope_mrip2_out(8)=selpar_slope_mrip2; selpar_slope_mrip2_out(1,7)=set_selpar_slope_mrip2;
       selpar_A50_mrip3_out(8)=selpar_A50_mrip3; selpar_A50_mrip3_out(1,7)=set_selpar_A50_mrip3;
       selpar_slope_mrip3_out(8)=selpar_slope_mrip3; selpar_slope_mrip3_out(1,7)=set_selpar_slope_mrip3;
	   selpar_A50_mrip4_out(8)=selpar_A50_mrip4; selpar_A50_mrip4_out(1,7)=set_selpar_A50_mrip4;
       selpar_slope_mrip4_out(8)=selpar_slope_mrip4; selpar_slope_mrip4_out(1,7)=set_selpar_slope_mrip4;
	   selpar_A50_mrip5_out(8)=selpar_A50_mrip5; selpar_A50_mrip5_out(1,7)=set_selpar_A50_mrip5;
       selpar_slope_mrip5_out(8)=selpar_slope_mrip5; selpar_slope_mrip5_out(1,7)=set_selpar_slope_mrip5;

	   selpar_Age0_HB_D_logit_out(8)=selpar_Age0_HB_D_logit; selpar_Age0_HB_D_logit_out(1,7)=set_selpar_Age0_HB_D_logit;
       selpar_Age1_HB_D_logit_out(8)=selpar_Age1_HB_D_logit; selpar_Age1_HB_D_logit_out(1,7)=set_selpar_Age1_HB_D_logit;
	   selpar_Age2_HB_D_logit_out(8)=selpar_Age2_HB_D_logit; selpar_Age2_HB_D_logit_out(1,7)=set_selpar_Age2_HB_D_logit;

	   selpar_A50_HBD4_out(8)=selpar_A50_HBD4; selpar_A50_HBD4_out(1,7)=set_selpar_A50_HBD4;
       selpar_slope_HBD4_out(8)=selpar_slope_HBD4; selpar_slope_HBD4_out(1,7)=set_selpar_slope_HBD4;
	   selpar_A502_HBD4_out(8)=selpar_A502_HBD4; selpar_A502_HBD4_out(1,7)=set_selpar_A502_HBD4;
       selpar_slope2_HBD4_out(8)=selpar_slope2_HBD4; selpar_slope2_HBD4_out(1,7)=set_selpar_slope2_HBD4;
	   selpar_A50_HBD5_out(8)=selpar_A50_HBD5; selpar_A50_HBD5_out(1,7)=set_selpar_A50_HBD5;
       selpar_slope_HBD5_out(8)=selpar_slope_HBD5; selpar_slope_HBD5_out(1,7)=set_selpar_slope_HBD5;
	   selpar_A502_HBD5_out(8)=selpar_A502_HBD5; selpar_A502_HBD5_out(1,7)=set_selpar_A502_HBD5;
       selpar_slope2_HBD5_out(8)=selpar_slope2_HBD5; selpar_slope2_HBD5_out(1,7)=set_selpar_slope2_HBD5;
	   
       log_q_Mbft_out(8)=log_q_Mbft; log_q_Mbft_out(1,7)=set_logq_Mbft;
	   log_q_Mcvt_out(8)=log_q_Mcvt; log_q_Mcvt_out(1,7)=set_logq_Mcvt;
	   //log_q_Vid_out(8)=log_q_Vid; log_q_Vid_out(1,7)=set_logq_Vid;
       log_q_cL_out(8)=log_q_cL; log_q_cL_out(1,7)=set_logq_cL;
	   log_q_HB_out(8)=log_q_HB; log_q_HB_out(1,7)=set_logq_HB;
                            
       log_avg_F_cL_out(8)=log_avg_F_cL; log_avg_F_cL_out(1,7)=set_log_avg_F_cL;
	   log_avg_F_cP_out(8)=log_avg_F_cP; log_avg_F_cP_out(1,7)=set_log_avg_F_cP;
       log_avg_F_HB_out(8)=log_avg_F_HB; log_avg_F_HB_out(1,7)=set_log_avg_F_HB;
       log_avg_F_cT_out(8)=log_avg_F_cT; log_avg_F_cT_out(1,7)=set_log_avg_F_cT;       
       log_avg_F_mrip_out(8)=log_avg_F_mrip; log_avg_F_mrip_out(1,7)=set_log_avg_F_mrip;
       log_avg_F_comm_D_out(8)=log_avg_F_comm_D; log_avg_F_comm_D_out(1,7)=set_log_avg_F_comm_D;
       log_avg_F_HB_D_out(8)=log_avg_F_HB_D; log_avg_F_HB_D_out(1,7)=set_log_avg_F_HB_D;
       log_avg_F_mrip_D_out(8)=log_avg_F_mrip_D; log_avg_F_mrip_D_out(1,7)=set_log_avg_F_mrip_D;
        
       log_rec_dev_out(styr_rec_dev, endyr_rec_phase2)=log_rec_dev;
       log_F_dev_cL_out(styr_cL_L,endyr_cL_L)=log_F_dev_cL;
	   log_F_dev_cP_out(styr_cP_L,endyr_cP_L)=log_F_dev_cP;
       log_F_dev_cT_out(styr_cT_L,endyr_cT_L)=log_F_dev_cT;
       log_F_dev_HB_out(styr_HB_L,endyr_HB_L)=log_F_dev_HB;
       log_F_dev_mrip_out(styr_mrip_L,endyr_mrip_L)=log_F_dev_mrip;
       log_F_dev_comm_D_out(styr_cL_D,endyr_cL_D)=log_F_dev_comm_D;
       log_F_dev_HB_D_out(styr_HB_D,endyr_HB_D)=log_F_dev_HB_D;
       log_F_dev_mrip_D_out(styr_mrip_D,endyr_mrip_D)=log_F_dev_mrip_D;
     #include "bsb_make_Robjectv14.cxx"   // write the S-compatible report

  } //endl last phase loop     
  
  
