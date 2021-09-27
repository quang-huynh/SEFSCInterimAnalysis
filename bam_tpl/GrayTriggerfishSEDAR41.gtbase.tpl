//#--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##
//##  SEDAR 41 benchmark  Gray Triggerfish assessment October 2015
//##  Converted from Gag update assessment (R:\PopDyn\Confidential\SEDAR\Updates2014\Gag\Assessment\BaseRun)
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

//Starting and ending year of the model (year data starts)
init_int styr;
init_int endyr;

// Starting year to estimate recruitment deviation from S-R curve
init_int styr_rec_dev;
// Ending year to estimate recruitment deviation from S-R curve
init_int endyr_rec_dev;
// Possible 3 phases of constraints on recruitment deviations
init_int endyr_rec_phase1;
init_int endyr_rec_phase2;

// Ending years for selectivity blocks
init_int endyr_selex_phase1;
//init_int endyr_selex_phase2; //gt has no selectivity blocks

// Number assessment years
number nyrs;
number nyrs_rec;
// This section MUST BE INDENTED!!!
 LOCAL_CALCS
   nyrs=endyr-styr+1.;
   nyrs_rec=endyr_rec_dev-styr_rec_dev+1.;
 END_CALCS
 
// Total number of ages in population model
init_int nages;
// Vector of ages for age bins in population model
init_vector agebins(1,nages);
 
// Total number of ages used to match age comps: plus group may differ from popn, first age must not
init_int nages_agec;
// Vector of ages for age bins in age comps
init_vector agebins_agec(1,nages_agec);

// Total number of length bins for each matrix and width of bins
init_int nlenbins;          //used to match data
init_number lenbins_width;  //width of length bins (mm)

// Vector of lengths for length bins (mm)(midpoint) 
init_vector lenbins(1,nlenbins);
 
// Max F used in spr and msy calcs
init_number max_F_spr_msy;
// Total number of iterations for spr calcs
init_int n_iter_spr;
// Total number of iterations for msy calcs
int n_iter_msy;
// Number years at end of time series over which to average sector F's, for weighted selectivities
 LOCAL_CALCS
		n_iter_msy=n_iter_spr; 
 END_CALCS

init_int selpar_n_yrs_wgted;
// Bias correction (set to 1.0 for no bias correction or a negative value to compute from rec variance)
init_number set_BiasCor;

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: observed data section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//################Commercial handline fleet #######################################
// Comm HL CPUE
init_int styr_cH_cpue;                                             
init_int endyr_cH_cpue;                                            
init_vector obs_cH_cpue(styr_cH_cpue,endyr_cH_cpue); //Observed CPUE
init_vector cH_cpue_cv(styr_cH_cpue,endyr_cH_cpue);  //CV of cpue

// Comm HL (Landings+Others+Discards; 1000 lb whole weight)
init_int styr_cH_L;
init_int endyr_cH_L;
init_vector obs_cH_L(styr_cH_L,endyr_cH_L);
init_vector cH_L_cv(styr_cH_L,endyr_cH_L);

// Comm HL length compositions (3 cm bins)
init_int nyr_cH_lenc;
init_ivector yrs_cH_lenc(1,nyr_cH_lenc);
init_vector nsamp_cH_lenc(1,nyr_cH_lenc);
init_vector nfish_cH_lenc(1,nyr_cH_lenc);
init_matrix obs_cH_lenc(1,nyr_cH_lenc,1,nlenbins);
 
// Comm HL age compositions
init_int nyr_cH_agec;
init_ivector yrs_cH_agec(1,nyr_cH_agec);
init_vector nsamp_cH_agec(1,nyr_cH_agec);
init_vector nfish_cH_agec(1,nyr_cH_agec);
init_matrix obs_cH_agec(1,nyr_cH_agec,1,nages_agec);

//################MARMAP Chevron trap #############################################
// SERFS chevron trap CPUE
init_int styr_mcvt_cpue;                                             
init_int endyr_mcvt_cpue;                                            
init_vector obs_mcvt_cpue(styr_mcvt_cpue,endyr_mcvt_cpue); //Observed CPUE
init_vector mcvt_cpue_cv(styr_mcvt_cpue,endyr_mcvt_cpue);  //CV of cpue

// SERFS chevron trap length compositions (3 cm bins)
init_int nyr_mcvt_lenc;
init_ivector yrs_mcvt_lenc(1,nyr_mcvt_lenc);
init_vector nsamp_mcvt_lenc(1,nyr_mcvt_lenc);
init_vector nfish_mcvt_lenc(1,nyr_mcvt_lenc);
init_matrix obs_mcvt_lenc(1,nyr_mcvt_lenc,1,nlenbins);
 
// SERFS chevron trap age compositions
init_int nyr_mcvt_agec;
init_ivector yrs_mcvt_agec(1,nyr_mcvt_agec);
init_vector nsamp_mcvt_agec(1,nyr_mcvt_agec);
init_vector nfish_mcvt_agec(1,nyr_mcvt_agec);
init_matrix obs_mcvt_agec(1,nyr_mcvt_agec,1,nages_agec);

//################SEFIS video index ###############################################
// SEFIS video CPUE
//init_int styr_vid_cpue;                                             
//init_int endyr_vid_cpue;                                            
//init_vector obs_vid_cpue(styr_vid_cpue,endyr_vid_cpue); //Observed CPUE
//init_vector vid_cpue_cv(styr_vid_cpue,endyr_vid_cpue);  //CV of cpue

//###################Headboat fleet ###############################################
// HB CPUE
init_int styr_HB_cpue;                                             
init_int endyr_HB_cpue;                                            
init_vector obs_HB_cpue(styr_HB_cpue,endyr_HB_cpue);//Observed CPUE
init_vector HB_cpue_cv(styr_HB_cpue,endyr_HB_cpue); //CV of cpue

// HB landings (1000s fish)
init_int styr_HB_L;
init_int endyr_HB_L;
init_vector obs_HB_L(styr_HB_L,endyr_HB_L);   //vector of observed landings by year 
init_vector HB_L_cv(styr_HB_L,endyr_HB_L);    //vector of CV of landings by year

// HB landings length compositions (3 cm bins)
init_int nyr_HB_lenc;
init_ivector yrs_HB_lenc(1,nyr_HB_lenc);
init_vector nsamp_HB_lenc(1,nyr_HB_lenc);
init_vector nfish_HB_lenc(1,nyr_HB_lenc);
init_matrix obs_HB_lenc(1,nyr_HB_lenc,1,nlenbins);
 
// HB landings age compositions
init_int nyr_HB_agec;
init_ivector yrs_HB_agec(1,nyr_HB_agec);
init_vector nsamp_HB_agec(1,nyr_HB_agec);
init_vector nfish_HB_agec(1,nyr_HB_agec);
init_matrix obs_HB_agec(1,nyr_HB_agec,1,nages_agec);

// HB  discards (1000 fish)
init_int styr_HB_D;
init_int endyr_HB_D;
init_vector obs_HB_released(styr_HB_D,endyr_HB_D);
init_vector HB_D_cv(styr_HB_D,endyr_HB_D);

// HB discard length compositions
init_int nyr_HB_D_lenc;
init_ivector yrs_HB_D_lenc(1,nyr_HB_D_lenc);
init_vector nsamp_HB_D_lenc(1,nyr_HB_D_lenc);
init_vector nfish_HB_D_lenc(1,nyr_HB_D_lenc);
init_matrix obs_HB_D_lenc(1,nyr_HB_D_lenc,1,nlenbins);

// !! cout << styr_HB_D << endl; 
//###################General Recreational fleet ###################################
// GR (MRIP) CPUE
init_int styr_GR_cpue;                                             
init_int endyr_GR_cpue;                                            
init_vector obs_GR_cpue(styr_GR_cpue,endyr_GR_cpue);//Observed CPUE
init_vector GR_cpue_cv(styr_GR_cpue,endyr_GR_cpue); //CV of cpue

// GR (MRIP) landings (1000s fish)
init_int styr_GR_L;
init_int endyr_GR_L;
init_vector obs_GR_L(styr_GR_L,endyr_GR_L);   //vector of observed landings by year 
init_vector GR_L_cv(styr_GR_L,endyr_GR_L);    //vector of CV of landings by year

// GR (MRIP) landings length compositions (3 cm bins)
init_int nyr_GR_lenc;
init_ivector yrs_GR_lenc(1,nyr_GR_lenc);
init_vector nsamp_GR_lenc(1,nyr_GR_lenc);
init_vector nfish_GR_lenc(1,nyr_GR_lenc);
init_matrix obs_GR_lenc(1,nyr_GR_lenc,1,nlenbins);

// // GR landings age compositions (single pooled age comp over 2004-2005 MRIP CB mode; weighted by number trips)
// init_int nyr_GR_agec;
// init_ivector yrs_GR_agec(1,nyr_GR_agec);
// init_vector nsamp_GR_agec(1,nyr_GR_agec);
// init_vector nfish_GR_agec(1,nyr_GR_agec);
// init_matrix obs_GR_agec(1,nyr_GR_agec,1,nages_agec);

// GR  discards (1000 fish)
init_int styr_GR_D;
init_int endyr_GR_D;
init_vector obs_GR_released(styr_GR_D,endyr_GR_D);
init_vector GR_D_cv(styr_GR_D,endyr_GR_D);

//  !! cout << obs_GR_released << endl; 
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: parameter section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//##################Single Parameter values and initial guesses ######################################
// Von Bert parameters in mm FL for the popn
init_vector set_Linf(1,7);
init_vector set_K(1,7);
init_vector set_t0(1,7);
//CV of length at age or its standard error for the popn
init_vector set_len_cv(1,7);

// FISHERY LANDINGS   Von Bert parameters in mm FL for all FD fish
init_vector set_Linf_L(1,7);
init_vector set_K_L(1,7);
init_vector set_t0_L(1,7);
//CV of length at age landings
init_vector set_len_cv_L(1,7);

// FISHERY INDEPENDENT  Von Bert parameters in mm FL mm all FI (SERFS chevron trap) fish
init_vector set_Linf_mcvt(1,7);
init_vector set_K_mcvt(1,7);
init_vector set_t0_mcvt(1,7);
//CV of length at age fishery-independent trap
init_vector set_len_cv_mcvt(1,7)

// Scalar used only for computing MSST. This might eventually be based on 0.25: MSST=(1-0.25)SSBmsy
init_vector set_M_constant(1,7);
    
// Spawner-recruit parameters (Initial guesses or fixed values)
init_vector set_steep(1,7);         //recruitment steepness
init_vector set_log_R0(1,7);        //recruitment R0
init_vector set_R_autocorr(1,7);    //recruitment autocorrelation
init_vector set_rec_sigma(1,7);     //recruitment standard deviation in log space

// Initial guesses or fixed values of estimated selectivity parameters
// Commercial handline selectivity
init_vector set_selpar_L50_cH1(1,7);
init_vector set_selpar_slope_cH1(1,7);

// init_vector set_selpar_L50_cH2(1,7);   //if using selectivity blocks
// init_vector set_selpar_slope_cH2(1,7); //if using selectivity blocks

// init_vector set_selpar_L50_cH3(1,7);   //if using selectivity blocks
// init_vector set_selpar_slope_cH3(1,7); //if using selectivity blocks

// SERFS chevron trap selectivity
init_vector set_selpar_L50_mcvt(1,7);
init_vector set_selpar_slope_mcvt(1,7);

// Headboat selectivity (set up as 2 blocks to allow for separate selectivities for specific sets of years)
init_vector set_selpar_L51_HB1(1,7);
init_vector set_selpar_slope1_HB1(1,7);
init_vector set_selpar_L52_HB1(1,7);
init_vector set_selpar_slope2_HB1(1,7);

init_vector set_selpar_L50_HB2(1,7);
init_vector set_selpar_slope_HB2(1,7);
init_vector set_selpar_afull_HB2(1,7);
init_vector set_selpar_sigma_HB2(1,7)
    
// init_vector set_selpar_L50_HB3(1,7); //if using selectivity blocks
// init_vector set_selpar_slope_HB3(1,7); //if using selectivity blocks

// Headboat discard selectivity
init_vector set_selpar_L50_HB_D(1,7);
init_vector set_selpar_slope_HB_D(1,7);
init_vector set_selpar_afull_HB_D(1,7);
init_vector set_selpar_sigma_HB_D(1,7);

// GR selectivity
init_vector set_selpar_L51_GR(1,7);
init_vector set_selpar_slope1_GR(1,7);
init_vector set_selpar_L52_GR(1,7);
init_vector set_selpar_slope2_GR(1,7);

//--index catchability-----------------------------------------------------------------------------------
init_vector set_log_q_cH(1,7);      //catchability coefficient (log) for comm handline index
init_vector set_log_q_mcvt(1,7);    //catchability coefficient (log) for SERFS chevron trap index
//init_vector set_log_q_vid(1,7);      //catchability coefficient (log) for SEFIS video index
init_vector set_log_q_HB(1,7);      //catchability coefficient (log) for headboat index
init_vector set_log_q_GR(1,7);      //catchability coefficient (log) for general rec index

// Initial F
init_vector set_F_init(1,7);  //scales initial F 
//--mean F's in log space -------------------------------------------------------------------------------
init_vector set_log_avg_F_cH(1,7); 
init_vector set_log_avg_F_HB(1,7);
init_vector set_log_avg_F_GR(1,7);
init_vector set_log_avg_F_HB_D(1,7);
init_vector set_log_avg_F_GR_D(1,7);

//##################Dev Vector Parameter values (vals) and bounds #######################################
//--F vectors--------------------------------------------------------------------------------------------
init_vector set_log_F_dev_cH(1,3);         
init_vector set_log_F_dev_HB(1,3);
init_vector set_log_F_dev_GR(1,3);
init_vector set_log_F_dev_HB_D(1,3);
init_vector set_log_F_dev_GR_D(1,3);
init_vector set_log_rec_dev(1,3);
init_vector set_log_Nage_dev(1,3);

init_vector set_log_F_dev_cH_vals(styr_cH_L,endyr_cH_L);
init_vector set_log_F_dev_HB_vals(styr_HB_L,endyr_HB_L);
init_vector set_log_F_dev_GR_vals(styr_GR_L,endyr_GR_L);
init_vector set_log_F_dev_HB_D_vals(styr_HB_D,endyr_HB_D);
init_vector set_log_F_dev_GR_D_vals(styr_GR_D,endyr_GR_D);
init_vector set_log_rec_dev_vals(styr_rec_dev,endyr_rec_dev);
init_vector set_log_Nage_dev_vals(2,nages);             

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: likelihood weights section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

init_number set_w_L;            //weight for landings
init_number set_w_D;            //weight for discards
init_number set_w_I_cH;         //weight for comm handline index
init_number set_w_I_mcvt;       //weight for SERFS chevron trap index
//init_number set_w_I_vid;        //weight for SEFIS video index
init_number set_w_I_HB;         //weight for headboat index
init_number set_w_I_GR;         //weight for MRFSS index
init_number set_w_lc_cH;        //weight for comm handline len comps
init_number set_w_lc_mcvt;      //weight for SERFS chevron trap len comps
init_number set_w_lc_HB;        //weight for headboat len comps
init_number set_w_lc_HB_D;      //weight for headboat discards len comps
init_number set_w_lc_GR;        //weight for GR (mrip) len comps 
init_number set_w_ac_cH;        //weight for comm handline age comps
init_number set_w_ac_mcvt;      //weight for SERFS chevron trap age comps
init_number set_w_ac_HB;        //weight for headboat age comps
//init_number set_w_ac_GR;        //weight for GR (MRIP) age comps
init_number set_w_Nage_init;    //for fitting initial abundance at age (excluding first age)
init_number set_w_rec;          //for fitting S-R curve
init_number set_w_rec_early;    //additional constraint on early years recruitment
init_number set_w_rec_end;      //additional constraint on ending years recruitment 
init_number set_w_fullF;        //penalty for any Fapex>3(removed in final phase of optimization)
init_number set_w_Ftune;        //weight applied to tuning F (removed in final phase of optimization)

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: miscellaneous stuff section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

// FL(mm)-weight(whole weight in kg) relationship: W=aL^b
init_number wgtpar_a;
init_number wgtpar_b;

// Weight(whole weight)- fecundity (units=kg) relationship: log(y)=a+bW
init_number fecpar_a;
init_number fecpar_b;
//init_number fecpar_batches;        //used if assuming age-invariant batch number
init_vector fecpar_batches(1,nages); //number of annual batches; may be important if age-dependent, otherwise just a scalar
init_number fecpar_scale;            //used for scaling annual egg production (10^X eggs)

// Maturity and proportion female at age
init_vector maturity_f_obs(1,nages);        //proportion females mature at age
init_vector prop_f_obs(1,nages);            //proportion female at age

init_number spawn_time_frac; //time of year of peak spawning, as a fraction of the year

// Natural mortality
init_vector set_M(1,nages);     //age-dependent: used in model
init_number max_obs_age;        //max observed age, used to scale M, if estimated

// Discard mortality constants
init_number set_Dmort_HB;
init_number set_Dmort_GR;

//Spawner-recruit parameters (initial guesses or fixed values)
init_int SR_switch;

// Rate of increase on q
init_int set_q_rate_phase;  //value sets estimation phase of rate increase, negative value turns it off
init_number set_q_rate;
// Density dependence on fishery q's 
init_int set_q_DD_phase;      //value sets estimation phase of random walk, negative value turns it off
init_number set_q_DD_beta;    //value of 0.0 is density indepenent
init_number set_q_DD_beta_se;
init_int set_q_DD_stage;      //age to begin counting biomass, should be near full exploitation

// Random walk on fishery q's 
init_int set_q_RW_phase;         //value sets estimation phase of random walk, negative value turns it off
init_number set_q_RW_rec_var;    //assumed variance of RW q

// Tune Fapex (tuning removed in final year of optimization)
init_number set_Ftune;
init_int set_Ftune_yr;

// Threshold sample sizes for length comps 
init_number minSS_cH_lenc;
init_number minSS_mcvt_lenc;
init_number minSS_HB_lenc;
init_number minSS_HB_D_lenc;  
init_number minSS_GR_lenc; 

//Threshold sample sizes for age comps
init_number minSS_cH_agec;
init_number minSS_mcvt_agec;
init_number minSS_HB_agec;
//init_number minSS_GR_agec;

//Ageing error matrix (columns are true ages, rows are ages as read for age comps: columns should sum to one)
init_matrix age_error(1,nages,1,nages);

// #######Indexing integers for year(iyear), age(iage),length(ilen) ##################################
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

init_number use_landings_growth;

init_number end_of_data_file;
//this section MUST BE INDENTED!!!
 LOCAL_CALCS
   if(end_of_data_file!=999)
   {
       cout << "*** WARNING: Data File NOT READ CORRECTLY ****" << endl;
       exit(0);  
   }
   else
   {cout << "Data File read correctly" << endl;} 
 END_CALCS   



PARAMETER_SECTION

 LOCAL_CALCS
  //POPULATION growth
  const double Linf_LO=set_Linf(2); const double Linf_HI=set_Linf(3); const double Linf_PH=set_Linf(4);
  const double K_LO=set_K(2); const double K_HI=set_K(3); const double K_PH=set_K(4);
  const double t0_LO=set_t0(2); const double t0_HI=set_t0(3); const double t0_PH=set_t0(4);  
  const double len_cv_LO=set_len_cv(2); const double len_cv_HI=set_len_cv(3); const double len_cv_PH=set_len_cv(4); 

  //FISHERY LANDINGS growth
  const double Linf_L_LO=set_Linf_L(2); const double Linf_L_HI=set_Linf_L(3); const double Linf_L_PH=set_Linf_L(4);
  const double K_L_LO=set_K_L(2); const double K_L_HI=set_K_L(3); const double K_L_PH=set_K_L(4);
  const double t0_L_LO=set_t0_L(2); const double t0_L_HI=set_t0_L(3); const double t0_L_PH=set_t0_L(4);  
  const double len_cv_L_LO=set_len_cv_L(2); const double len_cv_L_HI=set_len_cv_L(3); const double len_cv_L_PH=set_len_cv_L(4);

  //FISHERY INDEPENDENT growth (SERFS chevron trap)
  const double Linf_mcvt_LO=set_Linf_mcvt(2); const double Linf_mcvt_HI=set_Linf_mcvt(3); const double Linf_mcvt_PH=set_Linf_mcvt(4);
  const double K_mcvt_LO=set_K_mcvt(2); const double K_mcvt_HI=set_K_mcvt(3); const double K_mcvt_PH=set_K_mcvt(4);
  const double t0_mcvt_LO=set_t0_mcvt(2); const double t0_mcvt_HI=set_t0_mcvt(3); const double t0_mcvt_PH=set_t0_mcvt(4);  
  const double len_cv_mcvt_LO=set_len_cv_mcvt(2); const double len_cv_mcvt_HI=set_len_cv_mcvt(3); const double len_cv_mcvt_PH=set_len_cv_mcvt(4);
  
  const double M_constant_LO=set_M_constant(2); const double M_constant_HI=set_M_constant(3); const double M_constant_PH=set_M_constant(4);        

  const double steep_LO=set_steep(2); const double steep_HI=set_steep(3); const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); const double log_R0_HI=set_log_R0(3); const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); const double R_autocorr_HI=set_R_autocorr(3); const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); const double rec_sigma_HI=set_rec_sigma(3); const double rec_sigma_PH=set_rec_sigma(4);

  const double selpar_L50_cH1_LO=set_selpar_L50_cH1(2); const double selpar_L50_cH1_HI=set_selpar_L50_cH1(3); const double selpar_L50_cH1_PH=set_selpar_L50_cH1(4);
  const double selpar_slope_cH1_LO=set_selpar_slope_cH1(2); const double selpar_slope_cH1_HI=set_selpar_slope_cH1(3); const double selpar_slope_cH1_PH=set_selpar_slope_cH1(4);
//const double selpar_L50_cH2_LO=set_selpar_L50_cH2(2); const double selpar_L50_cH2_HI=set_selpar_L50_cH2(3); const double selpar_L50_cH2_PH=set_selpar_L50_cH2(4);
//const double selpar_slope_cH2_LO=set_selpar_slope_cH2(2); const double selpar_slope_cH2_HI=set_selpar_slope_cH2(3); const double selpar_slope_cH2_PH=set_selpar_slope_cH2(4);
//const double selpar_L50_cH3_LO=set_selpar_L50_cH3(2); const double selpar_L50_cH3_HI=set_selpar_L50_cH3(3); const double selpar_L50_cH3_PH=set_selpar_L50_cH3(4);
//const double selpar_slope_cH3_LO=set_selpar_slope_cH3(2); const double selpar_slope_cH3_HI=set_selpar_slope_cH3(3); const double selpar_slope_cH3_PH=set_selpar_slope_cH3(4);
  const double selpar_L50_mcvt_LO=set_selpar_L50_mcvt(2); const double selpar_L50_mcvt_HI=set_selpar_L50_mcvt(3); const double selpar_L50_mcvt_PH=set_selpar_L50_mcvt(4);
  const double selpar_slope_mcvt_LO=set_selpar_slope_mcvt(2); const double selpar_slope_mcvt_HI=set_selpar_slope_mcvt(3); const double selpar_slope_mcvt_PH=set_selpar_slope_mcvt(4); 

  const double selpar_L51_HB1_LO=set_selpar_L51_HB1(2); const double selpar_L51_HB1_HI=set_selpar_L51_HB1(3); const double selpar_L51_HB1_PH=set_selpar_L51_HB1(4);
  const double selpar_slope1_HB1_LO=set_selpar_slope1_HB1(2); const double selpar_slope1_HB1_HI=set_selpar_slope1_HB1(3); const double selpar_slope1_HB1_PH=set_selpar_slope1_HB1(4);
  const double selpar_L52_HB1_LO=set_selpar_L52_HB1(2); const double selpar_L52_HB1_HI=set_selpar_L52_HB1(3); const double selpar_L52_HB1_PH=set_selpar_L52_HB1(4);
  const double selpar_slope2_HB1_LO=set_selpar_slope2_HB1(2); const double selpar_slope2_HB1_HI=set_selpar_slope2_HB1(3); const double selpar_slope2_HB1_PH=set_selpar_slope2_HB1(4);

  const double selpar_L50_HB2_LO=set_selpar_L50_HB2(2); const double selpar_L50_HB2_HI=set_selpar_L50_HB2(3); const double selpar_L50_HB2_PH=set_selpar_L50_HB2(4);
  const double selpar_slope_HB2_LO=set_selpar_slope_HB2(2); const double selpar_slope_HB2_HI=set_selpar_slope_HB2(3); const double selpar_slope_HB2_PH=set_selpar_slope_HB2(4);
  const double selpar_afull_HB2_LO=set_selpar_afull_HB2(2); const double selpar_afull_HB2_HI=set_selpar_afull_HB2(3); const double selpar_afull_HB2_PH=set_selpar_afull_HB2(4);
  const double selpar_sigma_HB2_LO=set_selpar_sigma_HB2(2); const double selpar_sigma_HB2_HI=set_selpar_sigma_HB2(3); const double selpar_sigma_HB2_PH=set_selpar_sigma_HB2(4);

//const double selpar_L50_HB3_LO=set_selpar_L50_HB3(2); const double selpar_L50_HB3_HI=set_selpar_L50_HB3(3); const double selpar_L50_HB3_PH=set_selpar_L50_HB3(4);
//const double selpar_slope_HB3_LO=set_selpar_slope_HB3(2); const double selpar_slope_HB3_HI=set_selpar_slope_HB3(3); const double selpar_slope_HB3_PH=set_selpar_slope_HB3(4);

  const double selpar_L50_HB_D_LO=set_selpar_L50_HB_D(2); const double selpar_L50_HB_D_HI=set_selpar_L50_HB_D(3); const double selpar_L50_HB_D_PH=set_selpar_L50_HB_D(4);
  const double selpar_slope_HB_D_LO=set_selpar_slope_HB_D(2); const double selpar_slope_HB_D_HI=set_selpar_slope_HB_D(3); const double selpar_slope_HB_D_PH=set_selpar_slope_HB_D(4);
  const double selpar_afull_HB_D_LO=set_selpar_afull_HB_D(2); const double selpar_afull_HB_D_HI=set_selpar_afull_HB_D(3); const double selpar_afull_HB_D_PH=set_selpar_afull_HB_D(4);
  const double selpar_sigma_HB_D_LO=set_selpar_sigma_HB_D(2); const double selpar_sigma_HB_D_HI=set_selpar_sigma_HB_D(3); const double selpar_sigma_HB_D_PH=set_selpar_sigma_HB_D(4);

  const double selpar_L51_GR_LO=set_selpar_L51_GR(2); const double selpar_L51_GR_HI=set_selpar_L51_GR(3); const double selpar_L51_GR_PH=set_selpar_L51_GR(4); 
  const double selpar_slope1_GR_LO=set_selpar_slope1_GR(2); const double selpar_slope1_GR_HI=set_selpar_slope1_GR(3); const double selpar_slope1_GR_PH=set_selpar_slope1_GR(4);  
  const double selpar_L52_GR_LO=set_selpar_L52_GR(2); const double selpar_L52_GR_HI=set_selpar_L52_GR(3); const double selpar_L52_GR_PH=set_selpar_L52_GR(4);  
  const double selpar_slope2_GR_LO=set_selpar_slope2_GR(2); const double selpar_slope2_GR_HI=set_selpar_slope2_GR(3); const double selpar_slope2_GR_PH=set_selpar_slope2_GR(4); 

  const double log_q_cH_LO=set_log_q_cH(2); const double log_q_cH_HI=set_log_q_cH(3); const double log_q_cH_PH=set_log_q_cH(4);
  const double log_q_mcvt_LO=set_log_q_mcvt(2); const double log_q_mcvt_HI=set_log_q_mcvt(3); const double log_q_mcvt_PH=set_log_q_mcvt(4);  
//const double log_q_vid_LO=set_log_q_vid(2); const double log_q_vid_HI=set_log_q_vid(3); const double log_q_vid_PH=set_log_q_vid(4);  
  const double log_q_HB_LO=set_log_q_HB(2); const double log_q_HB_HI=set_log_q_HB(3); const double log_q_HB_PH=set_log_q_HB(4);
  const double log_q_GR_LO=set_log_q_GR(2); const double log_q_GR_HI=set_log_q_GR(3); const double log_q_GR_PH=set_log_q_GR(4);

  const double F_init_LO=set_F_init(2); const double F_init_HI=set_F_init(3); const double F_init_PH=set_F_init(4);
  const double log_avg_F_cH_LO=set_log_avg_F_cH(2); const double log_avg_F_cH_HI=set_log_avg_F_cH(3); const double log_avg_F_cH_PH=set_log_avg_F_cH(4);
  const double log_avg_F_HB_LO=set_log_avg_F_HB(2); const double log_avg_F_HB_HI=set_log_avg_F_HB(3); const double log_avg_F_HB_PH=set_log_avg_F_HB(4); 
  const double log_avg_F_GR_LO=set_log_avg_F_GR(2); const double log_avg_F_GR_HI=set_log_avg_F_GR(3); const double log_avg_F_GR_PH=set_log_avg_F_GR(4); 
  const double log_avg_F_HB_D_LO=set_log_avg_F_HB_D(2); const double log_avg_F_HB_D_HI=set_log_avg_F_HB_D(3); const double log_avg_F_HB_D_PH=set_log_avg_F_HB_D(4); 
  const double log_avg_F_GR_D_LO=set_log_avg_F_GR_D(2); const double log_avg_F_GR_D_HI=set_log_avg_F_GR_D(3); const double log_avg_F_GR_D_PH=set_log_avg_F_GR_D(4); 

  //-dev vectors-----------------------------------------------------------------------------------------------------------  
  const double log_F_dev_cH_LO=set_log_F_dev_cH(1); const double log_F_dev_cH_HI=set_log_F_dev_cH(2); const double log_F_dev_cH_PH=set_log_F_dev_cH(3);   
  const double log_F_dev_HB_LO=set_log_F_dev_HB(1); const double log_F_dev_HB_HI=set_log_F_dev_HB(2); const double log_F_dev_HB_PH=set_log_F_dev_HB(3);   
  const double log_F_dev_GR_LO=set_log_F_dev_GR(1); const double log_F_dev_GR_HI=set_log_F_dev_GR(2); const double log_F_dev_GR_PH=set_log_F_dev_GR(3);   

  const double log_F_dev_HB_D_LO=set_log_F_dev_HB_D(1); const double log_F_dev_HB_D_HI=set_log_F_dev_HB_D(2); const double log_F_dev_HB_D_PH=set_log_F_dev_HB_D(3);   
  const double log_F_dev_GR_D_LO=set_log_F_dev_GR_D(1); const double log_F_dev_GR_D_HI=set_log_F_dev_GR_D(2); const double log_F_dev_GR_D_PH=set_log_F_dev_GR_D(3);   

  const double log_rec_dev_LO=set_log_rec_dev(1); const double log_rec_dev_HI=set_log_rec_dev(2); const double log_rec_dev_PH=set_log_rec_dev(3);          
  const double log_Nage_dev_LO=set_log_Nage_dev(1); const double log_Nage_dev_HI=set_log_Nage_dev(2); const double log_Nage_dev_PH=set_log_Nage_dev(3);          

 END_CALCS
 
////--------------Growth-------------------------------------------------------------------------------- 
  //POPULATION GROWTH CURVE
  init_bounded_number Linf(Linf_LO,Linf_HI,Linf_PH);
  init_bounded_number K(K_LO,K_HI,K_PH);
  init_bounded_number t0(t0_LO,t0_HI,t0_PH);
  init_bounded_number len_cv_val(len_cv_LO,len_cv_HI,len_cv_PH); 

  //FISHERY LANDINGS GROWTH CURVE
  init_bounded_number Linf_L(Linf_L_LO,Linf_L_HI,Linf_L_PH);
  init_bounded_number K_L(K_L_LO,K_L_HI,K_L_PH);
  init_bounded_number t0_L(t0_L_LO,t0_L_HI,t0_L_PH);
  init_bounded_number len_cv_val_L(len_cv_L_LO,len_cv_L_HI,len_cv_L_PH);  
 
  //FISHERY INDEPENDENT (SERFS chevron trap) GROWTH CURVE
  init_bounded_number Linf_mcvt(Linf_mcvt_LO,Linf_mcvt_HI,Linf_mcvt_PH);
  init_bounded_number K_mcvt(K_mcvt_LO,K_mcvt_HI,K_mcvt_PH);
  init_bounded_number t0_mcvt(t0_mcvt_LO,t0_mcvt_HI,t0_mcvt_PH);
  init_bounded_number len_cv_val_mcvt(len_cv_mcvt_LO,len_cv_mcvt_HI,len_cv_mcvt_PH);
 
  vector Linf_out(1,8);
  vector K_out(1,8);
  vector t0_out(1,8);
  vector len_cv_val_out(1,8);

  vector Linf_L_out(1,8);
  vector K_L_out(1,8);
  vector t0_L_out(1,8);
  vector len_cv_val_L_out(1,8);

  vector Linf_mcvt_out(1,8);
  vector K_mcvt_out(1,8);
  vector t0_mcvt_out(1,8);
  vector len_cv_val_mcvt_out(1,8);

  //POPULATION GROWTH CURVE
  vector meanlen_FL(1,nages);   //mean fork length (mm) at age all fish
  vector wgt_g(1,nages);        //whole wgt in g
  vector wgt_kg(1,nages);       //whole wgt in kg
  vector wgt_mt(1,nages);       //whole wgt in mt
  vector wgt_klb(1,nages);      //whole wgt in 1000 lb
  vector wgt_lb(1,nages);       //whole wgt in lb  
  vector fecundity(1,nages);  //KC added

  //These  values used for intermediate and summary calcs as needed
  //FISHERY LANDINGS GROWTH CURVE
  vector meanlen_FL_L(1,nages);   //mean fork length (mm) at age all fish
  vector wgt_g_L(1,nages);        //whole wgt in g
  vector wgt_kg_L(1,nages);       //whole wgt in kg
  vector wgt_mt_L(1,nages);       //whole wgt in mt
  vector wgt_klb_L(1,nages);      //whole wgt in 1000 lb
  vector wgt_lb_L(1,nages);       //whole wgt in lb

  //FISHERY INDEPENDENT (SERFS chevron trap) GROWTH CURVE
  vector meanlen_FL_mcvt(1,nages);   //mean fork length (mm) at age all fish

  matrix len_cH_mm(styr,endyr,1,nages);          //mean length at age of commercial handline landings in mm (may differ from popn mean)
  matrix wholewgt_cH_klb(styr,endyr,1,nages);    //whole wgt of commercial handline landings in 1000 lb   
  matrix len_mcvt_mm(styr,endyr,1,nages);        //mean length at age of SERFS chevron trap in mm (may differ from popn mean)
  matrix wholewgt_mcvt_klb(styr,endyr,1,nages);  //whole wgt of SERFS chevron trap in 1000 lb   
  matrix len_HB_mm(styr,endyr,1,nages);          //mean length at age of HB landings in mm (may differ from popn mean)    
  matrix wholewgt_HB_klb(styr,endyr,1,nages);    //whole wgt of HB landings in 1000 lb  
  matrix len_GR_mm(styr,endyr,1,nages);          //mean length at age of GR landings in mm (may differ from popn mean)    
  matrix wholewgt_GR_klb(styr,endyr,1,nages);    //whole wgt of GR landings in 1000 lb  
  matrix len_HB_D_mm(styr,endyr,1,nages);        //mean length at age of HB discards in mm (may differ from popn mean)    
  matrix wholewgt_HB_D_klb(styr,endyr,1,nages);  //whole wgt of HB discards in 1000 lb  
  matrix len_GR_D_mm(styr,endyr,1,nages);        //mean length at age of GR discards in mm (may differ from popn mean)    
  matrix wholewgt_GR_D_klb(styr,endyr,1,nages);  //whole wgt of GR discards in 1000 lb    
  matrix lenprob(1,nages,1,nlenbins);            //distn of size at age (age-length key, 3 cm bins) in population
  matrix lenprob2(1,nages,1,nlenbins);            //distn of size at age (age-length key, 3 cm bins) in population
  number zscore_len;                             //standardized normal values used for computing lenprob
  number zscore_len2;                             //standardized normal values used for computing lenprob
  vector cprob_lenvec(1,nlenbins);               //cumulative probabilities used for computing lenprob
  vector cprob_lenvec2(1,nlenbins);               //cumulative probabilities used for computing lenprob
  number zscore_lzero;                           //standardized normal values for length = 0
  number cprob_lzero;                            //length probability mass below zero, used for computing lenprob
  
  number zscore_lzero2;                             //standardized normal values for length = 0
  number cprob_lzero2;                            //length probability mass below zero, used for computing lenprob

  //matrices below are used to match length comps
  matrix lenprob_cH(1,nages,1,nlenbins);     //distn of size at age in cH
  matrix lenprob_mcvt(1,nages,1,nlenbins);   //distn of size at age in SERFS chevron trap
  matrix lenprob_HB(1,nages,1,nlenbins);     //distn of size at age in HB
  matrix lenprob_HB_D(1,nages,1,nlenbins);   //distn of size at age in HB discards
  matrix lenprob_GR(1,nages,1,nlenbins);     //distn of size at age in GR

//  //init_bounded_dev_vector log_len_cv_dev(1,nages,-2,2,3)
//  number log_len_cv
  vector len_sd(1,nages);
  vector len_cv(1,nages); //for fishgraph 

  //LANDINGS
  vector len_sd_L(1,nages);
  vector len_cv_L(1,nages); //for fishgraph

  //SERFS chevron trap
  vector len_sd_mcvt(1,nages);
  vector len_cv_mcvt(1,nages); //for fishgraph

//----Predicted length and age compositions
  matrix pred_cH_lenc(1,nyr_cH_lenc,1,nlenbins);  
  matrix pred_mcvt_lenc(1,nyr_mcvt_lenc,1,nlenbins);  
  matrix pred_HB_lenc(1,nyr_HB_lenc,1,nlenbins); 
  matrix pred_HB_D_lenc(1,nyr_HB_D_lenc,1,nlenbins);  
  matrix pred_GR_lenc(1,nyr_GR_lenc,1,nlenbins); 

  matrix pred_cH_agec(1,nyr_cH_agec,1,nages_agec);
  matrix pred_cH_agec_allages(1,nyr_cH_agec,1,nages);
  matrix ErrorFree_cH_agec(1,nyr_cH_agec,1,nages);  
  matrix pred_mcvt_agec(1,nyr_mcvt_agec,1,nages_agec);
  matrix pred_mcvt_agec_allages(1,nyr_mcvt_agec,1,nages);  
  matrix ErrorFree_mcvt_agec(1,nyr_mcvt_agec,1,nages);
  matrix pred_HB_agec(1,nyr_HB_agec,1,nages_agec);
  matrix pred_HB_agec_allages(1,nyr_HB_agec,1,nages);  
  matrix ErrorFree_HB_agec(1,nyr_HB_agec,1,nages);
  //matrix pred_GR_agec(1,nyr_GR_agec,1,nages_agec);
  //matrix pred_GR_agec_allages(1,nyr_GR_agec,1,nages);  
  //matrix ErrorFree_GR_agec(1,nyr_GR_agec,1,nages);
  
//Effective sample size applied in multinomial distributions
  vector nsamp_cH_lenc_allyr(styr,endyr);
  vector nsamp_mcvt_lenc_allyr(styr,endyr);
  vector nsamp_HB_lenc_allyr(styr,endyr);
  vector nsamp_HB_D_lenc_allyr(styr,endyr); 
  vector nsamp_GR_lenc_allyr(styr,endyr);  
  vector nsamp_cH_agec_allyr(styr,endyr);
  vector nsamp_mcvt_agec_allyr(styr,endyr);
  vector nsamp_HB_agec_allyr(styr,endyr);
  //vector nsamp_GR_agec_allyr(styr,endyr);

//Nfish used in MCB analysis (not used in fitting)
  vector nfish_cH_lenc_allyr(styr,endyr);
  vector nfish_mcvt_lenc_allyr(styr,endyr);
  vector nfish_HB_lenc_allyr(styr,endyr);
  vector nfish_HB_D_lenc_allyr(styr,endyr);  
  vector nfish_GR_lenc_allyr(styr,endyr);  
  vector nfish_cH_agec_allyr(styr,endyr);
  vector nfish_mcvt_agec_allyr(styr,endyr);
  vector nfish_HB_agec_allyr(styr,endyr);
  //vector nfish_GR_agec_allyr(styr,endyr);

//Computed effective sample size for output (not used in fitting)
  vector neff_cH_lenc_allyr_out(styr,endyr);
  vector neff_mcvt_lenc_allyr_out(styr,endyr);
  vector neff_HB_lenc_allyr_out(styr,endyr);
  vector neff_GR_lenc_allyr_out(styr,endyr);  
  vector neff_HB_D_lenc_allyr_out(styr,endyr);  
  vector neff_cH_agec_allyr_out(styr,endyr);
  vector neff_mcvt_agec_allyr_out(styr,endyr);
  vector neff_HB_agec_allyr_out(styr,endyr);
  //vector neff_GR_agec_allyr_out(styr,endyr);

//-----Population----------------------------------------------------------------------------------------
  matrix N(styr,endyr+1,1,nages);           //Population numbers by year and age at start of yr
  matrix N_mdyr(styr,endyr,1,nages);        //Population numbers by year and age at mdpt of yr: used for comps and cpue
  matrix N_spawn(styr,endyr,1,nages);     //Population numbers by year and age at peaking spawning: used for SSB (proj yr ok bc of ssb on Jan1)  
  init_bounded_vector log_Nage_dev(2,nages,log_Nage_dev_LO,log_Nage_dev_HI,log_Nage_dev_PH);
  vector log_Nage_dev_output(1,nages);      //used in output. equals zero for first age
  matrix B(styr,endyr+1,1,nages);           //Population biomass by year and age at start of yr
  vector totB(styr,endyr+1);                //Total biomass by year
  vector totN(styr,endyr+1);                //Total abundance by year
  vector SSB(styr,endyr);                 //Total spawning biomass by year (female + male mature biomass) (proj yr ok bc of ssb on Jan1)
  vector MatFemB(styr,endyr);             //Total spawning biomass by year (mature female biomass) (proj yr ok bc of ssb on Jan1) 
  vector rec(styr,endyr+1);                 //Recruits by year
  vector prop_f(1,nages)  
  vector maturity_f(1,nages);               //Proportion of female mature at age
  vector reprod(1,nages);                   //vector used to compute spawning biomass (total mature female biomass)
  vector reprod2(1,nages); 

//---Stock-Recruit Function (Beverton-Holt, steepness parameterization)----------------------------------
  init_bounded_number log_R0(log_R0_LO,log_R0_HI,log_R0_PH);  //log(virgin Recruitment)
  vector log_R0_out(1,8);
  number R0;                                                  //virgin recruitment
  init_bounded_number steep(steep_LO,steep_HI,steep_PH);      //steepness
  vector steep_out(1,8);
  init_bounded_number rec_sigma(rec_sigma_LO,rec_sigma_HI,rec_sigma_PH);  //sd recruitment residuals
  vector rec_sigma_out(1,8);
  init_bounded_number R_autocorr(R_autocorr_LO,R_autocorr_HI,R_autocorr_PH);  //autocorrelation in SR 
  vector R_autocorr_out(1,8);

  number rec_sigma_sq;                        //square of rec_sigma      
  number rec_logL_add;                        //additive term in -logL term   
 
  init_bounded_dev_vector log_rec_dev(styr_rec_dev,endyr_rec_dev,log_rec_dev_LO,log_rec_dev_HI,log_rec_dev_PH);
  vector log_rec_dev_output(styr,endyr+1);             //used in t.series output. equals zero except for yrs in log_rec_dev
  vector log_rec_dev_out(styr_rec_dev,endyr_rec_dev);  //used in output for bound checking
  
  number var_rec_dev;                                //variance of log recruitment deviations, from yrs with unconstrainted S-R(XXXX-XXXX)
  number sigma_rec_dev;                              //sample SD of log residuals (may not equal rec_sigma 
                                             
  number BiasCor;                               //Bias correction in equilibrium recruits
  number S0;                                    //equal to spr_F0*R0 = virgin SSB
  number B0;                                    //equal to bpr_F0*R0 = virgin B  
  number R1;                                    //Recruits in styr
  number R_virgin;                              //unfished recruitment with bias correction
  vector SdS0(styr,endyr);                    //SSB / virgin SSB (projection yr possible bc of SSB on Jan1

//-------------------------------------------------------------------------------------------------------
////---Selectivity---------------------------------------------------- 

//Commercial handline-------------------------------------------------
  matrix sel_cH(styr,endyr,1,nages);
  init_bounded_number selpar_L50_cH1(selpar_L50_cH1_LO,selpar_L50_cH1_HI,selpar_L50_cH1_PH);
  init_bounded_number selpar_slope_cH1(selpar_slope_cH1_LO,selpar_slope_cH1_HI,selpar_slope_cH1_PH);
//  init_bounded_number selpar_L50_cH2(selpar_L50_cH2_LO,selpar_L50_cH2_HI,selpar_L50_cH2_PH);
//  init_bounded_number selpar_slope_cH2(selpar_slope_cH2_LO,selpar_slope_cH2_HI,selpar_slope_cH2_PH);
//  init_bounded_number selpar_L50_cH3(selpar_L50_cH3_LO,selpar_L50_cH3_HI,selpar_L50_cH3_PH);
//  init_bounded_number selpar_slope_cH3(selpar_slope_cH3_LO,selpar_slope_cH3_HI,selpar_slope_cH3_PH);
  vector selpar_L50_cH1_out(1,8);
  vector selpar_slope_cH1_out(1,8);
//  vector selpar_L50_cH2_out(1,8);
//  vector selpar_slope_cH2_out(1,8);
//  vector selpar_L50_cH3_out(1,8);
//  vector selpar_slope_cH3_out(1,8);

//SERFS chevron trap-------------------------------------------------
  matrix sel_mcvt(styr,endyr,1,nages);
  init_bounded_number selpar_L50_mcvt(selpar_L50_mcvt_LO,selpar_L50_mcvt_HI,selpar_L50_mcvt_PH);
  init_bounded_number selpar_slope_mcvt(selpar_slope_mcvt_LO,selpar_slope_mcvt_HI,selpar_slope_mcvt_PH);
  vector selpar_L50_mcvt_out(1,8);
  vector selpar_slope_mcvt_out(1,8);

//SEFIS video------------------------------------------------------
//  matrix sel_vid(styr,endyr,1,nages);
//  vector selpar_L50_vid_out(1,8);
//  vector selpar_slope_vid_out(1,8);

//Recreational (HB)-------------------------------------------------
  matrix sel_HB(styr,endyr,1,nages);
  init_bounded_number selpar_L51_HB1(selpar_L51_HB1_LO,selpar_L51_HB1_HI,selpar_L51_HB1_PH);
  init_bounded_number selpar_slope1_HB1(selpar_slope1_HB1_LO,selpar_slope1_HB1_HI,selpar_slope1_HB1_PH);
  init_bounded_number selpar_L52_HB1(selpar_L52_HB1_LO,selpar_L52_HB1_HI,selpar_L52_HB1_PH);
  init_bounded_number selpar_slope2_HB1(selpar_slope2_HB1_LO,selpar_slope2_HB1_HI,selpar_slope2_HB1_PH);

  init_bounded_number selpar_L50_HB2(selpar_L50_HB2_LO,selpar_L50_HB2_HI,selpar_L50_HB2_PH);
  init_bounded_number selpar_slope_HB2(selpar_slope_HB2_LO,selpar_slope_HB2_HI,selpar_slope_HB2_PH);
  init_bounded_number selpar_afull_HB2(selpar_afull_HB2_LO,selpar_afull_HB2_HI,selpar_afull_HB2_PH);
  init_bounded_number selpar_sigma_HB2(selpar_sigma_HB2_LO,selpar_sigma_HB2_HI,selpar_sigma_HB2_PH)

//  init_bounded_number selpar_L50_HB3(selpar_L50_HB3_LO,selpar_L50_HB3_HI,selpar_L50_HB3_PH);
//  init_bounded_number selpar_slope_HB3(selpar_slope_HB3_LO,selpar_slope_HB3_HI,selpar_slope_HB3_PH);
  
  vector selpar_L51_HB1_out(1,8);
  vector selpar_slope1_HB1_out(1,8);
  vector selpar_L52_HB1_out(1,8);
  vector selpar_slope2_HB1_out(1,8);

  vector selpar_L50_HB2_out(1,8);
  vector selpar_slope_HB2_out(1,8);
  vector selpar_afull_HB2_out(1,8);
  vector selpar_sigma_HB2_out(1,8);

//  vector selpar_L50_HB3_out(1,8);
//  vector selpar_slope_HB3_out(1,8);

//HB discard selectivities
  matrix sel_HB_D(styr,endyr,1,nages);
  init_bounded_number selpar_L50_HB_D(selpar_L50_HB_D_LO,selpar_L50_HB_D_HI,selpar_L50_HB_D_PH);
  init_bounded_number selpar_slope_HB_D(selpar_slope_HB_D_LO,selpar_slope_HB_D_HI,selpar_slope_HB_D_PH);
  init_bounded_number selpar_afull_HB_D(selpar_afull_HB_D_LO,selpar_afull_HB_D_HI,selpar_afull_HB_D_PH);
  init_bounded_number selpar_sigma_HB_D(selpar_sigma_HB_D_LO,selpar_sigma_HB_D_HI,selpar_sigma_HB_D_PH);

  vector selpar_L50_HB_D_out(1,8);
  vector selpar_slope_HB_D_out(1,8);
  vector selpar_afull_HB_D_out(1,8);
  vector selpar_sigma_HB_D_out(1,8);
 
//Recreational (GR)----------------------------------------
  matrix sel_GR(styr,endyr,1,nages);
  init_bounded_number selpar_L51_GR(selpar_L51_GR_LO,selpar_L51_GR_HI,selpar_L51_GR_PH);
  init_bounded_number selpar_slope1_GR(selpar_slope1_GR_LO,selpar_slope1_GR_HI,selpar_slope1_GR_PH);
  init_bounded_number selpar_L52_GR(selpar_L52_GR_LO,selpar_L52_GR_HI,selpar_L52_GR_PH);
  init_bounded_number selpar_slope2_GR(selpar_slope2_GR_LO,selpar_slope2_GR_HI,selpar_slope2_GR_PH);

  vector selpar_L51_GR_out(1,8);
  vector selpar_slope1_GR_out(1,8);
  vector selpar_L52_GR_out(1,8);
  vector selpar_slope2_GR_out(1,8);

//GR discard selectivity
  matrix sel_GR_D(styr,endyr,1,nages);
//  vector selpar_L50_GR_D_out(1,8); //set to HB_D selectivity
//  vector selpar_slope_GR_D_out(1,8); //set to HB_D selectivity
     
//Weighted total selectivity-----------------------------  
  //effort-weighted, recent selectivities
  vector sel_wgted_L(1,nages);  //toward landings 
  vector sel_wgted_D(1,nages);  //toward discards    
  vector sel_wgted_tot(1,nages);//toward Z, landings plus deads discards

//-------------------------------------------------------------------------------------------------------
//-------CPUE Predictions--------------------------------
  vector pred_cH_cpue(styr_cH_cpue,endyr_cH_cpue);        //predicted cH index (weight fish per effort)
  matrix N_cH(styr_cH_cpue,endyr_cH_cpue,1,nages);        //used to compute cH index
  vector pred_mcvt_cpue(styr_mcvt_cpue,endyr_mcvt_cpue);  //predicted SERFS chevron trap  index (weight fish per effort)
  matrix N_mcvt(styr_mcvt_cpue,endyr_mcvt_cpue,1,nages);  //used to compute SERFS chevron trap index
//  vector pred_vid_cpue(styr_vid_cpue,endyr_vid_cpue);   //predicted video index (weight fish per effort)
//  matrix N_vid(styr_vid_cpue,endyr_vid_cpue,1,nages);   //used to compute video index
  vector pred_HB_cpue(styr_HB_cpue,endyr_HB_cpue);        //predicted HB index (number fish per effort)
  matrix N_HB(styr_HB_cpue,endyr_HB_cpue,1,nages);        //used to compute HB index
  vector pred_GR_cpue(styr_GR_cpue,endyr_GR_cpue);        //predicted GR index (number fish per effort)
  matrix N_GR(styr_GR_cpue,endyr_GR_cpue,1,nages);        //used to compute GR index


//---Catchability (CPUE q's)----------------------------------------------------------
  init_bounded_number log_q_cH(log_q_cH_LO,log_q_cH_HI,log_q_cH_PH);
  init_bounded_number log_q_mcvt(log_q_mcvt_LO,log_q_mcvt_HI,log_q_mcvt_PH);  
//  init_bounded_number log_q_vid(log_q_vid_LO,log_q_vid_HI,log_q_vid_PH);  
  init_bounded_number log_q_HB(log_q_HB_LO,log_q_HB_HI,log_q_HB_PH);
  init_bounded_number log_q_GR(log_q_GR_LO,log_q_GR_HI,log_q_GR_PH);
  vector log_q_cH_out(1,8);
  vector log_q_mcvt_out(1,8);  
//  vector log_q_vid_out(1,8);   
  vector log_q_HB_out(1,8);
  vector log_q_GR_out(1,8);
  
//  init_bounded_number q_rate(0.001,0.1,set_q_rate_phase); //not estimated 
  number q_rate;
  vector q_rate_fcn_cH(styr_cH_cpue,endyr_cH_cpue);         //increase due to technology creep (saturates in 2003) 
  vector q_rate_fcn_mcvt(styr_mcvt_cpue,endyr_mcvt_cpue);   //increase due to technology creep (saturates in 2003)
//  vector q_rate_fcn_vid(styr_vid_cpue,endyr_vid_cpue);    //increase due to technology creep (saturates in 2003)  
  vector q_rate_fcn_HB(styr_HB_cpue,endyr_HB_cpue);         //increase due to technology creep (saturates in 2003) 
  vector q_rate_fcn_GR(styr_GR_cpue,endyr_GR_cpue);         //increase due to technology creep (saturates in 2003) 
  
//  init_bounded_number q_DD_beta(0.1,0.9,set_q_DD_phase);    //not estimated 
  number q_DD_beta;
  vector q_DD_fcn(styr,endyr);    //density dependent function as a multiple of q (scaled a la Katsukawa and Matsuda. 2003)
  number B0_q_DD;                 //B0 of ages q_DD_age plus
  vector B_q_DD(styr,endyr);      //annual biomass of ages q_DD_age plus

//Fishery dependent random walk catchability
//  init_bounded_vector q_RW_log_dev_HB(styr_HB_cpue,endyr_HB_cpue-1,-3.0,3.0,set_q_RW_phase); //not estimated in this model 
 vector q_RW_log_dev_cH(styr_cH_cpue,endyr_cH_cpue-1); 
 vector q_RW_log_dev_mcvt(styr_mcvt_cpue,endyr_mcvt_cpue-1); 
// vector q_RW_log_dev_vid(styr_vid_cpue,endyr_vid_cpue-1);  
 vector q_RW_log_dev_HB(styr_HB_cpue,endyr_HB_cpue-1); 
 vector q_RW_log_dev_GR(styr_GR_cpue,endyr_GR_cpue-1); 

//Fishery dependent catchability over time, may be constant
 vector q_cH(styr_cH_cpue,endyr_cH_cpue); 
 vector q_mcvt(styr_mcvt_cpue,endyr_mcvt_cpue) 
// vector q_vid(styr_vid_cpue,endyr_vid_cpue)  
 vector q_HB(styr_HB_cpue,endyr_HB_cpue); 
 vector q_GR(styr_GR_cpue,endyr_GR_cpue); 

//----------------------------------------------------------------------------------------------------------------------------------------------- 
//---Landings in numbers (total or 1000 fish) and in wgt (whole klb)--------------------------------------------------
  matrix L_cH_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cH_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  vector pred_cH_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages  
  vector pred_cH_L_klb(styr,endyr);      //yearly landings in 1000 lb whole summed over ages
  //vector obs_cH_L_wbias(styr,endyr);     //yearly landings observed, perhaps adjusted for multiplicitive bias

  matrix L_HB_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_HB_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age
  vector pred_HB_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages
  vector pred_HB_L_klb(styr,endyr);      //yearly landings in 1000 lb gutted summed over ages

  matrix L_GR_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_GR_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age
  vector pred_GR_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages
  vector pred_GR_L_klb(styr,endyr);      //yearly landings in 1000 lb gutted summed over ages

  matrix L_total_num(styr,endyr,1,nages);//total landings in number at age
  matrix L_total_klb(styr,endyr,1,nages);//landings in klb whole wgt at age 
  vector L_total_knum_yr(styr,endyr);    //total landings in 1000 fish by yr summed over ages  
  vector L_total_klb_yr(styr,endyr);     //total landings (klb whole wgt) by yr summed over ages
  
//---Discards (number dead fish) --------------------------------------------------  
  matrix D_HB_num(styr,endyr,1,nages);             //discards (numbers) at age
  matrix D_HB_klb(styr,endyr,1,nages);             //discards (1000 lb whole) at age
  vector pred_HB_D_knum(styr_HB_D,endyr_HB_D);     //yearly dead discards summed over ages
  vector obs_HB_D(styr_HB_D,endyr_HB_D);           //observed releases multiplied by discard mortality
  vector pred_HB_D_klb(styr_HB_D,endyr_HB_D);      //yearly dead discards in klb whole wgt summed over ages

  matrix D_GR_num(styr,endyr,1,nages);             //discards (numbers) at age
  matrix D_GR_klb(styr,endyr,1,nages);             //discards (1000 lb whole) at age
  vector pred_GR_D_knum(styr_GR_D,endyr_GR_D);     //yearly dead discards summed over ages
  vector obs_GR_D(styr_GR_D,endyr_GR_D);           //observed releases multiplied by discard mortality
  vector pred_GR_D_klb(styr_GR_D,endyr_GR_D);      //yearly dead discards in klb whole wgt summed over ages
  
  matrix D_total_num(styr,endyr,1,nages);          //total discards in number at age
  matrix D_total_klb(styr,endyr,1,nages);          //discards in klb whole wgt at age 
  vector D_total_knum_yr(styr,endyr);              //total discards in 1000 fish by yr summed over ages  
  vector D_total_klb_yr(styr,endyr);               //total discards (klb whole wgt) by yr summed over ages

  number Dmort_HB;
  number Dmort_GR;  

////---MSY calcs----------------------------------------------------------------------------
  number F_cH_prop;       //proportion of F_sum attributable to cH, last X=selpar_n_yrs_wgted yrs
  number F_HB_prop;       //proportion of F_sum attributable to HB, last X=selpar_n_yrs_wgted yrs 
  number F_GR_prop;       //proportion of F_sum attributable to GR, last X=selpar_n_yrs_wgted yrs
  number F_HB_D_prop;     //proportion of F_sum attributable to HB discards, last X=selpar_n_yrs_wgted yrs
  number F_GR_D_prop;     //proportion of F_sum attributable to GR discards, last X=selpar_n_yrs_wgted yrs
  
  number F_init_cH_prop;  //proportion of F_init attributable to cH, first X yrs, No discards in initial yrs
  number F_init_HB_prop;  //proportion of F_init attributable to HB, first X yrs
  number F_init_GR_prop;  //proportion of F_init attributable to GR, first X yrs
  number F_init_HB_D_prop;  //proportion of F_init attributable to HB discards, first X yrs
  number F_init_GR_D_prop;  //proportion of F_init attributable to GR discards, first X yrs

  number F_temp_sum;      //sum of geom mean Fsum's in last X yrs, used to compute F_fishery_prop

  vector F_end(1,nages);
  vector F_end_L(1,nages);  
  vector F_end_D(1,nages);    
  number F_end_apex;
  
  number SSB_msy_out;           //SSB (total mature biomass) at msy
  number F_msy_out;             //F at msy
  number msy_klb_out;           //max sustainable yield (1000 lb whole wgt)
  number msy_knum_out;          //max sustainable yield (1000 fish)  
  number D_msy_klb_out;         //discards associated with msy (1000 lb whole wgt)
  number D_msy_knum_out;        //discards associated with msy (1000 fish)  
  number B_msy_out;             //total biomass at MSY 
  number R_msy_out;             //equilibrium recruitment at F=Fmsy
  number spr_msy_out;           //spr at F=Fmsy

  number F20_dum;				//intermediate calculation for F20  added 1-19-16 here to...
  number F30_dum;				//intermediate calculation for F30
  number F40_dum;				//intermediate calculation for F40
  number F20_out;              	//F20
  number F30_out;              	//F30
  number F40_out;              	//F40
  number SSB_F30_out;  
  number B_F30_out;
  number R_F30_out;
  number L_F30_knum_out;
  number L_F30_klb_out;
  number D_F30_knum_out;
  number D_F30_klb_out;         //here
  
  vector N_age_msy(1,nages);         //numbers at age for MSY calculations: beginning of yr
  vector N_age_msy_spawn(1,nages);   //numbers at age for MSY calculations: time of peak spawning  
  vector L_age_msy(1,nages);         //landings at age for MSY calculations
  vector D_age_msy(1,nages);         //discard mortality (dead discards) at age for MSY calculations
  vector Z_age_msy(1,nages);         //total mortality at age for MSY calculations
  vector F_L_age_msy(1,nages);       //fishing mortality landings (not discards) at age for MSY calculations
  vector F_D_age_msy(1,nages);       //fishing mortality of discards at age for MSY calculations
  vector F_msy(1,n_iter_msy);        //values of full F to be used in equilibrium calculations
  vector spr_msy(1,n_iter_msy);      //reproductive capacity-per-recruit values corresponding to F values in F_msy
  vector R_eq(1,n_iter_msy);         //equilibrium recruitment values corresponding to F values in F_msy
  vector L_eq_klb(1,n_iter_msy);     //equilibrium landings(klb whole wgt) values corresponding to F values in F_msy
  vector L_eq_knum(1,n_iter_msy);    //equilibrium landings(1000 fish) values corresponding to F values in F_msy
  vector D_eq_klb(1,n_iter_msy);     //equilibrium discards(klb whole wgt) values corresponding to F values in F_msy
  vector D_eq_knum(1,n_iter_msy);    //equilibrium discards(1000 fish) values corresponding to F values in F_msy
  vector SSB_eq(1,n_iter_msy);       //equilibrium reproductive capacity values corresponding to F values in F_msy
  vector B_eq(1,n_iter_msy);         //equilibrium biomass values corresponding to F values in F_msy
  
  vector FdF_msy(styr,endyr);
  vector FdF30(styr,endyr);          //added 1-19-16
  vector SdSSB_msy(styr,endyr);    //(proj yr ok bc of ssb on Jan1)
  number SdSSB_msy_end;
  number FdF_msy_end;
  number FdF_msy_end_mean;           //geometric mean of last X yrs 
  vector SdSSB_F30(styr,endyr);	     //added 1-19-16 here to...
  vector Sdmsst_F30(styr,endyr);	 
  number SdSSB_F30_end;
  number Sdmsst_F30_end;
  number FdF30_end_mean;             //geometric mean of last selpar_n_yrs_wgted yrs    
  number Fend_mean_temp;	     //intermediate calc for geometric mean of last selpar_n_yrs_wgted yrs
  number Fend_mean;		     // here. geometric mean of last selpar_n_yrs_wgted yrs
 
  
  vector wgt_wgted_L_klb(1,nages);   //fishery-weighted average weight at age of landings in whole weight
  vector wgt_wgted_D_klb(1,nages);   //fishery-weighted average weight at age of discards in whole weight  
  number wgt_wgted_L_denom;          //used in intermediate calculations
  number wgt_wgted_D_denom;          //used in intermediate calculations

  number iter_inc_msy;               //increments used to compute msy, equals 1/(n_iter_msy-1)
  
////--------Mortality------------------------------------------------------------------

// Stuff immediately below used only if M is estimated
//  //init_bounded_number M_constant(0.1,0.2,1);                         //age-indpendent: used only for MSST
//  vector Mscale_ages(1,max_obs_age);
//  vector Mscale_len(1,max_obs_age);
//  vector Mscale_wgt_g(1,max_obs_age); 
//  vector M_lorenzen(1,max_obs_age);   
//  number cum_surv_1plus;  

  vector M(1,nages);                         //age-dependent natural mortality
  init_bounded_number M_constant(M_constant_LO,M_constant_HI,M_constant_PH);   //age-independent: used only for MSST
  vector M_constant_out(1,8);
  number smsy2msst;                           //scales Smsy to get msst using (1-M). Used only in output.
  number smsy2msst75;                         //scales Smsy to get msst using 75%. Used only in output.  
  
  matrix F(styr,endyr,1,nages);
  vector Fsum(styr,endyr);                   //Full fishing mortality rate by year
  vector Fapex(styr,endyr);                  //Max across ages, fishing mortality rate by year (may differ from Fsum bc of dome-shaped sel) 
  matrix Z(styr,endyr,1,nages);

  init_bounded_number log_avg_F_cH(log_avg_F_cH_LO,log_avg_F_cH_HI,log_avg_F_cH_PH);
  vector log_avg_F_cH_out(1,8);  
  init_bounded_dev_vector log_F_dev_cH(styr_cH_L,endyr_cH_L,log_F_dev_cH_LO,log_F_dev_cH_HI,log_F_dev_cH_PH);
  vector log_F_dev_cH_out(styr_cH_L,endyr_cH_L);
  matrix F_cH(styr,endyr,1,nages);
  vector F_cH_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cH;
  number log_F_dev_end_cH; 

  init_bounded_number log_avg_F_HB(log_avg_F_HB_LO,log_avg_F_HB_HI,log_avg_F_HB_PH);
  vector log_avg_F_HB_out(1,8); 
  init_bounded_dev_vector log_F_dev_HB(styr_HB_L,endyr_HB_L,log_F_dev_HB_LO,log_F_dev_HB_HI,log_F_dev_HB_PH);    
  vector log_F_dev_HB_out(styr_HB_L,endyr_HB_L);
  matrix F_HB(styr,endyr,1,nages);
  vector F_HB_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_HB;    
  number log_F_dev_end_HB;  

  init_bounded_number log_avg_F_GR(log_avg_F_GR_LO,log_avg_F_GR_HI,log_avg_F_GR_PH);
  vector log_avg_F_GR_out(1,8); 
  init_bounded_dev_vector log_F_dev_GR(styr_GR_L,endyr_GR_L,log_F_dev_GR_LO,log_F_dev_GR_HI,log_F_dev_GR_PH);    
  vector log_F_dev_GR_out(styr_GR_L,endyr_GR_L);
  matrix F_GR(styr,endyr,1,nages);
  vector F_GR_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_GR;    
  number log_F_dev_end_GR;  

  init_bounded_number log_avg_F_HB_D(log_avg_F_HB_D_LO,log_avg_F_HB_D_HI,log_avg_F_HB_D_PH);
  vector log_avg_F_HB_D_out(1,8);  
  init_bounded_dev_vector log_F_dev_HB_D(styr_HB_D,endyr_HB_D,log_F_dev_HB_D_LO,log_F_dev_HB_D_HI,log_F_dev_HB_D_PH);
  vector log_F_dev_HB_D_out(styr_HB_D,endyr_HB_D);
  matrix F_HB_D(styr,endyr,1,nages);
  vector F_HB_D_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_HB_D; 
  number log_F_dev_end_HB_D; 

  init_bounded_number log_avg_F_GR_D(log_avg_F_GR_D_LO,log_avg_F_GR_D_HI,log_avg_F_GR_D_PH);
  vector log_avg_F_GR_D_out(1,8);  
  init_bounded_dev_vector log_F_dev_GR_D(styr_GR_D,endyr_GR_D,log_F_dev_GR_D_LO,log_F_dev_GR_D_HI,log_F_dev_GR_D_PH);
  vector log_F_dev_GR_D_out(styr_GR_D,endyr_GR_D);
  matrix F_GR_D(styr,endyr,1,nages);
  vector F_GR_D_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_GR_D;  
  number log_F_dev_end_GR_D; 

  init_bounded_number F_init(F_init_LO,F_init_HI,F_init_PH); //scales early F for initialization
  vector F_init_out(1,8); 
  number F_init_denom;  //interim calculation
//  vector sel_initial(1,nages);       //initial selectivity (a combination of for-hire and commercial selectivities)

//---Per-recruit stuff----------------------------------------------------------------------------------
  vector N_age_spr(1,nages);         //numbers at age for SPR calculations: beginning of year
  vector N_age_spr_spawn(1,nages);   //numbers at age for SPR calculations: time of peak spawning  
  vector L_age_spr(1,nages);         //catch at age for SPR calculations
  vector Z_age_spr(1,nages);         //total mortality at age for SPR calculations
  vector spr_static(styr,endyr);     //vector of static SPR values by year
  vector F_L_age_spr(1,nages);       //fishing mortality of landings (not discards) at age for SPR calculations
  vector F_spr(1,n_iter_spr);        //values of full F to be used in per-recruit calculations
  vector spr_spr(1,n_iter_spr);      //reproductive capacity-per-recruit values corresponding to F values in F_spr
  vector spr_ratio(1,n_iter_spr);    // added 1-19-16  reproductive capacity-per-recruit relative to spr_F0 values corresponding to F values in F_spr
  vector L_spr(1,n_iter_spr);        //landings(lb whole)-per-recruit (ypr) values corresponding to F values in F_spr

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
 
  number sdnr_lc_cH;
  number sdnr_lc_mcvt;
  number sdnr_lc_HB;
  number sdnr_lc_HB_D;  
  number sdnr_lc_GR;  
 
  number sdnr_ac_cH;
  number sdnr_ac_mcvt;
  number sdnr_ac_HB;
  number sdnr_ac_GR;
    
  number sdnr_I_cH;
  number sdnr_I_mcvt;  
//  number sdnr_I_vid;  
  number sdnr_I_HB;
  number sdnr_I_GR;
        
////-------Objective function components-----------------------------------------------------------------------------
  number w_L;
  number w_D;
  
  number w_I_cH;
  number w_I_mcvt;  
//  number w_I_vid; 
  number w_I_HB;
  number w_I_GR;
   
  number w_lc_cH;
  number w_lc_mcvt;
  number w_lc_HB;
  number w_lc_HB_D;  
  number w_lc_GR;  

  number w_ac_cH; 
  number w_ac_mcvt;
  number w_ac_HB;
  //number w_ac_GR;
  
  number w_Nage_init;  
  number w_rec;
  number w_rec_early;
  number w_rec_end;
  number w_fullF;  
  number w_Ftune;

  number f_cH_L; 
  number f_HB_L; 
  number f_GR_L; 

  number f_HB_D; 
  number f_GR_D; 

  number f_cH_cpue;
  number f_mcvt_cpue;  
//  number f_vid_cpue;  
  number f_HB_cpue;
  number f_GR_cpue;
   
  number f_cH_lenc;
  number f_mcvt_lenc;
  number f_HB_lenc;
  number f_HB_D_lenc;  
  number f_GR_lenc;  

  number f_cH_agec;
  number f_mcvt_agec;
  number f_HB_agec;       
  //number f_GR_agec;

//  Penalties and constraints. Not all are used.
  number f_Nage_init;              //weight on log devs to estimate initial abundance (excluding first age)
  number f_rec_dev;                //weight on recruitment deviations to fit S-R curve
  number f_rec_dev_early;          //extra weight on deviations in first recruitment stanza
  number f_rec_dev_end;            //extra weight on deviations in ending recruitment stanza
  number f_fullF_constraint;       //penalty for Fapex>X
  number f_Ftune;                  //penalty for tuning F in Ftune yr.  Not applied in final optimization phase.
  number f_priors;                 //prior information on parameters
  
  //init_number xdum;
  objective_function_value fval;
  number fval_data;
  number grad_max;
  
//--Dummy variables ----
  number denom;                   //denominator used in some calculations
  number numer;                   //numerator used in some calculations
   
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//INITIALIZATION_SECTION


//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
GLOBALS_SECTION
  #include "admodel.h"     // Include AD class definitions
  #include "admb2r.cpp"    // Include S-compatible output functions (needs preceding)
  #include <time.h>
	time_t start,finish;
	long hour,minute,second;	
	double elapsed_time;
	
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
RUNTIME_SECTION
 maximum_function_evaluations 1000, 2000, 3000, 10000;
 convergence_criteria 1e-2, 1e-2, 1e-3, 1e-4;
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PRELIMINARY_CALCS_SECTION

// Set values of fixed parameters or set initial guess of estimated parameters
  Dmort_HB=set_Dmort_HB;
  obs_HB_D=Dmort_HB*obs_HB_released;
  Dmort_GR=set_Dmort_GR;
  obs_GR_D=Dmort_GR*obs_GR_released;

  Linf=set_Linf(1);
  K=set_K(1);
  t0=set_t0(1);
  len_cv_val=set_len_cv(1);

  // FISHERY LANDINGS Growth
  Linf_L=set_Linf_L(1);
  K_L=set_K_L(1);
  t0_L=set_t0_L(1);
  len_cv_val_L=set_len_cv_L(1);

  // FISHERY INDEPENDENT (SERFS chevron trap) Growth
  Linf_mcvt=set_Linf_mcvt(1);
  K_mcvt=set_K_mcvt(1);
  t0_mcvt=set_t0_mcvt(1);
  len_cv_val_mcvt=set_len_cv_mcvt(1);

  M=set_M; 
  M_constant=set_M_constant(1);
  smsy2msst=1.0-M_constant;
  smsy2msst75=0.75;  
//  for (iage=1;iage<=max_obs_age;iage++){Mscale_ages(iage)=iage;}
  
  log_R0=set_log_R0(1);
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);
  
  log_q_cH=set_log_q_cH(1);
  log_q_mcvt=set_log_q_mcvt(1);  
//  log_q_vid=set_log_q_vid(1);  
  log_q_HB=set_log_q_HB(1);
  log_q_GR=set_log_q_GR(1);
  
  q_rate=set_q_rate;
  q_rate_fcn_cH=1.0;   
  q_rate_fcn_mcvt=1.0;   
//  q_rate_fcn_vid=1.0;   
  q_rate_fcn_HB=1.0;   
  q_rate_fcn_GR=1.0; 
  q_DD_beta=set_q_DD_beta;
  q_DD_fcn=1.0;

  q_RW_log_dev_cH.initialize(); 
  q_RW_log_dev_mcvt.initialize();  
  q_RW_log_dev_HB.initialize(); 
  q_RW_log_dev_GR.initialize(); 
      
  if (set_q_rate_phase<0 & q_rate!=0.0)
  {
    for (iyear=styr_cH_cpue; iyear<=endyr_cH_cpue; iyear++)
      {   if (iyear>styr_cH_cpue & iyear <=2003) 
          {//q_rate_fcn_cH(iyear)=(1.0+q_rate)*q_rate_fcn_cH(iyear-1); //compound
             q_rate_fcn_cH(iyear)=(1.0+(iyear-styr_cH_cpue)*q_rate)*q_rate_fcn_cH(styr_cH_cpue);  //linear
          }
          if (iyear>2003) {q_rate_fcn_cH(iyear)=q_rate_fcn_cH(iyear-1);} 
      }   
    for (iyear=styr_HB_cpue; iyear<=endyr_HB_cpue; iyear++)
      {   if (iyear>styr_HB_cpue & iyear <=2003) 
          {//q_rate_fcn_HB(iyear)=(1.0+q_rate)*q_rate_fcn_HB(iyear-1); //compound
             q_rate_fcn_HB(iyear)=(1.0+(iyear-styr_HB_cpue)*q_rate)*q_rate_fcn_HB(styr_HB_cpue);  //linear
          }
          if (iyear>2003) {q_rate_fcn_HB(iyear)=q_rate_fcn_HB(iyear-1);} 
      }   
    for (iyear=styr_GR_cpue; iyear<=endyr_GR_cpue; iyear++)
      {   if (iyear>styr_GR_cpue & iyear <=2003) 
          {//q_rate_fcn_GR(iyear)=(1.0+q_rate)*q_rate_fcn_GR(iyear-1); //compound
             q_rate_fcn_GR(iyear)=(1.0+(iyear-styr_GR_cpue)*q_rate)*q_rate_fcn_GR(styr_GR_cpue);  //linear
          }
          if (iyear>2003) {q_rate_fcn_GR(iyear)=q_rate_fcn_GR(iyear-1);} 
      }   
  } //end q_rate conditional      

  w_L=set_w_L;
  w_D=set_w_D;
  
  w_I_cH=set_w_I_cH;
  w_I_mcvt=set_w_I_mcvt;
//  w_I_vid=set_w_I_vid;
  w_I_HB=set_w_I_HB;
  w_I_GR=set_w_I_GR;

  w_lc_cH=set_w_lc_cH;
  w_lc_mcvt=set_w_lc_mcvt;
  w_lc_HB=set_w_lc_HB;
  w_lc_HB_D=set_w_lc_HB_D;  
  w_lc_GR=set_w_lc_GR;   
  
  w_ac_cH=set_w_ac_cH;
  w_ac_mcvt=set_w_ac_mcvt;
  w_ac_HB=set_w_ac_HB;
  //w_ac_GR=set_w_ac_GR;

  w_Nage_init=set_w_Nage_init;
  w_rec=set_w_rec;
  w_rec_early=set_w_rec_early;
  w_rec_end=set_w_rec_end;
  w_fullF=set_w_fullF;
  w_Ftune=set_w_Ftune;

  F_init=set_F_init(1);

  log_avg_F_cH=set_log_avg_F_cH(1);
  log_avg_F_HB=set_log_avg_F_HB(1); 
  log_avg_F_GR=set_log_avg_F_GR(1); 
  log_avg_F_HB_D=set_log_avg_F_HB_D(1); 
  log_avg_F_GR_D=set_log_avg_F_GR_D(1); 
    
  log_F_dev_cH=set_log_F_dev_cH_vals;
  log_F_dev_HB=set_log_F_dev_HB_vals;
  log_F_dev_GR=set_log_F_dev_GR_vals;
  log_F_dev_HB_D=set_log_F_dev_HB_D_vals;
  log_F_dev_GR_D=set_log_F_dev_GR_D_vals;
 
  selpar_L50_cH1=set_selpar_L50_cH1(1);
  selpar_slope_cH1=set_selpar_slope_cH1(1);

  selpar_L50_mcvt=set_selpar_L50_mcvt(1);
  selpar_slope_mcvt=set_selpar_slope_mcvt(1);

  selpar_L51_HB1=set_selpar_L51_HB1(1);
  selpar_slope1_HB1=set_selpar_slope1_HB1(1);
  selpar_L52_HB1=set_selpar_L52_HB1(1);
  selpar_slope2_HB1=set_selpar_slope2_HB1(1);

  selpar_L50_HB2=set_selpar_L50_HB2(1);
  selpar_slope_HB2=set_selpar_slope_HB2(1);
  selpar_afull_HB2=set_selpar_afull_HB2(1);
  selpar_sigma_HB2=set_selpar_sigma_HB2(1);
  
// selpar_L50_HB3=set_selpar_L50_HB3(1);
// selpar_slope_HB3=set_selpar_slope_HB3(1);

  selpar_L50_HB_D=set_selpar_L50_HB_D(1);
  selpar_slope_HB_D=set_selpar_slope_HB_D(1);
  selpar_afull_HB_D=set_selpar_afull_HB_D(1);
  selpar_sigma_HB_D=set_selpar_sigma_HB_D(1);

  selpar_L51_GR=set_selpar_L51_GR(1);    
  selpar_slope1_GR=set_selpar_slope1_GR(1);  
  selpar_L52_GR=set_selpar_L52_GR(1);  
  selpar_slope2_GR=set_selpar_slope2_GR(1);  

 sqrt2pi=sqrt(2.*3.14159265);
 g2mt=0.000001;         //conversion of grams to metric tons
 g2kg=0.001;            //conversion of grams to kg 
 mt2klb=2.20462;        //conversion of metric tons to 1000 lb 
 mt2lb=mt2klb*1000.0;   //conversion of metric tons to lb
 g2klb=g2mt*mt2klb;     //conversion of grams to 1000 lb 
 dzero=0.00001;         
 huge_number=1.0e+10;   
 
 SSB_msy_out=0.0;

 iter_inc_msy=max_F_spr_msy/(n_iter_msy-1);
 iter_inc_spr=max_F_spr_msy/(n_iter_spr-1); 

 maturity_f=maturity_f_obs;
 prop_f=prop_f_obs;

 //lbins=lenbins; //NOT NEEDED
 
//Fill in sample sizes of comps, possibly sampled in nonconsec yrs 
//Used primarily for output in R object   

      nsamp_cH_lenc_allyr=missing;
      nsamp_mcvt_lenc_allyr=missing;
      nsamp_HB_lenc_allyr=missing;
      nsamp_HB_D_lenc_allyr=missing;  
      nsamp_GR_lenc_allyr=missing;  
      nsamp_cH_agec_allyr=missing;
      nsamp_mcvt_agec_allyr=missing;
      nsamp_HB_agec_allyr=missing;
      //nsamp_GR_agec_allyr=missing;

      nfish_cH_lenc_allyr=missing;
      nfish_mcvt_lenc_allyr=missing;
      nfish_HB_lenc_allyr=missing;
      nfish_HB_D_lenc_allyr=missing; 
      nfish_GR_lenc_allyr=missing;   
      nfish_cH_agec_allyr=missing;
      nfish_mcvt_agec_allyr=missing;
      nfish_HB_agec_allyr=missing;
      //nfish_GR_agec_allyr=missing;

      for (iyear=1; iyear<=nyr_cH_lenc; iyear++)
         {if (nsamp_cH_lenc(iyear)>=minSS_cH_lenc)
           {nsamp_cH_lenc_allyr(yrs_cH_lenc(iyear))=nsamp_cH_lenc(iyear);
            nfish_cH_lenc_allyr(yrs_cH_lenc(iyear))=nfish_cH_lenc(iyear);}}

      for (iyear=1; iyear<=nyr_mcvt_lenc; iyear++)
         {if (nsamp_mcvt_lenc(iyear)>=minSS_mcvt_lenc)
           {nsamp_mcvt_lenc_allyr(yrs_mcvt_lenc(iyear))=nsamp_mcvt_lenc(iyear);
            nfish_mcvt_lenc_allyr(yrs_mcvt_lenc(iyear))=nfish_mcvt_lenc(iyear);}}

      for (iyear=1; iyear<=nyr_HB_lenc; iyear++)
         {if (nsamp_HB_lenc(iyear)>=minSS_HB_lenc)
            {nsamp_HB_lenc_allyr(yrs_HB_lenc(iyear))=nsamp_HB_lenc(iyear);
             nfish_HB_lenc_allyr(yrs_HB_lenc(iyear))=nfish_HB_lenc(iyear);}}

      for (iyear=1; iyear<=nyr_HB_D_lenc; iyear++)
         {if (nsamp_HB_D_lenc(iyear)>=minSS_HB_D_lenc)  
            {nsamp_HB_D_lenc_allyr(yrs_HB_D_lenc(iyear))=nsamp_HB_D_lenc(iyear);  
             nfish_HB_D_lenc_allyr(yrs_HB_D_lenc(iyear))=nfish_HB_D_lenc(iyear);}}  

      for (iyear=1; iyear<=nyr_GR_lenc; iyear++)  
         {if (nsamp_GR_lenc(iyear)>=minSS_GR_lenc)  
            {nsamp_GR_lenc_allyr(yrs_GR_lenc(iyear))=nsamp_GR_lenc(iyear);  
             nfish_GR_lenc_allyr(yrs_GR_lenc(iyear))=nfish_GR_lenc(iyear);}}  

      for (iyear=1; iyear<=nyr_cH_agec; iyear++)
         {if (nsamp_cH_agec(iyear)>=minSS_cH_agec)
           {nsamp_cH_agec_allyr(yrs_cH_agec(iyear))=nsamp_cH_agec(iyear);
            nfish_cH_agec_allyr(yrs_cH_agec(iyear))=nfish_cH_agec(iyear);}}

      for (iyear=1; iyear<=nyr_mcvt_agec; iyear++)
         {if (nsamp_mcvt_agec(iyear)>=minSS_mcvt_agec)
           {nsamp_mcvt_agec_allyr(yrs_mcvt_agec(iyear))=nsamp_mcvt_agec(iyear);
            nfish_mcvt_agec_allyr(yrs_mcvt_agec(iyear))=nfish_mcvt_agec(iyear);}}

      for (iyear=1; iyear<=nyr_HB_agec; iyear++)
         {if (nsamp_HB_agec(iyear)>=minSS_HB_agec)
           {nsamp_HB_agec_allyr(yrs_HB_agec(iyear))=nsamp_HB_agec(iyear);
            nfish_HB_agec_allyr(yrs_HB_agec(iyear))=nfish_HB_agec(iyear);}} 

//      for (iyear=1; iyear<=nyr_GR_agec; iyear++)
//         {if (nsamp_GR_agec(iyear)>=minSS_GR_agec)
//           {nsamp_GR_agec_allyr(yrs_GR_agec(iyear))=nsamp_GR_agec(iyear);
//            nfish_GR_agec_allyr(yrs_GR_agec(iyear))=nfish_GR_agec(iyear);}}
             
//fill in Fs for msy and per-recruit analyses
  F_msy(1)=0.0;  
  for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
  for (ff=2;ff<=n_iter_spr;ff++) {F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}


//fill in F's, Catch matrices, and log rec dev with zero's
  F_cH.initialize(); L_cH_num.initialize();
  F_HB.initialize(); L_HB_num.initialize();
  F_GR.initialize(); L_GR_num.initialize();
  F_HB_D.initialize(); D_HB_num.initialize();
  F_GR_D.initialize(); D_GR_num.initialize();

  F_cH_out.initialize();
  F_HB_out.initialize();
  F_GR_out.initialize();
  F_HB_D_out.initialize();
  F_GR_D_out.initialize();

      
  sel_cH.initialize();
  sel_mcvt.initialize();
//  sel_vid.initialize();
  sel_HB.initialize();
  sel_HB_D.initialize();
  sel_HB.initialize();  
  
  log_rec_dev_output.initialize();  
  log_rec_dev=set_log_rec_dev_vals;
  log_Nage_dev_output.initialize();
  log_Nage_dev=set_log_Nage_dev_vals;
 
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
TOP_OF_MAIN_SECTION
  time(&start);
  arrmblsize=20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(1600);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(10000);

//>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PROCEDURE_SECTION
  
 //cout<<"start"<<endl;
 
 //get_M_at_age(); //Needed only if M is estimated
 
 get_length_weight_at_age(); 
 //cout << "got length, weight, fecundity transitions" << endl;
 get_reprod();
 //cout << "got repro stuff" << endl;
 get_length_at_age_dist(); 
 //cout << "got predicted length at age distribution" << endl;
 get_weight_at_age_landings();
 //cout << "got weight at age of landings" << endl; 
 get_spr_F0();
 //cout << "got F0 spr" << endl;
 get_selectivity(); 
 //cout << "got selectivity" << endl;
 get_mortality(); 
 //cout << "got mortalities" << endl;
 get_bias_corr(); 
 //cout << "got recruitment bias correction" << endl;
 get_numbers_at_age(); 
 //cout << "got numbers at age" << endl;
 get_landings_numbers();
 //cout << "got landings in numbers" << endl;
 get_landings_wgt();
 //cout << "got landings in wgt" << endl;
 get_dead_discards(); 
//   cout << "got dead discards in num and wgt" << endl;
//   cout << pred_GR_D_knum << endl;
 get_catchability_fcns(); 
 //cout << "got catchability_fcns" << endl;
 get_indices();
 //cout << "got indices" << endl;
 get_length_comps();
 //cout << "got length comps" << endl;
 get_age_comps();
 //cout << "got age comps" << endl;
 evaluate_objective_function();
 //cout << "objective function calculations complete" << endl;
 
 
FUNCTION get_length_weight_at_age

  //compute mean length (mm TL) and weight (whole) at age
    meanlen_FL=Linf*(1.0-mfexp(-K*(agebins-t0+0.5)));     //total length in mm

    wgt_kg=wgtpar_a*pow(meanlen_FL,wgtpar_b);             //whole wgt in kg 
    wgt_g=wgt_kg/g2kg;                                    //convert wgt in kg to weight in g    
    wgt_mt=wgt_g*g2mt;                                    //convert weight in g to weight in mt
    wgt_klb=mt2klb*wgt_mt;                                //1000 lb of whole wgt
    wgt_lb=mt2lb*wgt_mt;                                  //lb of whole wgt

    fecundity=(fecpar_a+fecpar_b*meanlen_FL)/fecpar_scale; //annual egg production of a mature female at age in units of fecpar_scale

  //THESE CALCULATIONS ASSUME THE FISHERY LANDINGS GROWTH CURVE
    meanlen_FL_L=Linf_L*(1.0-mfexp(-K_L*(agebins-t0_L+0.5))); //fork length in mm
    wgt_kg_L=wgtpar_a*pow(meanlen_FL_L,wgtpar_b);             //wgt in kg KC changed from wgt in metric tons (wgt_mt) because L-W relationship is mm FL and kg whole weight
    wgt_g_L=wgt_kg_L/g2kg;                                    //KC convert wgt in kg from L-W relationship to weight in g    
    wgt_mt_L=wgt_g_L*g2mt;                                    //KC convert weight in g to weight in mt
    wgt_klb_L=mt2klb*wgt_mt_L;                                //1000 lb of whole wgt
    wgt_lb_L=mt2lb*wgt_mt_L;                                  //1000 lb of whole wgt

  //THESE CALCULATIONS ASSUME THE FISHERY (SERFS chevron trap) GROWTH CURVE
    meanlen_FL_mcvt=Linf_mcvt*(1.0-mfexp(-K_mcvt*(agebins-t0_mcvt+0.5)));    //fork length in mm


FUNCTION get_reprod 

   //reprod is product of stuff going into reproductive capacity calcs
//  reprod=elem_prod(elem_prod(prop_f,maturity_f),fecundity);  
   reprod=elem_prod(elem_prod(elem_prod(prop_f,maturity_f),fecundity),fecpar_batches); 
   reprod2=elem_prod(elem_prod(prop_f,maturity_f),wgt_mt);

 
FUNCTION get_length_at_age_dist

  //compute matrix of length at age, based on the normal distribution
  
 dvar_vector length(1,nages);
  dvariable cvlen;
  
  if (use_landings_growth==1.0)														//LC added to allow for use of Landings Growth Curve
  	{
  	length=meanlen_FL_L;
  	cvlen=len_cv_val_L;
//        sdlen=len_sd_val_L;
        len_sd=len_sd_L;
  	}
   if (use_landings_growth==0.0) 	
  	{
  	length=meanlen_FL;
 	cvlen=len_cv_val;
 //     sdlen=len_sd_val;
        len_sd=len_sd;
	}
  
  for (iage=1;iage<=nages;iage++)
  {
    len_cv(iage)=cvlen; 
    len_sd(iage)=length(iage)*len_cv(iage); 
//      len_sd(iage)=sdlen;
//      len_cv(iage)=len_sd(iage)/length(iage);

    zscore_lzero=(0.0-length(iage))/len_sd(iage);  
    cprob_lzero=cumd_norm(zscore_lzero); 
    
    //first length bin
    zscore_len=((lenbins(1)+0.5*lenbins_width)-length(iage)) / len_sd(iage);  
    cprob_lenvec(1)=cumd_norm(zscore_len);          //includes any probability mass below zero
    lenprob(iage,1)=cprob_lenvec(1)-cprob_lzero;    //removes any probability mass below zero
         
    //most other length bins     
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len=((lenbins(ilen)+0.5*lenbins_width)-length(iage)) / len_sd(iage);
        cprob_lenvec(ilen)=cumd_norm(zscore_len);
        lenprob(iage,ilen)=cprob_lenvec(ilen)-cprob_lenvec(ilen-1);
      }
    //last length bin is a plus group
    zscore_len=((lenbins(nlenbins)-0.5*lenbins_width)-length(iage)) / len_sd(iage);
    lenprob(iage,nlenbins)=1.0-cumd_norm(zscore_len);
  
    lenprob(iage)=lenprob(iage)/(1.0-cprob_lzero);  //renormalize to account for any prob mass below size=0
  }
  //fleet and survey specific length probs, all assumed here to equal the popn
  lenprob_cH=lenprob;
//KC  lenprob_mcvt=lenprob;
  lenprob_HB=lenprob;
  lenprob_GR=lenprob;
  lenprob_HB_D=lenprob; 

    ///////////// THIS SECTION FOR SERFS chevron trap (mcvt)  ////////////////////////// 
  lenprob2.initialize();
  length=meanlen_FL_mcvt;
  cvlen=len_cv_val_mcvt;
//  sdlen=len_sd_val_mcvt;
  len_sd=len_sd_mcvt;
 
  for (iage=1;iage<=nages;iage++)
    {
    len_cv(iage)=cvlen;															
    len_sd(iage)=length(iage)*len_cv(iage);											
//      len_sd(iage)=sdlen;
//      len_cv(iage)=len_sd(iage)/length(iage);

      zscore_lzero2=(0.0-length(iage))/len_sd(iage);
      cprob_lzero2=cumd_norm(zscore_lzero2);
      
      zscore_len2=((lenbins(1)+0.5*lenbins_width)-length(iage)) / len_sd(iage);  					
      cprob_lenvec2(1)=cumd_norm(zscore_len2);         
      lenprob2(iage,1)=cprob_lenvec2(1)-cprob_lzero2;   
           
      //most other length bins     
      for (ilen=2;ilen<nlenbins;ilen++)  
        {
          zscore_len2=((lenbins(ilen)+0.5*lenbins_width)-length(iage)) / len_sd(iage);  				
          cprob_lenvec2(ilen)=cumd_norm(zscore_len2);
          lenprob2(iage,ilen)=cprob_lenvec2(ilen)-cprob_lenvec2(ilen-1);
        }
      //last length bin is a plus group
      zscore_len2=((lenbins(nlenbins)-0.5*lenbins_width)-length(iage)) / len_sd(iage); 
          lenprob2(iage,nlenbins)=1.0-cumd_norm(zscore_len2); 
      lenprob2(iage)=lenprob2(iage)/(1.0-cprob_lzero2);  //renormalize to account for any prob mass below size=0
    }
    //fleet and survey specific length probs
    lenprob_mcvt=lenprob2;
    
FUNCTION get_weight_at_age_landings
 
  //KC5--modified based on Lew code
   if(use_landings_growth==1.0){ 

  for (iyear=styr; iyear<=endyr; iyear++)
  {
    len_cH_mm(iyear)=meanlen_FL_L;
    wholewgt_cH_klb(iyear)=wgt_klb_L; //wholeweight used to match index 
    len_HB_mm(iyear)=meanlen_FL_L;
    wholewgt_HB_klb(iyear)=wgt_klb_L; //KC6  change _out to _Lc
    len_GR_mm(iyear)=meanlen_FL_L;
    wholewgt_GR_klb(iyear)=wgt_klb_L; //KC6   change _out to _L
 
    len_HB_D_mm(iyear)=meanlen_FL_L;
    wholewgt_HB_D_klb(iyear)=wgt_klb_L;  //KC6; change _out to _L
    len_GR_D_mm(iyear)=meanlen_FL_L;
    wholewgt_GR_D_klb(iyear)=wgt_klb_L;  //KC6  change _out to _L
  } 

  }
 
   if(use_landings_growth==0.0){ 

  for (iyear=styr; iyear<=endyr; iyear++)
  {
    len_cH_mm(iyear)=meanlen_FL; 
    wholewgt_cH_klb(iyear)=wgt_klb; //wholeweight used to match index  
    len_HB_mm(iyear)=meanlen_FL;
    wholewgt_HB_klb(iyear)=wgt_klb;  //KC6  KC6-got rid _out
    len_GR_mm(iyear)=meanlen_FL;
    wholewgt_GR_klb(iyear)=wgt_klb;  //KC6  KC6-got rid _out

    len_HB_D_mm(iyear)=meanlen_FL;
    wholewgt_HB_D_klb(iyear)=wgt_klb;  //KC6  KC6-got rid _out
    len_GR_D_mm(iyear)=meanlen_FL;
    wholewgt_GR_D_klb(iyear)=wgt_klb;  //KC6  KC6-got rid _out
  } 

  } 
 

FUNCTION get_spr_F0

  //at mdyr, apply half this yr's mortality, half next yr's
  N_spr_F0(1)=1.0*mfexp(-1.0*M(1)*spawn_time_frac); //at peak spawning time
  N_bpr_F0(1)=1.0;      //at start of year
  for (iage=2; iage<=nages; iage++)
  {
    //N_spr_F0(iage)=N_spr_F0(iage-1)*mfexp(-1.0*(M(iage-1));
    N_spr_F0(iage)=N_spr_F0(iage-1)*mfexp(-1.0*(M(iage-1)*(1.0-spawn_time_frac) + M(iage)*spawn_time_frac)); 
    N_bpr_F0(iage)=N_bpr_F0(iage-1)*mfexp(-1.0*(M(iage-1)));    
  }
  N_spr_F0(nages)=N_spr_F0(nages)/(1.0-mfexp(-1.0*M(nages))); //plus group (sum of geometric series)
  N_bpr_F0(nages)=N_bpr_F0(nages)/(1.0-mfexp(-1.0*M(nages)));
  
  spr_F0=sum(elem_prod(N_spr_F0,reprod));  
  bpr_F0=sum(elem_prod(N_bpr_F0,wgt_mt));    


FUNCTION get_selectivity

  //BLOCK 1 for selex. 
  for (iyear=styr; iyear<=endyr_selex_phase1; iyear++)
  {     
     sel_cH(iyear)=logistic(agebins, selpar_L50_cH1, selpar_slope_cH1);

     sel_mcvt(iyear)=logistic(agebins, selpar_L50_mcvt, selpar_slope_mcvt);
     //sel_mcvt(iyear)=logistic_exponential(agebins, selpar_L50_mcvt, selpar_slope_mcvt, selpar_sigma_mcvt, selpar_afull_mcvt);   

     //sel_vid(iyear)=logistic_exponential(agebins, selpar_L50_vid, selpar_slope_vid, selpar_sigma_vid, selpar_afull_vid);      
     //sel_vid(iyear)=logistic(agebins, selpar_L50_vid, selpar_slope_vid);
     //sel_vid(iyear)=sel_mcvt(iyear);

     sel_HB(iyear)=logistic_double(agebins, selpar_L51_HB1, selpar_slope1_HB1, selpar_L52_HB1, selpar_slope2_HB1);
    //sel_HB(iyear)=logistic(agebins, selpar_L50_HB1, selpar_slope_HB1);
    //sel_HB(iyear)=logistic_exponential(agebins, selpar_L50_HB1, selpar_slope_HB1, selpar_sigma_HB1, selpar_afull_HB1);
    
     sel_HB_D(iyear)=logistic_exponential(agebins, selpar_L50_HB_D, selpar_slope_HB_D, selpar_sigma_HB_D, selpar_afull_HB_D); 

     sel_GR(iyear)=logistic_double(agebins, selpar_L51_GR, selpar_slope1_GR, selpar_L52_GR, selpar_slope2_GR);
    //sel_GR(iyear)=logistic(agebins, selpar_L50_GR, selpar_slope_GR);

     sel_GR_D(iyear)=sel_HB_D(iyear); 
  }

//year-specific age at peak selectivity based on inspection of HB age comps (HB1=peak at age 3, HB2=peak at age 6
//if no annual age comp available then assume peak at age 3 (HB1); above loop
 
//     sel_HB(1991)=logistic_double(agebins, selpar_L50_HB2, selpar_slope_HB2, selpar_afull_HB2, selpar_sigma_HB2);
//     sel_HB(2007)=logistic_double(agebins, selpar_L50_HB2, selpar_slope_HB2, selpar_afull_HB2, selpar_sigma_HB2);
//     sel_HB(2008)=logistic_double(agebins, selpar_L50_HB2, selpar_slope_HB2, selpar_afull_HB2, selpar_sigma_HB2);
//     sel_HB(2009)=logistic_double(agebins, selpar_L50_HB2, selpar_slope_HB2, selpar_afull_HB2, selpar_sigma_HB2);
//     sel_HB(2011)=logistic_double(agebins, selpar_L50_HB2, selpar_slope_HB2, selpar_afull_HB2, selpar_sigma_HB2); 

//BLOCK 2 for selex. 
//  for (iyear=(endyr_selex_phase1+1); iyear<=endyr_selex_phase2; iyear++)
//   {     
//     sel_cH(iyear)=logistic(agebins, selpar_L50_cH2, selpar_slope_cH2);
//     sel_mcvt(iyear)=logistic_exponential(agebins, selpar_L50_mcvt, selpar_slope_mcvt, selpar_sigma_mcvt, selpar_afull_mcvt);      
//     sel_HB(iyear)=logistic(agebins, selpar_L50_HB2, selpar_slope_HB2);
     
//     sel_mcvt(iyear)(8,nages)=sel_mcvt(iyear)(7);

//   }
   
//  //BLOCK 3 for selex. 
//  for (iyear=(endyr_selex_phase2+1); iyear<=endyr; iyear++)
//   {
//     sel_cH(iyear)=logistic(agebins, selpar_L50_cH3, selpar_slope_cH3);
//     sel_mcvt(iyear)=logistic_exponential(agebins, selpar_L50_mcvt, selpar_slope_mcvt, selpar_sigma_mcvt, selpar_afull_mcvt);      
//     sel_HB(iyear)=logistic(agebins, selpar_L50_HB3, selpar_slope_HB3);
     
//     sel_mcvt(iyear)(8,nages)=sel_mcvt(iyear)(7);
//    }
   
//  sel_initial=sel_rec(styr);
 
  
FUNCTION get_mortality

  Fsum.initialize();
  Fapex.initialize();
  F.initialize();
  //initialization F is avg from first 3 yrs of observed landings
  log_F_dev_init_cH=sum(log_F_dev_cH(styr_cH_L,(styr_cH_L+2)))/3.0;         
  log_F_dev_init_HB=sum(log_F_dev_HB(styr_HB_L,(styr_HB_L+2)))/3.0;         
  log_F_dev_init_GR=sum(log_F_dev_GR(styr_GR_L,(styr_GR_L+2)))/3.0;         
  
  for (iyear=styr; iyear<=endyr; iyear++) 
  {
    if(iyear>=styr_cH_L & iyear<=endyr_cH_L)
    {  F_cH_out(iyear)=mfexp(log_avg_F_cH+log_F_dev_cH(iyear)); //}    
    // if (iyear<styr_cH_L){F_cH_out(iyear)=mfexp(log_avg_F_cH+log_F_dev_init_cH);}        
       F_cH(iyear)=sel_cH(iyear)*F_cH_out(iyear);
       Fsum(iyear)+=F_cH_out(iyear);
    }

    if(iyear>=styr_HB_L & iyear<=endyr_HB_L)
    {  F_HB_out(iyear)=mfexp(log_avg_F_HB+log_F_dev_HB(iyear)); //}    
    // if (iyear<styr_HB_L){F_HB_out(iyear)=mfexp(log_avg_F_HB+log_F_dev_init_HB);}        
       F_HB(iyear)=sel_HB(iyear)*F_HB_out(iyear);
       Fsum(iyear)+=F_HB_out(iyear);
    }
   
    if(iyear>=styr_GR_L & iyear<=endyr_GR_L) 
    {  F_GR_out(iyear)=mfexp(log_avg_F_GR+log_F_dev_GR(iyear)); //}    
    // if (iyear<styr_GR_L){F_GR_out(iyear)=mfexp(log_avg_F_GR+log_F_dev_init_GR);}        
       F_GR(iyear)=sel_GR(iyear)*F_GR_out(iyear); //general rec shares headboat selex 
       Fsum(iyear)+=F_GR_out(iyear);
    }

    if(iyear>=styr_HB_D & iyear<=endyr_HB_D)
    {  F_HB_D_out(iyear)=mfexp(log_avg_F_HB_D+log_F_dev_HB_D(iyear)); //}    
    // if (iyear<styr_HB_D){F_HB_D_out(iyear)=mfexp(log_avg_F_HB_D+log_F_dev_init_HB_D);}        
       F_HB_D(iyear)=sel_HB_D(iyear)*F_HB_D_out(iyear);
       Fsum(iyear)+=F_HB_D_out(iyear);
    }

    if(iyear>=styr_GR_D & iyear<=endyr_GR_D)
    {  F_GR_D_out(iyear)=mfexp(log_avg_F_GR_D+log_F_dev_GR_D(iyear)); //}    
    // if (iyear<styr_GR_D){F_GR_D_out(iyear)=mfexp(log_avg_F_GR_D+log_F_dev_init_GR_D);}        
       F_GR_D(iyear)=sel_HB_D(iyear)*F_GR_D_out(iyear); //general rec discards shares headboat discards selex
       Fsum(iyear)+=F_GR_D_out(iyear);
    }
 
    //Total F at age
    F(iyear)=F_cH(iyear);  //first in additive series (NO +=)
    F(iyear)+=F_HB(iyear);
    F(iyear)+=F_GR(iyear);
    F(iyear)+=F_HB_D(iyear);
    F(iyear)+=F_GR_D(iyear);
    
    Fapex(iyear)=max(F(iyear));
    Z(iyear)=M+F(iyear);
    
  }  //end iyear 
 

FUNCTION get_bias_corr

  var_rec_dev=norm2(log_rec_dev(styr_rec_dev,endyr_rec_dev)-
              sum(log_rec_dev(styr_rec_dev,endyr_rec_dev))/nyrs_rec)
              /(nyrs_rec-1.0);                           
  //if (set_BiasCor <= 0.0) {BiasCor=mfexp(var_rec_dev/2.0);}   //bias correction based on empirical residuals
  rec_sigma_sq=square(rec_sigma);
  if (set_BiasCor <= 0.0) {BiasCor=mfexp(rec_sigma_sq/2.0);}   //bias correction based on Rsigma               
  else {BiasCor=set_BiasCor;}


FUNCTION get_numbers_at_age

//Initialization
  R0=mfexp(log_R0);
  S0=spr_F0*R0;
  //R_virgin=(R0/((5.0*steep-1.0)*spr_F0))*
  //               (BiasCor*4.0*steep*spr_F0-spr_F0*(1.0-steep));
  //R_virgin=R0/(spr_F0/spr_F0)*BiasCor*(1.0+log(spr_F0/spr_F0)/steep); //Ricker
 
  R_virgin=SR_eq_func(R0, steep, spr_F0, spr_F0, BiasCor, SR_switch);
 
  B0=bpr_F0*R_virgin;   
  B0_q_DD=R_virgin*sum(elem_prod(N_bpr_F0(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages))); 
  
  F_init_denom=mfexp(log_avg_F_cH+log_F_dev_init_cH)+
               mfexp(log_avg_F_HB+log_F_dev_init_HB)+
               mfexp(log_avg_F_GR+log_F_dev_init_GR)+
               mfexp(log_avg_F_HB_D+log_F_dev_init_HB_D)+
               mfexp(log_avg_F_GR_D+log_F_dev_init_GR_D);

  F_init_cH_prop=mfexp(log_avg_F_cH+log_F_dev_init_cH)/F_init_denom;
  F_init_HB_prop=mfexp(log_avg_F_HB+log_F_dev_init_HB)/F_init_denom;
  F_init_GR_prop=mfexp(log_avg_F_GR+log_F_dev_init_GR)/F_init_denom;
  F_init_HB_D_prop=mfexp(log_avg_F_HB_D+log_F_dev_init_HB_D)/F_init_denom;
  F_init_GR_D_prop=mfexp(log_avg_F_GR_D+log_F_dev_init_GR_D)/F_init_denom;

  F_initial=sel_cH(styr)*F_init*F_init_cH_prop+
            sel_HB(styr)*F_init*F_init_HB_prop+
            sel_GR(styr)*F_init*F_init_GR_prop+
            sel_HB_D(styr)*F_init*F_init_HB_prop+
            sel_GR_D(styr)*F_init*F_init_GR_prop;

  Z_initial=M+F_initial;

//Initial equilibrium age structure
  N_spr_initial(1)=1.0*mfexp(-1.0*Z_initial(1)*spawn_time_frac); //at peak spawning time;
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   mfexp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  spr_initial=sum(elem_prod(N_spr_initial,reprod));  

  //if (styr==styr_rec_dev) {R1=(R0/((5.0*steep-1.0)*spr_initial))*
  //               (4.0*steep*spr_initial-spr_F0*(1.0-steep));} //without bias correction (deviation added later)
  //else {R1=(R0/((5.0*steep-1.0)*spr_initial))*
  //               (BiasCor*4.0*steep*spr_initial-spr_F0*(1.0-steep));} //with bias correction                 

  if (styr==styr_rec_dev) {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, 1.0, SR_switch);} //without bias correction (deviation added later)
  else {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, BiasCor, SR_switch);} //with bias correction
  if(R1<10.0) {R1=10.0;} //Avoid unrealistically low popn sizes during search algorithm

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

  SSB(styr)=sum(elem_prod(N_spawn(styr),reprod));  //KC4
  MatFemB(styr)=sum(elem_prod(N_spawn(styr),reprod2));
  B_q_DD(styr)=sum(elem_prod(N(styr)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
    
//Rest of years 
  for (iyear=styr; iyear<endyr; iyear++)
  {
    if(iyear<(styr_rec_dev-1)||iyear>(endyr_rec_dev-1)) //recruitment follows S-R curve (with bias correction) exactly
    {
        N(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch);
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod)); //KC
//KC        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod(iyear+1))); 
//KC        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod2(iyear+1))); 
        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod2));  //KC  
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));       
    }
    else   //recruitment follows S-R curve with lognormal deviation
    {
        N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch)*mfexp(log_rec_dev(iyear+1));
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));  //KC4
//KC        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod(iyear+1)));
//KC        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod2(iyear+1)));
        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod2));
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
    }
  }
  
  //last year (projection) has no recruitment variability
  N(endyr+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(endyr),SR_switch);
  N(endyr+1)(2,nages)=++elem_prod(N(endyr)(1,nages-1),(mfexp(-1.*Z(endyr)(1,nages-1))));
  N(endyr+1,nages)+=N(endyr,nages)*mfexp(-1.*Z(endyr,nages)); //plus group

  
FUNCTION get_landings_numbers //Baranov catch eqn

  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      L_cH_num(iyear,iage)=N(iyear,iage)*F_cH(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_HB_num(iyear,iage)=N(iyear,iage)*F_HB(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_GR_num(iyear,iage)=N(iyear,iage)*F_GR(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);        
    }          
    pred_cH_L_knum(iyear)=sum(L_cH_num(iyear))/1000.0;    
    pred_HB_L_knum(iyear)=sum(L_HB_num(iyear))/1000.0;
    pred_GR_L_knum(iyear)=sum(L_GR_num(iyear))/1000.0;
  }

 
FUNCTION get_landings_wgt

  for (iyear=styr; iyear<=endyr; iyear++)
  {    
    L_cH_klb(iyear)=elem_prod(L_cH_num(iyear),wholewgt_cH_klb(iyear));     //in 1000 lb whole weight  
    L_HB_klb(iyear)=elem_prod(L_HB_num(iyear),wholewgt_HB_klb(iyear));     //in 1000 lb whole weight 
    L_GR_klb(iyear)=elem_prod(L_GR_num(iyear),wholewgt_GR_klb(iyear));     //in 1000 lb whole weight 
    
    pred_cH_L_klb(iyear)=sum(L_cH_klb(iyear));
    pred_HB_L_klb(iyear)=sum(L_HB_klb(iyear));
    pred_GR_L_klb(iyear)=sum(L_GR_klb(iyear));    
  }

 
FUNCTION get_dead_discards

  //dead discards at age (number fish) 

  for (iyear=styr_HB_D; iyear<=endyr_HB_D; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_HB_num(iyear,iage)=N(iyear,iage)*F_HB_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }

    pred_HB_D_knum(iyear)=sum(D_HB_num(iyear))/1000.0;            //pred annual dead discards in 1000s (for matching data)
	
      pred_HB_D_klb(iyear)=sum(elem_prod(D_HB_num(iyear),wholewgt_HB_D_klb(iyear))); //annual dead discards in 1000 lb whole (for output only) 
  }

  for (iyear=styr_GR_D; iyear<=endyr_GR_D; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_GR_num(iyear,iage)=N(iyear,iage)*F_GR_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }
    pred_GR_D_knum(iyear)=sum(D_GR_num(iyear))/1000.0;            //pred annual dead discards in 1000s (for matching data)
    pred_GR_D_klb(iyear)=sum(elem_prod(D_GR_num(iyear),wholewgt_GR_D_klb(iyear))); //annual dead discards in 1000 lb whole (for output only) 
  }


FUNCTION get_catchability_fcns 
   
 //Get rate increase if estimated, otherwise fixed above
  if (set_q_rate_phase>0.0)
  {

      for (iyear=styr_cH_cpue; iyear<=endyr_cH_cpue; iyear++)
      {   if (iyear>styr_cH_cpue & iyear <=2003) 
          {//q_rate_fcn_cH(iyear)=(1.0+q_rate)*q_rate_fcn_cH(iyear-1); //compound
             q_rate_fcn_cH(iyear)=(1.0+(iyear-styr_cH_cpue)*q_rate)*q_rate_fcn_cH(styr_cH_cpue);  //linear
          }
          if (iyear>2003) {q_rate_fcn_cH(iyear)=q_rate_fcn_cH(iyear-1);} 
      }   
      for (iyear=styr_mcvt_cpue; iyear<=endyr_mcvt_cpue; iyear++) //KC added this section
      {   if (iyear>styr_mcvt_cpue & iyear <=2003) 
          {//q_rate_fcn_mcvt(iyear)=(1.0+q_rate)*q_rate_fcn_mcvt(iyear-1); //compound
             q_rate_fcn_mcvt(iyear)=(1.0+(iyear-styr_mcvt_cpue)*q_rate)*q_rate_fcn_mcvt(styr_mcvt_cpue);  //linear
          }
          if (iyear>2003) {q_rate_fcn_mcvt(iyear)=q_rate_fcn_mcvt(iyear-1);} 
      } 
//      for (iyear=styr_vid_cpue; iyear<=endyr_vid_cpue; iyear++) 
//      {   if (iyear>styr_vid_cpue & iyear <=2003) 
//          {//q_rate_fcn_vid(iyear)=(1.0+q_rate)*q_rate_fcn_vid(iyear-1); //compound
//             q_rate_fcn_vid(iyear)=(1.0+(iyear-styr_vid_cpue)*q_rate)*q_rate_fcn_vid(styr_vid_cpue);  //linear
//          }
//          if (iyear>2003) {q_rate_fcn_vid(iyear)=q_rate_fcn_vid(iyear-1);} 
//      }
      for (iyear=styr_HB_cpue; iyear<=endyr_HB_cpue; iyear++)
      {   if (iyear>styr_HB_cpue & iyear <=2003) 
          {//q_rate_fcn_HB(iyear)=(1.0+q_rate)*q_rate_fcn_HB(iyear-1); //compound
             q_rate_fcn_HB(iyear)=(1.0+(iyear-styr_HB_cpue)*q_rate)*q_rate_fcn_HB(styr_HB_cpue);  //linear
          }
          if (iyear>2003) {q_rate_fcn_HB(iyear)=q_rate_fcn_HB(iyear-1);} 
      }   
      for (iyear=styr_GR_cpue; iyear<=endyr_GR_cpue; iyear++)
      {   if (iyear>styr_GR_cpue & iyear <=2003) 
          {//q_rate_fcn_GR(iyear)=(1.0+q_rate)*q_rate_fcn_GR(iyear-1); //compound
             q_rate_fcn_GR(iyear)=(1.0+(iyear-styr_GR_cpue)*q_rate)*q_rate_fcn_GR(styr_GR_cpue);  //linear
          }
          if (iyear>2003) {q_rate_fcn_GR(iyear)=q_rate_fcn_GR(iyear-1);} 
      }   
      
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

 //cH  cpue
  q_cH(styr_cH_cpue)=mfexp(log_q_cH); 
  for (iyear=styr_cH_cpue; iyear<=endyr_cH_cpue; iyear++)
  {//index in weight units. original index in lb and re-scaled. predicted in klb whole weight, but difference in lb and klb is absorbed by q
      N_cH(iyear)=elem_prod(elem_prod(N_mdyr(iyear),sel_cH(iyear)),wholewgt_cH_klb(iyear));   
      //N_cH(iyear)=elem_prod(N_mdyr(iyear),sel_cH(iyear)); 
      pred_cH_cpue(iyear)=q_cH(iyear)*q_rate_fcn_cH(iyear)*q_DD_fcn(iyear)*sum(N_cH(iyear));
      if (iyear<endyr_cH_cpue){q_cH(iyear+1)=q_cH(iyear)*mfexp(q_RW_log_dev_cH(iyear));}
  }

 //mcvt  cpue 
  q_mcvt(styr_mcvt_cpue)=mfexp(log_q_mcvt); 
  for (iyear=styr_mcvt_cpue; iyear<=endyr_mcvt_cpue; iyear++)
  {//index in weight units. original index in lb and re-scaled. predicted in klb whole weight, but difference in lb and klb is absorbed by q
      //N_mcvt(iyear)=elem_prod(elem_prod(N_mdyr(iyear),sel_mcvt(iyear)),wholewgt_mcvt_klb(iyear));   
      N_mcvt(iyear)=elem_prod(N_mdyr(iyear),sel_mcvt(iyear));  
      pred_mcvt_cpue(iyear)=q_mcvt(iyear)*q_rate_fcn_mcvt(iyear)*q_DD_fcn(iyear)*sum(N_mcvt(iyear));
      if (iyear<endyr_mcvt_cpue){q_mcvt(iyear+1)=q_mcvt(iyear)*mfexp(q_RW_log_dev_mcvt(iyear));}
  }

//video  cpue 
//  q_vid(styr_vid_cpue)=mfexp(log_q_vid); 
//  for (iyear=styr_vid_cpue; iyear<=endyr_vid_cpue; iyear++)
//  {//index in weight units. original index in lb and re-scaled. predicted in klb whole weight, but difference in lb and klb is absorbed by q
//      N_vid(iyear)=elem_prod(N_mdyr(iyear),sel_vid(iyear));   
//      //N_vid(iyear)=elem_prod(N_mdyr(iyear),sel_vid(iyear)); 
//      pred_vid_cpue(iyear)=q_vid(iyear)*q_rate_fcn_vid(iyear)*q_DD_fcn(iyear)*sum(N_vid(iyear));
//      if (iyear<endyr_vid_cpue){q_vid(iyear+1)=q_vid(iyear)*mfexp(q_RW_log_dev_vid(iyear));}
//  }

 //HB  cpue
  q_HB(styr_HB_cpue)=mfexp(log_q_HB); 
  for (iyear=styr_HB_cpue; iyear<=endyr_HB_cpue; iyear++)
  {   
      N_HB(iyear)=elem_prod(N_mdyr(iyear),sel_HB(iyear)); 
      pred_HB_cpue(iyear)=q_HB(iyear)*q_rate_fcn_HB(iyear)*q_DD_fcn(iyear)*sum(N_HB(iyear));
      if (iyear<endyr_HB_cpue){q_HB(iyear+1)=q_HB(iyear)*mfexp(q_RW_log_dev_HB(iyear));}
  }

 //GR  cpue
  q_GR(styr_GR_cpue)=mfexp(log_q_GR); 
  for (iyear=styr_GR_cpue; iyear<=endyr_GR_cpue; iyear++)
  {   
      N_GR(iyear)=elem_prod(N_mdyr(iyear),sel_HB(iyear)); //GR uses HB selex
      pred_GR_cpue(iyear)=q_GR(iyear)*q_rate_fcn_GR(iyear)*q_DD_fcn(iyear)*sum(N_GR(iyear));
      if (iyear<endyr_GR_cpue){q_GR(iyear+1)=q_GR(iyear)*mfexp(q_RW_log_dev_GR(iyear));}
  } 


FUNCTION get_length_comps

 //comm handline
  for (iyear=1;iyear<=nyr_cH_lenc;iyear++) 
  {pred_cH_lenc(iyear)=(L_cH_num(yrs_cH_lenc(iyear))*lenprob_cH)/sum(L_cH_num(yrs_cH_lenc(iyear)));}

 //SERFS chevron trap
  for (iyear=1;iyear<=nyr_mcvt_lenc;iyear++) 
  {pred_mcvt_lenc(iyear)=(N_mcvt(yrs_mcvt_lenc(iyear))*lenprob_mcvt)/sum(N_mcvt(yrs_mcvt_lenc(iyear)));}
 
 //headboat
  for (iyear=1;iyear<=nyr_HB_lenc;iyear++) 
  {pred_HB_lenc(iyear)=(L_HB_num(yrs_HB_lenc(iyear))*lenprob_HB)/sum(L_HB_num(yrs_HB_lenc(iyear)));}
  
 //headboat discards
  for (iyear=1;iyear<=nyr_HB_D_lenc;iyear++) 
  {pred_HB_D_lenc(iyear)=(D_HB_num(yrs_HB_D_lenc(iyear))*lenprob_HB_D)/sum(D_HB_num(yrs_HB_D_lenc(iyear)));} 

 //GR (MRIP)
  for (iyear=1;iyear<=nyr_GR_lenc;iyear++) 
  {pred_GR_lenc(iyear)=(L_GR_num(yrs_GR_lenc(iyear))*lenprob_GR)/sum(L_GR_num(yrs_GR_lenc(iyear)));}
  	
       // cout << pred_HB_D_lenc <<endl;
  

FUNCTION get_age_comps 

  //Commercial handline
  for (iyear=1;iyear<=nyr_cH_agec;iyear++) 
  {
    ErrorFree_cH_agec(iyear)=L_cH_num(yrs_cH_agec(iyear))/sum(L_cH_num(yrs_cH_agec(iyear)));  
    //ErrorFree_cH_agec(iyear)=elem_prod(N(yrs_cH_agec(iyear)),sel_cH(yrs_cH_agec(iyear)));
    pred_cH_agec_allages(iyear)=age_error*(ErrorFree_cH_agec(iyear)/sum(ErrorFree_cH_agec(iyear)));   
    for (iage=1; iage<=nages_agec; iage++) {pred_cH_agec(iyear,iage)=pred_cH_agec_allages(iyear,iage);} 
    for (iage=(nages_agec+1); iage<=nages; iage++) {pred_cH_agec(iyear,nages_agec)+=pred_cH_agec_allages(iyear,iage);} //plus group                             
  }
 
  //SERFS chevron trap
  for (iyear=1;iyear<=nyr_mcvt_agec;iyear++) 
  {
    ErrorFree_mcvt_agec(iyear)=N_mcvt(yrs_mcvt_agec(iyear))/sum(N_mcvt(yrs_mcvt_agec(iyear)));
    //ErrorFree_mcvt_agec(iyear)=elem_prod(N(yrs_mcvt_agec(iyear)),sel_mcvt(yrs_mcvt_agec(iyear)));
    pred_mcvt_agec_allages(iyear)=age_error*(ErrorFree_mcvt_agec(iyear)/sum(ErrorFree_mcvt_agec(iyear)));     
    for (iage=1; iage<=nages_agec; iage++) {pred_mcvt_agec(iyear,iage)=pred_mcvt_agec_allages(iyear,iage);} 
    for (iage=(nages_agec+1); iage<=nages; iage++) {pred_mcvt_agec(iyear,nages_agec)+=pred_mcvt_agec_allages(iyear,iage);} //plus group                           
  }

 //Headboat
 for (iyear=1;iyear<=nyr_HB_agec;iyear++)
  {
    ErrorFree_HB_agec(iyear)=L_HB_num(yrs_HB_agec(iyear))/sum(L_HB_num(yrs_HB_agec(iyear)));
    pred_HB_agec_allages(iyear)=age_error*ErrorFree_HB_agec(iyear); 
    for (iage=1; iage<=nages_agec; iage++) {pred_HB_agec(iyear,iage)=pred_HB_agec_allages(iyear,iage);} 
    for (iage=(nages_agec+1); iage<=nages; iage++) {pred_HB_agec(iyear,nages_agec)+=pred_HB_agec_allages(iyear,iage);} //plus group                        
  }
 
 // //GR (MRIP)
 // for (iyear=1;iyear<=nyr_GR_agec;iyear++)
  // {
    // ErrorFree_GR_agec(iyear)=L_GR_num(yrs_GR_agec(iyear))/sum(L_GR_num(yrs_GR_agec(iyear)));
    // pred_GR_agec_allages(iyear)=age_error*ErrorFree_GR_agec(iyear); 
    // for (iage=1; iage<=nages_agec; iage++) {pred_GR_agec(iyear,iage)=pred_GR_agec_allages(iyear,iage);} 
    // for (iage=(nages_agec+1); iage<=nages; iage++) {pred_GR_agec(iyear,nages_agec)+=pred_GR_agec_allages(iyear,iage);} //plus group                        
  // }  


////--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_weighted_current 
  F_temp_sum=0.0;
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_cH+
        sum(log_F_dev_cH((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_HB+
        sum(log_F_dev_HB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_GR+
        sum(log_F_dev_GR((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_HB_D+
        sum(log_F_dev_HB_D((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_GR_D+
        sum(log_F_dev_GR_D((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);

  F_cH_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_cH+
        sum(log_F_dev_cH((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_HB_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_HB+
        sum(log_F_dev_HB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_GR_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_GR+
        sum(log_F_dev_GR((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_HB_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_HB_D+
        sum(log_F_dev_HB_D((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_GR_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_GR_D+
        sum(log_F_dev_GR_D((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  
  log_F_dev_end_cH=sum(log_F_dev_cH((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_HB=sum(log_F_dev_HB((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  
  log_F_dev_end_GR=sum(log_F_dev_GR((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  

  log_F_dev_end_HB_D=sum(log_F_dev_HB_D((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_GR_D=sum(log_F_dev_GR_D((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;

  F_end_L=sel_cH(endyr)*mfexp(log_avg_F_cH+log_F_dev_end_cH)+
          sel_HB(endyr)*mfexp(log_avg_F_HB+log_F_dev_end_HB)+
          sel_GR(endyr)*mfexp(log_avg_F_GR+log_F_dev_end_GR);     
 
  F_end_D=sel_HB_D(endyr)*mfexp(log_avg_F_HB_D+log_F_dev_end_HB_D)+
          sel_GR_D(endyr)*mfexp(log_avg_F_GR_D+log_F_dev_end_GR_D);  
    
  F_end=F_end_L+F_end_D;
  F_end_apex=max(F_end);
  
  sel_wgted_tot=F_end/F_end_apex;
  sel_wgted_L=elem_prod(sel_wgted_tot, elem_div(F_end_L,F_end));
  sel_wgted_D=elem_prod(sel_wgted_tot, elem_div(F_end_D,F_end));

  wgt_wgted_L_denom=F_cH_prop+F_HB_prop+F_GR_prop; 
  wgt_wgted_L_klb=F_cH_prop/wgt_wgted_L_denom*wholewgt_cH_klb(endyr)+      
                  F_HB_prop/wgt_wgted_L_denom*wholewgt_HB_klb(endyr)+  
                  F_GR_prop/wgt_wgted_L_denom*wholewgt_GR_klb(endyr);                       

  wgt_wgted_D_denom=F_HB_D_prop+F_GR_D_prop;
  wgt_wgted_D_klb=F_HB_D_prop/wgt_wgted_D_denom*wholewgt_HB_D_klb(endyr)+  
                 F_GR_D_prop/wgt_wgted_D_denom*wholewgt_GR_D_klb(endyr); 
 
  
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
      {N_age_msy(iage)=N_age_msy(iage-1)*mfexp(-1.*Z_age_msy(iage-1));}
    N_age_msy(nages)=N_age_msy(nages)/(1.0-mfexp(-1.*Z_age_msy(nages)));
    N_age_msy_spawn(1,(nages-1))=elem_prod(N_age_msy(1,(nages-1)),
                                   mfexp((-1.*Z_age_msy(1,(nages-1)))*spawn_time_frac));                 
    N_age_msy_spawn(nages)=(N_age_msy_spawn(nages-1)*(mfexp(-1.*(Z_age_msy(nages-1)*(1.0-spawn_time_frac) + 
                            Z_age_msy(nages)*spawn_time_frac) )))/(1.0-mfexp(-1.*Z_age_msy(nages)));
                     
    spr_msy(ff)=sum(elem_prod(N_age_msy_spawn,reprod)); 
        
    R_eq(ff)=SR_eq_func(R0, steep, spr_msy(1), spr_msy(ff), BiasCor, SR_switch);
    
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
    L_eq_klb(ff)=sum(elem_prod(L_age_msy,wgt_wgted_L_klb)); //in whole weight
    L_eq_knum(ff)=sum(L_age_msy)/1000.0;  
    D_eq_klb(ff)=sum(elem_prod(D_age_msy,wgt_wgted_D_klb)); //in whole weight   
    D_eq_knum(ff)=sum(D_age_msy)/1000.0;    
  }  
  
  msy_klb_out=max(L_eq_klb); //msy in whole weight
  
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
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z(iyear,iage-1));}
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
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));}
    N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
    N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
    N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));
    
    spr_spr(ff)=sum(elem_prod(N_age_spr_spawn,reprod));  //KC4
    L_spr(ff)=0.0;
    for (iage=1; iage<=nages; iage++)
    {
      L_age_spr(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      L_spr(ff)+=L_age_spr(iage)*wgt_wgted_L_klb(iage)*1000.0; //in lb gutted wgt
    }   
  }
  spr_ratio=spr_spr/spr_F0;	//added 1-19-16 here to...
  F20_dum=min(fabs(spr_ratio-0.2));
  F30_dum=min(fabs(spr_ratio-0.3));
  F40_dum=min(fabs(spr_ratio-0.4));
  for(ff=1; ff<=n_iter_spr; ff++)
  {   
      if (fabs(spr_ratio(ff)-0.2)==F20_dum) {F20_out=F_spr(ff);
	  
	  }
	  
	  if (fabs(spr_ratio(ff)-0.3)==F30_dum) {
		F30_out=F_spr(ff);
        SSB_F30_out=SSB_eq(ff);  //NOTE, this works bc F grid for msy calcs is the same as for spr calcs
		B_F30_out=B_eq(ff);
        R_F30_out=R_eq(ff);
        L_F30_knum_out=L_eq_knum(ff);
        L_F30_klb_out=L_eq_klb(ff);
		D_F30_knum_out=D_eq_knum(ff);
        D_F30_klb_out=D_eq_klb(ff);          	  
	  }
	  
	  if (fabs(spr_ratio(ff)-0.4)==F40_dum) {
		F40_out=F_spr(ff);
	  }
  }  //here.
  
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_miscellaneous_stuff

//switch here if var_rec_dev <=dzero 
  if(var_rec_dev>0.0)
   {sigma_rec_dev=sqrt(var_rec_dev);} //pow(var_rec_dev,0.5);  //sample SD of predicted residuals (may not equal rec_sigma)  
   else{sigma_rec_dev=0.0;}

//  len_cv=elem_div(len_sd,meanlen_FL);
  len_sd=elem_prod(len_cv,meanlen_FL);

 for(iage=1;iage<=nages;iage++)		//Added this loop to allow plotting of both Population and Landings length at age
    	{
//        len_sd(iage)=len_sd_val;
//        len_cv(iage)=len_sd(iage)/meanlen_FL(iage);
//        len_sd_L(iage)=len_sd_val_L;
//        len_cv_L(iage)=len_sd_L(iage)/meanlen_FL_L(iage);
//        len_sd_mcvt(iage)=len_sd_val_mcvt;
//        len_cv_mcvt(iage)=len_sd_mcvt(iage)/meanlen_FL_mcvt(iage);

    	len_cv(iage)=len_cv_val;
    	len_sd(iage)=meanlen_FL(iage)*len_cv(iage);
    	len_cv_L(iage)=len_cv_val_L;
	len_sd_L(iage)=meanlen_FL_L(iage)*len_cv_L(iage);
        len_cv_mcvt(iage)=len_cv_val_mcvt;
        len_sd_mcvt(iage)=meanlen_FL_mcvt(iage)*len_cv_mcvt(iage);
  	}

  //compute total landings- and discards-at-age in 1000 fish and klb whole weight
  L_total_num.initialize();
  L_total_klb.initialize();
  L_total_knum_yr.initialize();
  L_total_klb_yr.initialize();  
  D_total_num.initialize();
  D_total_klb.initialize();
  D_total_knum_yr.initialize();
  D_total_klb_yr.initialize();
 
  for(iyear=styr; iyear<=endyr; iyear++)
  {
             L_total_klb_yr(iyear)=pred_cH_L_klb(iyear)+pred_HB_L_klb(iyear)+pred_GR_L_klb(iyear);  
             L_total_knum_yr(iyear)=pred_cH_L_knum(iyear)+pred_HB_L_knum(iyear)+pred_GR_L_knum(iyear); 
                
        B(iyear)=elem_prod(N(iyear),wgt_mt);
        totN(iyear)=sum(N(iyear));
        totB(iyear)=sum(B(iyear));   
        
        if (iyear>=styr_HB_D && iyear<=endyr_HB_D)
        {
         D_total_knum_yr(iyear)+=pred_HB_D_knum(iyear);
         D_total_klb_yr(iyear)+=pred_HB_D_klb(iyear);
         D_HB_klb(iyear)=elem_prod(D_HB_num(iyear),wholewgt_HB_D_klb(iyear));     
        }    
    
        if (iyear>=styr_GR_D && iyear<=endyr_GR_D)
        {
         D_total_knum_yr(iyear)+=pred_GR_D_knum(iyear);
         D_total_klb_yr(iyear)+=pred_GR_D_klb(iyear);
         D_GR_klb(iyear)=elem_prod(D_GR_num(iyear),wholewgt_GR_D_klb(iyear));    
        }                          
       }

  L_total_num=L_cH_num+L_HB_num+L_GR_num;   //added landings at age in number fish
  L_total_klb=L_cH_klb+L_HB_klb+L_GR_klb;   //landings at age in klb whole weight

  D_total_num=(D_HB_num+D_GR_num);          //discards at age in number fish
  D_total_klb=D_HB_klb+D_GR_klb;            //discards at age in klb whole weight
 
 
  //Time series of interest

  B(endyr+1)=elem_prod(N(endyr+1),wgt_mt);
  totN(endyr+1)=sum(N(endyr+1));
  totB(endyr+1)=sum(B(endyr+1));  
  //N_spawn(endyr+1)=N(endyr+1);
  //SSB(endyr+1)=sum(elem_prod(N_spawn(endyr+1),reprod));  
  //MatFemB(endyr+1)=sum(elem_prod(N_spawn(endyr+1),reprod2)); 
  rec=column(N,1);
  SdS0=SSB/S0;

//  steep_sd=steep;
//  fullF_sd=Fsum;
  Fend_mean_temp=1.0;   //added 1-19-16 here to...
  for (iyear=1; iyear<=selpar_n_yrs_wgted; iyear++) {Fend_mean_temp*=Fapex(endyr-iyear+1);}
  Fend_mean=pow(Fend_mean_temp,(1.0/selpar_n_yrs_wgted));	 //here. 

  if(F_msy_out>0)
    {
      FdF_msy=Fapex/F_msy_out;
      FdF_msy_end=FdF_msy(endyr);
	  FdF_msy_end_mean=Fend_mean/F_msy_out;  //added 1-19-16
      //FdF_msy_end_mean=pow((FdF_msy(endyr)*FdF_msy(endyr-1)*FdF_msy(endyr-2)),(1.0/3.0));
    }
  if(SSB_msy_out>0)
    {
      SdSSB_msy=SSB/SSB_msy_out;
      SdSSB_msy_end=SdSSB_msy(endyr);
    }  
  if(F30_out>0)                  //added 1-19-16  here to...
    {
	  FdF30=Fapex/F30_out;
	  FdF30_end_mean=Fend_mean/F30_out;
	}
  if(SSB_F30_out>0)
    {
      SdSSB_F30=SSB/SSB_F30_out;
	  Sdmsst_F30=SSB/(smsy2msst*SSB_F30_out);  //KC-changed from smsy2msst75
      SdSSB_F30_end=SdSSB_F30(endyr);
	  Sdmsst_F30_end=Sdmsst_F30(endyr);
    }  							//here
   //fill in log recruitment deviations for yrs they are nonzero
   for(iyear=styr_rec_dev; iyear<=endyr_rec_dev; iyear++)
     {log_rec_dev_output(iyear)=log_rec_dev(iyear);}
   //fill in log Nage deviations for ages they are nonzero (ages2+)
   for(iage=2; iage<=nages; iage++)
     {log_Nage_dev_output(iage)=log_Nage_dev(iage);}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------     
FUNCTION get_effective_sample_sizes
      neff_cH_lenc_allyr_out=missing;
      neff_mcvt_lenc_allyr_out=missing;
      neff_HB_lenc_allyr_out=missing;
      neff_HB_D_lenc_allyr_out=missing;  
      neff_GR_lenc_allyr_out=missing;  
      
      neff_cH_agec_allyr_out=missing;
      neff_mcvt_agec_allyr_out=missing;
      neff_HB_agec_allyr_out=missing;
      //neff_GR_agec_allyr_out=missing;
   
      for (iyear=1; iyear<=nyr_cH_lenc; iyear++)
         {if (nsamp_cH_lenc(iyear)>=minSS_cH_lenc)
            {neff_cH_lenc_allyr_out(yrs_cH_lenc(iyear))=multinom_eff_N(pred_cH_lenc(iyear),obs_cH_lenc(iyear));} 
          else {neff_cH_lenc_allyr_out(yrs_cH_lenc(iyear))=-99;}
         }

      for (iyear=1; iyear<=nyr_mcvt_lenc; iyear++)
         {if (nsamp_mcvt_lenc(iyear)>=minSS_mcvt_lenc)
            {neff_mcvt_lenc_allyr_out(yrs_mcvt_lenc(iyear))=multinom_eff_N(pred_mcvt_lenc(iyear),obs_mcvt_lenc(iyear));} 
          else {neff_mcvt_lenc_allyr_out(yrs_mcvt_lenc(iyear))=-99;}
         }

      for (iyear=1; iyear<=nyr_HB_lenc; iyear++)
         {if (nsamp_HB_lenc(iyear)>=minSS_HB_lenc)
            {neff_HB_lenc_allyr_out(yrs_HB_lenc(iyear))=multinom_eff_N(pred_HB_lenc(iyear),obs_HB_lenc(iyear));} 
          else {neff_HB_lenc_allyr_out(yrs_HB_lenc(iyear))=-99;}
         }

      for (iyear=1; iyear<=nyr_HB_D_lenc; iyear++)  
         {if (nsamp_HB_D_lenc(iyear)>=minSS_HB_D_lenc)  
            {neff_HB_D_lenc_allyr_out(yrs_HB_D_lenc(iyear))=multinom_eff_N(pred_HB_D_lenc(iyear),obs_HB_D_lenc(iyear));}  
          else {neff_HB_D_lenc_allyr_out(yrs_HB_D_lenc(iyear))=-99;}  
         }  

      for (iyear=1; iyear<=nyr_GR_lenc; iyear++)
         {if (nsamp_GR_lenc(iyear)>=minSS_GR_lenc)
            {neff_GR_lenc_allyr_out(yrs_GR_lenc(iyear))=multinom_eff_N(pred_GR_lenc(iyear),obs_GR_lenc(iyear));} 
          else {neff_GR_lenc_allyr_out(yrs_GR_lenc(iyear))=-99;}
         }

      for (iyear=1; iyear<=nyr_cH_agec; iyear++)
         {if (nsamp_cH_agec(iyear)>=minSS_cH_agec)
            {neff_cH_agec_allyr_out(yrs_cH_agec(iyear))=multinom_eff_N(pred_cH_agec(iyear),obs_cH_agec(iyear));}                            
          else {neff_cH_agec_allyr_out(yrs_cH_agec(iyear))=-99;}
         }    
         
      for (iyear=1; iyear<=nyr_mcvt_lenc; iyear++)
         {if (nsamp_mcvt_lenc(iyear)>=minSS_mcvt_lenc)
            {neff_mcvt_lenc_allyr_out(yrs_mcvt_lenc(iyear))=multinom_eff_N(pred_mcvt_lenc(iyear),obs_mcvt_lenc(iyear));} 
          else {neff_mcvt_lenc_allyr_out(yrs_mcvt_lenc(iyear))=-99;}
         }

      for (iyear=1; iyear<=nyr_HB_agec; iyear++)
         {if (nsamp_HB_agec(iyear)>=minSS_HB_agec)
            {neff_HB_agec_allyr_out(yrs_HB_agec(iyear))=multinom_eff_N(pred_HB_agec(iyear),obs_HB_agec(iyear));}                            
          else {neff_HB_agec_allyr_out(yrs_HB_agec(iyear))=-99;}
         }

      // for (iyear=1; iyear<=nyr_GR_agec; iyear++)
         // {if (nsamp_GR_agec(iyear)>=minSS_GR_agec)
            // {neff_GR_agec_allyr_out(yrs_GR_agec(iyear))=multinom_eff_N(pred_GR_agec(iyear),obs_GR_agec(iyear));}                            
          // else {neff_GR_agec_allyr_out(yrs_GR_agec(iyear))=-99;}
         // }    
                           
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   


FUNCTION evaluate_objective_function
  //fval=square(xdum-9.0);
  
  fval=0.0;
  fval_data=0.0;  
//---likelihoods---------------------------

//---Indices-------------------------------

//  f_cH_cpue=0.0;
//  f_cH_cpue=lk_lognormal(pred_cH_cpue, obs_cH_cpue, cH_cpue_cv, w_I_cH);
//  fval+=f_cH_cpue;
//  fval_data+=f_cH_cpue;  

  f_mcvt_cpue=0.0;
  f_mcvt_cpue=lk_lognormal(pred_mcvt_cpue, obs_mcvt_cpue, mcvt_cpue_cv, w_I_mcvt);
  fval+=f_mcvt_cpue;
  fval_data+=f_mcvt_cpue;

//  f_vid_cpue=0.0;
//  f_vid_cpue=lk_lognormal(pred_vid_cpue, obs_vid_cpue, vid_cpue_cv, w_I_vid);
//  fval+=f_vid_cpue;
//  fval_data+=f_vid_cpue;

//  f_HB_cpue=0.0;
//  f_HB_cpue=lk_lognormal(pred_HB_cpue, obs_HB_cpue, HB_cpue_cv, w_I_HB);
//  fval+=f_HB_cpue;
//  fval_data+=f_HB_cpue;  

//  f_GR_cpue=0.0;
//  f_GR_cpue=lk_lognormal(pred_GR_cpue, obs_GR_cpue, GR_cpue_cv, w_I_GR);
//  fval+=f_GR_cpue;
//  fval_data+=f_GR_cpue;  

//---Landings-------------------------------
  
  //f_cH_L in 1000 lb whole wgt
  f_cH_L=lk_lognormal(pred_cH_L_klb(styr_cH_L,endyr_cH_L), obs_cH_L(styr_cH_L,endyr_cH_L),
                      cH_L_cv(styr_cH_L,endyr_cH_L), w_L);
  fval+=f_cH_L;
  fval_data+=f_cH_L;

  //f_HB_L in 1000 fish
  f_HB_L=lk_lognormal(pred_HB_L_knum(styr_HB_L,endyr_HB_L), obs_HB_L(styr_HB_L,endyr_HB_L), 
                      HB_L_cv(styr_HB_L,endyr_HB_L), w_L);
  fval+=f_HB_L;
  fval_data+=f_HB_L;  

  //f_GR_L in 1000 fish
  f_GR_L=lk_lognormal(pred_GR_L_knum(styr_GR_L,endyr_GR_L), obs_GR_L(styr_GR_L,endyr_GR_L), 
                      GR_L_cv(styr_GR_L,endyr_GR_L), w_L);
  fval+=f_GR_L;
  fval_data+=f_GR_L;  

//---Discards-------------------------------  

//  f_HB_D in 1000 fish
  f_HB_D=lk_lognormal(pred_HB_D_knum(styr_HB_D,endyr_HB_D), obs_HB_D(styr_HB_D,endyr_HB_D), 
                      HB_D_cv(styr_HB_D,endyr_HB_D), w_D);
  fval+=f_HB_D;
  fval_data+=f_HB_D;  

  //f_GR_D in 1000 fish
  f_GR_D=lk_lognormal(pred_GR_D_knum(styr_GR_D,endyr_GR_D), obs_GR_D(styr_GR_D,endyr_GR_D), 
                      GR_D_cv(styr_GR_D,endyr_GR_D), w_D);
  fval+=f_GR_D;
  fval_data+=f_GR_D;  
  
//---Length comps-------------------------------

//  f_cH_lenc
  f_cH_lenc=lk_robust_multinomial(nsamp_cH_lenc, pred_cH_lenc, obs_cH_lenc, nyr_cH_lenc, double(nlenbins), minSS_cH_lenc, w_lc_cH);
  //f_cH_lenc=lk_multinomial(nsamp_cH_lenc, pred_cH_lenc, obs_cH_lenc, nyr_cH_lenc, minSS_cH_lenc, w_lc_cH);
  fval+=f_cH_lenc;
  fval_data+=f_cH_lenc;
 
//  f_mcvt_lenc
  f_mcvt_lenc=lk_robust_multinomial(nsamp_mcvt_lenc, pred_mcvt_lenc, obs_mcvt_lenc, nyr_mcvt_lenc, double(nlenbins), minSS_mcvt_lenc, w_lc_mcvt);
  //f_mcvt_lenc=lk_multinomial(nsamp_mcvt_lenc, pred_mcvt_lenc, obs_mcvt_lenc, nyr_mcvt_lenc, minSS_mcvt_lenc, w_lc_mcvt);
  fval+=f_mcvt_lenc;
  fval_data+=f_mcvt_lenc;
  
//  f_HB_lenc
  f_HB_lenc=lk_robust_multinomial(nsamp_HB_lenc, pred_HB_lenc, obs_HB_lenc, nyr_HB_lenc, double(nlenbins), minSS_HB_lenc, w_lc_HB);
  //f_HB_lenc=lk_multinomial(nsamp_HB_lenc, pred_HB_lenc, obs_HB_lenc, nyr_HB_lenc, minSS_HB_lenc, w_lc_HB);
  fval+=f_HB_lenc;
  fval_data+=f_HB_lenc;

//  f_HB_lenc_D  
  f_HB_D_lenc=lk_robust_multinomial(nsamp_HB_D_lenc, pred_HB_D_lenc, obs_HB_D_lenc, nyr_HB_D_lenc, double(nlenbins), minSS_HB_D_lenc, w_lc_HB_D);  
  //f_HB_D_lenc=lk_multinomial(nsamp_HB_D_lenc, pred_HB_D_lenc, obs_HB_D_lenc, nyr_HB_D_lenc, minSS_HB_D_lenc, w_lc_HB_d);  
  fval+=f_HB_D_lenc;  
  fval_data+=f_HB_D_lenc;    

//  f_GR_lenc  
  f_GR_lenc=lk_robust_multinomial(nsamp_GR_lenc, pred_GR_lenc, obs_GR_lenc, nyr_GR_lenc, double(nlenbins), minSS_GR_lenc, w_lc_GR);
  //f_GR_lenc=lk_multinomial(nsamp_GR_lenc, pred_GR_lenc, obs_GR_lenc, nyr_GR_lenc, minSS_GR_lenc, w_lc_GR);
  fval+=f_GR_lenc;
  fval_data+=f_GR_lenc; 
//---Age comps-------------------------------

//  f_cH_agec
  f_cH_agec=lk_robust_multinomial(nsamp_cH_agec, pred_cH_agec, obs_cH_agec, nyr_cH_agec, double(nages_agec), minSS_cH_agec, w_ac_cH);
  //f_cH_agec=lk_multinomial(nsamp_cH_agec, pred_cH_agec, obs_cH_agec, nyr_cH_agec, minSS_cH_agec, w_ac_cH);
  fval+=f_cH_agec;
  fval_data+=f_cH_agec;

//  f_mcvt_agec
  f_mcvt_agec=lk_robust_multinomial(nsamp_mcvt_agec, pred_mcvt_agec, obs_mcvt_agec, nyr_mcvt_agec, double(nages_agec), minSS_mcvt_agec, w_ac_mcvt);
  //f_mcvt_agec=lk_multinomial(nsamp_mcvt_agec, pred_mcvt_agec, obs_mcvt_agec, nyr_mcvt_agec, minSS_mcvt_agec, w_ac_mcvt);
  fval+=f_mcvt_agec;
  fval_data+=f_mcvt_agec;

//  f_HB_agec
  f_HB_agec=lk_robust_multinomial(nsamp_HB_agec, pred_HB_agec, obs_HB_agec, nyr_HB_agec, double(nages_agec), minSS_HB_agec, w_ac_HB);
  //f_HB_agec=lk_multinomial(nsamp_HB_agec, pred_HB_agec, obs_HB_agec, nyr_HB_agec, minSS_HB_agec, w_ac_HB);
  fval+=f_HB_agec;
  fval_data+=f_HB_agec;

// //  f_GR_agec
  // f_GR_agec=lk_robust_multinomial(nsamp_GR_agec, pred_GR_agec, obs_GR_agec, nyr_GR_agec, double(nages_agec), minSS_GR_agec, w_ac_GR);
  // //f_GR_agec=lk_multinomial(nsamp_GR_agec, pred_GR_agec, obs_GR_agec, nyr_GR_agec, minSS_GR_agec, w_ac_GR);
  // fval+=f_GR_agec;
  // fval_data+=f_GR_agec;
//-----------Constraints and penalties--------------------------------
  
  //Light penalty applied to log_Nage_dev for deviation from zero. If not estimated, this penalty equals zero.
  f_Nage_init=norm2(log_Nage_dev);        
  fval+=w_Nage_init*f_Nage_init;
  
  f_rec_dev=0.0;
  //rec_sigma_sq=square(rec_sigma);
  rec_logL_add=nyrs_rec*log(rec_sigma);
  f_rec_dev=(square(log_rec_dev(styr_rec_dev) + rec_sigma_sq/2.0)/(2.0*rec_sigma_sq));
  for(iyear=(styr_rec_dev+1); iyear<=endyr_rec_dev; iyear++)
  {f_rec_dev+=(square(log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1) + rec_sigma_sq/2.0)/
               (2.0*rec_sigma_sq));}
  f_rec_dev+=rec_logL_add;            
  fval+=w_rec*f_rec_dev;
    

  f_rec_dev_early=0.0; //possible extra constraint on early rec deviations
  if (w_rec_early>0.0)
    { if (styr_rec_dev<endyr_rec_phase1)
        {  
          for(iyear=styr_rec_dev; iyear<=endyr_rec_phase1; iyear++)
          //{f_rec_dev_early+=(square(log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1) + rec_sigma_sq/2.0)/
          //                  (2.0*rec_sigma_sq)) + rec_logL_add;}
          {f_rec_dev_early+=square(log_rec_dev(iyear));}
        }
  fval+=w_rec_early*f_rec_dev_early;
  }
  
  f_rec_dev_end=0.0; //possible extra constraint on ending rec deviations
  if (w_rec_end>0.0)
  { if (endyr_rec_phase2<endyr_rec_dev)
        {  
          for(iyear=(endyr_rec_phase2+1); iyear<=endyr_rec_dev; iyear++)
          //{f_rec_dev_end+=(square(log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1) + rec_sigma_sq/2.0)/
          //                 (2.0*rec_sigma_sq)) + rec_logL_add;}
          {f_rec_dev_end+=square(log_rec_dev(iyear));}
        }
      fval+=w_rec_end*f_rec_dev_end;
   }  

  //Ftune penalty: does not apply in last phase
  f_Ftune=0.0; 
  if (w_Ftune>0.0)
  {if (set_Ftune>0.0 && !last_phase()) {f_Ftune=square(Fapex(set_Ftune_yr)-set_Ftune);}
   fval+=w_Ftune*f_Ftune;
  }

  //Penalty if apical F exceeds 3.0
  f_fullF_constraint=0.0;
  if (w_fullF>0.0)
  {for (iyear=styr; iyear<=endyr; iyear++)
        {if(Fapex(iyear)>3.0) {f_fullF_constraint+=(mfexp(Fapex(iyear)-3.0)-1.0);}}
   fval+=w_fullF*f_fullF_constraint;
  }
  
//  //Random walk components of fishery dependent indices
//  f_HB_RW_cpue=0.0;
//  for (iyear=styr_HB_cpue; iyear<endyr_HB_cpue; iyear++)
//      {f_HB_RW_cpue+=square(q_RW_log_dev_HB(iyear))/(2.0*set_q_RW_HB_var);}
//  fval+=f_HB_RW_cpue;   

  
//---Priors---------------------------------------------------
//neg_log_prior arguments: estimate, prior mean, prior var/-CV, pdf type
//Variance input as a negative value is considered to be CV in arithmetic space (CV=-1 implies loose prior) 
//pdf type 1=none, 2=lognormal, 3=normal, 4=beta 
  f_priors=0.0; 
  f_priors+=neg_log_prior(Linf,set_Linf(5),set_Linf(6),set_Linf(7));
  f_priors+=neg_log_prior(K,set_K(5),set_K(6),set_K(7));
  f_priors+=neg_log_prior(t0,set_t0(5),set_t0(6),set_t0(7));
  f_priors+=neg_log_prior(len_cv_val,set_len_cv(5),set_len_cv(6),set_len_cv(7));

  f_priors+=neg_log_prior(M_constant,set_M_constant(5),set_M_constant(6),set_M_constant(7));
  
  f_priors+=neg_log_prior(Linf_L,set_Linf_L(5),set_Linf_L(6),set_Linf_L(7));
  f_priors+=neg_log_prior(K_L,set_K_L(5),set_K_L(6),set_K_L(7));
  f_priors+=neg_log_prior(t0_L,set_t0_L(5),set_t0_L(6),set_t0_L(7));
  f_priors+=neg_log_prior(len_cv_val_L,set_len_cv_L(5),set_len_cv_L(6),set_len_cv_L(7));

  f_priors+=neg_log_prior(Linf_mcvt,set_Linf_mcvt(5),set_Linf_mcvt(6),set_Linf_mcvt(7));
  f_priors+=neg_log_prior(K_mcvt,set_K_mcvt(5),set_K_mcvt(6),set_K_mcvt(7));
  f_priors+=neg_log_prior(t0_mcvt,set_t0_mcvt(5),set_t0_mcvt(6),set_t0_mcvt(7));
  f_priors+=neg_log_prior(len_cv_val_mcvt,set_len_cv_mcvt(5),set_len_cv_mcvt(6),set_len_cv_mcvt(7));  

  f_priors+=neg_log_prior(M_constant,set_M_constant(5),set_M_constant(6),set_M_constant(7));
   
  f_priors+=neg_log_prior(steep,set_steep(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(log_R0,set_log_R0(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));
 
  f_priors+=neg_log_prior(selpar_L50_cH1,set_selpar_L50_cH1(5), set_selpar_L50_cH1(6), set_selpar_L50_cH1(7));
  f_priors+=neg_log_prior(selpar_slope_cH1,set_selpar_slope_cH1(5), set_selpar_slope_cH1(6), set_selpar_slope_cH1(7));
//  f_priors+=neg_log_prior(selpar_L50_cH2,set_selpar_L50_cH2(5), set_selpar_L50_cH2(6), set_selpar_L50_cH2(7));
//  f_priors+=neg_log_prior(selpar_slope_cH2,set_selpar_slope_cH2(5), set_selpar_slope_cH2(6), set_selpar_slope_cH2(7));
//  f_priors+=neg_log_prior(selpar_L50_cH3,set_selpar_L50_cH3(5), set_selpar_L50_cH3(6), set_selpar_L50_cH3(7));
//  f_priors+=neg_log_prior(selpar_slope_cH3,set_selpar_slope_cH3(5), set_selpar_slope_cH3(6), set_selpar_slope_cH3(7));

  f_priors+=neg_log_prior(selpar_L50_mcvt,set_selpar_L50_mcvt(5), set_selpar_L50_mcvt(6), set_selpar_L50_mcvt(7));
  f_priors+=neg_log_prior(selpar_slope_mcvt,set_selpar_slope_mcvt(5), set_selpar_slope_mcvt(6), set_selpar_slope_mcvt(7));
 
  f_priors+=neg_log_prior(selpar_L51_HB1,set_selpar_L51_HB1(5), set_selpar_L51_HB1(6), set_selpar_L51_HB1(7));
  f_priors+=neg_log_prior(selpar_slope1_HB1,set_selpar_slope1_HB1(5), set_selpar_slope1_HB1(6), set_selpar_slope1_HB1(7));
  f_priors+=neg_log_prior(selpar_L52_HB1,set_selpar_L52_HB1(5), set_selpar_L52_HB1(6), set_selpar_L52_HB1(7));
  f_priors+=neg_log_prior(selpar_slope2_HB1,set_selpar_slope2_HB1(5), set_selpar_slope2_HB1(6), set_selpar_slope2_HB1(7));

  f_priors+=neg_log_prior(selpar_L50_HB2,set_selpar_L50_HB2(5), set_selpar_L50_HB2(6), set_selpar_L50_HB2(7));
  f_priors+=neg_log_prior(selpar_slope_HB2,set_selpar_slope_HB2(5), set_selpar_slope_HB2(6), set_selpar_slope_HB2(7));
  f_priors+=neg_log_prior(selpar_afull_HB2,set_selpar_afull_HB2(5), set_selpar_afull_HB2(6), set_selpar_afull_HB2(7));
  f_priors+=neg_log_prior(selpar_sigma_HB2,set_selpar_sigma_HB2(5), set_selpar_sigma_HB2(6), set_selpar_sigma_HB2(7));

//  f_priors+=neg_log_prior(selpar_L50_HB3,set_selpar_L50_HB3(5), set_selpar_L50_HB3(6), set_selpar_L50_HB3(7));
//  f_priors+=neg_log_prior(selpar_slope_HB3,set_selpar_slope_HB3(5), set_selpar_slope_HB3(6), set_selpar_slope_HB3(7));

  f_priors+=neg_log_prior(selpar_L50_HB_D,set_selpar_L50_HB_D(5), set_selpar_L50_HB_D(6), set_selpar_L50_HB_D(7));
  f_priors+=neg_log_prior(selpar_slope_HB_D,set_selpar_slope_HB_D(5), set_selpar_slope_HB_D(6), set_selpar_slope_HB_D(7));
  f_priors+=neg_log_prior(selpar_afull_HB_D,set_selpar_afull_HB_D(5), set_selpar_afull_HB_D(6), set_selpar_afull_HB_D(7));
  f_priors+=neg_log_prior(selpar_sigma_HB_D,set_selpar_sigma_HB_D(5), set_selpar_sigma_HB_D(6), set_selpar_sigma_HB_D(7));

  f_priors+=neg_log_prior(selpar_L51_GR,set_selpar_L51_GR(5), set_selpar_L51_GR(6), set_selpar_L51_GR(7));
  f_priors+=neg_log_prior(selpar_slope1_GR,set_selpar_slope1_GR(5), set_selpar_slope1_GR(6), set_selpar_slope1_GR(7));
  f_priors+=neg_log_prior(selpar_L52_GR,set_selpar_L52_GR(5), set_selpar_L52_GR(6), set_selpar_L52_GR(7));
  f_priors+=neg_log_prior(selpar_slope2_GR,set_selpar_slope2_GR(5), set_selpar_slope2_GR(6), set_selpar_slope2_GR(7));

  f_priors+=neg_log_prior(log_q_cH,set_log_q_cH(5),set_log_q_cH(6),set_log_q_cH(7));
  f_priors+=neg_log_prior(log_q_mcvt,set_log_q_mcvt(5),set_log_q_mcvt(6),set_log_q_mcvt(7)); 
//  f_priors+=neg_log_prior(log_q_vid,set_log_q_vid(5),set_log_q_vid(6),set_log_q_vid(7));  
  f_priors+=neg_log_prior(log_q_HB,set_log_q_HB(5),set_log_q_HB(6),set_log_q_HB(7));
  f_priors+=neg_log_prior(log_q_GR,set_log_q_GR(5),set_log_q_GR(6),set_log_q_GR(7));
  
  f_priors+=neg_log_prior(F_init,set_F_init(5),set_F_init(6),set_F_init(7));
//  f_priors+=neg_log_prior(log_avg_F_cH,set_log_avg_F_cH(5),set_log_avg_F_cH(6),set_log_avg_F_cH(7));
//  f_priors+=neg_log_prior(log_avg_F_cL,set_log_avg_F_cL(5),set_log_avg_F_cL(6),set_log_avg_F_cL(7));
//  f_priors+=neg_log_prior(log_avg_F_HB,set_log_avg_F_HB(5),set_log_avg_F_HB(6),set_log_avg_F_HB(7));
//  f_priors+=neg_log_prior(log_avg_F_GR,set_log_avg_F_GR(5),set_log_avg_F_GR(6),set_log_avg_F_GR(7));
  
  fval+=f_priors;

  
  //cout << "fval = " << fval << "  fval_data = " << fval_data << endl;
  //cout << endl;

//----------------------------------------------------------------------------------
//Logistic function: 2 parameters
FUNCTION dvar_vector logistic(const dvar_vector& ages, const dvariable& L50, const dvariable& slope)
  //ages=vector of ages, L50=age at 50% selectivity, slope=rate of increase
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1./(1.+mfexp(-1.*slope*(ages-L50))); //logistic;  
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

//-----------------------------------------------------------------------------------                                                                           
//Logistic-exponential: 4 parameters (but 1 is fixed)
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
FUNCTION dvar_vector logistic_joint(const dvar_vector& ages, const dvariable& L501, const dvariable& slope1, const dvariable& L502, const dvariable& slope2, const dvariable& satval, const dvariable& joint)
  //ages=vector of ages, L501=age at 50% sel (ascending limb), slope1=rate of increase,L502=age at 50% sel (descending), slope1=rate of increase (ascending), 
  //satval=saturation value of descending limb, joint=location in age vector to join curves (may equal age or age + 1 if age-0 is included)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1.0; 
  for (iage=1; iage<=nages; iage++)
  {
   if (double(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope1*(ages(iage)-L501)));}  
   if (double(iage)>joint){Sel_Tmp(iage)=1.0-(1.0-satval)/(1.+mfexp(-1.*slope2*(ages(iage)-L502)));}  
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
//Spawner-recruit function (Beverton-Holt or Ricker)
FUNCTION dvariable SR_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& SSB, int func)
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, SSB=spawning biomass
  //func=1 for Beverton-Holt, 2 for Ricker
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  switch(func) {
    case 1: //Beverton-Holt
      Recruits_Tmp=((0.8*R0*h*SSB)/(0.2*R0*spr_F0*(1.0-h)+(h-0.2)*SSB));       
    break;
    case 2: //Ricker
      Recruits_Tmp=((SSB/spr_F0)*mfexp(h*(1-SSB/(R0*spr_F0))));       
    break;
  }
  RETURN_ARRAYS_DECREMENT();
  return Recruits_Tmp;
  
//-----------------------------------------------------------------------------------    
//Spawner-recruit equilibrium function (Beverton-Holt or Ricker)
FUNCTION dvariable SR_eq_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& spr_F, const dvariable& BC, int func)
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, spr_F=spawners per recruit @ F, BC=bias correction
  //func=1 for Beverton-Holt, 2 for Ricker
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  switch(func) {
    case 1: //Beverton-Holt
      Recruits_Tmp=(R0/((5.0*h-1.0)*spr_F))*(BC*4.0*h*spr_F-spr_F0*(1.0-h));    
    break;
    case 2: //Ricker
      Recruits_Tmp=R0/(spr_F/spr_F0)*(1.0+log(BC*spr_F/spr_F0)/h);      
    break;
  }
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
  //small_number is small value to avoid log(0) during search
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.00001;
  dvar_vector var(cv.indexmin(),cv.indexmax()); //variance in log space
  var=log(1.0+square(cv/wgt_dat));   // convert cv in arithmetic space to variance in log space
  LkvalTmp=sum(0.5*elem_div(square(log(elem_div((pred+small_number),(obs+small_number)))),var) );
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: multinomial
FUNCTION dvariable lk_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const double& minSS, const dvariable& wgt_dat)
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.00001;   //KC  
  LkvalTmp=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {LkvalTmp-=wgt_dat*nsamp(ii)*sum(elem_prod((obs_comp(ii)+small_number),
               log(elem_div((pred_comp(ii)+small_number), (obs_comp(ii)+small_number)))));
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;


//-----------------------------------------------------------------------------------
//Likelihood contribution: multinomial
FUNCTION dvariable lk_robust_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& wgt_dat)
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.00001;   //KC 
  LkvalTmp=0.0;
  dvar_matrix Eprime=elem_prod((1.0-obs_comp), obs_comp)+0.1/mbin; //E' of Francis 2011, p.1131  
  dvar_vector nsamp_wgt=nsamp*wgt_dat;
  //cout<<nsamp_wgt<<endl;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {LkvalTmp+= sum(0.5*log(Eprime(ii))-log(small_number+mfexp(elem_div((-square(obs_comp(ii)-pred_comp(ii))) , (Eprime(ii)*2.0/nsamp_wgt(ii)) ))) );
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;


//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//Likelihood contribution: priors
FUNCTION  dvariable neg_log_prior(dvariable pred, const double& prior, dvariable var, int pdf)
  //prior=prior point estimate, var=variance (if negative, treated as CV in arithmetic space), pred=predicted value, pdf=prior type (1=none, 2=lognormal, 3=normal, 4=beta)
    dvariable LkvalTmp;
    dvariable alpha, beta, ab_iq;
    dvariable big_number=1e10;
    LkvalTmp=0.0;
    // compute generic pdf's
    switch(pdf) {
        case 1: //option to turn off prior
          LkvalTmp=0.0;
          break;
        case 2: // lognormal
          if(prior<=0.0) cout << "YIKES: Don't use a lognormal distn for a negative prior" << endl;
          else if(pred<=0) LkvalTmp=big_number=1e10;
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
          else LkvalTmp=big_number;
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
REPORT_SECTION

  if (last_phase())  
  {
      cout<<"start report"<<endl;
      get_weighted_current();
      cout<<"got weighted"<<endl;
      get_msy();
      cout<<"got msy"<<endl;
 get_per_recruit_stuff();
      get_miscellaneous_stuff();
      cout<<"got misc stuff"<<endl;
      //cout<<"got per recruit"<<endl;  
      get_effective_sample_sizes();
     
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
      //cout << "BC Fmsy=" << F_msy_out<< "   BC SSBmsy=" << SSB_msy_out <<endl;
      cout <<"F status="<<FdF_msy_end<<endl;
      cout <<"Pop status="<<SdSSB_msy_end<<endl;
      cout << "h="<<steep<<"   R0="<<R0<<endl;
      //cout << "len_cv = "<<len_cv_val<<endl;
      //cout << "xdum " << xdum << endl;
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;  
     // cout << F_initial << endl;
      
      report << "TotalLikelihood " << fval << endl;
      report << "N" << endl;
      report << N<<endl;
      report << "F" << endl;
      report << F <<endl;
//        report <<"lenprob" <<endl;
//        report << lenprob<<endl;      

      sdnr_lc_cH=sdnr_multinomial(nyr_cH_lenc, lenbins, nsamp_cH_lenc, pred_cH_lenc, obs_cH_lenc, w_lc_cH);  
      sdnr_lc_mcvt=sdnr_multinomial(nyr_mcvt_lenc, lenbins, nsamp_mcvt_lenc, pred_mcvt_lenc, obs_mcvt_lenc, w_lc_mcvt);  
      sdnr_lc_HB=sdnr_multinomial(nyr_HB_lenc, lenbins, nsamp_HB_lenc, pred_HB_lenc, obs_HB_lenc, w_lc_HB); 
      sdnr_lc_HB_D=sdnr_multinomial(nyr_HB_D_lenc, lenbins, nsamp_HB_D_lenc, pred_HB_D_lenc, obs_HB_D_lenc, w_lc_HB_D);   
      sdnr_lc_GR=sdnr_multinomial(nyr_GR_lenc, lenbins, nsamp_GR_lenc, pred_GR_lenc, obs_GR_lenc, w_lc_GR);  

      sdnr_ac_cH=sdnr_multinomial(nyr_cH_agec, agebins_agec, nsamp_cH_agec, pred_cH_agec, obs_cH_agec, w_ac_cH);  
      sdnr_ac_mcvt=sdnr_multinomial(nyr_mcvt_agec, agebins_agec, nsamp_mcvt_agec, pred_mcvt_agec, obs_mcvt_agec, w_ac_mcvt);  
      sdnr_ac_HB=sdnr_multinomial(nyr_HB_agec, agebins_agec, nsamp_HB_agec, pred_HB_agec, obs_HB_agec, w_ac_HB);  
//      sdnr_ac_GR=sdnr_multinomial(nyr_GR_agec, agebins_agec, nsamp_GR_agec, pred_GR_agec, obs_GR_agec, w_ac_GR);

//      sdnr_I_cH=sdnr_lognormal(pred_cH_cpue, obs_cH_cpue, cH_cpue_cv, w_I_cH);
      sdnr_I_mcvt=sdnr_lognormal(pred_mcvt_cpue, obs_mcvt_cpue, mcvt_cpue_cv, w_I_mcvt);  
//      sdnr_I_vid=sdnr_lognormal(pred_vid_cpue, obs_vid_cpue, vid_cpue_cv, w_I_vid);  
//      sdnr_I_HB=sdnr_lognormal(pred_HB_cpue, obs_HB_cpue, HB_cpue_cv, w_I_HB);
//      sdnr_I_GR=sdnr_lognormal(pred_GR_cpue, obs_GR_cpue, GR_cpue_cv, w_I_GR);
      
      
      //#################################################################################################
      //##  Passing parameters to vector for bounds check plotting
      //################################################################################################# 
       Linf_out(8)=Linf; Linf_out(1,7)=set_Linf; 
       K_out(8)=K; K_out(1,7)=set_K;
       t0_out(8)=t0; t0_out(1,7)=set_t0;
       len_cv_val_out(8)=len_cv_val; len_cv_val_out(1,7)=set_len_cv;

       Linf_L_out(8)=Linf_L; Linf_L_out(1,7)=set_Linf_L; 
       K_L_out(8)=K_L; K_L_out(1,7)=set_K_L;
       t0_L_out(8)=t0_L; t0_L_out(1,7)=set_t0_L;
       len_cv_val_L_out(8)=len_cv_val_L; len_cv_val_L_out(1,7)=set_len_cv_L;       
       
       Linf_mcvt_out(8)=Linf_mcvt; Linf_mcvt_out(1,7)=set_Linf_mcvt; 
       K_mcvt_out(8)=K_mcvt; K_mcvt_out(1,7)=set_K_mcvt;
       t0_mcvt_out(8)=t0_mcvt; t0_mcvt_out(1,7)=set_t0_mcvt;
       len_cv_val_mcvt_out(8)=len_cv_val_mcvt; len_cv_val_mcvt_out(1,7)=set_len_cv_mcvt;

       log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
       M_constant_out(8)=M_constant; M_constant_out(1,7)=set_M_constant;
       steep_out(8)=steep; steep_out(1,7)=set_steep;
       rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
       R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;
       
       selpar_L50_cH1_out(8)=selpar_L50_cH1; selpar_L50_cH1_out(1,7)=set_selpar_L50_cH1;
       selpar_slope_cH1_out(8)=selpar_slope_cH1; selpar_slope_cH1_out(1,7)=set_selpar_slope_cH1;
//       selpar_L50_cH2_out(8)=selpar_L50_cH2; selpar_L50_cH2_out(1,7)=set_selpar_L50_cH2;
//       selpar_slope_cH2_out(8)=selpar_slope_cH2; selpar_slope_cH2_out(1,7)=set_selpar_slope_cH2;
//       selpar_L50_cH3_out(8)=selpar_L50_cH3; selpar_L50_cH3_out(1,7)=set_selpar_L50_cH3;
//       selpar_slope_cH3_out(8)=selpar_slope_cH3; selpar_slope_cH3_out(1,7)=set_selpar_slope_cH3;

       selpar_L50_mcvt_out(8)=selpar_L50_mcvt; selpar_L50_mcvt_out(1,7)=set_selpar_L50_mcvt;
       selpar_slope_mcvt_out(8)=selpar_slope_mcvt; selpar_slope_mcvt_out(1,7)=set_selpar_slope_mcvt;
 
       selpar_L51_HB1_out(8)=selpar_L51_HB1; selpar_L51_HB1_out(1,7)=set_selpar_L51_HB1;
       selpar_slope1_HB1_out(8)=selpar_slope1_HB1; selpar_slope1_HB1_out(1,7)=set_selpar_slope1_HB1;
       selpar_L52_HB1_out(8)=selpar_L52_HB1; selpar_L52_HB1_out(1,7)=set_selpar_L52_HB1;
       selpar_slope2_HB1_out(8)=selpar_slope2_HB1; selpar_slope2_HB1_out(1,7)=set_selpar_slope2_HB1;

       selpar_L50_HB2_out(8)=selpar_L50_HB2; selpar_L50_HB2_out(1,7)=set_selpar_L50_HB2;
       selpar_slope_HB2_out(8)=selpar_slope_HB2; selpar_slope_HB2_out(1,7)=set_selpar_slope_HB2;
       selpar_afull_HB2_out(8)=selpar_afull_HB2; selpar_afull_HB2_out(1,7)=set_selpar_afull_HB2;
       selpar_sigma_HB2_out(8)=selpar_sigma_HB2; selpar_sigma_HB2_out(1,7)=set_selpar_sigma_HB2;
       
//       selpar_L50_HB3_out(8)=selpar_L50_HB3; selpar_L50_HB3_out(1,7)=set_selpar_L50_HB3;
//       selpar_slope_HB3_out(8)=selpar_slope_HB3; selpar_slope_HB3_out(1,7)=set_selpar_slope_HB3;

       selpar_L50_HB_D_out(8)=selpar_L50_HB_D; selpar_L50_HB_D_out(1,7)=set_selpar_L50_HB_D;
       selpar_slope_HB_D_out(8)=selpar_slope_HB_D; selpar_slope_HB_D_out(1,7)=set_selpar_slope_HB_D;
       selpar_afull_HB_D_out(8)=selpar_afull_HB_D; selpar_afull_HB_D_out(1,7)=set_selpar_afull_HB_D;
       selpar_sigma_HB_D_out(8)=selpar_sigma_HB_D; selpar_sigma_HB_D_out(1,7)=set_selpar_sigma_HB_D;
//
       selpar_L51_GR_out(8)=selpar_L51_GR; selpar_L51_GR_out(1,7)=set_selpar_L51_GR;
       selpar_slope1_GR_out(8)=selpar_slope1_GR; selpar_slope1_GR_out(1,7)=set_selpar_slope1_GR;
       selpar_L52_GR_out(8)=selpar_L52_GR; selpar_L52_GR_out(1,7)=set_selpar_L52_GR;
       selpar_slope2_GR_out(8)=selpar_slope2_GR; selpar_slope2_GR_out(1,7)=set_selpar_slope2_GR;

       log_q_cH_out(8)=log_q_cH; log_q_cH_out(1,7)=set_log_q_cH;
       log_q_mcvt_out(8)=log_q_mcvt; log_q_mcvt_out(1,7)=set_log_q_mcvt;  
//       log_q_vid_out(8)=log_q_vid; log_q_vid_out(1,7)=set_log_q_vid;  
       log_q_HB_out(8)=log_q_HB; log_q_HB_out(1,7)=set_log_q_HB;
       log_q_GR_out(8)=log_q_GR; log_q_GR_out(1,7)=set_log_q_GR;
                     
       log_avg_F_cH_out(8)=log_avg_F_cH; log_avg_F_cH_out(1,7)=set_log_avg_F_cH;
       log_avg_F_HB_out(8)=log_avg_F_HB; log_avg_F_HB_out(1,7)=set_log_avg_F_HB;
       log_avg_F_GR_out(8)=log_avg_F_GR; log_avg_F_GR_out(1,7)=set_log_avg_F_GR;       
       log_avg_F_HB_D_out(8)=log_avg_F_HB_D; log_avg_F_HB_D_out(1,7)=set_log_avg_F_HB_D;
       log_avg_F_GR_D_out(8)=log_avg_F_GR_D; log_avg_F_GR_D_out(1,7)=set_log_avg_F_GR_D;
       F_init_out(8)=F_init; F_init_out(1,7)=set_F_init;
       
       log_rec_dev_out(styr_rec_dev, endyr_rec_dev)=log_rec_dev;
       log_F_dev_cH_out(styr_cH_L,endyr_cH_L)=log_F_dev_cH;
       log_F_dev_HB_out(styr_HB_L,endyr_HB_L)=log_F_dev_HB;
       log_F_dev_GR_out(styr_GR_L,endyr_GR_L)=log_F_dev_GR;
       log_F_dev_HB_D_out(styr_HB_D,endyr_HB_D)=log_F_dev_HB_D;
       log_F_dev_GR_D_out(styr_GR_D,endyr_GR_D)=log_F_dev_GR_D;
           
    #include "gtbase_make_Robject.cxx"   // write the R-compatible report

  } //endl last phase loop     
  
