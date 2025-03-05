//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##
//##  Beaufort Assessment Model (BAM) 
//##  NMFS, Beaufort Lab, Atlantic Fisheries Branch
//##  Application to Simulated Data
//##  Updated: Feb 2025. 
//##  Model has one survey and three fleets, first two are landings fleets, the third is a discard fleet.
//##  
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

// Starting and ending year of the model (year data start)
init_int styr;
init_int endyr;

//Starting year to estimate recruitment deviation from S-R curve
init_int styr_rec_dev;
//Ending year to estimate recruitment deviation from S-R curve
init_int endyr_rec_dev;
//possible 3 phases of constraints on recruitment deviations, turned on/off with model weights
init_int endyr_rec_phase1;
init_int endyr_rec_phase2;

// ending years for selectivity blocks
init_int endyr_selex_phase1; //code included, but not used in the simulated data

//number assessment years
number nyrs;
number nyrs_rec;
//this section MUST BE INDENTED!!!
 LOCAL_CALCS
    nyrs=endyr-styr+1.;
    nyrs_rec=endyr_rec_dev-styr_rec_dev+1.;
 END_CALCS

//Switch for spawner-recruit function (integer 1=Beverton-Holt, 2=Ricker, 3=Null (constant average recruitment))
init_int SR_switch;
//SPR proxy (as a ratio)
init_number spr_proxy; 
//Total number of ages in population model
init_int nages;
//Vector of ages for age bins in population model
init_vector agebins(1,nages); 

//Total number of ages used to match age comps: plus group may differ from popn, 
init_int nages_agec;
//Vector of ages for age bins in age comps LL
init_vector agebins_agec(1,nages_agec);

//Total number of length bins used to match length comps
init_int nlenbins;
init_int lenbins_width; // width of length bins (30mm) for use in age at length dist function
//Vector of lengths for length bins in length comps
init_vector lenbins(1,nlenbins); // 210-990 mm

//Max F used in spr and msy calcs
init_number max_F_spr_msy;
//Total number of iterations for spr calcs
init_int n_iter_spr;
//Total number of iterations for msy calcs
int n_iter_msy;
 LOCAL_CALCS
		n_iter_msy=n_iter_spr; 
 END_CALCS
//Number years at end of time series over which to average sector F's, for weighted selectivities; not used with a single fleet
init_int selpar_n_yrs_wgted;
//bias correction (set to 1.0 for no bias correction or a negative value to compute from rec variance)
init_number set_BiasCor;

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: observed data section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//###################Survey 1 ###########################################################
//CPUE: Survey1
init_int nyr_survey1_cpue;                                             
init_ivector yrs_survey1_cpue(1,nyr_survey1_cpue);                                            
init_vector obs_survey1_cpue(1,nyr_survey1_cpue);//Observed CPUE
init_vector survey1_cpue_cv(1,nyr_survey1_cpue); //CV of cpue

// Survey1 Age Compositions
init_int nyr_survey1_agec;
init_ivector yrs_survey1_agec(1,nyr_survey1_agec);
init_vector nsamp_survey1_agec(1,nyr_survey1_agec);
init_vector nfish_survey1_agec(1,nyr_survey1_agec);
init_matrix obs_survey1_agec(1,nyr_survey1_agec,1,nages_agec);

//################Fleet1 #######################################
// Landings 
init_int styr_fleet1_L;
init_int endyr_fleet1_L;
init_vector obs_fleet1_L(styr_fleet1_L,endyr_fleet1_L);
init_vector fleet1_L_cv(styr_fleet1_L,endyr_fleet1_L);
 
// Fleet1 landings age compositions
init_int nyr_fleet1_agec;
init_ivector yrs_fleet1_agec(1,nyr_fleet1_agec);
init_vector nsamp_fleet1_agec(1,nyr_fleet1_agec);
init_vector nfish_fleet1_agec(1,nyr_fleet1_agec);
init_matrix obs_fleet1_agec(1,nyr_fleet1_agec,1,nages_agec);

//################Fleet2 #######################################
// Landings 
init_int styr_fleet2_L;
init_int endyr_fleet2_L;
init_vector obs_fleet2_L(styr_fleet2_L,endyr_fleet2_L);
init_vector fleet2_L_cv(styr_fleet2_L,endyr_fleet2_L);
 
// Fleet2 landings age compositions
init_int nyr_fleet2_agec;
init_ivector yrs_fleet2_agec(1,nyr_fleet2_agec);
init_vector nsamp_fleet2_agec(1,nyr_fleet2_agec);
init_vector nfish_fleet2_agec(1,nyr_fleet2_agec);
init_matrix obs_fleet2_agec(1,nyr_fleet2_agec,1,nages_agec);

//################Fleet3  #######################################
// Dead discards (no age comps)
init_int styr_fleet3_D;
init_int endyr_fleet3_D;
init_vector obs_fleet3_D(styr_fleet3_D,endyr_fleet3_D);
init_vector fleet3_D_cv(styr_fleet3_D,endyr_fleet3_D);

// Fleet3 discards length compositions
init_int nyr_fleet3_lenc;
init_ivector yrs_fleet3_lenc(1,nyr_fleet3_lenc);
init_vector nsamp_fleet3_lenc(1,nyr_fleet3_lenc);
init_vector nfish_fleet3_lenc(1,nyr_fleet3_lenc);
init_matrix obs_fleet3_lenc(1,nyr_fleet3_lenc,1,nlenbins); // nyear x n length bins matrix of length comps

//################ Close kin mark recapture  #######################################
//init_int fyr_ckmr;                                    //first year of ckmr sampling
//init_int lyr_ckmr;                                    //last year of ckmr sampling
//init_int fage_ckmr;                                   //earliest age for ckmr
//init_int lage_ckmr;                                   //last age for ckmr
//init_vector nages_ckmr(fage_ckmr,lage_ckmr);          // vector of ages for ckmr
//init_int f_ckmr;                                      // index for first ckmr sample, i.e., 1
//init_int e_ckmr;                                      // final ckmr sample, i.e., however many pairs there are
//init_vector c1(f_ckmr,e_ckmr);                        // vector for cohort 1 (for hsps), e.g., vector of cohort IDs
//init_vector c2(f_ckmr,e_ckmr);                        // vector for cohort 2 (for hsps)
//init_vector s1(f_ckmr,e_ckmr);                        // vector of first-born sampling year (for POPs)
//init_vector HSPs(f_ckmr,e_ckmr);                      // vector of half-sibling pairs, e.g., vector of pairs within cohort
//init_vector POPs(f_ckmr,e_ckmr);                      // vector of parent-offspring pairs    
//init_vector not(f_ckmr,e_ckmr);                       // vector of pairs that are not hsps or pops (lots of these)

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: parameter section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//##################Single Parameter values and initial guesses #################################
// von Bert parms in TL mm population
init_vector set_Linf(1,7);
init_vector set_K(1,7);
init_vector set_t0(1,7);
init_vector set_len_cv(1,7);

//M scale used to scale age-dependent M. Logit transformed to stay on the space of (0,2) as a multiplier of age-dependent M.
init_vector set_M_logit(1,7);    

//Spawner-recruit parameters (Initial guesses or fixed values)
init_vector set_steep(1,7);         //recruitment steepness
init_vector set_log_R0(1,7);        //recruitment R0
init_vector set_R_autocorr(1,7);    //recruitment autocorrelation
init_vector set_rec_sigma(1,7);     //recruitment standard deviation in log space

init_vector set_log_dm_fleet1_ac(1,7);     //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_fleet2_ac(1,7);     //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_fleet3_lc(1,7);     //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_survey1_ac(1,7);    //Dirichlet-multinomial overdispersion parameter

//Initial guesses or fixed values of estimated selectivity parameters
init_vector set_selpar_A50_fleet1_B1(1,7);  //Fleet1 A50 (B1 for Block 1, multiple blocks not used in this simulation)
init_vector set_selpar_slope_fleet1_B1(1,7);//Fleet 1 slope
init_vector set_selpar_A50_fleet2_B1(1,7);  //Fleet 2 A50
init_vector set_selpar_slope_fleet2_B1(1,7);//Fleet 2 slope
init_vector set_selpar_A50_asc_fleet3_B1(1,7);//Fleet 3 A50 ascending limb
init_vector set_selpar_slope_asc_fleet3_B1(1,7);//Fleet 3 slope ascending limb
init_vector set_selpar_A50_dsc_fleet3_B1(1,7);//Fleet 3 A50 descending limb
init_vector set_selpar_slope_dsc_fleet3_B1(1,7);//Fleet 3 slope descending limb

init_vector set_selpar_A50_survey1(1,7);  //Survey A50
init_vector set_selpar_slope_survey1(1,7);//Survey slope

//--index catchability-----------------------------------------------------------------------------------
init_vector set_log_q_survey1(1,7);     //catchability coefficient (log) for survey index

//initial F
init_vector set_F_init(1,7);  //scales initial F; multiplicative scalar on early F 

//mean F in log space --------------------------------
init_vector set_log_avg_F_fleet1(1,7);
init_vector set_log_avg_F_fleet2(1,7);
init_vector set_log_avg_F_fleet3(1,7);

//close kin mark recapture probabilities and probability vectors-----------------------------------------
//init_number PHSP;
//init_number PPOP;
//init_vector PHSP_vec(f_ckmr,e_ckmr);
//init_vector PPOP_vec(f_ckmr,e_ckmr);

//##################Dev Vector constraints, bounds, and parameter values #################################
//F vectors---------------------------
init_vector set_log_F_dev_fleet1(1,3); 
init_vector set_log_F_dev_fleet2(1,3); 
init_vector set_log_F_dev_fleet3(1,3);
init_vector set_log_rec_dev(1,3);
init_vector set_log_Nage_dev(1,3);

init_vector set_log_F_dev_fleet1_vals(styr_fleet1_L,endyr_fleet1_L);
init_vector set_log_F_dev_fleet2_vals(styr_fleet2_L,endyr_fleet2_L);
init_vector set_log_F_dev_fleet3_vals(styr_fleet3_D,endyr_fleet3_D);
init_vector set_log_rec_dev_vals(styr_rec_dev,endyr_rec_dev);
init_vector set_log_Nage_dev_vals(2,nages);             

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: likelihood weights section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

init_number set_w_L;            //weight for fleet1, fleet2 landings
init_number set_w_D;            //weight for fleet 3 (discards)
init_number set_w_I_survey1;    //weight for survey1 index

init_number set_w_ac_fleet1;    //weight for fleet1 landings age comps
init_number set_w_ac_fleet2;    //weight for fleet2 landings age comps
init_number set_w_lc_fleet3;    //weight for fleet3 discard length comps
init_number set_w_ac_survey1;   //weight for survey1 age comps

init_number set_w_Nage_init;    //for fitting initial abundance at age (excluding first age)
init_number set_w_rec;          //for fitting S-R curve
init_number set_w_rec_early;    //additional constraint on early years recruitment; set to zero if not desired
init_number set_w_rec_end;      //additional constraint on ending years recruitment; set to zero if not desired
init_number set_w_Ftune;        //weight applied to tuning F (removed in final phase of optimization)

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: miscellaneous stuff section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//TL-weight relationship: W=aL^b
init_number wgtpar_a;  //Fixed value
init_number wgtpar_b;  //Fixed value 

//Proportion female at age
init_vector prop_f_obs(1,nages);           //proportion female at age

//Proportion mature female at age
init_vector maturity_f_obs(1,nages);       //proportion females mature at age

//Fecundity at length relationship
//init_vector set_fecpar_batches(1,nages); //leaving batches out for simplicity for the time being - not sure how/if we need/want to use it for ckmr
init_number set_fecpar_a;
init_number set_fecpar_b;
init_number set_fecpar_c;
init_number set_fecpar_thresh;
init_number set_fecpar_min;
//init_number fecpar_scale;

// Natural mortality
init_vector set_M(1,nages);     //age-independent

init_number spawn_time_frac; //time of year of peak spawning, as a fraction of the year

//Tune Fapex (tuning removed in final year of optimization)
init_int set_Ftune_yr;
init_number set_Ftune;
//threshold sample sizes for age comps
init_number minSS_agec;
//threshold sample sizes for length comps
init_number minSS_lenc;

//ageing error matrix (columns are true ages, rows are ages as read for age comps: columns should sum to one)
init_matrix age_error(1,nages,1,nages);

//input for deterministic F-based projections
init_int endyr_proj;   //last year of projections, by default, first year is endyr+1
init_int styr_regs;    //apply current F until styr_regs, then the projection F
init_int Fproj_switch; //Value to use for Fproj: 1=current F, 2=Fmsy, 3=F30, 4=F40
init_number Fproj_mult;//multiplier to the Fproj
int styr_proj;
 LOCAL_CALCS
   styr_proj=endyr+1;
 END_CALCS



// #######Indexing integers for year(iyear), age(iage), length(ilen) ###############
int iyear;
int iage;
int ilen;
int ff;

number pi;
number sqrt2pi;
number g2mt;                    //conversion of grams to metric tons 
number kg2mt;                   //conversion of kilograms to metric tons 
number mt2kg;                   //conversion of metric tons to kilogramps 
number g2kg;                    //conversion of grams to kg   
number g2klb;                   //conversion of grams to 1000 lb   
number mt2klb;                  //conversion of metric tons to 1000 lb
number mt2lb;                   //conversion of metric tons to lb
number dzero;                   //small additive constant to prevent division by zero
number huge_number;             //huge number, to avoid irregular parameter space
//number phsp;                  //intermediate probability of being a half-sibling pair
//number ppop;                  //intermediate probability of being a parent-offspring pair

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
  const double Linf_LO=set_Linf(2); const double Linf_HI=set_Linf(3); const double Linf_PH=set_Linf(4);
  const double K_LO=set_K(2); const double K_HI=set_K(3); const double K_PH=set_K(4);
  const double t0_LO=set_t0(2); const double t0_HI=set_t0(3); const double t0_PH=set_t0(4);  
  const double len_cv_LO=set_len_cv(2); const double len_cv_HI=set_len_cv(3); const double len_cv_PH=set_len_cv(4); 

  const double M_logit_LO=set_M_logit(2); const double M_logit_HI=set_M_logit(3); const double M_logit_PH=set_M_logit(4);

  const double steep_LO=set_steep(2); const double steep_HI=set_steep(3); const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); const double log_R0_HI=set_log_R0(3); const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); const double R_autocorr_HI=set_R_autocorr(3); const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); const double rec_sigma_HI=set_rec_sigma(3); const double rec_sigma_PH=set_rec_sigma(4);
  
  const double selpar_A50_survey1_LO=set_selpar_A50_survey1(2); const double selpar_A50_survey1_HI=set_selpar_A50_survey1(3); const double selpar_A50_survey1_PH=set_selpar_A50_survey1(4);
  const double selpar_slope_survey1_LO=set_selpar_slope_survey1(2); const double selpar_slope_survey1_HI=set_selpar_slope_survey1(3); const double selpar_slope_survey1_PH=set_selpar_slope_survey1(4);

  const double selpar_A50_fleet1_B1_LO=set_selpar_A50_fleet1_B1(2); const double selpar_A50_fleet1_B1_HI=set_selpar_A50_fleet1_B1(3); const double selpar_A50_fleet1_B1_PH=set_selpar_A50_fleet1_B1(4);
  const double selpar_slope_fleet1_B1_LO=set_selpar_slope_fleet1_B1(2); const double selpar_slope_fleet1_B1_HI=set_selpar_slope_fleet1_B1(3); const double selpar_slope_fleet1_B1_PH=set_selpar_slope_fleet1_B1(4);
  const double selpar_A50_fleet2_B1_LO=set_selpar_A50_fleet2_B1(2); const double selpar_A50_fleet2_B1_HI=set_selpar_A50_fleet2_B1(3); const double selpar_A50_fleet2_B1_PH=set_selpar_A50_fleet2_B1(4);
  const double selpar_slope_fleet2_B1_LO=set_selpar_slope_fleet2_B1(2); const double selpar_slope_fleet2_B1_HI=set_selpar_slope_fleet2_B1(3); const double selpar_slope_fleet2_B1_PH=set_selpar_slope_fleet2_B1(4);
  //fleet 3 double-logistic selectivity parameters
  const double selpar_A50_asc_fleet3_B1_LO=set_selpar_A50_asc_fleet3_B1(2); const double selpar_A50_asc_fleet3_B1_HI=set_selpar_A50_asc_fleet3_B1(3); const double selpar_A50_asc_fleet3_B1_PH=set_selpar_A50_asc_fleet3_B1(4);
  const double selpar_slope_asc_fleet3_B1_LO=set_selpar_slope_asc_fleet3_B1(2); const double selpar_slope_asc_fleet3_B1_HI=set_selpar_slope_asc_fleet3_B1(3); const double selpar_slope_asc_fleet3_B1_PH=set_selpar_slope_asc_fleet3_B1(4);
  const double selpar_A50_dsc_fleet3_B1_LO=set_selpar_A50_dsc_fleet3_B1(2); const double selpar_A50_dsc_fleet3_B1_HI=set_selpar_A50_dsc_fleet3_B1(3); const double selpar_A50_dsc_fleet3_B1_PH=set_selpar_A50_dsc_fleet3_B1(4);
  const double selpar_slope_dsc_fleet3_B1_LO=set_selpar_slope_dsc_fleet3_B1(2); const double selpar_slope_dsc_fleet3_B1_HI=set_selpar_slope_dsc_fleet3_B1(3); const double selpar_slope_dsc_fleet3_B1_PH=set_selpar_slope_dsc_fleet3_B1(4);

  const double log_q_survey1_LO=set_log_q_survey1(2); const double log_q_survey1_HI=set_log_q_survey1(3); const double log_q_survey1_PH=set_log_q_survey1(4);
  const double F_init_LO=set_F_init(2); const double F_init_HI=set_F_init(3); const double F_init_PH=set_F_init(4);
  const double log_avg_F_fleet1_LO=set_log_avg_F_fleet1(2); const double log_avg_F_fleet1_HI=set_log_avg_F_fleet1(3); const double log_avg_F_fleet1_PH=set_log_avg_F_fleet1(4);
  const double log_avg_F_fleet2_LO=set_log_avg_F_fleet2(2); const double log_avg_F_fleet2_HI=set_log_avg_F_fleet2(3); const double log_avg_F_fleet2_PH=set_log_avg_F_fleet2(4);
  const double log_avg_F_fleet3_LO=set_log_avg_F_fleet3(2); const double log_avg_F_fleet3_HI=set_log_avg_F_fleet3(3); const double log_avg_F_fleet3_PH=set_log_avg_F_fleet3(4);

  const double log_dm_fleet1_ac_LO=set_log_dm_fleet1_ac(2); const double log_dm_fleet1_ac_HI=set_log_dm_fleet1_ac(3); const double log_dm_fleet1_ac_PH=set_log_dm_fleet1_ac(4);
  const double log_dm_fleet2_ac_LO=set_log_dm_fleet2_ac(2); const double log_dm_fleet2_ac_HI=set_log_dm_fleet2_ac(3); const double log_dm_fleet2_ac_PH=set_log_dm_fleet2_ac(4);
  const double log_dm_fleet3_lc_LO=set_log_dm_fleet3_lc(2); const double log_dm_fleet3_lc_HI=set_log_dm_fleet3_lc(3); const double log_dm_fleet3_lc_PH=set_log_dm_fleet3_lc(4);
  const double log_dm_survey1_ac_LO=set_log_dm_survey1_ac(2); const double log_dm_survey1_ac_HI=set_log_dm_survey1_ac(3); const double log_dm_survey1_ac_PH=set_log_dm_survey1_ac(4);

  //-dev vectors-----------------------------------------------------------------------------------------------------------  
  const double log_F_dev_fleet1_LO=set_log_F_dev_fleet1(1); const double log_F_dev_fleet1_HI=set_log_F_dev_fleet1(2); const double log_F_dev_fleet1_PH=set_log_F_dev_fleet1(3);   
  const double log_F_dev_fleet2_LO=set_log_F_dev_fleet2(1); const double log_F_dev_fleet2_HI=set_log_F_dev_fleet2(2); const double log_F_dev_fleet2_PH=set_log_F_dev_fleet2(3);  
  const double log_F_dev_fleet3_LO=set_log_F_dev_fleet3(1); const double log_F_dev_fleet3_HI=set_log_F_dev_fleet3(2); const double log_F_dev_fleet3_PH=set_log_F_dev_fleet3(3);
  const double log_rec_dev_LO=set_log_rec_dev(1); const double log_rec_dev_HI=set_log_rec_dev(2); const double log_rec_dev_PH=set_log_rec_dev(3);          
  const double log_Nage_dev_LO=set_log_Nage_dev(1); const double log_Nage_dev_HI=set_log_Nage_dev(2); const double log_Nage_dev_PH=set_log_Nage_dev(3);          

 END_CALCS 

    
//--------------Growth--------------------------------------------------------------------------- 
  //Population growth parms and conversions
  init_bounded_number Linf(Linf_LO,Linf_HI,Linf_PH);
  init_bounded_number K(K_LO,K_HI,K_PH);
  init_bounded_number t0(t0_LO,t0_HI,t0_PH);
  init_bounded_number len_cv_val(len_cv_LO,len_cv_HI,len_cv_PH);  
  vector Linf_out(1,8);
  vector K_out(1,8);
  vector t0_out(1,8);
  vector len_cv_val_out(1,8);

  vector meanlen_TL(1,nages);        //mean total length (mm) at age all fish
  vector fecundity(1,nages);         // annual fec at age per mat female: fecpar_batches*(fecpar_a+(fecpar_b*TL^fecpar_c)), modified for age 1  
  //vector fecpar_batches(1,nages); //number of batches by age
  number fecpar_a;
  number fecpar_b;
  number fecpar_c;
  number fecpar_thresh; //smaller than this length, apply fecpar_min instead of the model-based value
  number fecpar_min;
 
  vector wgt_g(1,nages);        //whole wgt in g
  vector wgt_kg(1,nages);       //whole wgt in kg
  vector wgt_mt(1,nages);       //whole wgt in mt
  vector wgt_klb(1,nages);      //whole wgt in 1000 lb
  vector wgt_lb(1,nages);       //whole wgt in lb  

  //Natural mortality
  init_bounded_number M_logit(M_logit_LO,M_logit_HI,M_logit_PH);
  number M_scale; //Scalar on M in arithmetic space (back transformed from logit); constrained to (0,2)
  vector M_logit_out(1,8);

//----age-length key stuff
  matrix lenprob(1,nages,1,nlenbins);   //distribution of size at age (age-length key, 3 cm bins) in population
  number zscore_len;                    //standard normal values used for computing lenprob
  vector cprob_lenvec(1,nlenbins);      //cumulative probabilities used for computing lenprob
  number zscore_lzero;                  //standardized normal values for length = 0
  number cprob_lzero;                   //length probability mass below zero, used for computing lenprob

  matrix lenprob_fleet3(1,nages,1,nlenbins);  //dist of size at age in fleet 3

  vector len_sd(1,nages);
  vector len_cv(1,nages);
  
//----Predicted age compositions
  matrix pred_survey1_agec(1,nyr_survey1_agec,1,nages_agec);
  matrix pred_survey1_agec_allages(1,nyr_survey1_agec,1,nages);  
  matrix ErrorFree_survey1_agec(1,nyr_survey1_agec,1,nages);
  matrix pred_fleet1_agec(1,nyr_fleet1_agec,1,nages_agec);
  matrix pred_fleet1_agec_allages(1,nyr_fleet1_agec,1,nages);
  matrix ErrorFree_fleet1_agec(1,nyr_fleet1_agec,1,nages);  
  matrix pred_fleet2_agec(1,nyr_fleet2_agec,1,nages_agec);
  matrix pred_fleet2_agec_allages(1,nyr_fleet2_agec,1,nages);
  matrix ErrorFree_fleet2_agec(1,nyr_fleet2_agec,1,nages);  

//----Predicted length compositions
  matrix pred_fleet3_lenc(1,nyr_fleet3_lenc,1,nlenbins);
  
  vector obs_survey1_cpue_allyr(styr,endyr);
  vector pred_survey1_cpue_allyr(styr,endyr);
  vector survey1_cpue_cv_allyr(styr,endyr);
 
//sample size used to compute effective sample size for age comps
  vector nsamp_survey1_agec_allyr(styr,endyr);
  vector nsamp_fleet1_agec_allyr(styr,endyr);
  vector nsamp_fleet2_agec_allyr(styr,endyr);
  vector nsamp_fleet3_lenc_allyr(styr,endyr);
 
//Nfish used in uncertainty analysis (not used in fitting)
  vector nfish_survey1_agec_allyr(styr,endyr);
  vector nfish_fleet1_agec_allyr(styr,endyr);
  vector nfish_fleet2_agec_allyr(styr,endyr);
  vector nfish_fleet3_lenc_allyr(styr,endyr);
   
//Computed effective sample size for output (not used in fitting)
  vector neff_survey1_agec_allyr_out(styr,endyr);
  vector neff_fleet1_agec_allyr_out(styr,endyr);
  vector neff_fleet2_agec_allyr_out(styr,endyr);
  vector neff_fleet3_lenc_allyr_out(styr,endyr);

  
//-----Population-----------------------------------------------------------------------------------
  matrix N(styr,endyr+1,1,nages);           //Population numbers by year and age at start of yr
  matrix N_mdyr(styr,endyr,1,nages);        //Population numbers by year and age at mdpt of yr: used for fishery comps and cpue
  matrix N_spawn(styr,endyr,1,nages);       //Population numbers by year and age at peaking spawning: used for SSB  
  init_bounded_vector log_Nage_dev(2,nages,log_Nage_dev_LO,log_Nage_dev_HI,log_Nage_dev_PH);
  vector log_Nage_dev_output(1,nages);      //used in output. equals zero for first age
  matrix B(styr,endyr+1,1,nages);           //Population biomass by year and age at start of yr
  vector totB(styr,endyr+1);                //Total biomass by year
  vector totN(styr,endyr+1);                //Total abundance by year
  vector SSB(styr,endyr);                   //Total spawning capacity by year (may be based on fecundity)
  vector MatFemB(styr,endyr);               //Total spawning biomass by year (mature female biomass)  
  vector rec(styr,endyr+1);                 //Recruits by year
  vector prop_f(1,nages);					//Proportion female at age
  vector maturity_f(1,nages);               //Proportion of female mature at age
  vector reprod(1,nages);                   //vector used to compute spawning biomass 

//-----Close kin mark recapture-----------------------------------------------------------------------------------
//vector ro_denom(styr,endyr);              //denominator for reproductive output calcs
//matrix rel_ro(styr,endyr,1,nages);       //relative reproductive output

//---Stock-Recruit Function (Beverton-Holt, steepness parameterization)----------
  init_bounded_number log_R0(log_R0_LO,log_R0_HI,log_R0_PH);        //log(virgin Recruitment)
  vector log_R0_out(1,8);
  number R0;                                  //virgin recruitment
  init_bounded_number steep(steep_LO,steep_HI,steep_PH); //steepness
  vector steep_out(1,8);
  init_bounded_number rec_sigma(rec_sigma_LO,rec_sigma_HI,rec_sigma_PH);  //sd recruitment residuals  
  vector rec_sigma_out(1,8);
  init_bounded_number R_autocorr(R_autocorr_LO,R_autocorr_HI,R_autocorr_PH);  //autocorrelation in SR  
  vector R_autocorr_out(1,8);

  number rec_sigma_sq;                        //square of rec_sigma      
  number rec_logL_add;                        //additive term in -logL term   
 
  //vector of recruitment log devs
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
  vector SdS0(styr,endyr);                      //SSB / virgin SSB

  init_bounded_number log_dm_fleet1_ac(log_dm_fleet1_ac_LO,log_dm_fleet1_ac_HI,log_dm_fleet1_ac_PH);
  init_bounded_number log_dm_fleet2_ac(log_dm_fleet2_ac_LO,log_dm_fleet2_ac_HI,log_dm_fleet2_ac_PH);
  init_bounded_number log_dm_fleet3_lc(log_dm_fleet3_lc_LO,log_dm_fleet3_lc_HI,log_dm_fleet3_lc_PH);
  init_bounded_number log_dm_survey1_ac(log_dm_survey1_ac_LO,log_dm_survey1_ac_HI,log_dm_survey1_ac_PH);

  vector log_dm_fleet1_ac_out(1,8);
  vector log_dm_fleet2_ac_out(1,8);
  vector log_dm_fleet3_lc_out(1,8);
  vector log_dm_survey1_ac_out(1,8);

//-----------------------------------------------------------------------------------------------------------------------------------------------
////---Selectivity-------------------------------------------------------------------------

//survey1 selectivity-------------------------------------------------
  matrix sel_survey1(styr,endyr,1,nages);
  vector selvec_survey1(1,nages);
  init_bounded_number selpar_A50_survey1(selpar_A50_survey1_LO,selpar_A50_survey1_HI,selpar_A50_survey1_PH);
  init_bounded_number selpar_slope_survey1(selpar_slope_survey1_LO,selpar_slope_survey1_HI,selpar_slope_survey1_PH);
  vector selpar_A50_survey1_out(1,8);
  vector selpar_slope_survey1_out(1,8);
 
//fleet1 catch selectivity-------------------------------------------------
  matrix sel_fleet1(styr,endyr,1,nages);
  matrix sel_fleet2(styr,endyr,1,nages);
  vector selvec_fleet1_B1(1,nages);
  vector selvec_fleet2_B1(1,nages);
  init_bounded_number selpar_A50_fleet1_B1(selpar_A50_fleet1_B1_LO,selpar_A50_fleet1_B1_HI,selpar_A50_fleet1_B1_PH);
  init_bounded_number selpar_slope_fleet1_B1(selpar_slope_fleet1_B1_LO,selpar_slope_fleet1_B1_HI,selpar_slope_fleet1_B1_PH);
  init_bounded_number selpar_A50_fleet2_B1(selpar_A50_fleet2_B1_LO,selpar_A50_fleet2_B1_HI,selpar_A50_fleet2_B1_PH);
  init_bounded_number selpar_slope_fleet2_B1(selpar_slope_fleet2_B1_LO,selpar_slope_fleet2_B1_HI,selpar_slope_fleet2_B1_PH);
  //fleet 3 catch selectivity
  matrix sel_fleet3(styr,endyr,1,nages);
  vector selvec_fleet3_B1(1,nages);
  init_bounded_number selpar_A50_asc_fleet3_B1(selpar_A50_asc_fleet3_B1_LO,selpar_A50_asc_fleet3_B1_HI,selpar_A50_asc_fleet3_B1_PH);
  init_bounded_number selpar_slope_asc_fleet3_B1(selpar_slope_asc_fleet3_B1_LO,selpar_slope_asc_fleet3_B1_HI,selpar_slope_asc_fleet3_B1_PH);
  init_bounded_number selpar_A50_dsc_fleet3_B1(selpar_A50_dsc_fleet3_B1_LO,selpar_A50_dsc_fleet3_B1_HI,selpar_A50_dsc_fleet3_B1_PH);
  init_bounded_number selpar_slope_dsc_fleet3_B1(selpar_slope_dsc_fleet3_B1_LO,selpar_slope_dsc_fleet3_B1_HI,selpar_slope_dsc_fleet3_B1_PH);
  vector selpar_A50_fleet1_B1_out(1,8);
  vector selpar_slope_fleet1_B1_out(1,8);
  vector selpar_A50_fleet2_B1_out(1,8);
  vector selpar_slope_fleet2_B1_out(1,8);
  vector selpar_A50_asc_fleet3_B1_out(1,8);
  vector selpar_slope_asc_fleet3_B1_out(1,8);
  vector selpar_A50_dsc_fleet3_B1_out(1,8);
  vector selpar_slope_dsc_fleet3_B1_out(1,8);
	
//Weighted total selectivity--------------------------------------------  
  //F-weighted, recent selectivities
  vector sel_wgted_L(1,nages);  //toward landings 
  vector sel_wgted_D(1,nages);  //toward discards
  vector sel_wgted_tot(1,nages); //landings plus deads discards, used in computation of Z. Note, no discards in sim this study.

//-----------------------------------------------------------------------------------------------------------------------------------------------
//-------CPUE Predictions--------------------------------
  vector pred_survey1_cpue(1,nyr_survey1_cpue);                //predicted survey1 index (in numbers)
  matrix N_survey1(1,nyr_survey1_cpue,1,nages);                //used to compute survey1 index and age comps
 
//---Catchability (CPUE q's)----------------------------------------------------------
  init_bounded_number log_q_survey1(log_q_survey1_LO,log_q_survey1_HI,log_q_survey1_PH);
  vector log_q_survey1_out(1,8);
  
//----------------------------------------------------------------------------------------------------------------------------------------------- 

 vector log_F1_dev_output(styr,endyr);             //used in t.series output. equals zero except for yrs in log_F1_dev
 vector log_F2_dev_output(styr,endyr);             //used in t.series output. equals zero except for yrs in log_F2_dev
 vector log_F3_dev_output(styr,endyr);            //used in t.series output. equals zero except for yrs in log_F3_dev

//---Landings in numbers (total fish) and in wgt--------------------------------------------------
  matrix L_fleet1_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_fleet1_mt(styr,endyr,1,nages);    //landings (weight) at age    
  vector pred_fleet1_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages  
  vector pred_fleet1_L_mt(styr,endyr);       //yearly landings in mt summed over ages
 
  matrix L_fleet2_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_fleet2_mt(styr,endyr,1,nages);    //landings (weight) at age    
  vector pred_fleet2_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages  
  vector pred_fleet2_L_mt(styr,endyr);       //yearly landings in mt summed over ages
  
  matrix D_fleet3_num(styr,endyr,1,nages);   //dead discards (numbers) at age
  matrix D_fleet3_mt(styr,endyr,1,nages);    //dead discards (weight) at age    
  vector pred_fleet3_D_knum(styr,endyr);     //yearly dead discards in 1000 fish summed over ages  
  vector pred_fleet3_D_mt(styr,endyr);       //yearly dead discards in mt summed over ages
  
  matrix L_total_num(styr,endyr,1,nages);//total landings in number at age
  matrix L_total_mt(styr,endyr,1,nages); //landings in mt at age 
  vector L_total_knum_yr(styr,endyr);    //total landings in 1000 fish by yr summed over ages  
  vector L_total_mt_yr(styr,endyr);      //total landings (mt) by yr summed over ages
  matrix D_total_num(styr,endyr,1,nages);//total discards in number at age
  matrix D_total_mt(styr,endyr,1,nages);//total discards in mt at age
  vector D_total_knum_yr(styr,endyr);      //total discards in 1000 fish by yr summed over ages
  vector D_total_mt_yr(styr,endyr);        //total discards in mt by yr summed over ages
  
//---MSY calcs----------------------------------------------------------------------------
  number F_fleet1_prop;       //proportion of F_sum attributable to fleet1, last X yr
  number F_fleet2_prop;       //proportion of F_sum attributable to fleet2, last X yr
  number F_fleet3_prop;       //proportion of F_sum attributable to fleet3, last X yr
  
  number F_temp_sum;      //sum of geom mean Fsum's in last X yrs, used to compute F_fishery_prop

  vector F_end(1,nages);
  vector F_end_L(1,nages);   
  number F_end_apex;
  
  number SSB_msy_out;           //SSB (popn fecundity) at msy
  number F_msy_out;             //F at msy
  number msy_mt_out;           //max sustainable yield (mt)
  number msy_knum_out;          //max sustainable yield (1000 fish)  
  number B_msy_out;             //total biomass at MSY 
  number R_msy_out;             //equilibrium recruitment at F=Fmsy
  number spr_msy_out;           //spr at F=Fmsy

  number F30_dum;				//intermediate calculation for F30
  number F40_dum;				//intermediate calculation for F40
  number F50_dum;				//intermediate calculation for F50
  number Fproxy_dum;			//intermediate calculation for Fproxy  
  number F30_out;              	//F30
  number F40_out;              	//F40
  number F50_out;              	//F50
  number Fproxy_out;            //Fproxy
  number SSB_Fproxy_out;        //SSB, B, R, and L proxies assumes average recruitment taken from the assessment period 
  number B_Fproxy_out;
  number R_Fproxy_out;
  number L_Fproxy_knum_out;
  number L_Fproxy_mt_out;   
  number D_Fproxy_knum_out;    
  number D_Fproxy_mt_out;    

  vector wgt_wgted_L_mt(1,nages);   //fishery-weighted average weight at age of landings in whole weight
  vector wgt_wgted_D_mt(1,nages);   //fishery-weighted average weight at age of discards in whole weight

  vector N_age_msy(1,nages);         //numbers at age for MSY calculations: beginning of yr
  vector N_age_msy_spawn(1,nages);   //numbers at age for MSY calculations: time of peak spawning  
  vector L_age_msy(1,nages);         //landings at age for MSY calculations
  vector D_age_msy(1,nages);         //dead discards at age for MSY calculations
  vector Z_age_msy(1,nages);         //total mortality at age for MSY calculations
  vector F_L_age_msy(1,nages);       //fishing mortality landings (not discards) at age for MSY calculations
  vector F_D_age_msy(1,nages);       //fishing mortality discards at age for MSY calculations
  vector F_msy(1,n_iter_msy);        //values of full F to be used in equilibrium calculations
  vector spr_msy(1,n_iter_msy);      //reproductive capacity-per-recruit values corresponding to F values in F_msy
  vector R_eq(1,n_iter_msy);         //equilibrium recruitment values corresponding to F values in F_msy
  vector L_eq_mt(1,n_iter_msy);     //equilibrium landings(whole wgt) values corresponding to F values in F_msy
  vector L_eq_knum(1,n_iter_msy);    //equilibrium landings(1000 fish) values corresponding to F values in F_msy
  vector D_eq_mt(1,n_iter_msy);      //eq. discards (whole wgt) values corresponding to F values in F_msy
  vector D_eq_knum(1,n_iter_msy);    //eq. discards (1000 fish) values corresponding to F values in F_msy
  vector SSB_eq(1,n_iter_msy);       //equilibrium reproductive capacity values corresponding to F values in F_msy
  vector B_eq(1,n_iter_msy);         //equilibrium biomass values corresponding to F values in F_msy
  
  vector FdF_msy(styr,endyr);
  vector FdFproxy(styr,endyr);
  vector SdSSB_msy(styr,endyr);	 
  number SdSSB_msy_end;
  number FdF_msy_end;
  number FdF_msy_end_mean;           //geometric mean of last X yrs  
  vector SdSSB_Fproxy(styr,endyr);	 
  vector Sdmsst_Fproxy(styr,endyr);	 
  number SdSSB_Fproxy_end;
  number Sdmsst_Fproxy_end;
  number FdFproxy_end_mean;             //geometric mean of last selpar_n_yrs_wgted yrs  
  number Fend_mean_temp;			 //intermediate calc for geometric mean of last selpar_n_yrs_wgted yrs
  number Fend_mean;					 //geometric mean of last selpar_n_yrs_wgted yrs  
  vector L_age_Fproxy(1,nages);         //landings at age for Fproxy calculations
  vector D_age_Fproxy(1,nages);       //discards at age for Fproxy calculations
  number iter_inc_msy;               //increments used to compute msy, equals 1/(n_iter_msy-1)
  
//--------Mortality------------------------------------------------------------------

  vector M(1,nages);                         //age-dependent natural mortality
  number smsy2msst;                           //scales Smsy to get msst using 75%. Used only in output.
  
  matrix F(styr,endyr,1,nages);
  vector Fsum(styr,endyr);                   //Full fishing mortality rate by year
  vector Fapex(styr,endyr);                  //Max across ages, fishing mortality rate by year (may differ from Fsum bc of dome-shaped sel 
  matrix Z(styr,endyr,1,nages);
  
  init_bounded_number log_avg_F_fleet1(log_avg_F_fleet1_LO,log_avg_F_fleet1_HI,log_avg_F_fleet1_PH);
  vector log_avg_F_fleet1_out(1,8);  
  init_bounded_dev_vector log_F_dev_fleet1(styr_fleet1_L,endyr_fleet1_L,log_F_dev_fleet1_LO,log_F_dev_fleet1_HI,log_F_dev_fleet1_PH);
  vector log_F_dev_fleet1_out(styr_fleet1_L,endyr_fleet1_L);
  matrix F_fleet1(styr,endyr,1,nages);
  vector F_fleet1_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_fleet1;
  number log_F_dev_end_fleet1; 

  init_bounded_number log_avg_F_fleet2(log_avg_F_fleet2_LO,log_avg_F_fleet2_HI,log_avg_F_fleet2_PH);
  vector log_avg_F_fleet2_out(1,8);  
  init_bounded_dev_vector log_F_dev_fleet2(styr_fleet2_L,endyr_fleet2_L,log_F_dev_fleet2_LO,log_F_dev_fleet2_HI,log_F_dev_fleet2_PH);
  vector log_F_dev_fleet2_out(styr_fleet2_L,endyr_fleet2_L);
  matrix F_fleet2(styr,endyr,1,nages);
  vector F_fleet2_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_fleet2;
  number log_F_dev_end_fleet2; 
  
  init_bounded_number log_avg_F_fleet3(log_avg_F_fleet3_LO,log_avg_F_fleet3_HI,log_avg_F_fleet3_PH);
  vector log_avg_F_fleet3_out(1,8);
  init_bounded_dev_vector log_F_dev_fleet3(styr_fleet3_D,endyr_fleet3_D,log_F_dev_fleet3_LO,log_F_dev_fleet3_HI,log_F_dev_fleet3_PH);
  vector log_F_dev_fleet3_out(styr_fleet3_D,endyr_fleet3_D);
  matrix F_fleet3(styr,endyr,1,nages);
  vector F_fleet3_out(styr,endyr);  // use for intermediate calcs in fcn get_mortality
  number log_F_dev_init_fleet3;
  number log_F_dev_end_fleet3;

  init_bounded_number F_init(F_init_LO,F_init_HI,F_init_PH); //scales early F for initialization
  vector F_init_out(1,8); 
 
//---Per-recruit stuff----------------------------------------------------------------------------------
  vector N_age_spr(1,nages);         //numbers at age for SPR calculations: beginning of year
  vector N_age_spr_spawn(1,nages);   //numbers at age for SPR calculations: time of peak spawning  
  vector L_age_spr(1,nages);         //landings at age for SPR calculations
  vector D_age_spr(1,nages);         //discards at age for SPR calculations
  vector Z_age_spr(1,nages);         //total mortality at age for SPR calculations
  vector spr_static(styr,endyr);     //vector of static SPR values by year
  vector F_L_age_spr(1,nages);       //fishing mortality of landings (not discards) at age for SPR calculations
  vector F_D_age_spr(1,nages);       // ^ for discards in SPR calcs
  vector F_spr(1,n_iter_spr);        //values of full F to be used in per-recruit calculations
  vector D_spr(1,n_iter_spr);        // ^ for dicards in per-recruit calcs
  vector spr_spr(1,n_iter_spr);      //reproductive capacity-per-recruit values corresponding to F values in F_spr
  vector spr_ratio(1,n_iter_spr);    //reproductive capacity-per-recruit relative to spr_F0 values corresponding to F values in F_spr
  vector L_spr(1,n_iter_spr);        //landings(kg)-per-recruit (ypr) values corresponding to F values in F_spr
  
  vector N_spr_F0(1,nages);          //Used to compute spr at F=0: at time of peak spawning
  vector N_bpr_F0(1,nages);          //Used to compute bpr at F=0: at start of year  
  vector N_spr_initial(1,nages);     //Initial spawners per recruit at age given initial F
  vector N_initial_eq(1,nages);      //Initial equilibrium abundance at age
  vector F_initial(1,nages);         //initial F at age
  vector Z_initial(1,nages);         //initial Z at age
  number spr_initial;                //initial spawners per recruit
  number spr_F0;                     //Spawning biomass per recruit at F=0
  number bpr_F0;                     //Biomass per recruit at F=0
  number rec_mean;  				 //mean recruitment across years for computing spr-related biomass benchmarks	
  number iter_inc_spr;               //increments used to compute msy, equals max_F_spr_msy/(n_iter_spr-1)

////-------SDNR output-----------------------------------------------------------------------------
  number sdnr_I_survey1;
  number sdnr_ac_survey1;
  number sdnr_ac_fleet1;
  number sdnr_ac_fleet2;

  
////-------Objective function components-----------------------------------------------------------------------------
  number w_L;
  number w_D;
  number w_I_survey1;
  number w_ac_fleet1; 
  number w_ac_fleet2; 
  number w_lc_fleet3;
  number w_ac_survey1;

  number w_Nage_init;  
  number w_rec;
  number w_rec_early;
  number w_rec_end;
  number w_Ftune;

  number f_fleet1_L; 
  number f_fleet2_L; 
  number f_fleet3_D;
  number f_survey1_cpue;
  number f_fleet1_agec;
  number f_fleet2_agec;
  number f_fleet3_lenc;
  number f_survey1_agec;  

  //Close kin mark recapture
  //number nll_ckmr;

  //Penalties and constraints. Not all are used.
  number f_Nage_init;                //weight on log devs to estimate initial abundance (exfleet2uding first age)
  number f_rec_dev;                //weight on recruitment deviations to fit S-R curve
  number f_rec_dev_early;          //extra weight on deviations in first recruitment stanza
  number f_rec_dev_end;            //extra weight on deviations in ending recruitment stanza
  number f_Ftune;                  //penalty for tuning F in Ftune yr.  Not applied in final optimization phase.
  number f_priors;                 //prior information on parameters
  
  //init_number xdum;
  objective_function_value fval;
  number fval_data;
  number grad_max;
  
//--Dummy variables ----
  number denom;                   //denominator used in some calculations
  number numer;                   //numerator used in some calculations

  //---------- Projection quantities--------------------------------------------------------------
  number F_reg_proj;						   //value used to define the projections				
  vector F_proj(styr_proj,endyr_proj);         //F by yr for projections (=F_reg_proj after regulations start, current F till then)  
  vector L_knum_proj(styr_proj,endyr_proj);    //total landings in 1000 fish for projections  
  vector L_mt_proj(styr_proj,endyr_proj);      //total landings in weight (1000 lb) for projections
  vector D_knum_proj(styr_proj,endyr_proj);    //total dead discards in 1000 fish for projections  
  vector D_mt_proj(styr_proj,endyr_proj);      //total dead discards in weight (1000 lb) for projections
  
  vector B_proj(styr_proj,endyr_proj);         //Biomass for projections  
  vector SSB_proj(styr_proj,endyr_proj);       //SSB for projections  
  vector R_proj(styr_proj,endyr_proj);     	   //recruits for projections
  vector FL_age_proj(1,nages);      		      //F (landings) by age for projections   
  vector FD_age_proj(1,nages);                //F (dead discards) by age for projections
  
  matrix N_proj(styr_proj,endyr_proj,1,nages);           //Population numbers by year and age at start of yr
  matrix N_spawn_proj(styr_proj,endyr_proj,1,nages);     //Population numbers by year and age at peaking spawning: used for SSB in projections 
  matrix Z_proj(styr_proj,endyr_proj,1,nages);           //Z by year and age for projections 
  matrix L_age_proj(styr_proj,endyr_proj,1,nages);       //Projected landings at age in numbers 
  matrix D_age_proj(styr_proj,endyr_proj,1,nages);       //Projected dead discards at age in numbers
     
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
GLOBALS_SECTION
  #include "admodel.h"          // Include AD class definitions
  #include "admb2r.cpp"         // Include R-compatible output functions (needs preceding)
  #include <time.h>
	time_t start,finish;
	long hour,minute,second;	
	double elapsed_time;
	
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
RUNTIME_SECTION
 maximum_function_evaluations 1000, 2000,5000, 15000;
 convergence_criteria 1e-2, 1e-2,1e-3, 1e-4;
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PRELIMINARY_CALCS_SECTION
// Set values of fixed parameters or set initial guess of estimated parameters
  
  Linf=set_Linf(1);
  K=set_K(1);
  t0=set_t0(1);
  len_cv_val=set_len_cv(1);
  
  M_logit = set_M_logit(1);
  
  // M_scale = (2.0*exp(M_logit)/(1.0+exp(M_logit));
  // M=M_scale*set_M; 
  smsy2msst=0.75;
  
  log_R0=set_log_R0(1);
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);

  log_dm_fleet1_ac=set_log_dm_fleet1_ac(1);
  log_dm_fleet2_ac=set_log_dm_fleet2_ac(1);
  log_dm_fleet3_lc=set_log_dm_fleet3_lc(1);
  log_dm_survey1_ac=set_log_dm_survey1_ac(1);
   
  log_q_survey1=set_log_q_survey1(1);
   
  w_L=set_w_L;
  w_D=set_w_D;
  w_I_survey1=set_w_I_survey1;
  
  w_ac_fleet1=set_w_ac_fleet1;
  w_ac_fleet2=set_w_ac_fleet2;
  w_lc_fleet3=set_w_lc_fleet3;
  w_ac_survey1=set_w_ac_survey1;

  w_Nage_init=set_w_Nage_init;
  w_rec=set_w_rec;
  w_rec_early=set_w_rec_early;
  w_rec_end=set_w_rec_end;
  w_Ftune=set_w_Ftune;

  F_init=set_F_init(1);
  log_avg_F_fleet1=set_log_avg_F_fleet1(1);
  log_avg_F_fleet2=set_log_avg_F_fleet2(1);
  log_avg_F_fleet3=set_log_avg_F_fleet3(1);
  
  log_F_dev_fleet1=set_log_F_dev_fleet1_vals;
  log_F_dev_fleet2=set_log_F_dev_fleet2_vals;
  log_F_dev_fleet3=set_log_F_dev_fleet3_vals;

  selpar_A50_survey1=set_selpar_A50_survey1(1);
  selpar_slope_survey1=set_selpar_slope_survey1(1);
 
  selpar_A50_fleet1_B1=set_selpar_A50_fleet1_B1(1);
  selpar_slope_fleet1_B1=set_selpar_slope_fleet1_B1(1);
  selpar_A50_fleet2_B1=set_selpar_A50_fleet2_B1(1);
  selpar_slope_fleet2_B1=set_selpar_slope_fleet2_B1(1);
  selpar_A50_asc_fleet3_B1=set_selpar_A50_asc_fleet3_B1(1);
  selpar_slope_asc_fleet3_B1=set_selpar_slope_asc_fleet3_B1(1);
  selpar_A50_dsc_fleet3_B1=set_selpar_A50_dsc_fleet3_B1(1);
  selpar_slope_dsc_fleet3_B1=set_selpar_slope_dsc_fleet3_B1(1);
 
 pi=3.14159265;
 sqrt2pi=sqrt(2.*3.14159265);
 g2mt=0.000001;         //conversion of grams to metric tons
 g2kg=0.001;            //conversion of grams to kg 
 kg2mt=0.001;            //conversion of kilograms to metric tons  
 mt2kg=1000.0;            //conversion of metric tons to kilograms   
 mt2klb=2.20462;        //conversion of metric tons to 1000 lb 
 mt2lb=mt2klb*1000.0;   //conversion of metric tons to lb
 g2klb=g2mt*mt2klb;     //conversion of grams to 1000 lb 
 dzero=0.00001;         
 huge_number=1.0e+10;   
 
 SSB_msy_out=0.0;

 iter_inc_msy=max_F_spr_msy/(n_iter_msy-1);
 iter_inc_spr=max_F_spr_msy/(n_iter_spr-1); 

 prop_f=prop_f_obs;
 maturity_f=maturity_f_obs;
 fecpar_a=set_fecpar_a;
 fecpar_b=set_fecpar_b;
 fecpar_c=set_fecpar_c;
 fecpar_thresh=set_fecpar_thresh;
 fecpar_min=set_fecpar_min;
 //fecpar_batches=set_fecpar_batches;
 
//Fill in sample sizes of comps, possibly sampled in nonconsec yrs 
//Used primarily for output in R object   

  nsamp_survey1_agec_allyr=missing;
  nsamp_fleet1_agec_allyr=missing;
  nsamp_fleet2_agec_allyr=missing;
  nsamp_fleet3_lenc_allyr=missing;
      						
  for (iyear=1; iyear<=nyr_survey1_agec; iyear++)
         {if (nsamp_survey1_agec(iyear)>=minSS_agec)
           {nsamp_survey1_agec_allyr(yrs_survey1_agec(iyear))=nsamp_survey1_agec(iyear);
            nfish_survey1_agec_allyr(yrs_survey1_agec(iyear))=nfish_survey1_agec(iyear);}} 

  for (iyear=1; iyear<=nyr_fleet1_agec; iyear++)
         {if (nsamp_fleet1_agec(iyear)>=minSS_agec)
           {nsamp_fleet1_agec_allyr(yrs_fleet1_agec(iyear))=nsamp_fleet1_agec(iyear);
            nfish_fleet1_agec_allyr(yrs_fleet1_agec(iyear))=nfish_fleet1_agec(iyear);}}

 for (iyear=1; iyear<=nyr_fleet2_agec; iyear++)
         {if (nsamp_fleet2_agec(iyear)>=minSS_agec)
           {nsamp_fleet2_agec_allyr(yrs_fleet2_agec(iyear))=nsamp_fleet2_agec(iyear);
            nfish_fleet2_agec_allyr(yrs_fleet2_agec(iyear))=nfish_fleet2_agec(iyear);}}

 for (iyear=1; iyear<=nyr_fleet3_lenc; iyear++)
         {if (nsamp_fleet3_lenc(iyear)>=minSS_lenc)
           {nsamp_fleet3_lenc_allyr(yrs_fleet3_lenc(iyear))=nsamp_fleet3_lenc(iyear);
            nfish_fleet3_lenc_allyr(yrs_fleet3_lenc(iyear))=nfish_fleet3_lenc(iyear);}}           
			
  obs_survey1_cpue_allyr=missing; 
  pred_survey1_cpue_allyr=missing;
  survey1_cpue_cv_allyr=missing;	  
	  
//fill in Fs for msy and per-recruit analyses
  F_msy(1)=0.0;  
  for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
  for (ff=2;ff<=n_iter_spr;ff++) {F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}


//fill in F's, Catch matrices, and log rec dev with zero's
  F_fleet1.initialize();
  L_fleet1_num.initialize();
  F_fleet2.initialize();
  L_fleet2_num.initialize();
  F_fleet3.initialize();
  D_fleet3_num.initialize();
   
  F_fleet1_out.initialize();
  F_fleet2_out.initialize();
  F_fleet3_out.initialize();
  sel_survey1.initialize();
  sel_fleet1.initialize();
  sel_fleet2.initialize();
  sel_fleet3.initialize();
  
  log_F1_dev_output.initialize();
  log_F2_dev_output.initialize();
  log_F3_dev_output.initialize();
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
 get_natural_mortality();
 //cout << "got natural mortality" << endl;
 get_length_weight_at_age(); 
 //cout << "got length, weight, fecundity transitions" <<endl;
 get_reprod();
 //cout << "got repro stuff" << endl;
 get_len_at_age_dist();
 //cout << "got length at age dist" << end
 get_spr_F0();
 //cout << "got F0 spr" << endl;
 get_selectivity(); 
 //cout << "got selectivity" << endl;
 get_mortality(); 
 //cout << "got mortalities" << endl;
 get_bias_corr(); 
 //cout<< "got recruitment bias correction" << endl;
 get_numbers_at_age(); 
 //cout << "got numbers at age" << endl;
 get_landings_numbers();
 //cout << "got landings in numbers" << endl;
 get_landings_wgt();
 // cout << "got landings in wgt" << endl;
 get_indices();
 //cout << "got indices" << endl;
 get_age_comps();
 //cout<< "got age comps"<< endl;
 get_length_comps();
 //cout<< "got length comps"<< endl;
 evaluate_objective_function();
 //cout << "objective function calculations complete" << endl;
 
FUNCTION get_natural_mortality
  //scalar on natural mortality vector may be estimated
  M_scale = (2.0*exp(M_logit))/(1.0+exp(M_logit));
  M=M_scale*set_M;  
 
FUNCTION get_length_weight_at_age
	//population total length in mm
    //compute mean length (mm TL) and weight (whole) at age
    meanlen_TL=Linf*(1.0-exp(-K*(agebins-t0)));     
    wgt_kg=wgtpar_a*pow(meanlen_TL,wgtpar_b);             //whole wgt in kg 
    wgt_g=wgt_kg/g2kg;                                    //convert wgt in kg to weight in g    
    wgt_mt=wgt_g*g2mt;                                    //convert weight in g to weight in mt
    wgt_klb=mt2klb*wgt_mt;                                //1000 lb of whole wgt
    wgt_lb=mt2lb*wgt_mt;                                  //lb of whole wgt
    fecundity=fecpar_a+fecpar_b*pow(meanlen_TL,fecpar_c); //relationship provided by Wyanski et al.  
    for (iage=1; iage<=nages; iage++)
    { if(meanlen_TL(iage)<fecpar_thresh){fecundity(iage)=fecpar_min;} 
          else {break;}     
    }
 
FUNCTION get_reprod 
   //reprod is product of stuff going into reproductive capacity calcs
   //reprod=elem_prod((elem_prod(prop_f,maturity_f)),wgt_mt);
   reprod=elem_prod((elem_prod(prop_f,maturity_f)),fecundity); //Since this is in eggs, I think it can be used in RO calcs

FUNCTION get_len_at_age_dist
// SEDAR 73 uses population lenprob for gen rec discards
// code below borrowed directly for population lenprob in SEDAR 73

  //compute matrix of length at age, based on the normal distribution
    //population
  for (iage=1;iage<=nages;iage++)
   {len_cv(iage)=len_cv_val;
    len_sd(iage)=meanlen_TL(iage)*len_cv(iage);
    zscore_lzero=(0.0-meanlen_TL(iage))/len_sd(iage); 
    cprob_lzero=cumd_norm(zscore_lzero);

  //first length bin
  //population
    zscore_len=((lenbins(1)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage);
    cprob_lenvec(1)=cumd_norm(zscore_len);          //includes any probability mass below zero
    lenprob(iage,1)=cprob_lenvec(1)-cprob_lzero;    //removes any probability mass below zero

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

  //assign length probability to fleet
  lenprob_fleet3=lenprob;

FUNCTION get_spr_F0
  N_spr_F0(1)=1.0*exp(-1.0*M(1)*spawn_time_frac); //at peak spawning time
  N_bpr_F0(1)=1.0;      //at start of year
  for (iage=2; iage<=nages; iage++)
  {
    N_spr_F0(iage)=N_spr_F0(iage-1)*exp(-1.0*(M(iage-1)*(1.0-spawn_time_frac) + M(iage)*spawn_time_frac)); 
    N_bpr_F0(iage)=N_bpr_F0(iage-1)*exp(-1.0*(M(iage-1)));    
  }
  N_spr_F0(nages)=N_spr_F0(nages)/(1.0-exp(-1.0*M(nages))); //plus group (sum of geometric series)
  N_bpr_F0(nages)=N_bpr_F0(nages)/(1.0-exp(-1.0*M(nages)));
  
  spr_F0=sum(elem_prod(N_spr_F0,reprod)); 
  bpr_F0=sum(elem_prod(N_bpr_F0,wgt_mt));    

FUNCTION get_selectivity    
  selvec_survey1=logistic(agebins, selpar_A50_survey1, selpar_slope_survey1);	
  selvec_fleet1_B1=logistic(agebins, selpar_A50_fleet1_B1, selpar_slope_fleet1_B1);  
  selvec_fleet2_B1=logistic(agebins, selpar_A50_fleet2_B1, selpar_slope_fleet2_B1);  
  selvec_fleet3_B1=logistic_double(agebins, selpar_A50_asc_fleet3_B1, selpar_slope_asc_fleet3_B1, selpar_A50_dsc_fleet3_B1, selpar_slope_dsc_fleet3_B1);
 
   for (iyear=styr; iyear<=endyr_selex_phase1; iyear++)
   {     
 	   sel_survey1(iyear)(1,nages_agec)=selvec_survey1(1,nages_agec);
	   sel_fleet1(iyear)(1,nages_agec)=selvec_fleet1_B1(1,nages_agec);
	   sel_fleet2(iyear)(1,nages_agec)=selvec_fleet2_B1(1,nages_agec);
	   sel_fleet3(iyear)(1,nages_agec)=selvec_fleet3_B1(1,nages_agec);
   }
 for (iyear=(endyr_selex_phase1+1); iyear<=endyr; iyear++)  //architecture here in case second block is used
   {
	   sel_survey1(iyear)(1,nages_agec)=selvec_survey1(1,nages_agec);
	   sel_fleet1(iyear)(1,nages_agec)=selvec_fleet1_B1(1,nages_agec);
	   sel_fleet2(iyear)(1,nages_agec)=selvec_fleet2_B1(1,nages_agec);
	   sel_fleet3(iyear)(1,nages_agec)=selvec_fleet3_B1(1,nages_agec);
   }	   
  
FUNCTION get_mortality
  Fsum.initialize();
  Fapex.initialize();
  F.initialize();
  
  //log_F_dev_init_fleet1=sum(log_F_dev_fleet1(styr_fleet1_L,(styr_fleet1_L+2)))/3.0;  //initialization F is avg from first 3 yrs of observed landings                  
  log_F_dev_init_fleet1=log_F_dev_fleet1(styr_fleet1_L); //For sim study, only using first year to match operating model
  log_F_dev_init_fleet2=log_F_dev_fleet2(styr_fleet2_L); //For sim study, only using first year to match operating model
  log_F_dev_init_fleet3=log_F_dev_fleet3(styr_fleet3_D); //^
  
  for (iyear=styr; iyear<=endyr; iyear++) 
  {
    if(iyear>=styr_fleet1_L & iyear<=endyr_fleet1_L)
    {  F_fleet1_out(iyear)=exp(log_avg_F_fleet1+log_F_dev_fleet1(iyear));     
       F_fleet1(iyear)=sel_fleet1(iyear)*F_fleet1_out(iyear);
       Fsum(iyear)+=F_fleet1_out(iyear);
    }

    if(iyear>=styr_fleet2_L & iyear<=endyr_fleet2_L)
    {  F_fleet2_out(iyear)=exp(log_avg_F_fleet2+log_F_dev_fleet2(iyear));     
       F_fleet2(iyear)=sel_fleet2(iyear)*F_fleet2_out(iyear);
       Fsum(iyear)+=F_fleet2_out(iyear);
    }
    
    if(iyear>=styr_fleet3_D & iyear<=endyr_fleet3_D)
    {  F_fleet3_out(iyear)=exp(log_avg_F_fleet3+log_F_dev_fleet3(iyear));     
       F_fleet3(iyear)=sel_fleet3(iyear)*F_fleet3_out(iyear);
       Fsum(iyear)+=F_fleet3_out(iyear);
    }
	
    //Total F at age
    F(iyear)=F_fleet1(iyear); //first in additive series (NO +=)
    F(iyear)+=F_fleet2(iyear);
    F(iyear)+=F_fleet3(iyear);
	
    Fapex(iyear)=max(F(iyear));
    Z(iyear)=M+F(iyear);  
  }  //end iyear 
 

FUNCTION get_bias_corr
  var_rec_dev=norm2(log_rec_dev(styr_rec_dev,endyr_rec_dev)-
              sum(log_rec_dev(styr_rec_dev,endyr_rec_dev))/nyrs_rec)
              /(nyrs_rec-1.0);                           
  //if (set_BiasCor <= 0.0) {BiasCor=exp(var_rec_dev/2.0);}   //bias correction based on empirical residuals
  rec_sigma_sq=square(rec_sigma);
  if (set_BiasCor <= 0.0) {BiasCor=exp(rec_sigma_sq/2.0);}   //bias correction based on Rsigma               
  else {BiasCor=set_BiasCor;}

  
FUNCTION get_numbers_at_age
//Initialization
 
  R0=exp(log_R0);
  S0=spr_F0*R0;
  
  R_virgin=SR_eq_func(R0, steep, spr_F0, spr_F0, BiasCor, SR_switch);
  B0=bpr_F0*R_virgin;   
  
  F_initial=F_init*F(styr);  //initialization F is a scalar multiple (F_init) of F in year one
  Z_initial=M+F_initial;
  
//Initial equilibrium age structure
  N_spr_initial(1)=1.0*exp(-1.0*Z_initial(1)*spawn_time_frac); //at peak spawning time;
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   exp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-exp(-1.0*Z_initial(nages))); //plus group
  spr_initial=sum(elem_prod(N_spr_initial,reprod));

  if (styr==styr_rec_dev) {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, 1.0, SR_switch);} //without bias correction (deviation added later)
  else {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, BiasCor, SR_switch);} //with bias correction
  if(R1<10.0) {R1=10.0;} //Avoid unrealistically low popn sizes during search algorithm
  
//Compute equilibrium age structure for first year
  N_initial_eq(1)=R1;
  for (iage=2; iage<=nages; iage++)
  {
    N_initial_eq(iage)=N_initial_eq(iage-1)*
        exp(-1.0*(Z_initial(iage-1)));    
  }
  //plus group calculation
  N_initial_eq(nages)=N_initial_eq(nages)/(1.0-exp(-1.0*Z_initial(nages))); //plus group
  
//Add deviations, if estimated, to initial equilibrium N
  N(styr)(2,nages)=elem_prod(N_initial_eq(2,nages),exp(log_Nage_dev));
   
  if (styr==styr_rec_dev) {N(styr,1)=N_initial_eq(1)*exp(log_rec_dev(styr_rec_dev));}
  else {N(styr,1)=N_initial_eq(1);}
  
  N_mdyr(styr)(1,nages)=elem_prod(N(styr)(1,nages),(exp(-1.*(Z_initial(1,nages))*0.5))); //mid year 
  N_spawn(styr)(1,nages)=elem_prod(N(styr)(1,nages),(exp(-1.*(Z_initial(1,nages))*spawn_time_frac))); //peak spawning time 

  SSB(styr)=sum(elem_prod(N_spawn(styr),reprod));
  MatFemB(styr)=sum(elem_prod(N_spawn(styr),reprod));
    
//Rest of years 
  for (iyear=styr; iyear<endyr; iyear++)
  {
    if(iyear<(styr_rec_dev-1)||iyear>(endyr_rec_dev-1)) //recruitment follows S-R curve (with bias correction) exactly
    {
        N(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch);
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(exp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*exp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(exp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(exp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));   
        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));		
    }
    else   //recruitment follows S-R curve with lognormal deviation
    {
        N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch)*exp(log_rec_dev(iyear+1));
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(exp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*exp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(exp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(exp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));
		MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));		

    }
  }

  //last year (projection) has no recruitment variability
  N(endyr+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(endyr),SR_switch);
  N(endyr+1)(2,nages)=++elem_prod(N(endyr)(1,nages-1),(exp(-1.*Z(endyr)(1,nages-1))));
  N(endyr+1,nages)+=N(endyr,nages)*exp(-1.*Z(endyr,nages)); //plus group

//Time series of interest
  rec=column(N,1);
  SdS0=SSB/S0;

//---Predicted Landings------------------------
FUNCTION get_landings_numbers //Baranov catch eqn
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      L_fleet1_num(iyear,iage)=N(iyear,iage)*F_fleet1(iyear,iage)*
        (1.-exp(-1.*Z(iyear,iage)))/Z(iyear,iage);
	    L_fleet2_num(iyear,iage)=N(iyear,iage)*F_fleet2(iyear,iage)*
        (1.-exp(-1.*Z(iyear,iage)))/Z(iyear,iage);	
      D_fleet3_num(iyear,iage)=N(iyear,iage)*F_fleet3(iyear,iage)* // separate function in SEDAR 76, but does the same thing
        (1.-exp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }          
    pred_fleet1_L_knum(iyear)=sum(L_fleet1_num(iyear))/1000.0;  
	  pred_fleet2_L_knum(iyear)=sum(L_fleet2_num(iyear))/1000.0;  
	  pred_fleet3_D_knum(iyear)=sum(D_fleet3_num(iyear))/1000.0;
  }

 
FUNCTION get_landings_wgt
  for (iyear=styr; iyear<=endyr; iyear++)
  {    
    L_fleet1_mt(iyear)=elem_prod(L_fleet1_num(iyear),wgt_mt);     //in mt
    pred_fleet1_L_mt(iyear)=sum(L_fleet1_mt(iyear)); 
	
	  L_fleet2_mt(iyear)=elem_prod(L_fleet2_num(iyear),wgt_mt);     //in mt
    pred_fleet2_L_mt(iyear)=sum(L_fleet2_mt(iyear)); 
    
    D_fleet3_mt(iyear)=elem_prod(D_fleet3_num(iyear),wgt_mt);   //in mt
    pred_fleet3_D_mt(iyear)=sum(D_fleet3_mt(iyear));
  }

     
FUNCTION get_indices
//---Predicted CPUEs------------------------

 //Survey 1: cpue
  for (iyear=1; iyear<=nyr_survey1_cpue; iyear++)
  {   //index in number units, start of yr
      N_survey1(iyear)=elem_prod(N(yrs_survey1_cpue(iyear)),sel_survey1(yrs_survey1_cpue(iyear)));
	  pred_survey1_cpue(iyear)=exp(log_q_survey1)*sum(N_survey1(iyear));
      
	  obs_survey1_cpue_allyr(yrs_survey1_cpue(iyear))=obs_survey1_cpue(iyear);
	  pred_survey1_cpue_allyr(yrs_survey1_cpue(iyear))=pred_survey1_cpue(iyear);
	  survey1_cpue_cv_allyr(yrs_survey1_cpue(iyear))=survey1_cpue_cv(iyear);
  }



FUNCTION get_age_comps
   
 //survey1
 for (iyear=1;iyear<=nyr_survey1_agec;iyear++)
  {
    ErrorFree_survey1_agec(iyear)=N_survey1(iyear)/sum(N_survey1(iyear));
    pred_survey1_agec_allages(iyear)=age_error*ErrorFree_survey1_agec(iyear); 
    for (iage=1; iage<=nages_agec; iage++) {pred_survey1_agec(iyear,iage)=pred_survey1_agec_allages(iyear,iage);} 
    for (iage=(nages_agec+1); iage<=nages; iage++) {pred_survey1_agec(iyear,nages_agec)+=pred_survey1_agec_allages(iyear,iage);} //plus group                        
  }

  //fleet1
  for (iyear=1;iyear<=nyr_fleet1_agec;iyear++) 
  {
    ErrorFree_fleet1_agec(iyear)=L_fleet1_num(yrs_fleet1_agec(iyear))/sum(L_fleet1_num(yrs_fleet1_agec(iyear)));  
    pred_fleet1_agec_allages(iyear)=age_error*(ErrorFree_fleet1_agec(iyear)/sum(ErrorFree_fleet1_agec(iyear)));   
    for (iage=1; iage<=nages_agec; iage++) {pred_fleet1_agec(iyear,iage)=pred_fleet1_agec_allages(iyear,iage);} 
    for (iage=(nages_agec+1); iage<=nages; iage++) {pred_fleet1_agec(iyear,nages_agec)+=pred_fleet1_agec_allages(iyear,iage);} //plus group                             
  }
 
   //fleet2
  for (iyear=1;iyear<=nyr_fleet2_agec;iyear++) 
  {
    ErrorFree_fleet2_agec(iyear)=L_fleet2_num(yrs_fleet2_agec(iyear))/sum(L_fleet2_num(yrs_fleet2_agec(iyear)));  
    pred_fleet2_agec_allages(iyear)=age_error*(ErrorFree_fleet2_agec(iyear)/sum(ErrorFree_fleet2_agec(iyear)));   
    for (iage=1; iage<=nages_agec; iage++) {pred_fleet2_agec(iyear,iage)=pred_fleet2_agec_allages(iyear,iage);} 
    for (iage=(nages_agec+1); iage<=nages; iage++) {pred_fleet2_agec(iyear,nages_agec)+=pred_fleet2_agec_allages(iyear,iage);} //plus group                             
  }

FUNCTION get_length_comps
  //fleet3 (equivalent to gen rec in SEDAR 73)
  for (iyear=1;iyear<=nyr_fleet3_lenc;iyear++) 
  {pred_fleet3_lenc(iyear)=(D_fleet3_num(yrs_fleet3_lenc(iyear))*lenprob)/sum(D_fleet3_num(yrs_fleet3_lenc(iyear)));}

  
 //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_weighted_current 
  F_temp_sum=0.0;
  //Landings
  F_temp_sum+=exp((selpar_n_yrs_wgted*log_avg_F_fleet1+
        sum(log_F_dev_fleet1((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_temp_sum+=exp((selpar_n_yrs_wgted*log_avg_F_fleet2+
        sum(log_F_dev_fleet2((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  //Discards
  F_temp_sum+=exp((selpar_n_yrs_wgted*log_avg_F_fleet3+
        sum(log_F_dev_fleet3((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
	//Landings
  F_fleet1_prop=exp((selpar_n_yrs_wgted*log_avg_F_fleet1+
        sum(log_F_dev_fleet1((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_fleet2_prop=exp((selpar_n_yrs_wgted*log_avg_F_fleet2+
        sum(log_F_dev_fleet2((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  //Discards
  F_fleet3_prop=exp((selpar_n_yrs_wgted*log_avg_F_fleet3+
        sum(log_F_dev_fleet3((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
 
  
  
  F_end_apex=max(F_fleet1_prop*selvec_fleet1_B1  + F_fleet2_prop*selvec_fleet2_B1 + F_fleet3_prop*selvec_fleet3_B1);
  sel_wgted_L=(F_fleet1_prop*selvec_fleet1_B1  + F_fleet2_prop*selvec_fleet2_B1) / F_end_apex;
  sel_wgted_D=(F_fleet3_prop*selvec_fleet3_B1) / F_end_apex;
  sel_wgted_tot=sel_wgted_L+sel_wgted_D; //total weighted selex
 
  wgt_wgted_L_mt=wgt_mt; //weighted among fleets, only if multiple fleets are modeled 
  wgt_wgted_D_mt=wgt_mt; //
  
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_msy
  
  //compute values as functions of F
  for(ff=1; ff<=n_iter_msy; ff++)
  {
    //uses fleet-weighted F's; although only one fleet in this sim study
    Z_age_msy=0.0;
    F_L_age_msy=0.0;
    F_D_age_msy=0.0; 
          
    F_L_age_msy=F_msy(ff)*sel_wgted_L;
    F_D_age_msy=F_msy(ff)*sel_wgted_D;
    Z_age_msy=M+F_L_age_msy+F_D_age_msy;         
    
    N_age_msy(1)=1.0;
    for (iage=2; iage<=nages; iage++)
      {N_age_msy(iage)=N_age_msy(iage-1)*exp(-1.*Z_age_msy(iage-1));}
    N_age_msy(nages)=N_age_msy(nages)/(1.0-exp(-1.*Z_age_msy(nages)));
    N_age_msy_spawn(1,(nages-1))=elem_prod(N_age_msy(1,(nages-1)),
                                   exp((-1.*Z_age_msy(1,(nages-1)))*spawn_time_frac));                 
    N_age_msy_spawn(nages)=(N_age_msy_spawn(nages-1)*(exp(-1.*(Z_age_msy(nages-1)*(1.0-spawn_time_frac) + 
                            Z_age_msy(nages)*spawn_time_frac) )))/(1.0-exp(-1.*Z_age_msy(nages)));
                     
    spr_msy(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
	        
    R_eq(ff)=SR_eq_func(R0, steep, spr_msy(1), spr_msy(ff), BiasCor, SR_switch);
    
    if (R_eq(ff)<dzero) {R_eq(ff)=dzero;}    
    N_age_msy*=R_eq(ff);
    N_age_msy_spawn*=R_eq(ff);
    
    for (iage=1; iage<=nages; iage++)
    {
      L_age_msy(iage)=N_age_msy(iage)*(F_L_age_msy(iage)/Z_age_msy(iage))*
                      (1.-exp(-1.*Z_age_msy(iage)));
      D_age_msy(iage)=N_age_msy(iage)*(F_D_age_msy(iage)/Z_age_msy(iage))*
                      (1.-exp(-1.*Z_age_msy(iage)));
    }
    
    SSB_eq(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
	B_eq(ff)=sum(elem_prod(N_age_msy,wgt_mt));
    L_eq_mt(ff)=sum(elem_prod(L_age_msy,wgt_wgted_L_mt)); //in whole weight
    L_eq_knum(ff)=sum(L_age_msy)/1000.0;  
    D_eq_mt(ff)=sum(elem_prod(D_age_msy,wgt_wgted_D_mt)); 
    D_eq_knum(ff)=sum(D_age_msy)/1000.0;
  }  
  
  msy_mt_out=max(L_eq_mt); //msy in whole weight 
  
  for(ff=1; ff<=n_iter_msy; ff++)
  {
   if(L_eq_mt(ff) == msy_mt_out) 
      {    
        SSB_msy_out=SSB_eq(ff);
        B_msy_out=B_eq(ff);
        R_msy_out=R_eq(ff);
        msy_knum_out=L_eq_knum(ff);
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
     {N_age_spr(iage)=N_age_spr(iage-1)*exp(-1.*Z(iyear,iage-1));}
    N_age_spr(nages)=N_age_spr(nages)/(1.0-exp(-1.*Z(iyear,nages)));    
    N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                exp(-1.*Z(iyear)(1,(nages-1))*spawn_time_frac));
    N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (exp(-1.*(Z(iyear)(nages-1)*(1.0-spawn_time_frac) + Z(iyear)(nages)*spawn_time_frac) )))
                          /(1.0-exp(-1.*Z(iyear)(nages)));           
    spr_static(iyear)=sum(elem_prod(N_age_spr_spawn,reprod))/spr_F0;
  }
 
  //compute SSB/R and YPR as functions of F
  for(ff=1; ff<=n_iter_spr; ff++)
  {
    //uses fishery-weighted F's, same as in MSY calculations
    Z_age_spr=0.0;
    F_L_age_spr=0.0;
    F_D_age_spr=0.0;
    
	F_L_age_spr=F_spr(ff)*sel_wgted_L;
    F_D_age_spr=F_spr(ff)*sel_wgted_D;
    Z_age_spr=M+F_L_age_spr+F_D_age_spr;    

    N_age_spr(1)=1.0;
    for (iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*exp(-1.*Z_age_spr(iage-1));}
    N_age_spr(nages)=N_age_spr(nages)/(1-exp(-1.*Z_age_spr(nages)));
    N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   exp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
    N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (exp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-exp(-1.*Z_age_spr(nages)));
    spr_spr(ff)=sum(elem_prod(N_age_spr_spawn,reprod));
	L_spr(ff)=0.0;
	D_spr(ff)=0.0;
    for (iage=1; iage<=nages; iage++)
    {
      L_age_spr(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-exp(-1.*Z_age_spr(iage)));
      L_spr(ff)+=L_age_spr(iage)*wgt_wgted_L_mt(iage)*1000.0; //in kg whole wgt
      D_age_spr(iage)=N_age_spr(iage)*(F_D_age_spr(iage)/Z_age_spr(iage))*
                      (1.-exp(-1.*Z_age_spr(iage)));
      D_spr(ff)+=D_age_spr(iage)*wgt_wgted_D_mt(iage)*1000.0; 
    }   
  }
  spr_ratio=spr_spr/spr_F0;
  F30_dum=min(fabs(spr_ratio-0.3));
  F40_dum=min(fabs(spr_ratio-0.4));
  F50_dum=min(fabs(spr_ratio-0.5));
  Fproxy_dum=min(fabs(spr_ratio-spr_proxy));
  for(ff=1; ff<=n_iter_spr; ff++)
  {   
      if (fabs(spr_ratio(ff)-0.3)==F30_dum) {F30_out=F_spr(ff);}	  
	  if (fabs(spr_ratio(ff)-0.4)==F40_dum) {F40_out=F_spr(ff);}
	  if (fabs(spr_ratio(ff)-0.5)==F50_dum) {F50_out=F_spr(ff);}	
	  if (fabs(spr_ratio(ff)-spr_proxy)==Fproxy_dum) {Fproxy_out=F_spr(ff);}
  }
    
	rec=column(N,1);
    rec_mean=R0; //sum(rec(styr_rec_dev, endyr_rec_dev))/nyrs_rec; //avg recruitment across all years	
    // switch(SPR_rec_switch){
       // case 1: //Observed mean over years specified
          // rec_mean=sum(rec(styr_rec_spr, endyr_rec_spr))/nyrs_rec_spr;
          // break;
       // case 2: //Expected recruitment is R0 with bias correction (if that feature is turned on)
          // rec_mean=BiasCor*R0;
          // break;
       // default: // no such switch available
          // cout << "Error in input: SPR_rec_switch must be set to 1 or 2." << endl;
          // cout << "Presently it is set to " << SPR_rec_switch <<"."<< endl;
          // exit(0);          
	// }	  
   

  R_Fproxy_out=rec_mean;
  F_L_age_spr=Fproxy_out*sel_wgted_L;
  F_D_age_spr=Fproxy_out*sel_wgted_D;
  Z_age_spr=M+F_L_age_spr+F_D_age_spr;

  N_age_spr(1)=R_Fproxy_out;
  for (iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*exp(-1.*Z_age_spr(iage-1));}
  N_age_spr(nages)=N_age_spr(nages)/(1-exp(-1.*Z_age_spr(nages)));
  N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   exp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
  N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (exp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-exp(-1.*Z_age_spr(nages)));

  for (iage=1; iage<=nages; iage++)
    {
      L_age_Fproxy(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-exp(-1.*Z_age_spr(iage)));
      D_age_Fproxy(iage)=N_age_spr(iage)*(F_D_age_spr(iage)/Z_age_spr(iage))*
                      (1.-exp(-1.0*Z_age_spr(iage)));                      
    }
		
  SSB_Fproxy_out=sum(elem_prod(N_age_spr_spawn,reprod));
  B_Fproxy_out=sum(elem_prod(N_age_spr,wgt_mt));
  L_Fproxy_mt_out=sum(elem_prod(L_age_Fproxy,wgt_wgted_L_mt)); //in whole weight
  L_Fproxy_knum_out=sum(L_age_Fproxy)/1000.0;  
  D_Fproxy_mt_out=sum(elem_prod(D_age_Fproxy,wgt_wgted_D_mt)); //in whole weight   
  D_Fproxy_knum_out=sum(D_age_Fproxy)/1000.0;    



  


//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_miscellaneous_stuff

//switch here if var_rec_dev <=dzero 
  if(var_rec_dev>0.0)
   {sigma_rec_dev=sqrt(var_rec_dev);} //sample SD of predicted residuals (may not equal rec_sigma)  
   else{sigma_rec_dev=0.0;}

  //compute total landings- and discards-at-age in 1000 fish and mt whole weight
  L_total_num.initialize();
  L_total_mt.initialize();
  L_total_knum_yr.initialize();
  L_total_mt_yr.initialize(); 
  
  D_total_num.initialize();
  D_total_mt.initialize();
  D_total_knum_yr.initialize();
  D_total_mt_yr.initialize();
  
  
  for(iyear=styr; iyear<=endyr; iyear++)
  {
        L_total_mt_yr(iyear)=pred_fleet1_L_mt(iyear);
        L_total_knum_yr(iyear)=pred_fleet1_L_knum(iyear);
        
        D_total_mt_yr(iyear)=pred_fleet3_D_mt(iyear);
        D_total_knum_yr(iyear)=pred_fleet3_D_knum(iyear);
  
        B(iyear)=elem_prod(N(iyear),wgt_mt);
        totN(iyear)=sum(N(iyear));
        totB(iyear)=sum(B(iyear));   
  }
  
  L_total_num=L_fleet1_num;   //landings at age in number fish
  L_total_mt=L_fleet1_mt;     //landings at age in mt whole weight
  
  D_total_num=D_fleet3_num;   //discards at age in number of fish
  D_total_mt=D_fleet3_mt;     //discards at age in mt whole weight

  //Time series of interest  
  B(endyr+1)=elem_prod(N(endyr+1),wgt_mt);
  totN(endyr+1)=sum(N(endyr+1));
  totB(endyr+1)=sum(B(endyr+1));  
  SdS0=SSB/S0;
  
  Fend_mean_temp=1.0;
  for (iyear=1; iyear<=selpar_n_yrs_wgted; iyear++) {Fend_mean_temp*=Fapex(endyr-iyear+1);}
  Fend_mean=pow(Fend_mean_temp,(1.0/selpar_n_yrs_wgted)); //fishing status based on last selpar_n_yrs_wgted years of assessment	  
  if(F_msy_out>0)
    {
      FdF_msy=Fapex/F_msy_out;
      FdF_msy_end=FdF_msy(endyr);
      FdF_msy_end_mean=Fend_mean/F_msy_out;
    }
  if(SSB_msy_out>0)
    {
      SdSSB_msy=SSB/SSB_msy_out;
      SdSSB_msy_end=SdSSB_msy(endyr);
    }  

	if(Fproxy_out>0)
    {
	  FdFproxy=Fapex/Fproxy_out;
	  FdFproxy_end_mean=Fend_mean/Fproxy_out;
	}
  if(SSB_Fproxy_out>0)
    {
      SdSSB_Fproxy=SSB/SSB_Fproxy_out;
	  Sdmsst_Fproxy=SSB/(smsy2msst*SSB_Fproxy_out);
      SdSSB_Fproxy_end=SdSSB_Fproxy(endyr);
	  Sdmsst_Fproxy_end=Sdmsst_Fproxy(endyr);
    }  	
   //fill in log recruitment deviations for yrs they are nonzero
   for(iyear=styr_rec_dev; iyear<=endyr_rec_dev; iyear++)
     {log_rec_dev_output(iyear)=log_rec_dev(iyear);}
   //fill in log Nage deviations for ages they are nonzero (ages2+)
   for(iage=2; iage<=nages; iage++)
     {log_Nage_dev_output(iage)=log_Nage_dev(iage);}
   //fill in log F deviations for yrs they are nonzero
   for(iyear=styr_fleet1_L; iyear<=endyr_fleet1_L; iyear++)
     {log_F1_dev_output(iyear)=log_F_dev_fleet1(iyear);}
   for(iyear=styr_fleet2_L; iyear<=endyr_fleet2_L; iyear++)
     {log_F2_dev_output(iyear)=log_F_dev_fleet2(iyear);}	 
                           
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_projection
  
    switch(Fproj_switch){
       case 1: //F=Fcurrent
          F_reg_proj=Fend_mean;
          break;
       case 2: //F=Fmsy
          F_reg_proj=F_msy_out;
          break;
       case 3: //F=F30
          F_reg_proj=F30_out;
          break;     
       case 4: //F=F40
          F_reg_proj=F40_out;
          break;          		  
       default: // no such switch available
          cout << "Error in input: Projection switch Fproj_switch must be set to 1, 2, 3, or 4." << endl;
          cout << "Presently it is set to " << Fproj_switch <<"."<< endl;
          exit(0);          
   }

  N_proj(styr_proj)=N(endyr+1); //initial conditions computed previously
 
  for (iyear=styr_proj; iyear<=endyr_proj; iyear++) //recruitment follows S-R curve (with bias correction) exactly
  {     
        if (iyear<styr_regs) {F_proj(iyear)=Fend_mean;}
		else {F_proj(iyear)=Fproj_mult*F_reg_proj;}
		
		FL_age_proj=sel_wgted_L*F_proj(iyear);
		
        Z_proj(iyear)=M+FL_age_proj; //+FD_age_proj;
        N_spawn_proj(iyear)(1,nages)=elem_prod(N_proj(iyear)(1,nages),(exp(-1.*(Z_proj(iyear)(1,nages))*spawn_time_frac))); //peak spawning time
		SSB_proj(iyear)= sum(elem_prod(N_spawn_proj(iyear),reprod));
        B_proj(iyear)=sum(elem_prod(N_proj(iyear),wgt_mt)); //uses spawning weight
	     
		for (iage=1; iage<=nages; iage++)
			{L_age_proj(iyear,iage)=N_proj(iyear,iage)*FL_age_proj(iage)*(1.-exp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
			}          
        L_knum_proj(iyear)=sum(L_age_proj(iyear))/1000.0;
	    L_mt_proj(iyear)=sum(elem_prod(L_age_proj(iyear),wgt_wgted_L_mt));     //in mt
		
		if (iyear<endyr_proj) {
			N_proj(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB_proj(iyear),SR_switch);
			N_proj(iyear+1)(2,nages)=++elem_prod(N_proj(iyear)(1,nages-1),(exp(-1.*Z_proj(iyear)(1,nages-1))));
			N_proj(iyear+1,nages)+=N_proj(iyear,nages)*exp(-1.*Z_proj(iyear,nages)); //plus group		
		}
  }
   R_proj=column(N_proj,1);                          

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------     
// FUNCTION get_effective_sample_sizes  //only apply if using multinomial distn for comps; if dirichlet-multinomial, neff computed within cxx file
      // neff_survey1_agec_allyr_out=missing;
      // neff_fleet1_agec_allyr_out=missing;
        
      // for (iyear=1; iyear<=nyr_survey1_agec; iyear++)
         // {if (nsamp_survey1_agec(iyear)>=minSS_survey1_agec)
            // {neff_survey1_agec_allyr_out(yrs_survey1_agec(iyear))=multinom_eff_N(pred_survey1_agec(iyear),obs_survey1_agec(iyear));}                            
          // else {neff_survey1_agec_allyr_out(yrs_survey1_agec(iyear))=-99;}
         // }

      // for (iyear=1; iyear<=nyr_fleet1_agec; iyear++)
         // {if (nsamp_fleet1_agec(iyear)>=minSS_fleet1_agec)
            // {neff_fleet1_agec_allyr_out(yrs_fleet1_agec(iyear))=multinom_eff_N(pred_fleet1_agec(iyear),obs_fleet1_agec(iyear));}                            
          // else {neff_fleet1_agec_allyr_out(yrs_fleet1_agec(iyear))=-99;}
         // }      
   
 //---------------------------------------------------------
						   
FUNCTION evaluate_objective_function
 
  fval=0.0;
  fval_data=0.0;  
  
  //fval=square(xdum-2);
  
//---likelihoods---------------------------

//---Indices-------------------------------
  f_survey1_cpue=lk_lognormal(pred_survey1_cpue, obs_survey1_cpue, survey1_cpue_cv, w_I_survey1);
  fval+=f_survey1_cpue;
  fval_data+=f_survey1_cpue;

//---Landings-------------------------------
  
  //f_fleet1_L in mt
  f_fleet1_L=lk_lognormal(pred_fleet1_L_mt(styr_fleet1_L,endyr_fleet1_L), obs_fleet1_L(styr_fleet1_L,endyr_fleet1_L),
                      fleet1_L_cv(styr_fleet1_L,endyr_fleet1_L), w_L);
  fval+=f_fleet1_L;
  fval_data+=f_fleet1_L;

  //f_fleet2_L in mt
 f_fleet2_L=lk_lognormal(pred_fleet2_L_mt(styr_fleet2_L,endyr_fleet2_L), obs_fleet2_L(styr_fleet2_L,endyr_fleet2_L),
                      fleet2_L_cv(styr_fleet2_L,endyr_fleet2_L), w_L);
  fval+=f_fleet2_L;
  fval_data+=f_fleet2_L;
  
  //f_fleet3_D in mt
  f_fleet3_D=lk_lognormal(pred_fleet3_D_mt(styr_fleet3_D,endyr_fleet3_D), obs_fleet3_D(styr_fleet3_D,endyr_fleet3_D),
                      fleet3_D_cv(styr_fleet3_D,endyr_fleet3_D), w_D);
  fval+=f_fleet3_D;
  fval_data+=f_fleet3_D;
  
//---Age comps-------------------------------

  //f_survey1_agec
  //f_survey1_agec=lk_robust_multinomial(nsamp_survey1_agec, pred_survey1_agec, obs_survey1_agec, nyr_survey1_agec, double(nages_agec), minSS_agec, w_ac_survey1);
  //f_survey1_agec=lk_multinomial(nsamp_survey1_agec, pred_survey1_agec, obs_survey1_agec, nyr_survey1_agec, minSS_agec, w_ac_survey1);
  f_survey1_agec=lk_dirichlet_multinomial(nsamp_survey1_agec, pred_survey1_agec, obs_survey1_agec, nyr_survey1_agec, double(nages_agec), minSS_agec, log_dm_survey1_ac);
  fval+=f_survey1_agec;
  fval_data+=f_survey1_agec;
  
  //f_fleet1_agec
  //f_fleet1_agec=lk_robust_multinomial(nsamp_fleet1_agec, pred_fleet1_agec, obs_fleet1_agec, nyr_fleet1_agec, double(nages_agec), minSS_agec, w_ac_fleet1);
  //f_fleet1_agec=lk_multinomial(nsamp_fleet1_agec, pred_fleet1_agec, obs_fleet1_agec, nyr_fleet1_agec, minSS_agec, w_ac_fleet1);
  f_fleet1_agec=lk_dirichlet_multinomial(nsamp_fleet1_agec, pred_fleet1_agec, obs_fleet1_agec, nyr_fleet1_agec, double(nages_agec), minSS_agec, log_dm_fleet1_ac);
  fval+=f_fleet1_agec;
  fval_data+=f_fleet1_agec;

  //f_fleet2_agec 
  f_fleet2_agec=lk_dirichlet_multinomial(nsamp_fleet2_agec, pred_fleet2_agec, obs_fleet2_agec, nyr_fleet2_agec, double(nages_agec), minSS_agec, log_dm_fleet2_ac);
  fval+=f_fleet2_agec;
  fval_data+=f_fleet2_agec;

//---Length comps----------------------------
  f_fleet3_lenc=lk_dirichlet_multinomial(nsamp_fleet3_lenc, pred_fleet3_lenc, obs_fleet3_lenc, nyr_fleet3_lenc, double(nlenbins), minSS_lenc, log_dm_fleet3_lc);
  fval+=f_fleet3_lenc;
  fval_data+=f_fleet3_lenc;


//-----------Close kin mark recapture--------------------------------
//needs a for loop here!
// nll_ckmr -= log(PHSP + 1e-16) * HSPs(i) + log(PPOP + 1e-16) * POPs(i) + log(1.0 - PHSP + PPOP) * Not(i) 
// fval+=nll_ckmr;
// fval_data+=nll_ckmr;

//-----------Constraints and penalties--------------------------------
  
  //Light penalty applied to log_Nage_dev for deviation from zero. If not estimated, this penalty equals zero.
  f_Nage_init=norm2(log_Nage_dev);        
  fval+=w_Nage_init*f_Nage_init;
  
  // f_rec_dev=0.0;
  // rec_logL_add=nyrs_rec*log(rec_sigma);
  // f_rec_dev=(square(log_rec_dev(styr_rec_dev) + rec_sigma_sq/2.0)/(2.0*rec_sigma_sq));
  // for(iyear=(styr_rec_dev+1); iyear<=endyr; iyear++)
  // {f_rec_dev+=(square(log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1) + rec_sigma_sq/2.0)/
               // (2.0*rec_sigma_sq));}
  // f_rec_dev+=rec_logL_add;            
  // fval+=w_rec*f_rec_dev;
  f_rec_dev=0.0;
  rec_logL_add=0.5*nyrs_rec*log(2.0*pi);
  f_rec_dev=0.5*(square(log_rec_dev(styr_rec_dev))/rec_sigma_sq + log(rec_sigma_sq));
  for(iyear=(styr_rec_dev+1); iyear<=endyr_rec_dev; iyear++)
  {f_rec_dev+=0.5*(square(log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1))/rec_sigma_sq + log(rec_sigma_sq));}
  f_rec_dev+=rec_logL_add;            
  fval+=w_rec*f_rec_dev;    
 
 
  f_rec_dev_early=0.0; //possible extra constraint on early rec deviations
  if (w_rec_early>0.0)
    { if (styr_rec_dev<endyr_rec_phase1)
        {  
          for(iyear=styr_rec_dev; iyear<=endyr_rec_phase1; iyear++)
          {f_rec_dev_early+=square(log_rec_dev(iyear));}
        }
  fval+=w_rec_early*f_rec_dev_early;
  }
  
  f_rec_dev_end=0.0; //possible extra constraint on ending rec deviations
  if (w_rec_end>0.0)
  { if (endyr_rec_phase2<endyr_rec_dev)
        {  
          for(iyear=(endyr_rec_phase2+1); iyear<=endyr_rec_dev; iyear++)
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

  
//---Priors---------------------------------------------------
//neg_log_prior arguments: estimate, prior mean, prior var/-CV, pdf type
//Variance input as a negative value is considered to be CV in arithmetic space (CV=-1 implies loose prior) 
//pdf type 1=none, 2=lognormal, 3=normal, 4=beta 

  f_priors=0.0; 
  f_priors+=neg_log_prior(len_cv_val,set_len_cv(5),set_len_cv(6),set_len_cv(7));

  f_priors+=neg_log_prior(steep,set_steep(5),set_steep(6),set_steep(7)); 
  f_priors+=neg_log_prior(log_R0,set_log_R0(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));

  f_priors+=neg_log_prior(selpar_A50_survey1,set_selpar_A50_survey1(5), set_selpar_A50_survey1(6), set_selpar_A50_survey1(7));
  f_priors+=neg_log_prior(selpar_slope_survey1,set_selpar_slope_survey1(5), set_selpar_slope_survey1(6), set_selpar_slope_survey1(7));
    
  f_priors+=neg_log_prior(selpar_A50_fleet1_B1,set_selpar_A50_fleet1_B1(5), set_selpar_A50_fleet1_B1(6), set_selpar_A50_fleet1_B1(7));
  f_priors+=neg_log_prior(selpar_slope_fleet1_B1,set_selpar_slope_fleet1_B1(5), set_selpar_slope_fleet1_B1(6), set_selpar_slope_fleet1_B1(7));
  f_priors+=neg_log_prior(selpar_A50_fleet2_B1,set_selpar_A50_fleet2_B1(5), set_selpar_A50_fleet2_B1(6), set_selpar_A50_fleet2_B1(7));
  f_priors+=neg_log_prior(selpar_slope_fleet2_B1,set_selpar_slope_fleet2_B1(5), set_selpar_slope_fleet2_B1(6), set_selpar_slope_fleet2_B1(7));

  f_priors+=neg_log_prior(selpar_A50_asc_fleet3_B1,set_selpar_A50_asc_fleet3_B1(5), set_selpar_A50_asc_fleet3_B1(6), set_selpar_A50_asc_fleet3_B1(7));
  f_priors+=neg_log_prior(selpar_slope_asc_fleet3_B1,set_selpar_slope_asc_fleet3_B1(5), set_selpar_slope_asc_fleet3_B1(6), set_selpar_slope_asc_fleet3_B1(7));
  f_priors+=neg_log_prior(selpar_A50_dsc_fleet3_B1,set_selpar_A50_dsc_fleet3_B1(5), set_selpar_A50_dsc_fleet3_B1(6), set_selpar_A50_dsc_fleet3_B1(7));
  f_priors+=neg_log_prior(selpar_slope_dsc_fleet3_B1,set_selpar_slope_dsc_fleet3_B1(5), set_selpar_slope_dsc_fleet3_B1(6), set_selpar_slope_dsc_fleet3_B1(7));
 
  f_priors+=neg_log_prior(log_q_survey1,set_log_q_survey1(5),set_log_q_survey1(6),set_log_q_survey1(7));

  f_priors+=neg_log_prior(F_init,set_F_init(5),set_F_init(6),set_F_init(7));
  f_priors+=neg_log_prior(log_avg_F_fleet1,set_log_avg_F_fleet1(5),set_log_avg_F_fleet1(6),set_log_avg_F_fleet1(7));
  f_priors+=neg_log_prior(log_avg_F_fleet2,set_log_avg_F_fleet2(5),set_log_avg_F_fleet2(6),set_log_avg_F_fleet2(7));
  f_priors+=neg_log_prior(log_avg_F_fleet3,set_log_avg_F_fleet3(5),set_log_avg_F_fleet3(6),set_log_avg_F_fleet3(7));

  f_priors+=neg_log_prior(log_dm_fleet1_ac,set_log_dm_fleet1_ac(5),set_log_dm_fleet1_ac(6),set_log_dm_fleet1_ac(7));
  f_priors+=neg_log_prior(log_dm_fleet2_ac,set_log_dm_fleet2_ac(5),set_log_dm_fleet2_ac(6),set_log_dm_fleet2_ac(7));
  f_priors+=neg_log_prior(log_dm_survey1_ac,set_log_dm_survey1_ac(5),set_log_dm_survey1_ac(6),set_log_dm_survey1_ac(7));
  
  fval+=f_priors;

  //cout << "fval = " << fval << "  fval_data = " << fval_data << endl;

//----------------------------------------------------------------------------------
//Logistic function: 2 parameters
FUNCTION dvar_vector logistic(const dvar_vector& ages, const dvariable& L50, const dvariable& slope)
  //ages=vector of ages, L50=age at 50% selectivity, slope=rate of increase
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1./(1.+exp(-1.*slope*(ages-L50))); //logistic;  
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

  //----------------------------------------------------------------------------------
//Logistic function: 3 parameters (may peak <1)
FUNCTION dvar_vector logistic_peak(const dvar_vector& ages, const dvariable& L50, const dvariable& slope, const dvariable& peak)
  //ages=vector of ages, L50=age at 50% selectivity, slope=rate of increase, peak=asymptote
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=peak/(1.+exp(-1.*slope*(ages-L50))); //logistic;  
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
   if (ages(iage)<joint) {Sel_Tmp(iage)=1./(1.+exp(-1.*slope*(ages(iage)-L50)));}                                                                             
   if (ages(iage)>joint){Sel_Tmp(iage)=exp(-1.*square((ages(iage)-joint)/sigma));}                                                                            
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
  Sel_Tmp=elem_prod( (1./(1.+exp(-1.*slope1*(ages-L501)))),(1.-(1./(1.+exp(-1.*slope2*(ages-(L501+L502)))))) );     
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
   if (double(iage)<joint) {Sel_Tmp(iage)=1./(1.+exp(-1.*slope1*(ages(iage)-L501)));}  
   if (double(iage)>joint){Sel_Tmp(iage)=1.0-(1.0-satval)/(1.+exp(-1.*slope2*(ages(iage)-L502)));}  
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
  pars_tmp(2)=peak+1.0+(0.99*ages(nages)-peak-1.0)/(1.0+exp(-top));
  pars_tmp(3)=exp(ascwid);
  pars_tmp(4)=exp(deswid);
  pars_tmp(5)=1.0/(1.0+exp(-init));
  pars_tmp(6)=1.0/(1.0+exp(-final));
       
  sel_tmp_iq(1)=exp(-(square(ages(1)-pars_tmp(1))/pars_tmp(3)));
  sel_tmp_iq(2)=exp(-(square(ages(nages)-pars_tmp(2))/pars_tmp(4)));
  
  sel_step1=exp(-(square(ages-pars_tmp(1))/pars_tmp(3)));
  sel_step2=pars_tmp(5)+(1.0-pars_tmp(5))*(sel_step1-sel_tmp_iq(1))/(1.0-sel_tmp_iq(1));  
  sel_step3=exp(-(square(ages-pars_tmp(2))/pars_tmp(4)));
  sel_step4=1.0+(pars_tmp(6)-1.0)*(sel_step3-1.0)/(sel_tmp_iq(2)-1.0);
  sel_step5=1.0/ (1.0+exp(-(20.0* elem_div((ages-pars_tmp(1)), (1.0+sfabs(ages-pars_tmp(1)))) )));
  sel_step6=1.0/(1.0+exp(-(20.0*elem_div((ages-pars_tmp(2)),(1.0+sfabs(ages-pars_tmp(2)))) )));  

  Sel_Tmp=elem_prod(sel_step2,(1.0-sel_step5))+ 
          elem_prod(sel_step5,((1.0-sel_step6)+ elem_prod(sel_step4,sel_step6)) ); 
 
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
    
///-----------------------------------------------------------------------------------    
//Spawner-recruit function (Beverton-Holt, Ricker, or Null)
FUNCTION dvariable SR_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& SSB, int func)
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, SSB=spawning biomass
  //func=1 for Beverton-Holt, 2 for Ricker, 3 for Null
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  switch(func) {
    case 1: //Beverton-Holt
      Recruits_Tmp=((0.8*R0*h*SSB)/(0.2*R0*spr_F0*(1.0-h)+(h-0.2)*SSB));       
    break;
    case 2: //Ricker
      Recruits_Tmp=((SSB/spr_F0)*exp(h*(1-SSB/(R0*spr_F0))));       
    break;
	case 3: //Null
      Recruits_Tmp=R0;       
    break;
  }
  RETURN_ARRAYS_DECREMENT();
  return Recruits_Tmp;
  
  
//-----------------------------------------------------------------------------------    
//Spawner-recruit equilibrium function (Beverton-Holt, Ricker, or Null)
FUNCTION dvariable SR_eq_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& spr_F, const dvariable& BC, int func)
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, spr_F=spawners per recruit @ F, BC=bias correction
  //func=1 for Beverton-Holt, 2 for Ricker, 3 for Null
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  switch(func) {
    case 1: //Beverton-Holt
      Recruits_Tmp=(R0/((5.0*h-1.0)*spr_F))*(BC*4.0*h*spr_F-spr_F0*(1.0-h));    
    break;
    case 2: //Ricker
      Recruits_Tmp=R0/(spr_F/spr_F0)*(1.0+log(BC*spr_F/spr_F0)/h);      
    break;
    case 3: //None
      Recruits_Tmp=BC*R0;      
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
  //small_number is small value to avoid log(0) during searfleet1
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
  dvariable small_number=0.00001;
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
  dvariable small_number=0.00001;
  LkvalTmp=0.0;
  dvar_matrix Eprime=elem_prod((1.0-obs_comp), obs_comp)+0.1/mbin; //E' of Francis 2011, p.1131  
  dvar_vector nsamp_wgt=nsamp*wgt_dat;
 
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {LkvalTmp+= sum(0.5*log(Eprime(ii))-log(small_number+exp(elem_div((-square(obs_comp(ii)-pred_comp(ii))) , (Eprime(ii)*2.0/nsamp_wgt(ii)) ))) );
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
  LkvalTmp=0.0; 
  dvar_vector nsamp_adjust=nsamp*exp(log_dir_par);
  
  for (int ii=1; ii<=ncomp; ii++)
  {
	if (nsamp(ii)>=minSS)
    {
		LkvalTmp-=gammln(nsamp_adjust(ii))-gammln(nsamp(ii)+nsamp_adjust(ii));
		LkvalTmp-=sum(gammln(nsamp(ii)*obs_comp(ii)+nsamp_adjust(ii)*pred_comp(ii)));
        LkvalTmp+=sum(gammln(nsamp_adjust(ii)*pred_comp(ii)));		
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: logistic normal (aka multivariate logistic in iSCAM; logistic normal in Francis' terminology)
//NOTE, lk_logisitic_normal function is relatively untested
FUNCTION dvariable lk_logistic_normal(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS)
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.0001;
  LkvalTmp=0.0;
  dvar_matrix nu=pred_comp+0.0;
  dvar_matrix pred_plus=pred_comp+small_number;
  dvar_matrix obs_plus=obs_comp+small_number;

  dvariable nu_mean;
  dvariable nu_sum_sq;
  dvariable tau_hat_sq;
  dvariable year_count; //keeps track of years included in likelihood (i.e., that meet the sample size requirement)
   
  LkvalTmp=0.0;
  nu_sum_sq=0.0;
  year_count=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {
		year_count+=1.0;
		nu_mean=sum( log(obs_plus(ii))-log(pred_plus(ii))  )/mbin;	//year-specific mean log residual
		for (int jj=1; jj<=mbin;jj++)
		{
			nu(ii,jj) = log(obs_plus(ii,jj)) - log(pred_plus(ii,jj)) - nu_mean;
			nu_sum_sq += square(nu(ii,jj));
		}
    }
  }  
  if (year_count>0.0)
  {
	  tau_hat_sq = nu_sum_sq/((mbin-1.0)*year_count);
	  LkvalTmp = (mbin-1.0)*year_count*log(tau_hat_sq);
  }
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

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
          if(ab_iq<=0) {
			cout << "Parameter input error: For beta priors, mu*(1-mu)/var must be greater than one. Try decreasing var." << endl;
            exit(0);
		  }
		  if(pred>=0 && pred<=1) LkvalTmp= (1.0-alpha)*log(pred)+(1.0-beta)*log(1.0-pred)-gammln(alpha+beta)+gammln(alpha)+gammln(beta);
          else LkvalTmp=big_number;
          break;
        default: // no such prior pdf currently available
          cout << "Parameter input error: Prior must be either 1(none), 2(lognormal), 3(normal), or 4(beta)." << endl;
          cout << "Presently at least one is " << pdf << endl;
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
      cout << "TotalLikelihood " << fval << endl;
	  //cout << "xdum = " << xdum << endl;
    if (last_phase())  
    {

      //cout<<"start report"<<endl;
       get_weighted_current();
      //cout<<"got weighted"<<endl;
       get_msy();
      //cout<<"got msy"<<endl;
	  get_per_recruit_stuff();
      //cout<<"got per recruit"<<endl;  
      get_miscellaneous_stuff();
      //cout<<"got misc stuff"<<endl;
	  get_projection();
      //cout<<"got projection"<<endl;
      // get_effective_sample_sizes();
      //cout<<"got effective_sample_sizes"<<endl;
	  
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
	  cout<<"*******************************************"<<endl;
      
      cout <<endl;     
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;
      cout << "Fmsy=" << F_msy_out<< "   SSBmsy=" << SSB_msy_out <<endl;
      cout << "F status="<<FdF_msy_end<<endl;
      cout << "Pop status="<<SdSSB_msy_end<<endl;
      switch(SR_switch) {
	    case 1: //Beverton-Holt
			cout << "Recruitment model = Beverton-Holt" <<endl;
			cout << "h="<<steep<<"   R0="<<R0<<endl;
			break;
		case 2: //Ricker
			cout << "Recruitment model = Ricker" <<endl;
			cout << "h="<<steep<<"   R0="<<R0<<endl; 
			break;
		case 3: //Null
			cout << "Recruitment model = Null" <<endl;
			cout << "R0="<<R0<<endl;
			break;	
  }

	  
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;  
      
      //report << "TotalLikelihood " << fval << endl;
      //report << "N" << endl;
      //report << N<<endl;
      //report << "F" << endl;
      //report << F <<endl;   

      sdnr_ac_survey1=sdnr_multinomial(nyr_survey1_agec, agebins_agec, nsamp_survey1_agec, pred_survey1_agec, obs_survey1_agec, w_ac_survey1);
      sdnr_ac_fleet1=sdnr_multinomial(nyr_fleet1_agec, agebins_agec, nsamp_fleet1_agec, pred_fleet1_agec, obs_fleet1_agec, w_ac_fleet1);
	  sdnr_ac_fleet2=sdnr_multinomial(nyr_fleet2_agec, agebins_agec, nsamp_fleet2_agec, pred_fleet2_agec, obs_fleet2_agec, w_ac_fleet2);
      sdnr_I_survey1=sdnr_lognormal(pred_survey1_cpue, obs_survey1_cpue, survey1_cpue_cv, w_I_survey1);
       
      //#################################################################################################
      //##  Pass parameters to vector for bounds check plotting
      //################################################################################################# 

       Linf_out(8)=Linf; Linf_out(1,7)=set_Linf; 
       K_out(8)=K; K_out(1,7)=set_K;
       t0_out(8)=t0; t0_out(1,7)=set_t0;
       len_cv_val_out(8)=len_cv_val; len_cv_val_out(1,7)=set_len_cv;
       
	   M_logit_out(8)=M_logit; M_logit_out(1,7)=set_M_logit;
	   
       log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
       steep_out(8)=steep; steep_out(1,7)=set_steep;
       rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
       R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;
 
	   log_dm_fleet1_ac_out(8)=log_dm_fleet1_ac; log_dm_fleet1_ac_out(1,7)=set_log_dm_fleet1_ac;
	   log_dm_fleet2_ac_out(8)=log_dm_fleet2_ac; log_dm_fleet2_ac_out(1,7)=set_log_dm_fleet2_ac;
       log_dm_fleet3_lc_out(8)=log_dm_fleet3_lc; log_dm_fleet3_lc_out(1,7)=set_log_dm_fleet3_lc;
	   log_dm_survey1_ac_out(8)=log_dm_survey1_ac; log_dm_survey1_ac_out(1,7)=set_log_dm_survey1_ac;
 
       selpar_A50_survey1_out(8)=selpar_A50_survey1; selpar_A50_survey1_out(1,7)=set_selpar_A50_survey1;
       selpar_slope_survey1_out(8)=selpar_slope_survey1; selpar_slope_survey1_out(1,7)=set_selpar_slope_survey1;
	   
       selpar_A50_fleet1_B1_out(8)=selpar_A50_fleet1_B1; selpar_A50_fleet1_B1_out(1,7)=set_selpar_A50_fleet1_B1;
       selpar_slope_fleet1_B1_out(8)=selpar_slope_fleet1_B1; selpar_slope_fleet1_B1_out(1,7)=set_selpar_slope_fleet1_B1;
       selpar_A50_fleet2_B1_out(8)=selpar_A50_fleet2_B1; selpar_A50_fleet2_B1_out(1,7)=set_selpar_A50_fleet2_B1;
       selpar_slope_fleet2_B1_out(8)=selpar_slope_fleet2_B1; selpar_slope_fleet2_B1_out(1,7)=set_selpar_slope_fleet2_B1;
       selpar_A50_asc_fleet3_B1_out(8)=selpar_A50_asc_fleet3_B1; selpar_A50_asc_fleet3_B1_out(1,7)=set_selpar_A50_asc_fleet3_B1;
       selpar_slope_asc_fleet3_B1_out(8)=selpar_slope_asc_fleet3_B1; selpar_slope_asc_fleet3_B1_out(1,7)=set_selpar_slope_asc_fleet3_B1;
       selpar_A50_dsc_fleet3_B1_out(8)=selpar_A50_dsc_fleet3_B1; selpar_A50_dsc_fleet3_B1_out(1,7)=set_selpar_A50_dsc_fleet3_B1;
       selpar_slope_dsc_fleet3_B1_out(8)=selpar_slope_dsc_fleet3_B1; selpar_slope_dsc_fleet3_B1_out(1,7)=set_selpar_slope_dsc_fleet3_B1;
	   
       log_q_survey1_out(8)=log_q_survey1; log_q_survey1_out(1,7)=set_log_q_survey1;
                    
       log_avg_F_fleet1_out(8)=log_avg_F_fleet1; log_avg_F_fleet1_out(1,7)=set_log_avg_F_fleet1;
	     log_avg_F_fleet2_out(8)=log_avg_F_fleet2; log_avg_F_fleet2_out(1,7)=set_log_avg_F_fleet2;
	     log_avg_F_fleet3_out(8)=log_avg_F_fleet3; log_avg_F_fleet3_out(1,7)=set_log_avg_F_fleet3;
       F_init_out(8)=F_init; F_init_out(1,7)=set_F_init;
       
       log_rec_dev_out(styr_rec_dev, endyr_rec_dev)=log_rec_dev;
       log_F_dev_fleet1_out(styr_fleet1_L,endyr_fleet1_L)=log_F_dev_fleet1;
       log_F_dev_fleet2_out(styr_fleet2_L,endyr_fleet2_L)=log_F_dev_fleet2;
       log_F_dev_fleet3_out(styr_fleet3_D,endyr_fleet3_D)=log_F_dev_fleet3;
       #include "BAM-Sim.cxx"   // write the R-compatible report

    } //endl last phase loop     
  
