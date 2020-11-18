{
    // Diego Venegas Vargas
    // Prospect
    // Neutrino Directionality Study
    // #include "TMath.h"
    #include <fstream>
    #include <iostream>
    #include <vector>
    #include <string>
    #include <cmath>
    #include <cctype>
    #include <ctime>   
    #include <algorithm>
    #include <TSystem.h>
    using namespace std;
    // SetupProspectStyle();
    // First, we define variables and parameters
    TCanvas* c1 = new TCanvas("c1","x"); 
    TCanvas* c2 = new TCanvas("c2","y"); 
    TCanvas* c3 = new TCanvas("c3","z"); 
    TCanvas* c4 = new TCanvas("c4","r"); 
    // TCanvas* c5 = new TCanvas("c5","z_dist"); 
    // TCanvas* c6 = new TCanvas("c6","no_zcut"); 
    // TCanvas* c7 = new TCanvas("c7","3D IBD"); 
    // string a = "/home/prospect-collab/converted_data/Analyzed/Analyzed_2020A_IBD_v23.1/"; // path to file
    string a = "/p/lustre2/psptexp/converted/analyzed/Analyzed_2020A_IBD_v23.1/"; // path to file
    string b = "/AD1_IBD_2020.root"; // name of each root file for all files
    double Emin = 0.0; // minum value of energy
    double Emax = 12.0; // maximum value of energy
    int E_bins = 60; // number of bins in Energy histogram
    double E_cut_min = 0.8;
    double E_cut_max = 7.2;
    double dE_min = 0.0; // min value of Esmear-Emax_seg
    double dE_max = 1.2; // min value of Esmear-Emax_seg
    int dE_bins = 60;// number of bins for Esmear-Emax_seg histogram
    double bin_width= E_bins/(Emax-Emin);
    double atm_scaling=1.00025443769309;//atmospheric scaling factor
    int x_bins=14;
    double x_max=13.500000;
    double x_min=-0.50000000;
    double x_bw = (x_max-x_min)/x_bins;
    int y_bins=11;
    double y_max=10.500000;
    double y_min=-0.50000000;
    double y_bw = (y_max-y_min)/y_bins;
    double z_max_mm=400.00; //default cell_length is [mm]
    double z_min_mm=-400.00;
    int z_bins=18;
    const int s_array = 18;
    double x_pos_true = 5.969; // true value for x distance
    double y_pos_true = 5.086; // true value for y distance
    double z_pos_true = 1.19; // true value for z distance
    double r_true = TMath::Sqrt(TMath::Power((x_pos_true),2) + TMath::Power((y_pos_true),2) + TMath::Power((z_pos_true),2)); 
    double x_pos [s_array]; // array to hold calculated x positions from fit
    double y_pos [s_array]; // array to hold calculated y positions from fit
    double z_pos [s_array]; // array to hold calculated z positions from fit
    double r [s_array]; // array to hold calculated r positions from fit
    double x_pos_err [s_array]; // array to hold calculated x positions errors from fit
    double y_pos_err [s_array]; // array to hold calculated y positions errors from fit
    double z_pos_err [s_array]; // array to hold calculated z positions errors from fit
    double r_err [s_array]; // array to hold calculated z positions errors from fit
    double x_sim_pos [s_array]; // array to hold calculated x simulation positions from fit
    double y_sim_pos [s_array]; // array to hold calculated y simulation positions from fit
    double z_sim_pos [s_array]; // array to hold calculated z simulation positions from fit
    double r_sim [s_array]; // array to hold calculated r simulation positions from fit
    double x_sim_pos_err [s_array]; // array to hold calculated x simulation positions errors from fit
    double y_sim_pos_err [s_array]; // array to hold calculated y simulation positions errors from fit
    double z_sim_pos_err [s_array]; // array to hold calculated z simulation positions errors from fit
    double r_sim_err [s_array]; // array to hold calculated z simulation positions errors from fit
    double nbinsz [s_array]; // array which will hold the x values of all graphs (bins in z)
    double nbinsz_err[s_array]; // array for error in x for each graph, all zeros
    double nbinsz_sim [s_array]; // array which will hold the x values of all graphs (bins in z)
    double nbinsz_sim_err[s_array]; // array for error in x for each graph, all zeros
    double x_const = 6.5;
    double y_const = 5.0;
    double length_const = 0.14;
    double z_max=z_max_mm/800.0; //default cell_length is [mm]
    double z_min=z_min_mm/800.0; 

    TH1F* x_hist = new TH1F("x_hist", "x_hist", z_bins, 0.5, 18.5);
    TH1F* y_hist = new TH1F("y_hist", "y_hist", z_bins, 0.5, 18.5);
    TH1F* z_hist = new TH1F("z_hist", "z_hist", z_bins, 0.5, 18.5);
    TH1F* r_hist = new TH1F("r_hist", "r_hist", z_bins, 0.5, 18.5);
    TH1F* x_hist_sim = new TH1F("x_hist_sim", "x_hist_sim", z_bins, 0.5, 18.5);
    TH1F* y_hist_sim = new TH1F("y_hist_sim", "y_hist_sim", z_bins, 0.5, 18.5);
    TH1F* z_hist_sim = new TH1F("z_hist_sim", "z_hist_sim", z_bins, 0.5, 18.5);
    TH1F* r_hist_sim = new TH1F("r_hist_sim", "r_hist_sim", z_bins, 0.5, 18.5);



    for (int cur_zbin = 1;cur_zbin<=z_bins;cur_zbin++){
        double cell_length = (z_max_mm-z_min_mm)/cur_zbin;
        // double z_max=z_max_mm/cell_length; //default cell_length is [mm]
        // double z_min=z_min_mm/cell_length;
        double z_bw = (z_max-z_min)/cur_zbin;
        std::cout<<"cell_length:"<<cell_length<<std::endl;
        std::cout<<"z_max:"<<z_max<<std::endl;
        std::cout<<"z_min:"<<z_min<<std::endl;
        std::cout<<"z_bw:"<<z_bw<<std::endl;
        Float_t Esmear; // Energy 
        Int_t maxseg; // Max Energy deposited in cell
        Int_t segmult; // Max Energy deposited in cell
        Float_t ncapt_dt; //neutron acpture time
        Float_t prompt_psd; //prompt signal PSD
        Float_t xyz[1][3]; // position of prompt signal (mm)
        Float_t n_xyz[1][3]; // position of prompt signal (mm)
        double on_live_time=0.0; //live time for Rxon periods
        double on_run_time=0.0; //live time for Rxon periods
        double off_run_time=0.0; //live time for  Rxoff periods
        double off_live_time=0.0; //live time for  Rxoff periods
        double sim_live_time=0.0; //live time for simulation
        double sim_run_time=0.0; //live time for simulation
        nbinsz[cur_zbin-1]= cur_zbin; 
        nbinsz_err[cur_zbin-1]= 0.0; 
        ////////////////////////////////////////////////////////////////////////////////
        // We begin by obtaining the segment map
        TH1F* seg_map = new TH1F("seg_map","seg_map",154,-0.5,153.5);//seg_map histogram
        // TFile* base_file = new TFile("/home/dcvenega/baseline.root"); // open baseline.root file which contains the PositionTree101 (Flat Tree of the Segment Baselines)
        TFile* base_file = new TFile("/g/g90/dvargas/nu_dir/baseline.root"); // open baseline.root file which contains the PositionTree101 (Flat Tree of the Segment Baselines)
        TTree* ttree = (TTree*) base_file->Get("PositionTree101");// get "PositionTree101" tree
        Long64_t NEntries = ttree->GetEntries(); // number of entreies in tree
        Int_t Segment;
        Double_t Baseline;
        ttree->SetBranchAddress("Segment",&Segment); //grab segment branch
        ttree->SetBranchAddress("Baseline",&Baseline); //grab Baseline branch
        for(Long_t i = 0; i< NEntries ; i ++)
        {
            ttree->GetEntry(i);
            seg_map->SetBinContent(i+1,Baseline); // fill
        }
        base_file->Close();
        //////////////////////////////////////////////////////////////////////////////////////
        
        //////////////////////////////////////////////////////////////////////////////////////
        // Now we define histograms that will hold the total RxOn(Off) data, as well as the IBD scpectrum
        // This histograms are 2D histograms since we will determine the number of IBD detected in each cell
        // currently only using combined_Huber-Standard_map.root updated with p2x from AndrewM
        // All histograms are cloned from the segment_efficiency histogram, just to keep an uniform format
        // TFile* f2 = new TFile("/home/dcvenega/efficiency.root");
        TFile* f2 = new TFile("/g/g90/dvargas/nu_dir/efficiency.root");
        if ( f2->IsZombie() ){
            std::cout << "Error Opening file : " << "files/eff_map.root"  << "\n" ;
        }
        TH2D* efficiency = ( TH2D* ) f2->Get("segment_efficiency");
        TH2F* seg_baseline =  ( TH2F* ) efficiency->Clone("seg_baseline");
        TH2F* cell_map =  ( TH2F* ) efficiency->Clone("cell_map");
        cell_map->Scale(0.0);
        seg_baseline->Scale(0.0);
        // 3D histograms
        TH3D* RxOn_3D = new TH3D("RxOn_3D","RxOn_3D",x_bins, x_min, x_max,y_bins, y_min, y_max,cur_zbin, z_min, z_max);
        TH3D* RxOff_3D = new TH3D("RxOff_3D","RxOff_3D",x_bins, x_min, x_max,y_bins, y_min, y_max,cur_zbin, z_min, z_max); 
        TH3D* RxIBD_3D = new TH3D("RxIBD_3D","RxIBD_3D",x_bins, x_min, x_max,y_bins, y_min, y_max,cur_zbin, z_min, z_max); 
        TH3D* IBD_sim_3D = new TH3D("IBD_sim_3D","IBD_sim_3D",x_bins, x_min, x_max,y_bins, y_min, y_max,cur_zbin, z_min, z_max); 
        //////////////////////////////////////////////////////////////////////////////////////

        // //////////////////////////////////////////////////////////////////////////////////////
        // //Here we loop through the files on the 2019XList, which contains both RxOn and RxOff files
        // // file names with a 0(1) at the end indicate a RxOff(On) period
        // int count=0;
        ifstream list("2019Xlist.txt");
        for( std::string line; getline( list, line ); ){ // get first root file
            // REACTOR ON FILES
            if(line[line.size()-1] == '1'){
                string file_name = a + line.substr(0,line.size()-2)+b;
                if(!gSystem->AccessPathName(file_name.c_str())){//gSystem->AccessPathName returns true if file_name doesn't exist
                    TFile* temp_file= new TFile (file_name.c_str()); // open it
                    // TFile* temp_file= new TFile ("/home/prospect-collab/converted_data/Analyzed/Analyzed_2020A_IBD_v23.1/WetCommissioning/series015/s015_f00000_ts1520293010/AD1_IBD_2020.root"); // open it
                    TTree* Tibd=(TTree*) temp_file->Get("P2kIBDPlugin/Tibd"); // get the Tibd tree
                    TVectorT<double>* temp_time = (TVectorT<double>*) temp_file->Get("runtime") ; // obtain the runtime for each file
                    TVectorT<double>* prompt_time = (TVectorT<double>*) temp_file->Get("accumulated/P2kIBDPlugin.tveto_prompt") ; // obtain the prompt time for each file
                    TVectorT<double>* delayed_time = (TVectorT<double>*) temp_file->Get("accumulated/P2kIBDPlugin.tveto_delayed") ; // obtain the delayed time for each file
                    double a = (temp_time->Max());
                    double b = prompt_time->Max();
                    double c = delayed_time->Max();
                    double temp_dt_on = (a/(a-b)) * (a/(a-c));// calculated deadtime for file
                    on_run_time = temp_time->Max() + on_run_time; // add runtime's file to RxOn's total runtime
                    on_live_time =  (temp_time->Max()/temp_dt_on) + on_live_time ; // add livetime's file to total livetime
                    // get branches of Tibd tree
                    Tibd->SetMakeClass(1); 
                    Tibd->SetBranchStatus("*", 0); // clear branch
                    Tibd->SetBranchStatus("Esmear", 1); 
                    Tibd->SetBranchAddress("Esmear",&Esmear);
                    Tibd->SetBranchStatus("maxseg", 1); 
                    Tibd->SetBranchAddress("maxseg",&maxseg);
                    Tibd->SetBranchStatus("ncapt_dt", 1); 
                    Tibd->SetBranchAddress("ncapt_dt",&ncapt_dt);
                    Tibd->SetBranchStatus("xyz", 1); 
                    Tibd->SetBranchAddress("xyz",&xyz);
                    Tibd->SetBranchStatus("n_xyz", 1); 
                    Tibd->SetBranchAddress("n_xyz",&n_xyz);
                    // define histograms for accidental and correlated+accidental events
                    // Temporary for 3D histograms
                    TH3D* RxOn_Cor_3D = ( TH3D* ) RxOn_3D->Clone("RxOn_Cor_3D");
                    TH3D* RxOn_Acc_3D = ( TH3D* ) RxOn_3D->Clone("RxOn_Acc_3D");
                    TH3D* On_temp_3D = ( TH3D* ) RxOn_3D->Clone("On_temp_3D");
                    On_temp_3D->Scale(0.0);
                    RxOn_Cor_3D->Scale(0.0);
                    RxOn_Acc_3D->Scale(0.0);
                    Long64_t Nentries = Tibd->GetEntries(); // get number of entries of tree
                    for(Long_t i = 0; i< Nentries ; i ++){ // loop through all entries
                        Tibd->GetEntry(i);// grab ith entry
                        if ( Esmear <= 7.2 && Esmear >= 0.8 ) { // IBD energy cut
                            Int_t x = maxseg % 14; // x position of segment
                            Int_t y = maxseg / 14; // y position of segment
                            Float_t z = (xyz[0][2])/(cell_length); // z position along the cell in [mm]
                            if (ncapt_dt/1000. > 1. && ncapt_dt/1000. < 120. ) { // correlated+accidental events condition
                                RxOn_Cor_3D->Fill(x,y,z);
                            } else { // accidental only events
                                RxOn_Acc_3D->Fill(x,y,z);
                            }
                        }
                    }
                    //Now that we have accidental and correlated histograms filled, we scale accidental events by 
                    // 0.01 times the deadtime correction and subract them from the correlated+accidental events
                    // Then we add the resulting histogram to the final histogram for reactor On
                    On_temp_3D->Add(RxOn_Cor_3D);
                    On_temp_3D->Add(RxOn_Acc_3D,-0.01*temp_dt_on);
                    RxOn_3D->Add(On_temp_3D);
                    //delete temporary histograms
                    delete RxOn_Cor_3D;
                    delete RxOn_Acc_3D;
                    delete On_temp_3D;
                    delete temp_file;
                }
            }
            else{
                // REACTOR OFF FILES
                if(line[line.size()-1] == '0'){
                string file_name = a + line.substr(0,line.size()-2)+b;
                    if(!gSystem->AccessPathName(file_name.c_str())){//gSystem->AccessPathName returns true if file_name doesn't exist
                        TFile* temp_file= new TFile (file_name.c_str()); // open it
                        TTree* Tibd=(TTree*) temp_file->Get("P2kIBDPlugin/Tibd"); // get the Tibd tree
                        TVectorT<double>* temp_time = (TVectorT<double>*) temp_file->Get("runtime") ; // obtain the runtime for each file
                        TVectorT<double>* prompt_time = (TVectorT<double>*) temp_file->Get("accumulated/P2kIBDPlugin.tveto_prompt") ; // obtain the prompt time for each file
                        TVectorT<double>* delayed_time = (TVectorT<double>*) temp_file->Get("accumulated/P2kIBDPlugin.tveto_delayed") ; // obtain the delayed time for each file
                        double a = (temp_time->Max());
                        double b = prompt_time->Max();
                        double c = delayed_time->Max();
                        double temp_dt_off = (a/(a-b)) * (a/(a-c));// calculated deadtime for file
                        off_run_time = temp_time->Max() + off_run_time; // add runtime's file to RxOff's total runtime
                        off_live_time =  (temp_time->Max()/temp_dt_off) + off_live_time ; // add livetime's file to total livetime
                        // get branches of Tibd tree
                        Tibd->SetMakeClass(1); 
                        Tibd->SetBranchStatus("*", 0); // clear branch
                        Tibd->SetBranchStatus("Esmear", 1); 
                        Tibd->SetBranchAddress("Esmear",&Esmear);
                        Tibd->SetBranchStatus("maxseg", 1); 
                        Tibd->SetBranchAddress("maxseg",&maxseg);
                        Tibd->SetBranchStatus("ncapt_dt", 1); 
                        Tibd->SetBranchAddress("ncapt_dt",&ncapt_dt);
                        Tibd->SetBranchStatus("xyz", 1); 
                        Tibd->SetBranchAddress("xyz",&xyz);
                        Tibd->SetBranchStatus("n_xyz", 1); 
                        Tibd->SetBranchAddress("n_xyz",&n_xyz);
                        // define histograms for accidental and correlated+accidental events
                        // Temporary for 3D histograms
                        TH3D* RxOff_Cor_3D = ( TH3D* ) RxOff_3D->Clone("RxOff_Cor_3D");
                        TH3D* RxOff_Acc_3D = ( TH3D* ) RxOff_3D->Clone("RxOff_Acc_3D");
                        TH3D* Off_temp_3D = ( TH3D* ) RxOff_3D->Clone("Off_temp_3D");
                        Off_temp_3D->Scale(0.0);
                        RxOff_Cor_3D->Scale(0.0);
                        RxOff_Acc_3D->Scale(0.0);
                        Long64_t Nentries = Tibd->GetEntries(); // get number of entries of tree
                        for(Long_t i = 0; i< Nentries ; i ++){ // Loop through all entries
                            Tibd->GetEntry(i); // get ith entry of tree
                            if ( Esmear <= 7.2 && Esmear >= 0.8 ) { //IBD energy cut
                                Int_t x = maxseg % 14; // x position of cell
                                Int_t y = maxseg / 14; // y position of cell
                                Float_t z = (xyz[0][2])/(cell_length); // z position along the cell
                                if (ncapt_dt/1000. > 1. && ncapt_dt/1000. < 120. ) { // correlated+accidental condition
                                    RxOff_Cor_3D->Fill(x,y,z);
                                } else {
                                    RxOff_Acc_3D->Fill(x,y,z);
                                }
                            }
                        }
                        //Now that we have accidental and correlated histograms filled, we scale accidental events by 
                        // 0.01 times the deadtime correction and subract them from the correlated+accidental events
                        // Then we add the resulting histogram to the final histogram for reactor Off
                        Off_temp_3D->Add(RxOff_Cor_3D);
                        Off_temp_3D->Add(RxOff_Acc_3D,-0.01*temp_dt_off);
                        RxOff_3D->Add(Off_temp_3D);
                        //delete temporary histograms
                        delete RxOff_Cor_3D;
                        delete RxOff_Acc_3D;
                        delete Off_temp_3D;
                        delete temp_file;
                    }
                }
            }
            // count++;
        }

        // SIMULATION HISTOGRAMS/// 
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // First we fill up the histograms for reactor on fromt the full 2019Blist.txt
        ifstream list_sim("sim_files.txt");
        for( std::string line; getline( list_sim, line ); ){ // get first root file
            // numfiles++;
            if(!gSystem->AccessPathName(line.c_str())){//gSystem->AccessPathName returns true if line doesn't exist
                ///////////////////////////////////////////////////////////////////////////////////
                TFile* temp_file= new TFile(line.c_str()); // open it
                TTree* Tibd=(TTree*) temp_file->Get("P2kIBDPlugin/Tibd"); // get the Tibd tree
                TVectorT<double>* temp_time = (TVectorT<double>*) temp_file->Get("runtime") ; // obtain the runtime for each file
                TVectorT<double>* prompt_time = (TVectorT<double>*) temp_file->Get("accumulated/P2kIBDPlugin.tveto_prompt") ; // obtain the deadtime for each file
                TVectorT<double>* delayed_time = (TVectorT<double>*) temp_file->Get("accumulated/P2kIBDPlugin.tveto_delayed") ; // obtain the deadtime for each file
                double a = (temp_time->Max());
                double b = prompt_time->Max();
                double c = delayed_time->Max();
                double temp_dt_sim = (a/(a-b)) * (a/(a-c));
                sim_run_time = temp_time->Max() + sim_run_time; // total runtime
                sim_live_time =  (temp_time->Max()/temp_dt_sim) + sim_live_time ; // total livetime
                Tibd->SetMakeClass(1); 
                Tibd->SetBranchStatus("*", 0); // clear branch
                Tibd->SetBranchStatus("Esmear", 1); 
                Tibd->SetBranchAddress("Esmear",&Esmear);
                Tibd->SetBranchStatus("maxseg", 1); 
                Tibd->SetBranchAddress("maxseg",&maxseg);
                Tibd->SetBranchStatus("ncapt_dt", 1); 
                Tibd->SetBranchAddress("ncapt_dt",&ncapt_dt);
                Tibd->SetBranchStatus("xyz", 1); 
                Tibd->SetBranchAddress("xyz",&xyz);
                Tibd->SetBranchStatus("n_xyz", 1); 
                Tibd->SetBranchAddress("n_xyz",&n_xyz);
                TH3D* Sim_Cor_3D = ( TH3D* ) IBD_sim_3D->Clone("Sim_Cor_3D");
                TH3D* Sim_Acc_3D = ( TH3D* ) IBD_sim_3D->Clone("Sim_Acc_3D");
                TH3D* Sim_temp_3D = ( TH3D* ) IBD_sim_3D->Clone("Sim_temp_3D");
                Sim_temp_3D->Scale(0.0);
                Sim_Cor_3D->Scale(0.0);
                Sim_Acc_3D->Scale(0.0);
                Long64_t Nentries = Tibd->GetEntries(); // get number of entries of tree
                for(Long_t i = 0; i< Nentries ; i ++){ // Loop through all entries
                    Tibd->GetEntry(i); // get ith entry of tree
                    if ( Esmear <= 7.2 && Esmear >= 0.8 ) { //IBD energy cut
                        Int_t x = maxseg % 14; // x position of cell
                        Int_t y = maxseg / 14; // y position of cell
                        Float_t z = (xyz[0][2])/(cell_length); // z position along the cell
                        if (ncapt_dt/1000. > 1. && ncapt_dt/1000. < 120. ) { // correlated+accidental condition
                            Sim_Cor_3D->Fill(x,y,z);
                        } else {
                            Sim_Acc_3D->Fill(x,y,z);
                        }
                    }
                }
                //Now that we have accidental and correlated histograms filled, we scale accidental events by 
                // 0.01 times the deadtime correction and subract them from the correlated+accidental events
                // Then we add the resulting histogram to the final histogram for reactor Sim
                Sim_temp_3D->Add(Sim_Cor_3D);
                Sim_temp_3D->Add(Sim_Acc_3D,-0.01*temp_dt_sim);
                IBD_sim_3D->Add(Sim_temp_3D);
                //delete temporary histograms
                delete Sim_Cor_3D;
                delete Sim_Acc_3D;
                delete Sim_temp_3D;
                delete temp_file;
            }
        }    
        // //////////////////////////////////////////////////////////////////////////////////////
        // // Now we generate the final IBD histograms 
        RxIBD_3D->Add(RxOn_3D);
        RxIBD_3D->Add(RxOff_3D,-on_live_time/off_live_time); // final 3D histogram for IBD 
        // //////////////////////////////////////////////////////////////////////////////////////

        // //////////////////////////////////////////////////////////////////////////////////////
        // // Correction for efficiency 
        for(Int_t z=1; z<=cur_zbin;z++){
            for ( Int_t lx = 0 ; lx <= 13 ; lx ++  ){
                for ( Int_t ly = 0 ; ly <= 10 ; ly ++  ){
                    Int_t x = lx + 1;
                    Int_t y = ly + 1;
                    Double_t OnContent = RxIBD_3D->GetBinContent(x,y,z);
                    Double_t SimContent = IBD_sim_3D->GetBinContent(x,y,z);

                    Double_t OnError = RxIBD_3D->GetBinError(x,y,z) ;
                    Double_t SimError = IBD_sim_3D->GetBinError(x,y,z) ;

                    Double_t eff = efficiency->GetBinContent(x,y,z);

                    if ( eff == 0 ) {
                        RxIBD_3D->SetBinContent( x, y, z, 0 );
                        RxIBD_3D->SetBinError( x, y, z, 0) ;
                        IBD_sim_3D->SetBinContent( x, y, z, 0 );
                        IBD_sim_3D->SetBinError( x, y, z, 0) ;
                    } else {
                        // eff = 1.0;
                        RxIBD_3D->SetBinContent( x, y, z, OnContent/eff );
                        RxIBD_3D->SetBinError( x, y, z, OnError/eff ) ;
                        IBD_sim_3D->SetBinContent( x, y, z, SimContent/eff );
                        IBD_sim_3D->SetBinError( x, y, z, SimError/eff ) ;
                    }
                }
            }
        }    
        ////////////////////////////////////////////////////////////////////////////////////////
        
        ////////////////////////////////////////////////////////////////////////////////////////
        // // Make Plots
        // c1->cd();
        // RxIBD_3D->SetTitle("3D IBD");
        // RxIBD_3D->GetXaxis()->SetTitle("X Segment");
        // RxIBD_3D->GetYaxis()->SetTitle("Y Segment");
        // RxIBD_3D->GetZaxis()->SetTitle("Z_pos [mm]");
        // RxIBD_3D->Draw("LEGO2");
        ////////////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////////////
        // // Write Histograms
        // c1->cd();
        TString his_file;
        TString his_name_data;
        TString his_name_sim;
        his_file.Form("histograms.root");
        his_name_data.Form("IBD_3D_%i_zdiv_data",cur_zbin);
        his_name_sim.Form("IBD_3D_%i_zdiv_sim",cur_zbin);
        TH3D* his_data = ( TH3D* ) RxIBD_3D->Clone(his_name_data);
        TH3D* his_sim = ( TH3D* ) IBD_sim_3D->Clone(his_name_sim);
        TFile* file_J = new TFile(his_file,"UPDATE");
        file_J->cd();
        his_data->Write();
        his_sim->Write();
        file_J->Write();
        file_J->Close();
        // RxIBD_3D->SetTitle("3D IBD");
        // RxIBD_3D->GetXaxis()->SetTitle("X Segment");
        // RxIBD_3D->GetYaxis()->SetTitle("Y Segment");
        // RxIBD_3D->GetZaxis()->SetTitle("Z_pos [mm]");
        // RxIBD_3D->Draw("LEGO2");
        ////////////////////////////////////////////////////////////////////////////////////////

        // ////////////////////////////////////////////////////////////////////////////////////////
        // // Here we make the fit for the data
        // //---------------------------------------//
        // //----Normal Plot Data----//
        // //---------------------------------------//
        TF3 *f = new TF3("f","[0]/(TMath::Power((x-[1]),2) + TMath::Power((y-[2]),2) + TMath::Power((z-[3]),2))",x_min,x_max,y_min,y_max,z_min,z_max); // define new fit function
        TFitResultPtr r1 = RxIBD_3D->Fit("f","SR");
        //grab parameters from fit
        double p1 = (-1)*(f->GetParameter(1)); 
        double p2 = (-1)*(f->GetParameter(2));
        double p3 = (f->GetParameter(3));
        double p1_err = (f->GetParError(1)); 
        double p2_err = (f->GetParError(2));
        double p3_err = (f->GetParError(3));
        // append to x, y, z and r arrays 
        x_pos[cur_zbin-1] = (p1+x_const)*length_const;
        y_pos[cur_zbin-1] = (p2+y_const)*length_const;
        z_pos[cur_zbin-1] = (p3)*length_const;
        x_hist->SetBinContent(cur_zbin,(p1+x_const)*length_const);
        y_hist->SetBinContent(cur_zbin,(p2+y_const)*length_const);
        z_hist->SetBinContent(cur_zbin,(p3)*length_const);
        x_hist->SetBinError(cur_zbin,TMath::Sqrt(TMath::Power((p1_err*length_const),2)+TMath::Power((x_bw/2),2)));
        y_hist->SetBinError(cur_zbin,TMath::Sqrt(TMath::Power((p2_err*length_const),2)+TMath::Power((y_bw/2),2)));
        z_hist->SetBinError(cur_zbin,TMath::Sqrt(TMath::Power((p3_err*length_const),2)+TMath::Power((z_bw/2),2)));

        r[cur_zbin-1] = TMath::Sqrt(TMath::Power((x_pos[cur_zbin-1]),2) + TMath::Power((y_pos[cur_zbin-1]),2) + TMath::Power((z_pos[cur_zbin-1]),2)); 
        x_pos_err[cur_zbin-1] = TMath::Sqrt(TMath::Power((p1_err*length_const),2)+TMath::Power((x_bw/2),2));
        y_pos_err[cur_zbin-1] = TMath::Sqrt(TMath::Power((p2_err*length_const),2)+TMath::Power((y_bw/2),2));
        z_pos_err[cur_zbin-1] = TMath::Sqrt(TMath::Power((p3_err*length_const),2)+TMath::Power((z_bw/2),2));
        // now the error propagation for r is a bit more involved
        double x_sq_err = (TMath::Power((x_pos[cur_zbin-1]),2)*(2)*x_pos_err[cur_zbin-1])/x_pos[cur_zbin-1];
        double y_sq_err = (TMath::Power((y_pos[cur_zbin-1]),2)*(2)*y_pos_err[cur_zbin-1])/y_pos[cur_zbin-1];
        double z_sq_err = (TMath::Power((z_pos[cur_zbin-1]),2)*(2)*z_pos_err[cur_zbin-1])/z_pos[cur_zbin-1];
        double r_sq_err = TMath::Sqrt(TMath::Power((x_sq_err),2) + TMath::Power((y_sq_err),2) + TMath::Power((z_sq_err),2));
        r_err[cur_zbin-1] = ((0.5)*(r[cur_zbin-1])*r_sq_err)/(TMath::Power((r[cur_zbin-1]),2));
        r_hist->SetBinContent(cur_zbin,TMath::Sqrt(TMath::Power((x_pos[cur_zbin-1]),2) + TMath::Power((y_pos[cur_zbin-1]),2) + TMath::Power((z_pos[cur_zbin-1]),2)));
        r_hist->SetBinError(cur_zbin,((0.5)*(r[cur_zbin-1])*r_sq_err)/(TMath::Power((r[cur_zbin-1]),2)));

        //////////////////////////////////////////////////////////////////////////////////////
        // // Here we make the fit for the simulation
        // //---------------------------------------//
        // //----Normal Plot Simulation----//
        // //---------------------------------------//
        TF3 *f_sim = new TF3("f_sim","[0]/(TMath::Power((x-[1]),2) + TMath::Power((y-[2]),2) + TMath::Power((z-[3]),2))",x_min,x_max,y_min,y_max,z_min,z_max); // define new fit function
        TFitResultPtr r2 = IBD_sim_3D->Fit("f_sim","SR");
        //grab parameters from fit
        double p1_sim = (-1)*(f_sim->GetParameter(1)); 
        double p2_sim = (-1)*(f_sim->GetParameter(2));
        double p3_sim = (f_sim->GetParameter(3));
        double p1_sim_err = (f_sim->GetParError(1)); 
        double p2_sim_err = (f_sim->GetParError(2));
        double p3_sim_err = (f_sim->GetParError(3));
        // append to x, y, z and r arrays 
        x_sim_pos[cur_zbin-1] = (p1_sim+x_const)*length_const;
        y_sim_pos[cur_zbin-1] = (p2_sim+y_const)*length_const;
        z_sim_pos[cur_zbin-1] = (p3_sim)*length_const;
        x_hist_sim->SetBinContent(cur_zbin,(p1_sim+x_const)*length_const);
        y_hist_sim->SetBinContent(cur_zbin,(p2_sim+y_const)*length_const);
        z_hist_sim->SetBinContent(cur_zbin,(p3_sim)*length_const);
        x_hist_sim->SetBinError(cur_zbin,TMath::Sqrt(TMath::Power((p1_sim_err*length_const),2)+TMath::Power((x_bw/2),2)));
        y_hist_sim->SetBinError(cur_zbin,TMath::Sqrt(TMath::Power((p2_sim_err*length_const),2)+TMath::Power((y_bw/2),2)));
        z_hist_sim->SetBinError(cur_zbin,TMath::Sqrt(TMath::Power((p3_sim_err*length_const),2)+TMath::Power((z_bw/2),2)));
        r_sim[cur_zbin-1] = TMath::Sqrt(TMath::Power((x_sim_pos[cur_zbin-1]),2) + TMath::Power((y_sim_pos[cur_zbin-1]),2) + TMath::Power((z_sim_pos[cur_zbin-1]),2)); 
        x_sim_pos_err[cur_zbin-1] = TMath::Sqrt(TMath::Power((p1_sim_err*length_const),2)+TMath::Power((x_bw/2),2));
        y_sim_pos_err[cur_zbin-1] = TMath::Sqrt(TMath::Power((p2_sim_err*length_const),2)+TMath::Power((y_bw/2),2));
        z_sim_pos_err[cur_zbin-1] = TMath::Sqrt(TMath::Power((p3_sim_err*length_const),2)+TMath::Power((z_bw/2),2));
        std::cout<<"NDivisions:"<<cur_zbin<<std::endl;
        std::cout<<"statistical error x:"<<p1_sim_err*length_const<<" , systematic error x:"<<x_bw/2<<std::endl;
        std::cout<<"statistical error y:"<<p2_sim_err*length_const<<" , systematic error y:"<<y_bw/2<<std::endl;
        std::cout<<"statistical error z:"<<p3_sim_err*length_const<<" , systematic error z:"<<z_bw/2<<std::endl;
        // now the error propagation for r is a bit more involved
        double x_sim_sq_err = (TMath::Power((x_sim_pos[cur_zbin-1]),2)*(2)*x_sim_pos_err[cur_zbin-1])/x_sim_pos[cur_zbin-1];
        double y_sim_sq_err = (TMath::Power((y_sim_pos[cur_zbin-1]),2)*(2)*y_sim_pos_err[cur_zbin-1])/y_sim_pos[cur_zbin-1];
        double z_sim_sq_err = (TMath::Power((z_sim_pos[cur_zbin-1]),2)*(2)*z_sim_pos_err[cur_zbin-1])/z_sim_pos[cur_zbin-1];
        double r_sim_sq_err = TMath::Sqrt(TMath::Power((x_sim_sq_err),2) + TMath::Power((y_sim_sq_err),2) + TMath::Power((z_sim_sq_err),2));
        r_sim_err[cur_zbin-1] = ((0.5)*(r_sim[cur_zbin-1])*r_sim_sq_err)/(TMath::Power((r_sim[cur_zbin-1]),2));
        r_hist_sim->SetBinContent(cur_zbin,TMath::Sqrt(TMath::Power((x_sim_pos[cur_zbin-1]),2) + TMath::Power((y_sim_pos[cur_zbin-1]),2) + TMath::Power((z_sim_pos[cur_zbin-1]),2)));
        r_hist_sim->SetBinError(cur_zbin,((0.5)*(r_sim[cur_zbin-1])*r_sim_sq_err)/(TMath::Power((r_sim[cur_zbin-1]),2)));
        ////////////////////////////////////////////////////////////////////////////////////////////
        

        //---------------------------------------//
        //----Angular Info Plot Data----//
        //---------------------------------------//
        // TF3 *f = new TF3("f","[0]/(TMath::Power((x-[1]*TMath::Sin(-8.59)*TMath::Cos([2])),2) + TMath::Power((y-[1]*TMath::Sin(-8.59)*TMath::Sin([2])),2) + TMath::Power((z-[1]*TMath::Cos(-8.59)),2))",x_min,x_max,y_min,y_max,z_min,z_max); // define new fit function
        // TFitResultPtr r1 = RxIBD_3D->Fit("f","SR");
        //grab parameters from fit
        // double p1 = (-1)*(f->GetParameter(1)); 
        // double p2 = (-1)*(f->GetParameter(2));
        // double p3 = (f->GetParameter(3));
        // double p1_err = (f->GetParError(1)); 
        // double p2_err = (f->GetParError(2));
        // double p3_err = (f->GetParError(3));
        // // append to x, y, z and r arrays 
        // x_pos[cur_zbin-1] = (p1+x_const)*length_const;
        // y_pos[cur_zbin-1] = (p2+y_const)*length_const;
        // z_pos[cur_zbin-1] = (p3)*length_const;
        // r[cur_zbin-1] = TMath::Sqrt(TMath::Power((x_pos[cur_zbin-1]),2) + TMath::Power((y_pos[cur_zbin-1]),2) + TMath::Power((z_pos[cur_zbin-1]),2)); 
        // x_pos_err[cur_zbin-1] = TMath::Sqrt(TMath::Power((p1_err*length_const),2)+TMath::Power((x_bw/2),2));
        // y_pos_err[cur_zbin-1] = TMath::Sqrt(TMath::Power((p2_err*length_const),2)+TMath::Power((y_bw/2),2));
        // z_pos_err[cur_zbin-1] = TMath::Sqrt(TMath::Power((p3_err*length_const),2)+TMath::Power((z_bw/2),2));
        // // now the error propagation for r is a bit more involved
        // double x_sq_err = (TMath::Power((x_pos[cur_zbin-1]),2)*(2)*x_pos_err[cur_zbin-1])/x_pos[cur_zbin-1];
        // double y_sq_err = (TMath::Power((y_pos[cur_zbin-1]),2)*(2)*y_pos_err[cur_zbin-1])/y_pos[cur_zbin-1];
        // double z_sq_err = (TMath::Power((z_pos[cur_zbin-1]),2)*(2)*z_pos_err[cur_zbin-1])/z_pos[cur_zbin-1];
        // double r_sq_err = TMath::Sqrt(TMath::Power((x_sq_err),2) + TMath::Power((y_sq_err),2) + TMath::Power((z_sq_err),2));
        // r_err[cur_zbin-1] = ((0.5)*(r[cur_zbin-1])*r_sq_err)/(TMath::Power((r[cur_zbin-1]),2));
        // //////////////////////////////////////////////////////////////////////////////////////

        //////////////////////////////////////////////////////////////////////////////////////
        // // Here we make the fit for the simulation
        // //---------------------------------------//
        // //----Angular Info Plot Data----//
        // //---------------------------------------//
        // TF3 *f_sim = new TF3("f_sim","[0]/(TMath::Power((x-[1]*TMath::Sin(-8.59)*TMath::Cos([2])),2) + TMath::Power((y-[1]*TMath::Sin(-8.59)*TMath::Sin([2])),2) + TMath::Power((z-[1]*TMath::Cos(-8.59)),2))",x_min,x_max,y_min,y_max,z_min,z_max); // define new fit function
        // TFitResultPtr r2 = IBD_sim_3D->Fit("f_sim","SR");
        // //grab parameters from fit
        // double p1_sim = (-1)*(f_sim->GetParameter(1)); 
        // double p2_sim = (-1)*(f_sim->GetParameter(2));
        // double p3_sim = (f_sim->GetParameter(3));
        // double p1_sim_err = (f_sim->GetParError(1)); 
        // double p2_sim_err = (f_sim->GetParError(2));
        // double p3_sim_err = (f_sim->GetParError(3));
        // // append to x, y, z and r arrays 
        // x_sim_pos[cur_zbin-1] = (p1_sim+x_const)*length_const;
        // y_sim_pos[cur_zbin-1] = (p2_sim+y_const)*length_const;
        // z_sim_pos[cur_zbin-1] = (p3_sim)*length_const;
        // r_sim[cur_zbin-1] = TMath::Sqrt(TMath::Power((x_sim_pos[cur_zbin-1]),2) + TMath::Power((y_sim_pos[cur_zbin-1]),2) + TMath::Power((z_sim_pos[cur_zbin-1]),2)); 
        // x_sim_pos_err[cur_zbin-1] = TMath::Sqrt(TMath::Power((p1_sim_err*length_const),2)+TMath::Power((x_bw/2),2));
        // y_sim_pos_err[cur_zbin-1] = TMath::Sqrt(TMath::Power((p2_sim_err*length_const),2)+TMath::Power((y_bw/2),2));
        // z_sim_pos_err[cur_zbin-1] = TMath::Sqrt(TMath::Power((p3_sim_err*length_const),2)+TMath::Power((z_bw/2),2));
        // std::cout<<"NDivisions:"<<cur_zbin<<std::endl;
        // std::cout<<"statistical error x:"<<p1_sim_err*length_const<<" , systematic error x:"<<x_bw/2<<std::endl;
        // std::cout<<"statistical error y:"<<p2_sim_err*length_const<<" , systematic error y:"<<y_bw/2<<std::endl;
        // std::cout<<"statistical error z:"<<p3_sim_err*length_const<<" , systematic error z:"<<z_bw/2<<std::endl;
        // // now the error propagation for r is a bit more involved
        // double x_sim_sq_err = (TMath::Power((x_sim_pos[cur_zbin-1]),2)*(2)*x_sim_pos_err[cur_zbin-1])/x_sim_pos[cur_zbin-1];
        // double y_sim_sq_err = (TMath::Power((y_sim_pos[cur_zbin-1]),2)*(2)*y_sim_pos_err[cur_zbin-1])/y_sim_pos[cur_zbin-1];
        // double z_sim_sq_err = (TMath::Power((z_sim_pos[cur_zbin-1]),2)*(2)*z_sim_pos_err[cur_zbin-1])/z_sim_pos[cur_zbin-1];
        // double r_sim_sq_err = TMath::Sqrt(TMath::Power((x_sim_sq_err),2) + TMath::Power((y_sim_sq_err),2) + TMath::Power((z_sim_sq_err),2));
        // r_sim_err[cur_zbin-1] = ((0.5)*(r_sim[cur_zbin-1])*r_sim_sq_err)/(TMath::Power((r_sim[cur_zbin-1]),2));
        ////////////////////////////////////////////////////////////////////////////////////////////


        ////////////////////////////////////////////////////////////////////////////////////////
        // //Write histograms to output root file
        // TString output_name;
        // output_name.Form("nu_dir_study.root");
        // TFile* file1 = new TFile(output_name,"RECREATE");
        // file1->cd();
        // RxIBD[0]->Write("IBDCell_cut1");
        // RxIBD[1]->Write("IBDCell_cut2");
        // RxIBD[2]->Write("IBDCell_cut3");
        // RxIBD[3]->Write("IBDCell_cut4");
        // RxIBD_All->Write("IBDCell_allZ");
        // RxIBD_3D->Write("IBDCell_3D");
        // z_IBD->Write("z_pos_IBD");
        // file1->Write();
        // file1->Close();
        // delete file1;
        // // // // c1->SaveAs("PRD_Plots/IBD.pdf");
        ////////////////////////////////////////////////////////////////////////////////////////
        
        ////////////////////////////////////////////////////////////////////////////////////////
        //get rid of all pointers
        // delete f;
        // delete f2;
        delete RxOn_3D;
        delete RxOff_3D;
        delete RxIBD_3D;
        delete IBD_sim_3D;
        delete seg_map;
        // delete base_file;
        // delete ttree;
        // delete efficiency;
        // delete seg_baseline;
        // delete cell_map; 
        
        ////////////////////////////////////////////////////////////////////////////////////////
    }
    // Now that we are done with histograms and fitting we can create the graphs
    //---------------------------------------//
    //----create graph for x----//
    //---------------------------------------//
    c1->cd();
    x_hist->SetMarkerStyle(kFullCircle);
    x_hist->SetMarkerColor(kBlack);
    x_hist->SetMarkerSize(.7);
    x_hist->SetLineColor(kBlack);
    x_hist_sim->SetMarkerStyle(kFullCircle);
    x_hist_sim->SetMarkerColor(kBlue);
    x_hist_sim->SetMarkerSize(.7);
    x_hist_sim->SetLineColor(kBlue);
    x_hist->GetYaxis()->SetTitle("x_pos [m]");
    x_hist->GetXaxis()->SetTitle("Number of cell divisions");
    x_hist->Draw("E0");
    x_hist_sim->Draw("same E0");

    // TGraphErrors *x_gr = new TGraphErrors(z_bins,nbinsz,x_pos,nbinsz_err,x_pos_err);
    // x_gr->SetName("x_gr");
    // x_gr->SetTitle("Data ");
    // x_gr->SetMarkerColor(4);
    // x_gr->SetLineColor(4);
    // x_gr->SetMarkerStyle(21);
    // x_gr->GetXaxis()->SetTitle("Z bins");
    // x_gr->GetYaxis()->SetTitle("x_pos [m]");
    // // x_gr->Draw("AP");
    // TGraphErrors *x_sim_gr = new TGraphErrors(z_bins,nbinsz,x_sim_pos,nbinsz_err,x_sim_pos_err);
    // x_sim_gr->SetName("x_sim_gr");
    // x_sim_gr->SetTitle("Simulation");
    // x_sim_gr->SetMarkerColor(3);
    // x_sim_gr->SetLineColor(3);
    // x_sim_gr->SetMarkerStyle(20);
    // // x_sim_gr->Draw("AP same");
    // auto mg_x = new TMultiGraph("mg_x","mg_x");
    // mg_x->SetTitle("x pos for different z binnings;Z bins; x_pos [m]");
    // mg_x->Add(x_gr);
    // mg_x->Add(x_sim_gr);
    // mg_x->Draw("apl");
    // c1->BuildLegend();
    TLine *l_x = new TLine(nbinsz[0],x_pos_true,nbinsz[z_bins-1],x_pos_true);
    l_x->SetLineStyle(kSolid);
    l_x->SetLineColor(kRed);
    l_x->SetLineWidth(4);
    l_x->Draw("same");
    
    //---------------------------------------//
    //----create graph for y----//
    //---------------------------------------//
    c2->cd();
    y_hist->SetMarkerStyle(kFullCircle);
    y_hist->SetMarkerColor(kBlack);
    y_hist->SetMarkerSize(.7);
    y_hist->SetLineColor(kBlack);
    y_hist_sim->SetMarkerStyle(kFullCircle);
    y_hist_sim->SetMarkerColor(kBlue);
    y_hist_sim->SetMarkerSize(.7);
    y_hist_sim->SetLineColor(kBlue);
    y_hist->GetYaxis()->SetTitle("y_pos [m]");
    y_hist->GetXaxis()->SetTitle("Number of cell divisions");
    y_hist->Draw("E0");
    y_hist_sim->Draw("same E0");

    // TGraphErrors *y_gr = new TGraphErrors(z_bins,nbinsz,y_pos,nbinsz_err,y_pos_err);
    // y_gr->SetName("y_gr");
    // y_gr->SetTitle("Data");
    // y_gr->SetMarkerColor(4);
    // y_gr->SetLineColor(4);
    // y_gr->SetMarkerStyle(21);
    // TGraphErrors *y_sim_gr = new TGraphErrors(z_bins,nbinsz,y_sim_pos,nbinsz_err,y_sim_pos_err);
    // y_sim_gr->SetName("y_sim_gr");
    // y_sim_gr->SetTitle("Simulation");
    // y_sim_gr->SetMarkerColor(3);
    // y_sim_gr->SetLineColor(3);
    // y_sim_gr->SetMarkerStyle(20);
    // auto mg_y = new TMultiGraph("mg_y","mg_y");
    // mg_y->SetTitle("y pos for different z binnings;Z bins; y_pos [m]");
    // mg_y->Add(y_gr);
    // mg_y->Add(y_sim_gr);
    // mg_y->Draw("apl");
    // c2->BuildLegend();
    TLine *l_y = new TLine(nbinsz[0],y_pos_true,nbinsz[z_bins-1],y_pos_true);
    l_y->SetLineStyle(kSolid);
    l_y->SetLineColor(kRed);
    l_y->SetLineWidth(4);
    l_y->Draw();
    // y_gr->GetXaxis()->SetTitle("Z bins");
    // y_gr->GetYaxis()->SetTitle("y_pos [m]");
    //---------------------------------------//
    //----create graph for z----//
    //---------------------------------------//
    c3->cd();
    z_hist->SetMarkerStyle(kFullCircle);
    z_hist->SetMarkerColor(kBlack);
    z_hist->SetMarkerSize(.7);
    z_hist->SetLineColor(kBlack);
    z_hist_sim->SetMarkerStyle(kFullCircle);
    z_hist_sim->SetMarkerColor(kBlue);
    z_hist_sim->SetMarkerSize(.7);
    z_hist_sim->SetLineColor(kBlue);
    z_hist->GetYaxis()->SetTitle("z_pos [m]");
    z_hist->GetXaxis()->SetTitle("Number of cell divisions");
    z_hist->Draw("E0");
    z_hist_sim->Draw("same E0");
    // TGraphErrors *z_gr = new TGraphErrors(z_bins,nbinsz,z_pos,nbinsz_err,z_pos_err);
    // z_gr->SetName("z_gr");
    // z_gr->SetTitle("Data");
    // z_gr->SetMarkerColor(4);
    // z_gr->SetLineColor(4);
    // z_gr->SetMarkerStyle(21);
    // // z_gr->Draw("AP");
    // TGraphErrors *z_sim_gr = new TGraphErrors(z_bins,nbinsz,z_sim_pos,nbinsz_err,z_sim_pos_err);
    // z_sim_gr->SetName("z_sim_gr");
    // z_sim_gr->SetTitle("Simulation ");
    // z_sim_gr->SetMarkerColor(3);
    // z_sim_gr->SetLineColor(3);
    // z_sim_gr->SetMarkerStyle(20);
    // // z_sim_gr->Draw("AP same");
    // auto mg_z = new TMultiGraph("mg_z","mg_z");
    // mg_z->SetTitle("z pos for different z binnings;Z bins; z_pos [m]");
    // mg_z->Add(z_gr);
    // mg_z->Add(z_sim_gr);
    // mg_z->Draw("apl");
    // c3->BuildLegend();
    TLine *l_z = new TLine(nbinsz[0],z_pos_true,nbinsz[z_bins-1],z_pos_true);
    l_z->SetLineStyle(kSolid);
    l_z->SetLineColor(kRed);
    l_z->SetLineWidth(4);
    l_z->Draw("same");
    // z_gr->GetXaxis()->SetTitle("Z bins");
    // z_gr->GetYaxis()->SetTitle("z_pos [m]");
    //---------------------------------------//
    //----create graph for r----//
    //---------------------------------------//
    c4->cd();
    r_hist->SetMarkerStyle(kFullCircle);
    r_hist->SetMarkerColor(kBlack);
    r_hist->SetMarkerSize(.7);
    r_hist->SetLineColor(kBlack);
    r_hist_sim->SetMarkerStyle(kFullCircle);
    r_hist_sim->SetMarkerColor(kBlue);
    r_hist_sim->SetMarkerSize(.7);
    r_hist_sim->SetLineColor(kBlue);
    r_hist->GetYaxis()->SetTitle("r_pos [m]");
    r_hist->GetXaxis()->SetTitle("Number of cell divisions");
    r_hist->Draw("E0");
    r_hist_sim->Draw("same E0");
    // TGraphErrors *r_gr = new TGraphErrors(z_bins,nbinsz,r,nbinsz_err,r_err);
    // r_gr->SetName("r_gr");
    // r_gr->SetTitle("Data");
    // r_gr->SetMarkerColor(4);
    // r_gr->SetLineColor(4);
    // r_gr->SetMarkerStyle(21);
    // // r_gr->Draw("AP");
    // TGraphErrors *r_sim_gr = new TGraphErrors(z_bins,nbinsz,r_sim,nbinsz_err,r_sim_err);
    // r_sim_gr->SetName("r_sim_gr");
    // r_sim_gr->SetTitle("Simulation");
    // r_sim_gr->SetMarkerColor(3);
    // r_sim_gr->SetLineColor(3);
    // r_sim_gr->SetMarkerStyle(20);
    // // r_sim_gr->Draw("AP same");
    // auto mg_r = new TMultiGraph("mg_r","mg_r");
    // mg_r->SetTitle("r pos for different z binnings;Z bins; r_pos [m]");
    // mg_r->Add(r_gr);
    // mg_r->Add(r_sim_gr);
    // mg_r->Draw("apl");
    // c4->BuildLegend();
    TLine *l_r = new TLine(nbinsz[0],r_true,nbinsz[z_bins-1],r_true);
    l_r->SetLineStyle(kSolid);
    l_r->SetLineColor(kRed);
    l_r->SetLineWidth(3);
    l_r->Draw("same");
    // r_gr->GetXaxis()->SetTitle("Z bins");
    // r_gr->GetYaxis()->SetTitle("r_pos [m]");
    //---------------------------------------//
    //----Write file----//
    //---------------------------------------//
    TString output_name;
    // output_name.Form("graps.root");
    // TFile* file = new TFile(output_name,"RECREATE");
    // file->cd();
    // x_gr->Write();
    // y_gr->Write();
    // z_gr->Write();
    // r_gr->Write();
    // x_sim_gr->Write();
    // y_sim_gr->Write();
    // z_sim_gr->Write();
    // r_sim_gr->Write();
    // file->Write();
    // file->Close();
    output_name.Form("pos_histograms.root");
    TFile* file2 = new TFile(output_name,"UPDATE");
    file2->cd();
    x_hist->Write();
    y_hist->Write();
    z_hist->Write();
    r_hist->Write();
    x_hist_sim->Write();
    y_hist_sim->Write();
    z_hist_sim->Write();
    r_hist_sim->Write();
    file2->Write();
    file2->Close();

}
