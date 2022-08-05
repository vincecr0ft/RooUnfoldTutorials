// Pim Jordi Verschuuren
// Writen December 2020

#include <time.h>

Double_t Smear(Double_t xtrue, Double_t NP3, Double_t NP2){

  Double_t xeff= NP3 - abs(xtrue)/600;
  Double_t x= gRandom->Rndm();
  if (x>xeff) return -99999;
  Double_t xsmear= gRandom->Gaus(0,NP2*sqrt(xtrue/30));     // bias and smear
  return xtrue+xsmear;
}


void Generate_exp(){

  // Variable binning.
  Double_t binning[13] = {0,2,4,6,8,10,12,14,18,25,35,45,60};

  Int_t n_bins = 12;

  // Nuisance parameter 1 affecting a scale factor.
  Double_t NP1_vals[3] = {0.95,1.,1.05};

  // Nuisance parameter 2 affecting the smearing.
  Double_t NP2_vals[3] = {3,4,5};

  // Nuisance parameter 3 affecting the efficiency.
  Double_t NP3_vals[3] = {0.94,0.95,0.96};

  // Nuisance parameter 4 affecting the bkg curve.
  Double_t NP4_vals[3] = {-0.0895,-0.09,-0.0905};

  // Nuisance parameter 5 affecting the sig curve.
  Double_t NP5_vals[3] = {-0.099,-0.1,-0.11};

  //Double_t signal_exp = -0.14;

  Double_t total_scale = 0.01;

  // 7 different nuisance parameter sets.
  // 1 nominal, 4 up and 4 down variations.
  TMatrixD NP_vals(11,5);
  std::vector<std::string> names;

  // Nominal
  NP_vals(0,0) = NP1_vals[1];
  NP_vals(0,1) = NP2_vals[1];
  NP_vals(0,2) = NP3_vals[1];
  NP_vals(0,3) = NP4_vals[1];
  NP_vals(0,4) = NP5_vals[1];
  names.push_back("nom");

  // NP1
  NP_vals(1,0) = NP1_vals[0];
  NP_vals(1,1) = NP2_vals[1];
  NP_vals(1,2) = NP3_vals[1];
  NP_vals(1,3) = NP4_vals[1];
  NP_vals(1,4) = NP5_vals[1];
  NP_vals(2,0) = NP1_vals[2];
  NP_vals(2,1) = NP2_vals[1];
  NP_vals(2,2) = NP3_vals[1];
  NP_vals(2,3) = NP4_vals[1];
  NP_vals(2,4) = NP5_vals[1];
  names.push_back("NP1up");
  names.push_back("NP1down");
  

  // NP2
  NP_vals(3,0) = NP1_vals[1];
  NP_vals(3,1) = NP2_vals[0];
  NP_vals(3,2) = NP3_vals[1];
  NP_vals(3,3) = NP4_vals[1];
  NP_vals(3,4) = NP5_vals[1];
  NP_vals(4,0) = NP1_vals[1];
  NP_vals(4,1) = NP2_vals[2];
  NP_vals(4,2) = NP3_vals[1];
  NP_vals(4,3) = NP4_vals[1];
  NP_vals(4,4) = NP5_vals[1];
  names.push_back("NP2up");
  names.push_back("NP2down");


  // NP3
  NP_vals(5,0) = NP1_vals[1];
  NP_vals(5,1) = NP2_vals[1];
  NP_vals(5,2) = NP3_vals[0];
  NP_vals(5,3) = NP4_vals[1];
  NP_vals(5,4) = NP5_vals[1];
  NP_vals(6,0) = NP1_vals[1];
  NP_vals(6,1) = NP2_vals[1];
  NP_vals(6,2) = NP3_vals[2];
  NP_vals(6,3) = NP4_vals[1];
  NP_vals(6,4) = NP5_vals[1];
  names.push_back("NP3up");
  names.push_back("NP3down");


  // NP4
  NP_vals(7,0) = NP1_vals[1];
  NP_vals(7,1) = NP2_vals[1];
  NP_vals(7,2) = NP3_vals[1];
  NP_vals(7,3) = NP4_vals[0];
  NP_vals(7,4) = NP5_vals[1];
  NP_vals(8,0) = NP1_vals[1];
  NP_vals(8,1) = NP2_vals[1];
  NP_vals(8,2) = NP3_vals[1];
  NP_vals(8,3) = NP4_vals[2];
  NP_vals(8,4) = NP5_vals[1];
  names.push_back("NP4up");
  names.push_back("NP4down");

  // NP4
  NP_vals(9,0) = NP1_vals[1];
  NP_vals(9,1) = NP2_vals[1];
  NP_vals(9,2) = NP3_vals[1];
  NP_vals(9,3) = NP4_vals[1];
  NP_vals(9,4) = NP5_vals[0];
  NP_vals(10,0) = NP1_vals[1];
  NP_vals(10,1) = NP2_vals[1];
  NP_vals(10,2) = NP3_vals[1];
  NP_vals(10,3) = NP4_vals[1];
  NP_vals(10,4) = NP5_vals[2];
  names.push_back("NP5up");
  names.push_back("NP5down");

  Int_t n_MC_sig_events = 10000000;
  Int_t n_MC_bkg_events = 10000000;

  Int_t n_data_sig_events = 100000;
  Int_t n_data_bkg_events = 300000;

  Double_t sig_lumiScale = (double)n_data_sig_events/n_MC_sig_events;
  Double_t bkg_lumiScale = (double)n_data_bkg_events/n_MC_bkg_events;

  RooRealVar xtrue_sig("xtrue_sig","xtrue_sig",0,60);
  RooRealVar xtrue_bkg("xtrue_bkg","xtrue_bkg",0,60);

  //RooRandom::randomGenerator()->SetSeed(time(NULL));

  TFile output("input_file.root","RECREATE");
      
  for (int i = 0; i < names.size(); i++){

    Double_t NP1_val = NP_vals(i,0);
    Double_t NP2_val = NP_vals(i,1);
    Double_t NP3_val = NP_vals(i,2);
    Double_t NP4_val = NP_vals(i,3);
    Double_t NP5_val = NP_vals(i,4);

    std::string name = names.at(i);

    std::string truth_name("truth");
    std::string reco_name("reco");
    std::string reco_bkg_name("reco_bkg");
    std::string response_name("response");

    truth_name.append(name);
    reco_name.append(name);
    reco_bkg_name.append(name);
    response_name.append(name);

    RooRealVar lambda_sig("lambda_sig","lambda_sig",NP5_val);
    RooRealVar lambda_bkg("lambda_bkg","lambda_bkg",NP4_val);

    RooExponential sigPDF("sigPDF","sigPDF", xtrue_sig, lambda_sig);
    RooExponential bkgPDF("bkgPDF","bkgPDF", xtrue_bkg, lambda_bkg);

    RooDataSet* MC_sig_Events = sigPDF.generate(xtrue_sig,n_MC_sig_events);
    RooDataSet* MC_sig_test_Events = sigPDF.generate(xtrue_sig,n_MC_sig_events);
    RooDataSet* data_sig_Events = sigPDF.generate(xtrue_sig,n_data_sig_events);

    RooDataSet* MC_bkg_Events = bkgPDF.generate(xtrue_bkg,n_MC_bkg_events);
    RooDataSet* data_bkg_Events = bkgPDF.generate(xtrue_bkg,n_data_bkg_events);

    TH1D truthHist(truth_name.c_str(),truth_name.c_str(),n_bins,binning);
    TH1D recoHist(reco_name.c_str(),reco_name.c_str(),n_bins,binning);
    TH1D recoBkgHist(reco_bkg_name.c_str(),reco_bkg_name.c_str(),n_bins,binning);
    TH2D responseHist(response_name.c_str(),response_name.c_str(),n_bins, binning, n_bins, binning);  
    
    TH1D truthTestHist("truthnom_test","truthnom_test",n_bins,binning);
    TH1D dataHist("reconom_test","reconom_test",n_bins,binning);

    for (int j = 0; j < n_MC_sig_events; j++){
      
      Double_t xtrue_sig_val = ((RooRealVar*)(MC_sig_Events->get(j))->find(xtrue_sig.GetName()))->getVal();
      
      Double_t xreco_sig_val = Smear(xtrue_sig_val, NP3_val, NP2_val);
      
      truthHist.Fill(xtrue_sig_val, sig_lumiScale*total_scale);
      
      if (xreco_sig_val != -99999){
	recoHist.Fill(xreco_sig_val, sig_lumiScale*NP1_val*total_scale);
	responseHist.Fill(xreco_sig_val, xtrue_sig_val, sig_lumiScale*NP1_val*total_scale);
      }
    }

    for (int j = 0; j < n_MC_bkg_events; j++){
      
      Double_t xtrue_bkg_val = ((RooRealVar*)(MC_bkg_Events->get(j))->find(xtrue_bkg.GetName()))->getVal();
      
      Double_t xreco_bkg_val = Smear(xtrue_bkg_val, NP3_val, NP2_val);
      
      if (xreco_bkg_val != -99999){
	recoBkgHist.Fill(xreco_bkg_val, bkg_lumiScale*NP1_val*total_scale);
      }
    }

    recoHist.Write();
    recoBkgHist.Write();
    truthHist.Write();
    responseHist.Write();

    if (i > 0) continue;

    for (int j = 0; j < n_data_sig_events; j++){
      
      Double_t xtrue_val = ((RooRealVar*)(data_sig_Events->get(j))->find(xtrue_sig.GetName()))->getVal();
      
      Double_t xreco_val = Smear(xtrue_val, NP3_val, NP2_val);
      
      if (xreco_val != -99999){
	dataHist.Fill(xreco_val,total_scale);
      }    
    }

    for (int j = 0; j < n_data_bkg_events; j++){
      
      Double_t xtrue_val = ((RooRealVar*)(data_bkg_Events->get(j))->find(xtrue_bkg.GetName()))->getVal();
      
      Double_t xreco_val = Smear(xtrue_val, NP3_val, NP2_val);

      if (xreco_val != -99999){
	dataHist.Fill(xreco_val,total_scale);
      }    
    }

    for (int j = 0; j < n_MC_sig_events; j++){
      
      Double_t xtrue_sig_val = ((RooRealVar*)(MC_sig_test_Events->get(j))->find(xtrue_sig.GetName()))->getVal();
      
      truthTestHist.Fill(xtrue_sig_val, sig_lumiScale*total_scale);
      
    }

    truthTestHist.Write();
    dataHist.Write();

    delete MC_sig_Events;
    delete data_sig_Events;
    delete MC_bkg_Events;
    delete data_bkg_Events;
  }
  
  output.Close();

  
}
