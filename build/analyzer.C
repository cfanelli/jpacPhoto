#include <stdio.h>
#include <time.h>
#include "TGraphErrors.h"
#include <TH2F.h>
#include <TH1F.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TH3F.h>
#include <TF1.h>
#include <TChain.h>
#include <TString.h>
#include <vector>
#include <string>
#include "TROOT.h"
#include "TFrame.h"
#include "TLatex.h"
#include "TFile.h"
#include "TMath.h"
#include "lhcbStyle.h"
#include "TGenPhaseSpace.h"
#include "TRandom.h"
#include <math.h>
#include "TStyle.h"

#include "analyzer.h"

#define Mp 0.938


//----------------------------------------------------------------------------//

vector <double> calc_pol(double , double , TGraphErrors*, TGraph*, TH1F*);
int calc_bin(double val, double min, double width);


vector <double> calc_Asy(double Yperp, double Ypara);

//----------------------------------------------------------------------------//




//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



void analyzer(double Emin=-99., double Emax=99., double tmin=-99., double tmax=99., string big_label="t_channel"){

  cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<endl;
  cout<<"Emin: "<< Emin << ", Emax: " << Emax << ", tmin: " << tmin << ", tmax: "<< tmax << endl;
  cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<endl;


  char filename[256];
  if(big_label == "t_channel") sprintf(filename, "./outputfiles/5q_mc_t-chan__HorM_large.root");
  if(big_label == "s_channel") sprintf(filename, "./outputfiles/5q_mc_Pc__HorM_large.root");

  //--------------------------------------------------------------------------//

  TH1F* h_beamE_jpsi_para = new TH1F("h_beamE_jpsi_para",";beam E [GeV];", 3200, 8.3, 11.5);
  TH1F* h_beamE_jpsi_perp = new TH1F("h_beamE_jpsi_perp",";beam E [GeV];", 3200, 8.3, 11.5);

  TH1F* h_phi_para = new TH1F("h_phi_para",";#phi_{p}(rad);", 24, -180, 180); //-M_PI, M_PI
  TH1F* h_phi_perp = new TH1F("h_phi_perp",";#phi_{p}(rad);", 24, -180, 180);

  TH1F* h_M_para = new TH1F("h_M_para",";M_{e+e-} [GeV/c^{2}];", 500, 2.5, 3.5);
  TH1F* h_M_perp = new TH1F("h_M_perp",";M_{e+e-} [GeV/c^{2}];", 500, 2.5, 3.5);

  TH1F* h_Sigma = new TH1F("h_Sigma",";#phi [deg];#Sigma;", 24, -180, 180);

  TH1F* h_t = new TH1F("h_t","Mandelstam t in bin of energy;t [GeV^{2}];;", 500, 0, 5);


  //--------------------------------------------------------------------------//

  const char* rootfile_P = "./inputfiles/beamP_v3.root";// = "./inputfiles/beamP.root";

  TFile*polfile = TFile::Open(rootfile_P);
  TGraphErrors*pol = (TGraphErrors*) polfile->Get("pol");
  TGraph*uppol = (TGraph*) polfile->Get("uppol");

  const char* treename = "kinematics";

  double beam_E, beam_P, beam_h;
  double ep_px, ep_py, ep_pz, ep_E;
  double em_px, em_py, em_pz, em_E;
  double prec_px, prec_py, prec_pz, prec_E;
  double phi_psi, theta_psi, phi_ee, theta_ee;


  TFile* origin = new TFile(filename);
  TTree*tree;
  origin->GetObject(treename,tree);

  //--------------------------------------------------------------------------//


  // beam
  tree->SetBranchAddress("beam_E", &beam_E);
  tree->SetBranchAddress("beam_P", &beam_P);
  tree->SetBranchAddress("beam_h", &beam_h);
  // ep
  tree->SetBranchAddress("ep_px", &ep_px);
  tree->SetBranchAddress("ep_py", &ep_py);
  tree->SetBranchAddress("ep_pz", &ep_pz);
  tree->SetBranchAddress("ep_E", &ep_E);
  // em
  tree->SetBranchAddress("em_px", &em_px);
  tree->SetBranchAddress("em_py", &em_py);
  tree->SetBranchAddress("em_pz", &em_pz);
  tree->SetBranchAddress("em_E", &em_E);
  // prec
  tree->SetBranchAddress("prec_px", &prec_px);
  tree->SetBranchAddress("prec_py", &prec_py);
  tree->SetBranchAddress("prec_pz", &prec_pz);
  tree->SetBranchAddress("prec_E", &prec_E);
  // angles
  tree->SetBranchAddress("phi_psi", &phi_psi);
  tree->SetBranchAddress("theta_psi", &theta_psi);
  tree->SetBranchAddress("phi_ee", &phi_ee);
  tree->SetBranchAddress("theta_ee", &theta_ee);

  //TLorentzVector ep_P4, em_P4, prec_P4, jpsi_P4, ptar_P4, beam_P4, Mand_t, Mand_s;

  event phys;

  phys.ptar_P4.SetPxPyPzE(0.,0.,0.,Mp);

  phys.Emin = Emin;
  phys.Emax = Emax;
  phys.tmin = tmin;
  phys.tmax = tmax;


  double t_mand, s_mand;

  bool accept = false;

  vector <double> tmp_pol;

  TRandom3 r;
  double ran;

  //--------------------------------------------------------------------------//
  //                            (1)  TREE LOOP
  //--------------------------------------------------------------------------//
  for(int i=0;i<tree->GetEntries();i++)
  {
      tree->GetEntry(i);
      accept = false;

      phys.ep_P4.SetPxPyPzE(ep_px,ep_py,ep_pz,ep_E);
      phys.em_P4.SetPxPyPzE(em_px,em_py,em_pz,em_E);
      phys.prec_P4.SetPxPyPzE(prec_px,prec_py,prec_pz,prec_E);
      phys.beam_P4.SetPxPyPzE(0.,0.,beam_E,beam_E);

      phys.jpsi_P4 = phys.ep_P4 + phys.em_P4;

      phys.Mand_t = phys.ptar_P4 - phys.prec_P4;
      phys.Mand_s = phys.beam_P4 + phys.ptar_P4;

      phys.t_mand = phys.Mand_t * phys.Mand_t;
      phys.s_mand = phys.Mand_s * phys.Mand_s;


      //cout<<"t: "<< phys.t_mand <<", s: "<< phys.s_mand <<endl;


      if(phys.beam_P4.E()>Emin && phys.beam_P4.E()<Emax) h_t->Fill(abs(phys.t_mand));


      //---------------------------//
      ran = r.Uniform(0.,1.);
      accept = phys.selection() && phys.detector_acc(ran);

      if(accept==false) continue;

      //cout<< "selected energy: "<< phys.beam_P4.E() <<endl;
      //---------------------------//

      // only accepted events from here


      // two polarization states
      if(beam_h== 1.){
         h_beamE_jpsi_para->Fill(phys.beam_P4.E());
         h_phi_para->Fill(phys.prec_P4.Phi()*TMath::RadToDeg());
         h_M_para->Fill(phys.jpsi_P4.M());
       }
      if(beam_h==-1.){
        h_beamE_jpsi_perp->Fill(phys.beam_P4.E());
        h_phi_perp->Fill(phys.prec_P4.Phi()*TMath::RadToDeg());
        h_M_perp->Fill(phys.jpsi_P4.M());
      }


  }
  //--------------------------------------------------------------------------//


  //--------------------------------------------------------------------------//
  //                        (2) ANALYZE HISTOGRAMS
  //--------------------------------------------------------------------------//


  //calculate polarization
  tmp_pol.clear();
  tmp_pol = calc_pol(Emin,Emax,pol,uppol,h_beamE_jpsi_para);


  //lhcbStyle();

  //TCanvas* c = new TCanvas();

  //h-para->SetLi
  //h_phi_para->Draw("EP");
  //h_phi_perp->Draw("EP same");


  char out_file[256];
  sprintf(out_file,"./outputfiles/out_ana_Emin_%1.3f_Emax_%1.3f_tmin_%1.3f_tmax_%1.3f.root",Emin, Emax, tmin, tmax);




  //--------------------------------------------------------------------------//


  //--------------------------------------------------------------------------//
  //                       (3)  BEAM ASYMMETRY SURVEY
  //--------------------------------------------------------------------------//

  vector <double> Asy;

  double tmp_Ypara, tmp_Yperp;
  for(int i=0;i<h_phi_para->GetNbinsX();i++){
    Asy.clear();

    tmp_Ypara = h_phi_para->GetBinContent(i);
    tmp_Yperp = h_phi_perp->GetBinContent(i);

    Asy = calc_Asy(tmp_Yperp,tmp_Ypara);

    //cout<<"phi: "<< h_phi_para->GetBinCenter(i) << ", A: "<< Asy.at(0) << ", sigma(A): " << Asy.at(1) << endl;

    h_Sigma->SetBinContent(i,Asy.at(0));
    h_Sigma->SetBinError(i,Asy.at(1));

  }


  //--------------------------------------------------------------------------//
  //                            CUSTOM PLOTS
  //--------------------------------------------------------------------------//

  char Emin_label[256], Emax_label[256], tmin_label[256], tmax_label[256];

  sprintf(Emin_label,"%1.3f",Emin);
  sprintf(Emax_label,"%1.3f",Emax);
  sprintf(tmin_label,"%1.3f",tmin);
  sprintf(tmax_label,"%1.3f",tmax);

  if(Emin==-99.) sprintf(Emin_label,"8.3");
  if(Emax==99.) sprintf(Emax_label,"12");

  if(tmin==-99.) sprintf(tmin_label,"0.");
  if(tmax==99.) sprintf(tmax_label,"+#infty");


  TCanvas*c = new TCanvas();
  c->SetTitle(big_label.c_str());


  lhcbStyle();

  c->Divide(1,3);


  c->cd(1);

  gROOT->SetStyle("Plain");
  lhcbStyle();
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.10);

  h_phi_para->SetLineWidth(2);
  h_phi_para->SetLineColor(kBlue);
  h_phi_para->SetMarkerColor(kBlue);
  h_phi_para->GetYaxis()->SetRangeUser(0.,h_phi_para->GetMaximum()*1.5);
  h_phi_para->GetXaxis()->SetTitle("#phi_{p} [deg]");
  h_phi_para->GetXaxis()->GetCenterTitle();
  h_phi_para->GetXaxis()->SetTitleSize(0.07);
  h_phi_para->GetYaxis()->SetTitle("entries");
  h_phi_para->GetYaxis()->SetTitleSize(0.07);
  h_phi_para->GetYaxis()->GetCenterTitle();

  h_phi_para->Draw("EP");

  TLatex *chlabel1 = new TLatex();
  chlabel1-> SetNDC();
  chlabel1 -> SetTextFont(62);
  chlabel1 -> SetTextColor(1);
  chlabel1 -> SetTextSize(0.045);
  chlabel1 -> SetTextAlign(22);
  chlabel1 -> SetTextAngle(0);
  char put_text1[256];
  sprintf(put_text1,"PARA, E_{#gamma}#in(%s,%s)",Emin_label, Emax_label);
  chlabel1 -> DrawLatex(0.7, 0.75, put_text1);
  sprintf(put_text1,"-t #in(%s,%s)",tmin_label, tmax_label);
  chlabel1 -> DrawLatex(0.7, 0.65, put_text1);

  //TCanvas*c2 = new TCanvas();

  c->cd(2);

  gROOT->SetStyle("Plain");
  lhcbStyle();
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.10);


  h_phi_perp->SetLineWidth(2);
  h_phi_perp->SetLineColor(kBlue);
  h_phi_perp->SetMarkerColor(kBlue);
  h_phi_perp->GetYaxis()->SetRangeUser(0.,h_phi_para->GetMaximum());
  h_phi_perp->GetXaxis()->SetTitle("#phi_{p} [deg]");
  h_phi_perp->GetXaxis()->GetCenterTitle();
  h_phi_perp->GetXaxis()->SetTitleSize(0.07);
  h_phi_perp->GetYaxis()->SetTitle("entries");
  h_phi_perp->GetYaxis()->SetTitleSize(0.07);
  h_phi_perp->GetYaxis()->GetCenterTitle();
  h_phi_perp->Draw("EP");

  //TLatex *chlabel1 = new TLatex();
  chlabel1-> SetNDC();
  chlabel1 -> SetTextFont(62);
  chlabel1 -> SetTextColor(1);
  chlabel1 -> SetTextSize(0.045);
  chlabel1 -> SetTextAlign(22);
  chlabel1 -> SetTextAngle(0);
  //char put_text1[256];
  sprintf(put_text1,"PERP, E_{#gamma}#in(%s,%s)",Emin_label, Emax_label);
  chlabel1 -> DrawLatex(0.7, 0.75, put_text1);
  sprintf(put_text1,"-t #in(%s,%s)",tmin_label, tmax_label);
  chlabel1 -> DrawLatex(0.7, 0.65, put_text1);

  //TCanvas*c3 = new TCanvas();

  c->cd(3);

  lhcbStyle();

  gROOT->SetStyle("Plain");
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.10);

  lhcbStyle();


  h_Sigma->SetLineWidth(2);
  h_Sigma->SetLineColor(kRed);
  h_Sigma->SetMarkerColor(kRed);

  h_Sigma->GetXaxis()->SetTitle("#phi_{p} [deg]");
  h_Sigma->GetXaxis()->SetTitleSize(0.07);
  h_Sigma->GetYaxis()->SetTitle("A(#phi)");
  h_Sigma->GetYaxis()->SetTitleSize(0.07);
  h_Sigma->GetXaxis()->GetCenterTitle();
  h_Sigma->GetYaxis()->GetCenterTitle();

  h_Sigma->Draw("EP");

  TLatex *chlabelf = new TLatex();
  chlabelf-> SetNDC();
  chlabelf -> SetTextFont(62);
  chlabelf -> SetTextColor(1);
  chlabelf -> SetTextSize(0.045);
  chlabelf -> SetTextAlign(22);
  chlabelf -> SetTextAngle(0);
  char put_textf[256];
  sprintf(put_textf,"Asym, E_{#gamma}#in(%s,%s)",Emin_label, Emax_label);
  chlabelf -> DrawLatex(0.7, 0.75, put_textf);
  sprintf(put_textf,"-t #in(%s,%s)",tmin_label, tmax_label);
  chlabelf -> DrawLatex(0.7, 0.65, put_textf);

  TCanvas *ct = new TCanvas();

  gROOT->SetStyle("Plain");
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.10);

  lhcbStyle();


  h_t->SetLineWidth(2);
  h_t->SetLineColor(kRed);
  h_t->SetMarkerColor(kRed);

  h_t->GetXaxis()->SetTitle("-t [GeV^{2}]");
  h_t->GetXaxis()->SetTitleSize(0.07);
  h_t->GetYaxis()->SetTitle("entries");
  h_t->GetYaxis()->SetTitleSize(0.07);
  h_t->GetXaxis()->GetCenterTitle();
  h_t->GetYaxis()->GetCenterTitle();
  sprintf(put_textf,"t-Mand, E_{#gamma}#in(%s,%s)",Emin_label, Emax_label);
  h_t->Draw("EP");

  chlabelf -> DrawLatex(0.6, 0.75, put_textf);




  //--------------------------------------------------------------------------//
  //                              SAVE PLOTS
  //--------------------------------------------------------------------------//

  TFile *ofile = new TFile(out_file, "RECREATE");

  h_M_para->Write();
  h_M_perp->Write();
  h_phi_para->Write();
  h_phi_perp->Write();
  h_beamE_jpsi_para->Write();
  h_beamE_jpsi_perp->Write();
  h_Sigma->Write();


}



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//







vector <double> calc_pol(double Emin, double Emax, TGraphErrors* gr, TGraph* gr2, TH1F* hm){

  vector <double> res;
  double tmpE, tmpP, tmpSigma;
  int tmpW;

  double pol(0.), sigma(0.);

  double xmin = hm->GetXaxis()->GetXmin();
  double xmax = hm->GetXaxis()->GetXmax();
  double xwidth = hm->GetBinWidth(1);


  int binEmin = calc_bin(Emin, xmin, xwidth);
  int binEmax = calc_bin(Emax, xmin, xwidth);

  double Pavg(0.), dPavg(0.), den(0.);

  for(int i=binEmin; i<=binEmax; i++){
    tmpE = hm->GetBinCenter(i);
    tmpW = hm->GetBinContent(i);
    tmpP = gr->Eval(tmpE);
    tmpSigma = gr2->Eval(tmpE) - tmpP;

    Pavg+=tmpW*tmpP;
    den+=tmpW;

    dPavg+=tmpW*tmpW*tmpSigma*tmpSigma;

  }

  dPavg = sqrt(dPavg);

  if(den>0.){
    Pavg = Pavg/den;
    dPavg = dPavg/den;
  }
  else{
    Pavg = 0.;
    dPavg = 0.;
  }

  //cout<< Emin << ", "<< gr->Eval(Emin) << ", "<< gr2->Eval(Emin) << endl;


  res.push_back(Pavg);
  res.push_back(dPavg);


  //cout<<"den: "<< den << endl;
  //cout<<"avg P: "<< Pavg << " +/- " << dPavg << endl;

  return res;
}

int calc_bin(double val, double min, double width){

  int bin;

  bin = (int) ((val - min) / width);

  return bin;
}

vector <double> calc_Asy(double Yperp, double Ypara){

  vector <double> res;

  double tmp;
  double Ytot = Yperp + Ypara;

  tmp = (Yperp - Ypara)/(Yperp + Ypara);

  res.push_back(tmp);


  //cout<< "Ypara: "<< Ypara << ", Yperp: "<< Yperp << ", Ytot: "<< Ytot << endl;

  tmp = 2./(Ytot*Ytot) * sqrt(Ypara*Ypara*Yperp + Yperp*Yperp*Ypara);

  res.push_back(tmp);

  return res;
}
