/// Orignal developed by Yuanjing, Ji
/// Modified by Yu Hu for the d-lambda correlation
//#include "style.h"
//#include "draw.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TLine.h"

void getKinematics(TLorentzVector& b, double const mass)
{
    double const pt = gRandom->Uniform(0.,4.);
    double const y = gRandom->Uniform(-1., 0.);
    double const phi = TMath::TwoPi() * gRandom->Rndm();
    
    double const mT = sqrt(mass * mass + pt * pt);
    double const pz = mT * sinh(y);
    double const E = mT * cosh(y);
    
    b.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

double getptweight(TLorentzVector b)
{
  //let's get the same pt weight as embedding data
  TF1* fH3Ldydpt[4]; 
  double par[4][3]={{0.250571, 101, 2.99089}, {0.234906, 21.0216*1e1, 2.99089}, {0.174251, 253*1e2, 2.99089}, { 0.133087, 6264.09*1e3, 2.99089}};
  for (int irap=0;irap<4;irap++) {
    fH3Ldydpt[irap] = new TF1(Form("fH3Ldydpt[%d]",irap), "2*TMath::Pi()*[1]*x*exp(-(sqrt([2]*[2]+x*x))/[0])", 0,5);
    fH3Ldydpt[irap]->SetParameters(par[irap]);
  }

  double ycm=-1.045;
  double pt=b.Pt();
  double y=-1*(b.Rapidity() - ycm); 
  double ptweight=1;
  if (y>=0 && y<=0.25) ptweight = fH3Ldydpt[0]->Eval(pt);
  if (y<0  && y>=-0.25) ptweight = fH3Ldydpt[0]->Eval(pt);
  if (y<-0.25 && y>=-0.5) ptweight = fH3Ldydpt[1]->Eval(pt);
  if (y<-0.5 && y>=-0.75) ptweight = fH3Ldydpt[2]->Eval(pt);
  if (y<-0.75) ptweight = fH3Ldydpt[3]->Eval(pt);

  //cout << "weight_pt " << ptweight << endl;
  return double(ptweight);
}

double getPfromM(double M, double m1, double m2)
{
  //give M, m1 m2, calculate the input p
  return sqrt(pow(M*M+m1*m1-m2*m2,2)-4*M*M*m1*m1)/(2*M);
}

double getm1fromP(double M, double m2, double p)
{
  //give the input M, m2, p; calculate the m1
  return sqrt(-2*sqrt(M*M*(m2*m2 + p*p)) + M*M + m2*m2);
}

void pseudoPtEtaPhi(TLorentzVector& l, TH1F *hpt, TH1F *heta, TH1F *hphi)
{
  // based on the given sample, give pseudo Pt, Eta, and Phi
  double pt=l.Pt();
  double eta=l.Eta();
  double phi=l.Phi();
  double m=l.M();
  //  cout << "O: " << pt << " " << eta << " " << phi << endl;
  double pdpt=hpt->GetRandom();
  double pdeta=heta->GetRandom();
  double pdphi=hphi->GetRandom();
  //  cout << "I: " << pdpt << " " << pdeta << " " << pdphi << endl;
  l.SetPtEtaPhiM(pt+pdpt, eta+pdeta, phi+pdphi, m);
  //  cout << "F: " << l.Pt() << " " << l.Eta() << " " << l.Phi() << endl;
}

void pseudoPtEtaPhiMass(TLorentzVector& l, TH1F *hpt, TH1F *heta, TH1F *hphi, TH1F *hmass)
{
  // based on the given sample, give pseudo Pt, Eta, and Phi
  double pt=l.Pt();
  double eta=l.Eta();
  double phi=l.Phi();
  double m=l.M();
  //  cout << "O: " << pt << " " << eta << " " << phi << " " << m << endl;
  double pdpt=hpt->GetRandom();
  double pdeta=heta->GetRandom();
  double pdphi=hphi->GetRandom();
  double pdmass=hmass->GetRandom();
  //  cout << "I: " << pdpt << " " << pdeta << " " << pdphi << " " <<  pdmass<< endl;
  l.SetPtEtaPhiM(pt+pdpt, eta+pdeta, phi+pdphi, m+pdmass);
  //  cout << "F: " << l.Pt() << " " << l.Eta() << " " << l.Phi() <<  " " << l.M() << endl;
}

void H3Ldecay_quasi2body()
{
  TStopwatch time;
  time.Start();
  
  //  style();
  // const Double_t M_H = 2.230;
  // const Double_t M_H = 2.9913;
  const Double_t M_d = 1.875613;
  // const Double_t M_Xi = 1.32171;
  const Double_t M_Xi = 1.115683;
  const Double_t M_pi = 0.13957;
  const Double_t M_p = 0.93827;
  const Double_t M_pi0 = 0.134975;
  const Double_t M_Ld = 1.115683;
  const Double_t M_vLd = 1.115;
  const Double_t M_vd = 1.875613;
  const Double_t M_H3L= 2.99089;
  const Double_t M_H3L_x=2.99130;
  Double_t masses2[2] = { M_vLd, M_vd} ;
  Double_t masses3[3] = { M_p, M_pi, M_d} ;
  Double_t masses22[2] = { M_p, M_pi} ;

  int const nevents = 1000000;

  TGenPhaseSpace event1, event2,event3;
  TRandom3 *gRandom = new TRandom3();
  TRandom* vLRand = new TRandom();

  char * histname=new char[100];
  
  TH2F* hpd_dpi = new TH2F("hpd_dpi","hpd_dpi;M(pd) (GeV/c^{2});M(dpi) (GeV/c^{2})", 100, 2.8, 2.86, 100, 2.01, 2.06 );
  TH2F* hppi_dpi= new TH2F("hppi_dpi","hppi_dpi;M(p#pi) (GeV/c^{2});M(d#pi) (GeV/c^{2})", 100,1.05,1.15, 100, 2.01,2.06);
  TH2F* hppi_pd = new TH2F("hppi_pd","hppi_pd;M(p#pi) (GeV/c^{2});M(pd) (GeV/c^{2})", 100, 1.05,1.15, 100, 2.8, 2.86 );
  TH2F* hTpi_Td= new TH2F("hTpi_Td","hTpi_Td;Tpi (MeV);Td (MeV)", 100, 0, 36, 100, 0, 36);

  TH2F * hdphideta_MC;
  TH1F * hkstar_MC;
  TH2F * hdphideta_EBD;
  TH1F * hkstar_EBD;

  // let's check the  pt eta phi differences for d and vL samples
  TH1F * hptdiff_d = new TH1F("hptdiff_d","hptdiff_d",120,-.06,0.06);
  TH1F * hetadiff_d = new TH1F("hetadiff_d","hetadiff_d",100,-.03,0.03);
  TH1F * hphidiff_d = new TH1F("hphidiff_d","hphidiff_d",200,-.1,0.1);

  TH1F * hptdiff_l = new TH1F("hptdiff_l","hptdiff_l",120,-.06,0.06);
  TH1F * hetadiff_l = new TH1F("hetadiff_l","hetadiff_l",100,-.03,0.03);
  TH1F * hphidiff_l = new TH1F("hphidiff_l","hphidiff_l",200,-.1,0.1);  
  TH1F * hmassdiff_l = new TH1F("hmassdiff_l","hmassdiff_l",200,-.05,0.05);


  // load the root file for pt,eta,phi sampling
  TFile * f = new TFile("fout_H3L_MC_0050_015pt_sys_new.root","READ");
  // TH1F * hptmcd_p = (TH1F *)f->Get("hptmcdiff_p"); 
  // TH1F * hetamcd_p = (TH1F *)f->Get("hetamcdiff_p"); 
  // TH1F * hphimcd_p = (TH1F *)f->Get("hphimcdiff_p"); 

  // TH1F * hptmcd_pi = (TH1F *)f->Get("hptmcdiff_pi");
  // TH1F * hetamcd_pi = (TH1F *)f->Get("hetamcdiff_pi");
  // TH1F * hetamcd_pi = (TH1F *)f->Get("hetamcdiff_pi");

  TH1F * hptmcd_d = (TH1F *)f->Get("hptmcdiff_d");
  TH1F * hetamcd_d = (TH1F *)f->Get("hetamcdiff_d");
  TH1F * hphimcd_d = (TH1F *)f->Get("hphimcdiff_d");

  TH1F * hptmcd_l = (TH1F *)f->Get("hptmcdiff_l"); 
  TH1F * hetamcd_l = (TH1F *)f->Get("hetamcdiff_l"); 
  TH1F * hphimcd_l = (TH1F *)f->Get("hphimcdiff_l");
  TH1F * hmassmcd_l = (TH1F *)f->Get("hmassmcdiff_l"); 

  // creat a file to save the result hists
  TFile* file = new TFile("H3L_Dz2b_double.root","recreate");

  sprintf(histname,"hdphideta_MC");
  hdphideta_MC= new TH2F(histname,";dphi;deta",600,-3.,3.,600,-3.,3.);
  sprintf(histname,"hkstar_MC");
  hkstar_MC = new TH1F(histname,histname, 800, 0, 200);
  sprintf(histname,"hdphideta_EBD");
  hdphideta_EBD= new TH2F(histname,";dphi;deta",600,-3.,3.,600,-3.,3.);
  sprintf(histname,"hkstar_EBD");
  hkstar_EBD = new TH1F(histname,histname, 800, 0, 200);

  cout << "The corresponding p: " <<  getPfromM(M_H3L,M_vLd,M_vd) << endl;

  for (int ie=0;ie<nevents;ie++){
    if(ie%(nevents/10)==0){cout << "Event " << ie << "; " <<double(ie)/double(nevents)*100. << "%" << endl;}

    //generate quasi-2 body decay 
    TLorentzVector H3Lf3b; 
    TLorentzVector H3Lq2b; 

    //give the M_H3L a momentum distribution by doing the ptweight
    getKinematics(H3Lq2b,M_H3L);
    double ptweight=getptweight(H3Lq2b);

    event1.SetDecay(H3Lq2b, 2, masses2); //lambda & d
    Double_t weight1 = event1.Generate();
    TLorentzVector *vLd= event1.GetDecay(0);
    TLorentzVector *d = event1.GetDecay(1);

    //ckstar
    TLorentzVector dMC;
    TLorentzVector vLdMC;
    TLorentzVector dEBD;
    TLorentzVector vLdEBD;

    // TLV for d and vLd, MC 
    dMC.SetXYZM(d->X(),d->Y(),d->Z(),M_d);
    vLdMC.SetXYZM(vLd->X(),vLd->Y(),vLd->Z(),M_Ld);
    // TLV for d and vLd, embedding 
    dEBD.SetXYZM(d->X(),d->Y(),d->Z(),M_d);
    pseudoPtEtaPhi(dEBD, hptmcd_d, hetamcd_d, hphimcd_d);
    vLdEBD.SetXYZM(vLd->X(),vLd->Y(),vLd->Z(),M_Ld);
    pseudoPtEtaPhiMass(vLdEBD, hptmcd_l, hetamcd_l, hphimcd_l, hmassmcd_l);


    //QA, check pt, eta, phi differ after pseudo embedding
    double pdiff;
    pdiff=dEBD.Pt()-dMC.Pt();
    hptdiff_d->Fill(pdiff);
    pdiff=dEBD.Eta()-dMC.Eta();
    hetadiff_d->Fill(pdiff);
    pdiff=dEBD.Phi()-dMC.Phi();
    hphidiff_d->Fill(pdiff);

    pdiff=vLdEBD.Pt()-vLdMC.Pt();
    hptdiff_l->Fill(pdiff);
    pdiff=vLdEBD.Eta()-vLdMC.Eta();
    hetadiff_l->Fill(pdiff);
    pdiff=vLdEBD.Phi()-vLdMC.Phi();
    hphidiff_l->Fill(pdiff);
    pdiff=vLdEBD.M()-vLdMC.M();
    hmassdiff_l->Fill(pdiff);

    
    // check the deta dphi, MC
    double dphi4pi, dphiMC, detaMC, dphiEBD, detaEBD;
    dphi4pi=dMC.Phi()-vLdMC.Phi();
    dphiMC=dphi4pi;
    if(dphi4pi<-TMath::Pi()){dphiMC=dphi4pi+2*TMath::Pi();}
    if(dphi4pi>TMath::Pi()){dphiMC=dphi4pi-2*TMath::Pi();}
    detaMC=dMC.Eta()-vLdMC.Eta();
    // check the deta dphi, embedding
    dphi4pi=dEBD.Phi()-vLdEBD.Phi();
    dphiEBD=dphi4pi;
    if(dphi4pi<-TMath::Pi()){dphiEBD=dphi4pi+2*TMath::Pi();}
    if(dphi4pi>TMath::Pi()){dphiEBD=dphi4pi-2*TMath::Pi();}
    detaEBD=dEBD.Eta()-vLdEBD.Eta();

    //kstar, MC
    TLorentzVector HMC = dMC+ vLdMC;
    TLorentzVector QvectMC = (vLdMC-dMC);
    double PinvMC = HMC.Mag();
    double Q1MC = (M_Ld*M_Ld-M_d*M_d)/PinvMC;
    double QMC=sqrt(Q1MC*Q1MC-QvectMC.Mag2());
    double kstarMC = QMC/2.0 * 1000; // convert to MeV

    //kstar, EBD
    TLorentzVector HEBD = dEBD+ vLdEBD;
    TLorentzVector QvectEBD = (vLdEBD-dEBD);
    double PinvEBD = HEBD.Mag();
    double Q1EBD = (M_Ld*M_Ld-M_d*M_d)/PinvEBD;
    double QEBD=sqrt(Q1EBD*Q1EBD-QvectEBD.Mag2());
    double kstarEBD = QEBD/2.0 * 1000; // convert to MeV

    //Add ckstar
    hkstar_MC->Fill(kstarMC,weight1*ptweight);
    hdphideta_MC->Fill(dphiMC,detaMC,weight1*ptweight);
    //Add ckstar
    hkstar_EBD->Fill(kstarEBD,weight1*ptweight);
    hdphideta_EBD->Fill(dphiEBD,detaEBD,weight1*ptweight);

    
    //free memory
    //delete d, vLd;
      
  } // entry loop

  hptdiff_d->Write();
  hetadiff_d->Write();
  hphidiff_d->Write();

  hptdiff_l->Write();
  hetadiff_l->Write();
  hphidiff_l->Write();
  hmassdiff_l->Write();
  
  hdphideta_MC->Write();
  hkstar_MC->Write();
  hdphideta_EBD->Write();
  hkstar_EBD->Write();

  
  time.Stop();
  time.Print();
}

// void H3Ldecay_free3body()
// {
//   TStopwatch time;
//   time.Start();
  
//   //  style();
//   // const Double_t M_H = 2.230;
//   // const Double_t M_H = 2.9913;
//   const Double_t M_d = 1.8756;
//   // const Double_t M_Xi = 1.32171;
//   const Double_t M_Xi = 1.115683;
//   const Double_t M_pi = 0.13957;
//   const Double_t M_p = 0.93827;
//   const Double_t M_pi0 = 0.134975;
//   const Double_t M_Ld = 1.115683;
//   const Double_t M_vLd = 1.115;
//   const Double_t M_vd = 1.875613;
//   const Double_t M_H3L= 2.99089+(0.00041-0.00013); //binding energy = 0.13MeV
//   Double_t masses2[2] = { M_vLd, M_vd} ;
//   Double_t masses3[3] = { M_p, M_pi, M_d} ;
//   Double_t masses22[2] = { M_p, M_pi} ;

//   int const nevents = 1000;

//   TGraph* g = new TGraph("xsection.csv");

//   TGenPhaseSpace event1, event2,event3;
//   TRandom3 *gRandom = new TRandom3();

//   TH2F* hpd_dpi = new TH2F("hpd_dpi","hpd_dpi;M(pd) (GeV/c^{2});M(dpi) (GeV/c^{2})", 100, 2.8, 2.86, 100, 2.01, 2.06 );
//   TH2F* hppi_dpi= new TH2F("hppi_dpi","hppi_dpi;M(p#pi) (GeV/c^{2});M(d#pi) (GeV/c^{2})", 100,1.05,1.15, 100, 2.01,2.06);
//   TH2F* hppi_pd = new TH2F("hppi_pd","hppi_pd;M(p#pi) (GeV/c^{2});M(pd) (GeV/c^{2})", 100, 1.05,1.15, 100, 2.8, 2.86 );
//   TH2F* hTpi_Td= new TH2F("hTpi_Td","hTpi_Td;Tpi (MeV);Td (MeV)", 100, 0, 36, 100, 0, 36);
//   TH2F * hdphideta = new TH2F("hdphideta","hdphideta",200,-5,5,200,-5,5);

  
//   for (int ie=0;ie<nevents;ie++){
//       TLorentzVector H3Lf3b; 
//       TLorentzVector H3Lq2b; 

//       //generate free 3  body
//       H3Lf3b.SetXYZM( 0,0,0,M_H3L);
//       // bool isallow = event3.SetDecay(H3Lf3b, 3, masses3,"Fermi");
//       bool isallow = event3.SetDecay(H3Lf3b, 3, masses3);
//       if (isallow){
//       double weight = event3.Generate();
//       TLorentzVector *p= event3.GetDecay(0);
//       TLorentzVector *pi = event3.GetDecay(1);
//       TLorentzVector *d = event3.GetDecay(2);
//       // // cout << (*d).Pt()<< " "<<(*d).M() << endl;

//       TLorentzVector ppi = *p + *pi;
//       TLorentzVector pd = *p + *d;
//       TLorentzVector dpi = *pi+*d;

//       double Tpi = ( sqrt(M_pi*M_pi+(*pi).P()*(*pi).P())-M_pi )*1000;
//       double Tp = ( (*p).P()*(*p).P()/2./M_p)*1000;
//       double Td = ( (*p).P()*(*p).P()/2./M_d)* 1000;
//       // cout << Td+Tpi+Tp << endl;
//       // cout << ( (*d).Vect() + (*pi).Vect() + (*p).Vect() ).Mag() << endl;

//       TLorentzVector LA1 = ppi;
//       TLorentzVector LA2 = *d;
//       double dphi4pi, dphi, deta;
//       dphi4pi=LA1.Phi()-LA2.Phi();
//       dphi=dphi4pi;
//       if(dphi4pi<=-TMath::Pi()){dphi=dphi4pi+2*TMath::Pi();}
//       if(dphi4pi>=TMath::Pi()){dphi=dphi4pi-2*TMath::Pi();}
//       //	deta=rapidity1-rapidity2;
//       deta=LA1.Eta()-LA2.Eta();
      
//       hdphideta->Fill(dphi,deta);

      
//       // cout << ppi.M()<<" "<<pd.M()<<endl;
//       double weight_gaus = TMath::Gaus(Tpi, 32, 1, true);
//       if (Td<0.5) weight_gaus = weight_gaus* (g->Eval(Td));
//       else weight_gaus = 0;
//       // double weight_gaus = 1;
//       // if ( Td<0.3 && Tpi<33 && Tpi>31 ){
//         hppi_pd->Fill( ppi.M(), pd.M(),weight*weight_gaus); 
//         hppi_dpi->Fill( ppi.M(), dpi.M(),weight*weight_gaus); 
//         hpd_dpi->Fill( pd.M(), dpi.M(),weight*weight_gaus); 
//         hTpi_Td->Fill( Tpi, Td, weight*weight_gaus);
//       // }
//     }
//     else {cout << "decay not allowed" << endl;}
//   }

//   TFile* file = new TFile("H3L_Dz3b.root","recreate");
//   hppi_pd->Write();
//   hpd_dpi->Write();
//   hppi_dpi->Write();
//   hTpi_Td->Write();
//   hdphideta->Write();
  
//   time.Stop();
//   time.Print();
// }

void H3Ldecay()
{
  //  H3Ldecay_free3body();
  H3Ldecay_quasi2body();
}


