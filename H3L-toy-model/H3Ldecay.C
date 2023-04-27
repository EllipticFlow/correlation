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
  bool f3body=true;

  TGenPhaseSpace event1, event2,event3;
  double weight1, weight2, weight3;
  double pdiff;
  TRandom3 *gRandom = new TRandom3();
  TRandom* vLRand = new TRandom();

  char * histname=new char[100];
  
  TH2F * hdphideta_MC;
  TH1F * hkstar_MC;
  TH2F * hdphideta_EBD;
  TH1F * hkstar_EBD;
  TH2F * hdphideta_3body;
  TH1F * hkstar_3body;

  
  // let's check the  pt eta phi differences for d and vL samples
  TH1F * hptdiff_p = new TH1F("hptdiff_p","hptdiff_p",120,-.06,0.06);
  TH1F * hetadiff_p = new TH1F("hetadiff_p","hetadiff_p",100,-.03,0.03);
  TH1F * hphidiff_p = new TH1F("hphidiff_p","hphidiff_p",200,-.1,0.1);

  TH1F * hptdiff_pi = new TH1F("hptdiff_pi","hptdiff_pi",120,-.06,0.06);
  TH1F * hetadiff_pi = new TH1F("hetadiff_pi","hetadiff_pi",100,-.03,0.03);
  TH1F * hphidiff_pi = new TH1F("hphidiff_pi","hphidiff_pi",200,-.1,0.6);

  TH1F * hptdiff_d = new TH1F("hptdiff_d","hptdiff_d",120,-.06,0.06);
  TH1F * hetadiff_d = new TH1F("hetadiff_d","hetadiff_d",100,-.03,0.03);
  TH1F * hphidiff_d = new TH1F("hphidiff_d","hphidiff_d",200,-.1,0.1);

  TH1F * hptdiff_l = new TH1F("hptdiff_l","hptdiff_l",120,-.06,0.06);
  TH1F * hetadiff_l = new TH1F("hetadiff_l","hetadiff_l",100,-.03,0.03);
  TH1F * hphidiff_l = new TH1F("hphidiff_l","hphidiff_l",200,-.1,0.1);  
  TH1F * hmassdiff_l = new TH1F("hmassdiff_l","hmassdiff_l",200,-.05,0.05);
  TH1F * hmassdiff_l3body = new TH1F("hmassdiff_l3body","hmassdiff_l3body",200,-.05,0.05);


  // load the root file for pt,eta,phi sampling
  TFile * f = new TFile("fout_H3L_MC_0050_015pt_sys_new.root","READ");
  TH1F * hptmcd_p = (TH1F *)f->Get("hptmcdiff_p"); 
  TH1F * hetamcd_p = (TH1F *)f->Get("hetamcdiff_p"); 
  TH1F * hphimcd_p = (TH1F *)f->Get("hphimcdiff_p"); 

  TH1F * hptmcd_pi = (TH1F *)f->Get("hptmcdiff_pi");
  TH1F * hetamcd_pi = (TH1F *)f->Get("hetamcdiff_pi");
  TH1F * hphimcd_pi = (TH1F *)f->Get("hphimcdiff_pi");

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

  TH2F* hTpi_Td;
  
  if(f3body){
    sprintf(histname,"hdphideta_3body");
    hdphideta_3body= new TH2F(histname,";dphi;deta",600,-3.,3.,600,-3.,3.);
    sprintf(histname,"hkstar_3body");
    hkstar_3body = new TH1F(histname,histname, 800, 0, 200);
    
    hTpi_Td = new TH2F("hTpi_Td","hTpi_Td;Tpi (MeV);Td (MeV)", 4000, 0, 40, 4000, 0, 40);
  }
  
  cout << "The corresponding p: " <<  getPfromM(M_H3L,M_vLd,M_vd) << endl;
  
  for (int ie=0;ie<nevents;ie++){
       if(ie%(nevents/10)==0){cout << "Event " << ie << "; " <<double(ie)/double(nevents)*100. << "%" << endl;}

    //generate quasi-2 body decay 
    TLorentzVector H3Lf3b; 
    TLorentzVector H3Lq2b; 

    //give the M_H3L a momentum distribution by doing the ptweight
    //    H3Lq2b.SetXYZM(0.,0.,0.,M_H3L);
    getKinematics(H3Lq2b,M_H3L);
    double ptweight=getptweight(H3Lq2b);
    event1.SetDecay(H3Lq2b, 2, masses2); //lambda & d
    weight1 = event1.Generate();
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
    // TLV for d and vLd, embedding, smearing 
    dEBD.SetXYZM(d->X(),d->Y(),d->Z(),M_d);
    pseudoPtEtaPhi(dEBD, hptmcd_d, hetamcd_d, hphimcd_d);
    vLdEBD.SetXYZM(vLd->X(),vLd->Y(),vLd->Z(),M_Ld);
    pseudoPtEtaPhiMass(vLdEBD, hptmcd_l, hetamcd_l, hphimcd_l, hmassmcd_l);
    //pseudoPtEtaPhi(vLdEBD, hptmcd_l, hetamcd_l, hphimcd_l);
 
    
    TLorentzVector pMC, piMC, pEBD, piEBD;
    if(f3body){
      //let's try to do a further decay for the MC vLd 
      event2.SetDecay(vLdMC, 2, masses22);
      weight2 = event2.Generate();
      TLorentzVector *p= event2.GetDecay(0);
      TLorentzVector *pi = event2.GetDecay(1);

      //check the Dalitz Plots
      //boost to the rest frame for Kamada's selection
      TVector3 boost_vec = -H3Lq2b.BoostVector(); // Boost vector
      TLorentzVector piboost, pboost, dboost;
      piboost.SetXYZM( pi->X(), pi->Y(), pi->Z(), pi->M());
      pboost.SetXYZM( p->X(), p->Y(), p->Z(), p->M());
      dboost.SetXYZM( d->X(), d->Y(), d->Z(), d->M());
      piboost.Boost(boost_vec); // Boost b
      pboost.Boost(boost_vec); // Boost b
      dboost.Boost(boost_vec); // Boost b

      double Tpi = ( sqrt(M_pi*M_pi+(piboost).P()*(piboost).P())-M_pi )*1000;
      double Tp = ( (pboost).P()*(pboost).P()/2./M_p)*1000;
      double Td = ( (dboost).P()*(dboost).P()/2./M_d)* 1000;
      // cout << (*pi).X() << " " << (*p).X() << " " << (*d).X() << " "  << (*pi).X()+(*p).X()+(*d).X() << endl;
      // cout << (piboost).X() << " " << (pboost).X() << " " << (dboost).X() << " "  << (piboost).X()+(pboost).X()+(dboost).X() << endl;
      //  cout << (piboost).P() << " " << (pboost).P() << " " << (dboost).P() << endl;

      //without decay selection, for all possible decay kinematics
      hTpi_Td->Fill( Tpi, Td, weight1*weight2*ptweight);
      
      pMC.SetXYZM(p->X(),p->Y(),p->Z(),M_p);
      piMC.SetXYZM(pi->X(),pi->Y(),pi->Z(),M_pi);
      /// smearing step
      pEBD.SetXYZM(p->X(),p->Y(),p->Z(),M_p);
      pseudoPtEtaPhi(pEBD, hptmcd_p, hetamcd_p, hphimcd_p);
      piEBD.SetXYZM(pi->X(),pi->Y(),pi->Z(),M_pi);
      pseudoPtEtaPhi(piEBD, hptmcd_pi, hetamcd_pi, hphimcd_pi);

      pdiff=pEBD.Pt()-pMC.Pt();
      hptdiff_p->Fill(pdiff);
      pdiff=pEBD.Eta()-pMC.Eta();
      hetadiff_p->Fill(pdiff);
      pdiff=pEBD.Phi()-pMC.Phi();
      hphidiff_p->Fill(pdiff);

      pdiff=piEBD.Pt()-piMC.Pt();
      hptdiff_pi->Fill(pdiff);
      pdiff=piEBD.Eta()-piMC.Eta();
      hetadiff_pi->Fill(pdiff);
      pdiff=piEBD.Phi()-piMC.Phi();
      hphidiff_pi->Fill(pdiff);      
    }
    
    

    //QA, check pt, eta, phi differ after pseudo embedding
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


    if(f3body){
      //after smearing, add pi+p to vld
      TLorentzVector vLd3body=piEBD+pEBD;
      pdiff=vLd3body.M()-vLdMC.M();
      hmassdiff_l3body->Fill(pdiff);
      
      double dphi3body, deta3body;
      dphi4pi=dEBD.Phi()-vLd3body.Phi();
      dphi3body=dphi4pi;
      if(dphi4pi<-TMath::Pi()){dphi3body=dphi4pi+2*TMath::Pi();}
      if(dphi4pi>TMath::Pi()){dphi3body=dphi4pi-2*TMath::Pi();}
      deta3body=dEBD.Eta()-vLd3body.Eta();

      //kstar, 3body
      TLorentzVector H3body = dEBD+ vLd3body;
      TLorentzVector Qvect3body = (vLd3body-dEBD);
      double Pinv3body = H3body.Mag();
      double Q13body = (M_Ld*M_Ld-M_d*M_d)/Pinv3body;
      double Q3body=sqrt(Q13body*Q13body-Qvect3body.Mag2());
      double kstar3body = Q3body/2.0 * 1000; // convert to MeV

      //Add ckstar
      hkstar_3body->Fill(kstar3body,weight1*weight2*ptweight);
      hdphideta_3body->Fill(dphi3body,deta3body,weight1*weight2*ptweight);      
    }
    
    
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

  if(f3body){
    hptdiff_p->Write();
    hetadiff_p->Write();
    hphidiff_p->Write();

    hptdiff_pi->Write();
    hetadiff_pi->Write();
    hphidiff_pi->Write();

    hmassdiff_l3body->Write();  

    hdphideta_3body->Write();
    hkstar_3body->Write();

    hTpi_Td->Write();    
  }
  
  time.Stop();
  time.Print();
}


void H3Ldecay_free3body()
{
  TStopwatch time;
  time.Start();
  
  //  style();
  // const Double_t M_H = 2.230;
  // const Double_t M_H = 2.9913;
  const Double_t M_d = 1.8756;
  // const Double_t M_Xi = 1.32171;
  const Double_t M_Xi = 1.115683;
  const Double_t M_pi = 0.13957;
  const Double_t M_p = 0.93827;
  const Double_t M_pi0 = 0.134975;
  const Double_t M_Ld = 1.115683;
  const Double_t M_vLd = 1.115;
  const Double_t M_vd = 1.875613;
  //  const Double_t M_H3L= 2.99089+(0.00041-0.00013); //binding energy = 0.13MeV
  //  const Double_t M_H3L= 2.99089; //STAR 0.41
  const Double_t M_H3L=2.99089+0.00041;
  Double_t masses2[2] = { M_vLd, M_vd} ;
  Double_t masses3[3] = { M_p, M_pi, M_d} ;
  Double_t masses22[2] = { M_p, M_pi} ;

  int const nevents = 1000000;

  //load the paper data
  TGraph* g = new TGraph("Td_Kamada_32MeV.csv");
  //  TGraph* gham = new TGraph("H3L_phasespace_Hammer.csv");
  
  TGenPhaseSpace event1, event2,event3;
  TRandom3 *gRandom = new TRandom3();

  char * histname=new char[100];
  

  TH2F * hdphideta_MC;
  TH1F * hkstar_MC;
  TH2F * hdphideta_EBD;
  TH1F * hkstar_EBD;

  TH2F * hdphideta_MC_Kamada;
  TH1F * hkstar_MC_Kamada;
  TH2F * hdphideta_EBD_Kamada;
  TH1F * hkstar_EBD_Kamada;

  sprintf(histname,"hdphideta_MC");
  hdphideta_MC= new TH2F(histname,";dphi;deta",600,-3.,3.,600,-3.,3.);
  sprintf(histname,"hkstar_MC");
  hkstar_MC = new TH1F(histname,histname, 800, 0, 200);
  sprintf(histname,"hdphideta_EBD");
  hdphideta_EBD= new TH2F(histname,";dphi;deta",600,-3.,3.,600,-3.,3.);
  sprintf(histname,"hkstar_EBD");
  hkstar_EBD = new TH1F(histname,histname, 800, 0, 200);

  sprintf(histname,"hdphideta_MC_Kamada");
  hdphideta_MC_Kamada= new TH2F(histname,";dphi;deta",600,-3.,3.,600,-3.,3.);
  sprintf(histname,"hkstar_MC_Kamada");
  hkstar_MC_Kamada = new TH1F(histname,histname, 800, 0, 200);
  sprintf(histname,"hdphideta_EBD_Kamada");
  hdphideta_EBD_Kamada= new TH2F(histname,";dphi;deta",600,-3.,3.,600,-3.,3.);
  sprintf(histname,"hkstar_EBD_Kamada");
  hkstar_EBD_Kamada = new TH1F(histname,histname, 800, 0, 200);

  
  // let's check the  pt eta phi differences for d and vL samples
  TH1F * hptdiff_p = new TH1F("hptdiff_p","hptdiff_p",120,-.06,0.06);
  TH1F * hetadiff_p = new TH1F("hetadiff_p","hetadiff_p",100,-.03,0.03);
  TH1F * hphidiff_p = new TH1F("hphidiff_p","hphidiff_p",200,-.1,0.1);

  TH1F * hptdiff_pi = new TH1F("hptdiff_pi","hptdiff_pi",120,-.06,0.06);
  TH1F * hetadiff_pi = new TH1F("hetadiff_pi","hetadiff_pi",100,-.03,0.03);
  TH1F * hphidiff_pi = new TH1F("hphidiff_pi","hphidiff_pi",200,-.1,0.6);

  TH1F * hptdiff_d = new TH1F("hptdiff_d","hptdiff_d",120,-.06,0.06);
  TH1F * hetadiff_d = new TH1F("hetadiff_d","hetadiff_d",100,-.03,0.03);
  TH1F * hphidiff_d = new TH1F("hphidiff_d","hphidiff_d",200,-.1,0.1);

  TH1F * hptdiff_l = new TH1F("hptdiff_l","hptdiff_l",120,-.06,0.06);
  TH1F * hetadiff_l = new TH1F("hetadiff_l","hetadiff_l",100,-.03,0.03);
  TH1F * hphidiff_l = new TH1F("hphidiff_l","hphidiff_l",200,-.1,0.1);  
  TH1F * hmassdiff_l = new TH1F("hmassdiff_l","hmassdiff_l",200,-.05,0.05);
  TH1F * hmassdiff_l3body = new TH1F("hmassdiff_l3body","hmassdiff_l3body",200,-.05,0.05);

  // load the root file for pt,eta,phi sampling
  TFile * f = new TFile("fout_H3L_MC_0050_015pt_sys_new.root","READ");
  TH1F * hptmcd_p = (TH1F *)f->Get("hptmcdiff_p"); 
  TH1F * hetamcd_p = (TH1F *)f->Get("hetamcdiff_p"); 
  TH1F * hphimcd_p = (TH1F *)f->Get("hphimcdiff_p"); 

  TH1F * hptmcd_pi = (TH1F *)f->Get("hptmcdiff_pi");
  TH1F * hetamcd_pi = (TH1F *)f->Get("hetamcdiff_pi");
  TH1F * hphimcd_pi = (TH1F *)f->Get("hphimcdiff_pi");

  TH1F * hptmcd_d = (TH1F *)f->Get("hptmcdiff_d");
  TH1F * hetamcd_d = (TH1F *)f->Get("hetamcdiff_d");
  TH1F * hphimcd_d = (TH1F *)f->Get("hphimcdiff_d");

  TH1F * hptmcd_l = (TH1F *)f->Get("hptmcdiff_l"); 
  TH1F * hetamcd_l = (TH1F *)f->Get("hetamcdiff_l"); 
  TH1F * hphimcd_l = (TH1F *)f->Get("hphimcdiff_l");
  TH1F * hmassmcd_l = (TH1F *)f->Get("hmassmcdiff_l"); 

  
  TH2F* hpd_dpi = new TH2F("hpd_dpi","hpd_dpi;M(pd) (GeV/c^{2});M(dpi) (GeV/c^{2})", 100, 2.8, 2.86, 100, 2.01, 2.06 );
  TH2F* hppi_dpi= new TH2F("hppi_dpi","hppi_dpi;M(p#pi) (GeV/c^{2});M(d#pi) (GeV/c^{2})", 100,1.05,1.15, 100, 2.01,2.06);
  TH2F* hppi_pd = new TH2F("hppi_pd","hppi_pd;M(p#pi) (GeV/c^{2});M(pd) (GeV/c^{2})", 100, 1.05,1.15, 100, 2.8, 2.86 );
  TH2F* hTpi_Td= new TH2F("hTpi_Td","hTpi_Td;Tpi (MeV);Td (MeV)", 350, 0, 35, 350, 0, 35);

  TH2F* hpd_dpi_Kamada = new TH2F("hpd_dpi_Kamada","hpd_dpi;M(pd) (GeV/c^{2});M(dpi) (GeV/c^{2})", 100, 2.8, 2.86, 100, 2.01, 2.06 );
  TH2F* hppi_dpi_Kamada = new TH2F("hppi_dpi_Kamada","hppi_dpi;M(p#pi) (GeV/c^{2});M(d#pi) (GeV/c^{2})", 100,1.05,1.15, 100, 2.01,2.06);
  TH2F* hppi_pd_Kamada = new TH2F("hppi_pd_Kamada","hppi_pd;M(p#pi) (GeV/c^{2});M(pd) (GeV/c^{2})", 100, 1.05,1.15, 100, 2.8, 2.86 );
  TH2F* hTpi_Td_Kamada = new TH2F("hTpi_Td_Kamada","hTpi_Td;Tpi (MeV);Td (MeV)", 5000, 30, 35, 5000, 0, 5);

  TH1F* hTd_Kamada = new TH1F("hTd_Kamada","hTd_Kamada;Td;",5000,0,5);
  TH1F* hkd_Kamada = new TH1F("hkd_Kamada","hkd_Kamada;kd;",5000,0,5);
    
  for (int ie=0;ie<nevents;ie++){
       if(ie%(nevents/10)==0){cout << "Event " << ie << "; " <<double(ie)/double(nevents)*100. << "%" << endl;}

    TLorentzVector H3Lf3b; 
    TLorentzVector H3Lq2b; 

    //generate free 3  body
    getKinematics(H3Lf3b,M_H3L);
    double ptweight=getptweight(H3Lf3b);
    //    H3Lf3b.SetXYZM( 0,0,0,M_H3L);
    // bool isallow = event3.SetDecay(H3Lf3b, 3, masses3,"Fermi");
    bool isallow = event3.SetDecay(H3Lf3b, 3, masses3);
    if (isallow){
      double weight = event3.Generate();
      TLorentzVector *p= event3.GetDecay(0);
      TLorentzVector *pi = event3.GetDecay(1);
      TLorentzVector *d = event3.GetDecay(2);

      TLorentzVector ppi = *p + *pi;
      TLorentzVector pd = *p + *d;
      TLorentzVector dpi = *pi+*d;

      //boost to the rest frame for Kamada's selection
      TVector3 boost_vec = -H3Lf3b.BoostVector(); // Boost vector
      TLorentzVector piboost, pboost, dboost;
      piboost.SetXYZM( pi->X(), pi->Y(), pi->Z(), pi->M());
      pboost.SetXYZM( p->X(), p->Y(), p->Z(), p->M());
      dboost.SetXYZM( d->X(), d->Y(), d->Z(), d->M());
      piboost.Boost(boost_vec); // Boost b
      pboost.Boost(boost_vec); // Boost b
      dboost.Boost(boost_vec); // Boost b
      
      double Tpi = ( sqrt(M_pi*M_pi+(piboost).P()*(piboost).P())-M_pi )*1000;
      double Tp = ( (pboost).P()*(pboost).P()/2./M_p)*1000;
      double Td = ( (dboost).P()*(dboost).P()/2./M_d)* 1000;
      // cout << (*pi).X() << " " << (*p).X() << " " << (*d).X() << " "  << (*pi).X()+(*p).X()+(*d).X() << endl;
      // cout << (piboost).X() << " " << (pboost).X() << " " << (dboost).X() << " "  << (piboost).X()+(pboost).X()+(dboost).X() << endl;
      // cout << (piboost).P() << " " << (pboost).P() << " " << (dboost).P() << endl;

      /// without decay selection, for all possible decay kinematics
      hppi_pd->Fill( ppi.M(), pd.M(),weight*ptweight); 
      hppi_dpi->Fill( ppi.M(), dpi.M(),weight*ptweight); 
      hpd_dpi->Fill( pd.M(), dpi.M(),weight*ptweight); 
      hTpi_Td->Fill( Tpi, Td, weight*ptweight);

      
      //weight factor based on Fig7, Gaus on Tpi, Kamada, arXiv:nucl-th/9709035
      double weight_Kamada = TMath::Gaus(Tpi, 32, 1, true);
      if (Td<0.5) weight_Kamada = weight_Kamada* (g->Eval(Td));
      else weight_Kamada = 0;

      // Dalitz plot from Kamada's study 
      hppi_pd_Kamada->Fill( ppi.M(), pd.M(),weight*ptweight*weight_Kamada); 
      hppi_dpi_Kamada->Fill( ppi.M(), dpi.M(),weight*ptweight*weight_Kamada); 
      hpd_dpi_Kamada->Fill( pd.M(), dpi.M(),weight*ptweight*weight_Kamada); 
      hTpi_Td_Kamada->Fill( Tpi, Td, weight*ptweight*weight_Kamada);
      hTd_Kamada->Fill(Td,weight*ptweight*weight_Kamada);
      hkd_Kamada->Fill((dboost).P(),weight*ptweight*weight_Kamada);

      //ckstar
      TLorentzVector dMC, pMC, piMC, vLdMC;
      TLorentzVector dEBD, pEBD, piEBD, vLdEBD;

      // TLV for d and vLd, MC 
      dMC.SetXYZM(d->X(),d->Y(),d->Z(),M_d);
      //      vLdMC.SetXYZM(ppi->X(),ppi->Y(),ppi->Z(),M_Ld);
      vLdMC= *p + *pi;
      
      // TLV for d and vLd, embedding, smearing 
      dEBD.SetXYZM(d->X(),d->Y(),d->Z(),M_d);
      pseudoPtEtaPhi(dEBD, hptmcd_d, hetamcd_d, hphimcd_d);
      pEBD.SetXYZM(p->X(),p->Y(),p->Z(),M_p);
      pseudoPtEtaPhi(pEBD, hptmcd_p, hetamcd_p, hphimcd_p);
      piEBD.SetXYZM(pi->X(),pi->Y(),pi->Z(),M_pi);
      pseudoPtEtaPhi(piEBD, hptmcd_pi, hetamcd_pi, hphimcd_pi);
      vLdEBD=piEBD+pEBD;

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
      hkstar_MC->Fill(kstarMC,weight*ptweight);
      hdphideta_MC->Fill(dphiMC,detaMC,weight*ptweight);
      //Add ckstar
      hkstar_EBD->Fill(kstarEBD,weight*ptweight);
      hdphideta_EBD->Fill(dphiEBD,detaEBD,weight*ptweight);

      //Add ckstar
      hkstar_MC_Kamada->Fill(kstarMC,weight*ptweight*weight_Kamada);
      hdphideta_MC_Kamada->Fill(dphiMC,detaMC,weight*ptweight*weight_Kamada);
      //Add ckstar
      hkstar_EBD_Kamada->Fill(kstarEBD,weight*ptweight*weight_Kamada);
      hdphideta_EBD_Kamada->Fill(dphiEBD,detaEBD,weight*ptweight*weight_Kamada);


    }
    else {cout << "decay not allowed" << endl;}
  }

  TFile* file = new TFile("H3L_Dz3b.root","recreate");

  hppi_pd->Write();
  hpd_dpi->Write();
  hppi_dpi->Write();  
  hTpi_Td->Write();
  hkstar_MC->Write();
  hdphideta_MC->Write();
  hkstar_EBD->Write();
  hdphideta_EBD->Write();

  hppi_pd_Kamada->Write();
  hpd_dpi_Kamada->Write();
  hppi_dpi_Kamada->Write();  
  hTpi_Td_Kamada->Write();
  hkstar_MC_Kamada->Write();
  hdphideta_MC_Kamada->Write();
  hkstar_EBD_Kamada->Write();
  hdphideta_EBD_Kamada->Write();

  hTd_Kamada->Write();
  hkd_Kamada->Write();
  
  time.Stop();
  time.Print();
}

void H3Ldecay()
{
  H3Ldecay_free3body();
  //  H3Ldecay_quasi2body();
}


