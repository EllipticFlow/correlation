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
    float const pt = gRandom->Uniform(0.,4.);
    float const y = gRandom->Uniform(-1., 0.);
    float const phi = TMath::TwoPi() * gRandom->Rndm();
    
    float const mT = sqrt(mass * mass + pt * pt);
    float const pz = mT * sinh(y);
    float const E = mT * cosh(y);
    
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
  return ptweight;
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

  int const nevents = 200000;

  TGenPhaseSpace event1, event2,event3;
  TRandom3 *gRandom = new TRandom3();
  TRandom* vLRand = new TRandom();

  char * histname=new char[100];
  
  TH2F* hpd_dpi = new TH2F("hpd_dpi","hpd_dpi;M(pd) (GeV/c^{2});M(dpi) (GeV/c^{2})", 100, 2.8, 2.86, 100, 2.01, 2.06 );
  TH2F* hppi_dpi= new TH2F("hppi_dpi","hppi_dpi;M(p#pi) (GeV/c^{2});M(d#pi) (GeV/c^{2})", 100,1.05,1.15, 100, 2.01,2.06);
  TH2F* hppi_pd = new TH2F("hppi_pd","hppi_pd;M(p#pi) (GeV/c^{2});M(pd) (GeV/c^{2})", 100, 1.05,1.15, 100, 2.8, 2.86 );
  TH2F* hTpi_Td= new TH2F("hTpi_Td","hTpi_Td;Tpi (MeV);Td (MeV)", 100, 0, 36, 100, 0, 36);

  TH2F * hdphideta[10];
  TH1F * hkstar[10];
  TH1F * hvldM[10];

  TFile* file = new TFile("H3L_Dz2b.root","recreate");

  for(int i=0;i<10;i++)
    {
      sprintf(histname,"hdphideta_%d",i);
      hdphideta[i]= new TH2F(histname,";dphi;deta",600,-3.,3.,600,-3.,3.);
      sprintf(histname,"hkstar_%d",i);
      hkstar[i] = new TH1F(histname,histname, 800, 0, 200);
      sprintf(histname,"hvldM_%d",i);
      hvldM[i] = new TH1F(histname,histname,200,1.1,1.13);
    }
  
  for(int imvld=0; imvld<9; imvld++){
    // a search of v-lambda mass
    //double M_vLd_y=1.11+0.001*imvld;
    double M_vLd_y= M_vLd;

    //a search on the v-lambda mass width
    double M_vLd_width=0.+0.0005*imvld;
    
    cout << "Case " << imvld << "; v-Lambda Mass: " << M_vLd_y << "; Width is " << M_vLd_width <<endl;
    for (int ie=0;ie<nevents;ie++){
      if(ie%100000==0){cout << ie << endl;}
    
      //give a distribution to lambda mass
      //double mvld=vLRand->Gaus(M_vLd_y,0.0005);
      double mvld=vLRand->Gaus(M_vLd_y,M_vLd_width);
      double masses2_g[2] = {mvld, M_vd};

      // the mvld must be smaller than 1.115276
      if(mvld>1.115276)continue;
      hvldM[imvld]->Fill(mvld);

      //generate quasi-2 body decay 
      TLorentzVector H3Lf3b; 
      TLorentzVector H3Lq2b; 

      //give the M_H3L a momentum distribution by doing the ptweight
      getKinematics(H3Lq2b,M_H3L);
      double ptweight=getptweight(H3Lq2b);
      //    event1.SetDecay(H3Lq2b, 2, masses2); //lambda & d
      event1.SetDecay(H3Lq2b, 2, masses2_g); //lambda & d
      Double_t weight1 = event1.Generate();
      TLorentzVector *vLd= event1.GetDecay(0);
      TLorentzVector *d = event1.GetDecay(1);

      // cout << (*d).Pt()<< " "<<(*d).M() << endl;
      // event2.SetDecay(*vLd, 2, masses22);
      // double weight2 = event2.Generate();
      // TLorentzVector *p= event2.GetDecay(0);
      // TLorentzVector *pi = event2.GetDecay(1);
      // TLorentzVector ppi = *p + *pi;
      // TLorentzVector pd = *p + *d;
      // TLorentzVector dpi = *pi+*d;

      // double Tpi = ( sqrt(M_pi*M_pi+(*pi).P()*(*pi).P())-M_pi )*1000;
      // double Tp = ( (*p).P()*(*p).P()/2./M_p)*1000;
      // double Td = ( (*p).P()*(*p).P()/2./M_vd)* 1000; 

      //ckstar
      TLorentzVector* dr = new TLorentzVector();
      TLorentzVector* vLdr = new TLorentzVector();
      // in our data analysis we use real deuteron mass
      dr->SetXYZM(d->X(),d->Y(),d->Z(),M_d);
      vLdr->SetXYZM(vLd->X(),vLd->Y(),vLd->Z(),M_Ld);

      double dphi4pi, dphi, deta;
      dphi4pi=dr->Phi()-vLdr->Phi();
      dphi=dphi4pi;
      if(dphi4pi<-TMath::Pi()){dphi=dphi4pi+2*TMath::Pi();}
      if(dphi4pi>TMath::Pi()){dphi=dphi4pi-2*TMath::Pi();}
      deta=dr->Eta()-vLdr->Eta();
    
      TLorentzVector H = *dr+ *vLdr;
      TLorentzVector Qvect = (*vLdr-*dr);

      double corr_weight;
      double Pinv = H.Mag();
      double Q1 = (M_Ld*M_Ld-M_d*M_d)/Pinv;
      double Q=sqrt(Q1*Q1-Qvect.Mag2());
      double kstar = Q/2.0 * 1000; // convert to MeV

      // cout << Td+Tpi+Tp << endl;
      // cout << ( (*d).Vect() + (*pi).Vect() + (*p).Vect() ).Mag() << endl;
      // cout << ppi.M()<<" "<<pd.M()<<endl;
      // hppi_pd->Fill( ppi.M(), pd.M(),weight1*weight2); 
      // hppi_dpi->Fill( ppi.M(), dpi.M(),weight1*weight2); 
      // hpd_dpi->Fill( pd.M(), dpi.M(),weight1*weight2); 
      // hTpi_Td->Fill( Tpi, Td,weight1*weight2);

      //Add ckstar
      hkstar[imvld]->Fill(kstar,weight1*ptweight);
      hdphideta[imvld]->Fill(dphi,deta,weight1*ptweight);

      //free memory
      delete dr, d, vLdr, vLd, H, Qvect;
      
    } // entry loop

    cout << "The corresponding p: " <<  getPfromM(M_H3L,M_vLd_y,M_vd) << endl;
    hvldM[imvld]->Write();    
    hdphideta[imvld]->Write();
    hkstar[imvld]->Write();
  } //lambda mass loop

  // hppi_pd->Write();
  // hpd_dpi->Write();
  // hppi_dpi->Write();
  // hTpi_Td->Write();
  
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
  const Double_t M_H3L= 2.99089+(0.00041-0.00013); //binding energy = 0.13MeV
  Double_t masses2[2] = { M_vLd, M_vd} ;
  Double_t masses3[3] = { M_p, M_pi, M_d} ;
  Double_t masses22[2] = { M_p, M_pi} ;

  int const nevents = 1000;

  TGraph* g = new TGraph("xsection.csv");

  TGenPhaseSpace event1, event2,event3;
  TRandom3 *gRandom = new TRandom3();

  TH2F* hpd_dpi = new TH2F("hpd_dpi","hpd_dpi;M(pd) (GeV/c^{2});M(dpi) (GeV/c^{2})", 100, 2.8, 2.86, 100, 2.01, 2.06 );
  TH2F* hppi_dpi= new TH2F("hppi_dpi","hppi_dpi;M(p#pi) (GeV/c^{2});M(d#pi) (GeV/c^{2})", 100,1.05,1.15, 100, 2.01,2.06);
  TH2F* hppi_pd = new TH2F("hppi_pd","hppi_pd;M(p#pi) (GeV/c^{2});M(pd) (GeV/c^{2})", 100, 1.05,1.15, 100, 2.8, 2.86 );
  TH2F* hTpi_Td= new TH2F("hTpi_Td","hTpi_Td;Tpi (MeV);Td (MeV)", 100, 0, 36, 100, 0, 36);
  TH2F * hdphideta = new TH2F("hdphideta","hdphideta",200,-5,5,200,-5,5);

  
  for (int ie=0;ie<nevents;ie++){
      TLorentzVector H3Lf3b; 
      TLorentzVector H3Lq2b; 

      //generate free 3  body
      H3Lf3b.SetXYZM( 0,0,0,M_H3L);
      // bool isallow = event3.SetDecay(H3Lf3b, 3, masses3,"Fermi");
      bool isallow = event3.SetDecay(H3Lf3b, 3, masses3);
      if (isallow){
      double weight = event3.Generate();
      TLorentzVector *p= event3.GetDecay(0);
      TLorentzVector *pi = event3.GetDecay(1);
      TLorentzVector *d = event3.GetDecay(2);
      // // cout << (*d).Pt()<< " "<<(*d).M() << endl;

      TLorentzVector ppi = *p + *pi;
      TLorentzVector pd = *p + *d;
      TLorentzVector dpi = *pi+*d;

      double Tpi = ( sqrt(M_pi*M_pi+(*pi).P()*(*pi).P())-M_pi )*1000;
      double Tp = ( (*p).P()*(*p).P()/2./M_p)*1000;
      double Td = ( (*p).P()*(*p).P()/2./M_d)* 1000;
      // cout << Td+Tpi+Tp << endl;
      // cout << ( (*d).Vect() + (*pi).Vect() + (*p).Vect() ).Mag() << endl;

      TLorentzVector LA1 = ppi;
      TLorentzVector LA2 = *d;
      double dphi4pi, dphi, deta;
      dphi4pi=LA1.Phi()-LA2.Phi();
      dphi=dphi4pi;
      if(dphi4pi<=-TMath::Pi()){dphi=dphi4pi+2*TMath::Pi();}
      if(dphi4pi>=TMath::Pi()){dphi=dphi4pi-2*TMath::Pi();}
      //	deta=rapidity1-rapidity2;
      deta=LA1.Eta()-LA2.Eta();
      
      hdphideta->Fill(dphi,deta);

      
      // cout << ppi.M()<<" "<<pd.M()<<endl;
      double weight_gaus = TMath::Gaus(Tpi, 32, 1, true);
      if (Td<0.5) weight_gaus = weight_gaus* (g->Eval(Td));
      else weight_gaus = 0;
      // double weight_gaus = 1;
      // if ( Td<0.3 && Tpi<33 && Tpi>31 ){
        hppi_pd->Fill( ppi.M(), pd.M(),weight*weight_gaus); 
        hppi_dpi->Fill( ppi.M(), dpi.M(),weight*weight_gaus); 
        hpd_dpi->Fill( pd.M(), dpi.M(),weight*weight_gaus); 
        hTpi_Td->Fill( Tpi, Td, weight*weight_gaus);
      // }
    }
    else {cout << "decay not allowed" << endl;}
  }

  TFile* file = new TFile("H3L_Dz3b.root","recreate");
  hppi_pd->Write();
  hpd_dpi->Write();
  hppi_dpi->Write();
  hTpi_Td->Write();
  hdphideta->Write();
  
  time.Stop();
  time.Print();
}

void H3Ldecay()
{
  H3Ldecay_free3body();
  H3Ldecay_quasi2body();
}


