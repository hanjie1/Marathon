#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <cassert>
#include "GetFile.h"

using namespace std;
void Multi_tracks()
{
  int getTree=0;
  int arm=-1; //LHRS: arm=0; RHRS: arm=1;
  TString kin;
  TChain* T=GetTree(arm,kin);

  if(arm==-1){cout<<"Don't know which arm"<<endl;exit(0);}
  int LHRS=0,RHRS=0;
  if(arm==0)LHRS=1;
  if(arm==1)RHRS=1;

  Double_t th_max=0.06;
  Double_t th_min=-0.05;
  Double_t ph_max=0.03;
  Double_t ph_min=-0.03;
  Double_t dp_max=0.04;
  Double_t dp_min=-0.04;
  Double_t vz_max=0.1;
  Double_t vz_min=-0.1;

  TCut TRK1,TRK2,Ep,CK,ACC,VZ,trigger2,trigger1,beta;
  if(LHRS){
     beta = "L.tr.beta>0";
     Ep = "(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>0.75";
     ACC = Form("L.gold.th>%f && L.gold.th<%f && L.gold.ph>%f && L.gold.ph<%f && L.gold.dp>%f && L.gold.dp<%f",th_min,th_max,ph_min,ph_max,dp_min,dp_max);
     VZ = Form("L.tr.vz[0]>%f && L.tr.vz[0]<%f",vz_min,vz_max);
     CK = "L.cer.asum_c>2000";
     trigger2 = "(DL.evtypebits>>2)&1";
     trigger1 = "(DL.evtypebits>>1)&1";
     TRK1 = "L.tr.n==1";
     TRK2 = "L.tr.n>1";
   }
  if(RHRS){
     beta = "R.tr.beta>0";
     Ep = "(R.ps.e+R.sh.e)/(1000*R.gold.p)>0.75";
     ACC = Form("R.gold.th>%f && R.gold.th<%f && R.gold.ph>%f && R.gold.ph<%f && R.gold.dp>%f && R.gold.dp<%f",th_min,th_max,ph_min,ph_max,dp_min,dp_max);
     VZ = Form("R.tr.vz[0]>%f && R.tr.vz[0]<%f",vz_min,vz_max);
     CK = "R.cer.asum_c>2000";
     trigger2 = "(DR.evtypebits>>5)&1";
     trigger1 = "(DR.evtypebits>>4)&1";
     TRK1 = "R.tr.n==1";
     TRK2 = "R.tr.n>1";
   }

   Double_t ngood_ele = T->GetEntries(beta+Ep+ACC+CK+VZ+trigger2);
   Double_t ntrack1 = T->GetEntries(beta+Ep+ACC+CK+VZ+trigger2+TRK1);
   Double_t ntrack2 = T->GetEntries(beta+Ep+ACC+CK+VZ+trigger2+TRK2);

   cout<<"-----------------------------"<<endl;
   cout<<"one track events:  "<<ntrack1/ngood_ele<<endl;
   cout<<"two or more track events:  "<<ntrack2/ngood_ele<<endl;
   cout<<"-----------------------------"<<endl;

   TCanvas *c1=new TCanvas("c1","c1",1500,1500);
   TH1F *htrk = new TH1F("htrk","htrk",);


}
