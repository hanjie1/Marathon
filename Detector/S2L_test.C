#include "GetTrees.h"

void S2L_test(int nrun)
{
     TChain *T;
     T=GetTree(nrun,"T");
     if(T==NULL)return;

     Int_t nentries=T->GetEntries();
     if(nentries<30000)return;

     Double_t ped_val[16]={0.0},ped_wid[16]={0.0},peak[16]={0.0},peak_wid[16]={0.0};
     
    for(int ii=0;ii<16;ii++){
        delete gROOT->FindObject(Form("hs2_lc[%d]",ii));
    }
     TH1F *hs2_lc[16];
     TCanvas *c2=new TCanvas("c2");
     c2->Divide(4,4);
     int flag=0;
     for(int ii=0;ii<16;ii++){
         c2->cd(ii+1);
         hs2_lc[ii]=new TH1F(Form("hs2_lc[%d]",ii),"s2L a_c",300,-200,1000);
         T->Draw(Form("L.s2.la_c[%d]>>hs2_lc[%d]",ii,ii));
         Int_t tmp1=hs2_lc[ii]->FindBin(150);
         Double_t ninte=hs2_lc[ii]->Integral(tmp1,300);
         if(ninte<100)return;

         TF1 *g1=new TF1("g1","landau",200,500);

         hs2_lc[ii]->Fit(g1,"QR");
         Double_t mean1 = g1->GetParameter(1);
         Double_t sigma1 = g1->GetParameter(2);
         if(mean1<200||mean1>500){
                peak[ii]=mean1;
                peak_wid[ii]=sigma1;
                flag=1;
                gPad->SetLogy();
                continue;
         }

         hs2_lc[ii]->Fit(g1,"NQ","",mean1-100,mean1+150);
         Double_t mean2 = g1->GetParameter(1);
         Double_t sigma2 = g1->GetParameter(2);
         if(mean2<200||mean2>500){
                peak[ii]=mean1;
                peak_wid[ii]=sigma1;
                flag=1;
                hs2_lc[ii]->Fit(g1,"Q","",200,500);
                gPad->SetLogy();
                continue;
         }

         hs2_lc[ii]->Fit(g1,"Q","",mean2-100,mean2+150);
         Double_t mean3 = g1->GetParameter(1);
         Double_t sigma3 = g1->GetParameter(2);
         if(mean3<200||mean3>500){
                peak[ii]=mean2;
                peak_wid[ii]=sigma2;
                flag=1;
                hs2_lc[ii]->Fit(g1,"Q","",mean1-100,mean1+150);
                gPad->SetLogy();
                continue;
         }

         peak[ii] = mean3;
         peak_wid[ii] = sigma3;
         gPad->SetLogy();

    cout<<peak[ii]<<"  "<<peak_wid[ii]<<"  "<<flag<<endl;

    }



}

