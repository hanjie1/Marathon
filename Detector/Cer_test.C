#include "GetTrees.h"

void Cer_test(int nrun)
{
     TChain *T;
     T=GetTree(nrun,"T");
     if(T==NULL)return;

     Int_t nentries=T->GetEntries();
     if(nentries<30000)return;

     Double_t ped_val[10]={0.0},ped_wid[10]={0.0},peak[10]={0.0},peak_wid[10]={0.0};
     
    for(int ii=0;ii<10;ii++){
        delete gROOT->FindObject(Form("hcer_p[%d]",ii));
        delete gROOT->FindObject(Form("hcer_c[%d]",ii));
    }
     TH1F *hcer_p[10];
/*     TCanvas *c1=new TCanvas("c1");
     c1->Divide(5,2);
     for(int ii=0;ii<10;ii++){
         c1->cd(ii+1);
         hcer_p[ii]=new TH1F(Form("hcer_p[%d]",ii),"cer a_p",300,-200,1000);
	 T->Draw(Form("L.cer.a_p[%d]>>hcer_p[%d]",ii,ii));
         Int_t max_bin = hcer_p[ii]->GetMaximumBin();
         Double_t bin_center = hcer_p[ii]->GetBinCenter(max_bin);

         TF1 *g1=new TF1("g1","gaus",bin_center-20,bin_center+20);  

         hcer_p[ii]->Fit(g1,"NQR");
         Double_t min = g1->GetParameter(1) - 2*g1->GetParameter(2);
         Double_t max = g1->GetParameter(1) + 2*g1->GetParameter(2);

         hcer_p[ii]->Fit(g1,"Q","",min,max);
         ped_val[ii] = g1->GetParameter(1);
         ped_wid[ii] = g1->GetParameter(2);
         cout<<ped_val[ii]<<"  "<<ped_wid[ii]<<endl;
         gPad->SetLogy();
    }  
*/

     TH1F *hcer_c[10];
     TCanvas *c2=new TCanvas("c2");
     c2->Divide(5,2);
     int flag=0;
     for(int ii=0;ii<10;ii++){
         c2->cd(ii+1);
         hcer_c[ii]=new TH1F(Form("hcer_c[%d]",ii),"cer a_c",300,-200,1000);
	 T->Draw(Form("L.cer.a_c[%d]>>hcer_c[%d]",ii,ii));
         Int_t tmp1=hcer_c[ii]->FindBin(150);
         Double_t ninte=hcer_c[ii]->Integral(tmp1,300);
         cout<<"!!!:  "<<ninte<<endl;
         if(ninte<200)return;

         TF1 *g1=new TF1("g1","gaus",150,450);  

         hcer_c[ii]->Fit(g1,"NQR");
         Double_t mean1 = g1->GetParameter(1);
         Double_t sigma1 = g1->GetParameter(2);
         if(mean1<150||mean1>450){
		peak[ii]=mean1;
                peak_wid[ii]=sigma1;
                flag=1;
                gPad->SetLogy();
                continue;
         }

         hcer_c[ii]->Fit(g1,"Q","",mean1-sigma1,mean1+sigma1);
         Double_t mean2 = g1->GetParameter(1);
         Double_t sigma2 = g1->GetParameter(2);
         if(mean2<150||mean2>450){
		peak[ii]=mean1;
                peak_wid[ii]=sigma1;
                flag=1;
                hcer_c[ii]->Fit(g1,"Q","",150,450);
                gPad->SetLogy();
                continue;
         }

         hcer_c[ii]->Fit(g1,"NQ","",mean2-sigma2,mean2+sigma2);
         Double_t mean3 = g1->GetParameter(1);
         Double_t sigma3 = g1->GetParameter(2);
         if(mean3<150||mean3>450){
		peak[ii]=mean2;
                peak_wid[ii]=sigma2;
                flag=1;
                hcer_c[ii]->Fit(g1,"Q","",mean1-sigma1,mean1+sigma1);
                gPad->SetLogy();
                continue;
         }
         
         peak[ii] = mean3;
         peak_wid[ii] = sigma3;
         cout<<peak[ii]<<"  "<<peak_wid[ii]<<"  "<<flag<<endl;

      //   cout<<mean2<<"  "<<mean3<<endl;
         gPad->SetLogy();
    }  
}

