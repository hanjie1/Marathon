#include "GetTrees.h"

void Check(int nrun)
{
     TChain *T;
     T=GetTree(nrun,"T");
     if(T==NULL)return;

     Int_t nentries=T->GetEntries();
     if(nentries<30000)return;

     Double_t ped_val[16]={0.0},ped_wid[16]={0.0},peak[16]={0.0},peak_wid[16]={0.0};
     
    for(int ii=0;ii<16;ii++){
        delete gROOT->FindObject(Form("hs2_rp[%d]",ii));
        delete gROOT->FindObject(Form("hs2_rc[%d]",ii));
    }
     TH1F *hs2_lp[16];
     TCanvas *c1=new TCanvas("c1");
     c1->Divide(4,4);
     for(int ii=0;ii<16;ii++){
         c1->cd(ii+1);
         hs2_lp[ii]=new TH1F(Form("hs2_rp[%d]",ii),"s2 ra_p",300,-200,1000);
	 T->Draw(Form("R.s2.ra_p[%d]>>hs2_rp[%d]",ii,ii));
         Int_t max_bin = hs2_lp[ii]->GetMaximumBin();
         Double_t bin_center = hs2_lp[ii]->GetBinCenter(max_bin);

         TF1 *g1=new TF1("g1","gaus",bin_center-20,bin_center+20);  

         hs2_lp[ii]->Fit(g1,"NQR");
         Double_t min = g1->GetParameter(1) - g1->GetParameter(2);
         Double_t max = g1->GetParameter(1) + g1->GetParameter(2);

         hs2_lp[ii]->Fit(g1,"Q","",min,max);
         ped_val[ii] = g1->GetParameter(1);
         ped_wid[ii] = g1->GetParameter(2);
         gPad->SetLogy();
    }  

     TH1F *hs2_lc[16];
     TCanvas *c2=new TCanvas("c2");
     c2->Divide(4,4);
     int flag=0;
     for(int ii=0;ii<16;ii++){
         c2->cd(ii+1);
         hs2_lc[ii]=new TH1F(Form("hs2_rc[%d]",ii),"s2R a_c",300,-200,800);
         T->Draw(Form("R.s2.ra_c[%d]>>hs2_rc[%d]",ii,ii));
         Int_t tmp1=hs2_lc[ii]->FindBin(130);
         Double_t ninte=hs2_lc[ii]->Integral(tmp1,300);
         if(ninte<200)return;

         TF1 *g1=new TF1("g1","landau",130,350);

         hs2_lc[ii]->Fit(g1,"QR");
         Double_t mean1 = g1->GetParameter(1);
         Double_t sigma1 = g1->GetParameter(2);
         if(mean1<120||mean1>350){
                peak[ii]=mean1;
                peak_wid[ii]=sigma1;
                flag=1;
                gPad->SetLogy();
                continue;
         }

         hs2_lc[ii]->Fit(g1,"NQ","",mean1-40,mean1+100);
         Double_t mean2 = g1->GetParameter(1);
         Double_t sigma2 = g1->GetParameter(2);
         if(mean2<120||mean2>350){
                peak[ii]=mean1;
                peak_wid[ii]=sigma1;
                flag=1;
                hs2_lc[ii]->Fit(g1,"Q","",120,350);
                gPad->SetLogy();
                continue;
         }

         hs2_lc[ii]->Fit(g1,"Q","",mean2-40,mean2+100);
         Double_t mean3 = g1->GetParameter(1);
         Double_t sigma3 = g1->GetParameter(2);
         if(mean3<120||mean3>350){
                peak[ii]=mean2;
                peak_wid[ii]=sigma2;
                flag=1;
                hs2_lc[ii]->Fit(g1,"Q","",mean1-40,mean1+100);
                gPad->SetLogy();
                continue;
         }

         peak[ii] = mean3;
         peak_wid[ii] = sigma3;
         gPad->SetLogy();
    }

    ofstream outfile1,outfile2,outfile3,outfile4;
    outfile1.open("OUT/Rs2R_ped.dat",fstream::app);
    outfile2.open("OUT/Rs2R_ped_wid.dat",fstream::app);
    outfile3.open("OUT/Rs2R_peak.dat",fstream::app);
    outfile4.open("OUT/Rs2R_peak_wid.dat",fstream::app);
    outfile1<<nrun<<"  ";
    outfile2<<nrun<<"  ";
    outfile3<<nrun<<"  ";
    outfile4<<nrun<<"  ";
    for(int ii=0;ii<16;ii++){
      outfile1<<fixed<<setprecision(2)<<ped_val[ii]<<"  ";
      outfile2<<fixed<<setprecision(2)<<ped_wid[ii]<<"  ";
      outfile3<<fixed<<setprecision(2)<<peak[ii]<<"  ";
      outfile4<<fixed<<setprecision(2)<<peak_wid[ii]<<"  ";
    }
    outfile1<<endl;
    outfile2<<endl;
    outfile3<<flag<<endl;
    outfile4<<flag<<endl;

    outfile1.close();
    outfile2.close();
    outfile3.close();
    outfile4.close();
    c1->Print(Form("PDF/Rs2R_%d.pdf[",nrun));
    c1->Print(Form("PDF/Rs2R_%d.pdf",nrun));
    c2->Print(Form("PDF/Rs2R_%d.pdf",nrun));
    c2->Print(Form("PDF/Rs2R_%d.pdf]",nrun));

    return;

}

void RS2R_check(){
     for(int ii=90837;ii<92020;ii++){
	 Check(ii);
     }
     exit(0);
}

