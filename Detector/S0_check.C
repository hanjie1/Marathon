#include "GetTrees.h"

void Check(int nrun)
{
     TChain *T;
     T=GetTree(nrun,"T");
     if(T==NULL)return;

     Int_t nentries=T->GetEntries();
     if(nentries<30000)return;

     Double_t ped_val[2]={0.0},ped_wid[2]={0.0},peak[2]={0.0},peak_wid[2]={0.0};
     
     delete gROOT->FindObject("hs0_lp");
     delete gROOT->FindObject("hs0_lc");
     delete gROOT->FindObject("hs0_rp");
     delete gROOT->FindObject("hs0_rc");

     TCanvas *c1=new TCanvas("c1");
     c1->Divide(2,1);
     c1->cd(1);
     TH1F *hs0_lp=new TH1F("hs0_lp","s0 la_p",550,-200,2000);
     T->Draw("L.s0.la_p>>hs0_lp");
     TF1 *f1=new TF1("f1","gaus",-100,100);
     hs0_lp->Fit(f1,"NQR");
     Double_t min = f1->GetParameter(1) - f1->GetParameter(2);
     Double_t max = f1->GetParameter(1) + f1->GetParameter(2);

     hs0_lp->Fit(f1,"Q","",min,max);
     ped_val[0] = f1->GetParameter(1);
     ped_wid[0] = f1->GetParameter(2);
     gPad->SetLogy();

     c1->cd(2);
     TH1F *hs0_rp=new TH1F("hs0_rp","s0 ra_p",550,-200,2000);
     T->Draw("L.s0.ra_p>>hs0_rp");
     TF1 *f2=new TF1("f2","gaus",-100,100);
     hs0_rp->Fit(f2,"NQR");
     min = f2->GetParameter(1) - f2->GetParameter(2);
     max = f2->GetParameter(1) + f2->GetParameter(2);

     hs0_rp->Fit(f2,"Q","",min,max);
     ped_val[1] = f2->GetParameter(1);
     ped_wid[1] = f2->GetParameter(2);
     gPad->SetLogy();

     int flag=0;
     Double_t mean1=0,mean2=0,mean3=0;
     Double_t sigma1=0,sigma2=0,sigma3=0;
     TCanvas *c2=new TCanvas("c2");
     c2->Divide(2,1);
     c2->cd(1);
     TH1F *hs0_lc=new TH1F("hs0_lc","s0 la_c",550,-200,2000);
     T->Draw("L.s0.la_c>>hs0_lc");
     Int_t tmp1=hs0_lc->FindBin(150);
     Double_t ninte=hs0_lc->Integral(tmp1,550);
     if(ninte<1000)return;
     TF1 *f3=new TF1("f3","landau",200,700);
     hs0_lc->Fit(f3,"QR");

     mean1 = f3->GetParameter(1);
     sigma1 = f3->GetParameter(2);
     if(mean1<200||mean1>700){
        peak[0]=mean1;
        peak_wid[0]=sigma1;
        flag=1;
        gPad->SetLogy();
        goto end1;
     }

     hs0_lc->Fit(f3,"NQ","",mean1-150,mean1+200);
     mean2 = f3->GetParameter(1);
     sigma2 = f3->GetParameter(2);
     if(mean2<200||mean2>700){
        peak[0]=mean1;
        peak_wid[0]=sigma1;
        flag=1;
        hs0_lc->Fit(f3,"Q","",mean1-150,mean1+200);
        gPad->SetLogy();
        goto end1;
     }

     hs0_lc->Fit(f3,"NQ","",mean2-150,mean2+200);
     mean3 = f3->GetParameter(1);
     sigma3 = f3->GetParameter(2);
     if(mean3<200||mean3>700){
        peak[0]=mean2;
        peak_wid[0]=sigma2;
        flag=1;
        hs0_lc->Fit(f3,"Q","",mean2-150,mean2+200);
        gPad->SetLogy();
        goto end1;
     }

     peak[0] = mean3;
     peak_wid[0] = sigma3;
     gPad->SetLogy();

     end1:
     c2->cd(2);
     TH1F *hs0_rc=new TH1F("hs0_rc","s0 ra_c",550,-200,2000);
     T->Draw("L.s0.ra_c>>hs0_rc");
     TF1 *f4=new TF1("f4","landau",200,700);
     hs0_rc->Fit(f4,"QR");

     mean1 = f4->GetParameter(1);
     sigma1 = f4->GetParameter(2);
     if(mean1<200||mean1>700){
        peak[1]=mean1;
        peak_wid[1]=sigma1;
        flag=1;
        gPad->SetLogy();
        goto end;
     }

     hs0_rc->Fit(f4,"NQ","",mean1-150,mean1+200);
     mean2 = f4->GetParameter(1);
     sigma2 = f4->GetParameter(2);
     if(mean2<200||mean2>700){
        peak[1]=mean1;
        peak_wid[1]=sigma1;
        flag=1;
        hs0_rc->Fit(f4,"Q","",mean1-150,mean1+200);
        gPad->SetLogy();
        goto end;
     }

     hs0_rc->Fit(f4,"NQ","",mean2-150,mean2+200);
     mean3 = f4->GetParameter(1);
     sigma3 = f4->GetParameter(2);
     if(mean3<200||mean3>700){
        peak[1]=mean2;
        peak_wid[1]=sigma2;
        flag=1;
        hs0_rc->Fit(f4,"Q","",mean2-150,mean2+200);
        gPad->SetLogy();
        goto end;
     }

     peak[1] = mean3;
     peak_wid[1] = sigma3;
     gPad->SetLogy();

    end:

    ofstream outfile1,outfile2,outfile3,outfile4;
    outfile1.open("OUT/Ls0_ped.dat",fstream::app);
    outfile2.open("OUT/Ls0_ped_wid.dat",fstream::app);
    outfile3.open("OUT/Ls0_peak.dat",fstream::app);
    outfile4.open("OUT/Ls0_peak_wid.dat",fstream::app);
    outfile1<<nrun<<"  ";
    outfile2<<nrun<<"  ";
    outfile3<<nrun<<"  ";
    outfile4<<nrun<<"  ";
    for(int ii=0;ii<2;ii++){
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
    c1->Print(Form("PDF/s0_%d.pdf[",nrun));
    c1->Print(Form("PDF/s0_%d.pdf",nrun));
    c2->Print(Form("PDF/s0_%d.pdf",nrun));
    c2->Print(Form("PDF/s0_%d.pdf]",nrun));

    return;

}

void S0_check(){
     for(int ii=1206;ii<2826;ii++){
	 Check(ii);
     }
     exit(0);
}

