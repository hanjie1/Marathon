void plot_allphys()
{
     TFile *f1=new TFile("Phys_all.root");
     TString target[4]={"H1","D2","He3","H3"};
     int kin[11]={0,1,2,3,4,5,7,9,11,13,15};

     TH1F *hH1[5];
     TH1F *hD2[11];
     TH1F *hHe3[11];
     TH1F *hH3[11];
     for(int ii=0;ii<11;ii++){
         if(ii<5){
            hH1[ii]=(TH1F *)f1->Get(Form("%s_kin%d_dp","H1",kin[ii]));
            hH1[ii]->Rebin(10);     
         }
         hD2[ii]=(TH1F *)f1->Get(Form("%s_kin%d_dp","D2",kin[ii]));
         hHe3[ii]=(TH1F *)f1->Get(Form("%s_kin%d_dp","He3",kin[ii]));
         hH3[ii]=(TH1F *)f1->Get(Form("%s_kin%d_dp","H3",kin[ii]));
     }

    int color[11]={1,2,3,4,6,7,8,9,11,38,46};

    TCanvas *c1=new TCanvas("c1");
    THStack *hsH1=new THStack("hsH1","Theta distribution");
    for(int ii=0;ii<5;ii++){
        hsH1->Add(hH1[ii]);
        hH1[ii]->SetLineColor(color[ii]);
    }
    hsH1->Draw("nostack");

    auto leg=new TLegend(0.7,0.6,0.85,0.85);
    for(int ii=0;ii<5;ii++)
       leg->AddEntry(hH1[ii],Form("%s kin%d","H1",kin[ii]),"L");
    leg->Draw();

    TCanvas *c2=new TCanvas("c2");
    THStack *hsD2=new THStack("hsD2","dp distribution");
    for(int ii=0;ii<11;ii++){
        hsD2->Add(hD2[ii]);
        hD2[ii]->SetLineColor(color[ii]);
    }
    hsD2->Draw("nostack");

    auto leg1=new TLegend(0.7,0.6,0.85,0.85);
    for(int ii=0;ii<11;ii++)
       leg1->AddEntry(hD2[ii],Form("%s kin%d","D2",kin[ii]),"L");
    leg1->Draw();

    TCanvas *c3=new TCanvas("c3");
    THStack *hsH3=new THStack("hsH3","dp distribution");
    for(int ii=0;ii<11;ii++){
        hsH3->Add(hH3[ii]);
        hH3[ii]->SetLineColor(color[ii]);
    }
    hsH3->Draw("nostack");

    auto leg2=new TLegend(0.7,0.6,0.85,0.85);
    for(int ii=0;ii<11;ii++)
       leg2->AddEntry(hH3[ii],Form("%s kin%d","H3",kin[ii]),"L");
    leg2->Draw();

    TCanvas *c4=new TCanvas("c4");
    THStack *hsHe3=new THStack("hsHe3","dp distribution");
    for(int ii=0;ii<11;ii++){
        hsHe3->Add(hHe3[ii]);
        hHe3[ii]->SetLineColor(color[ii]);
    }
    hsHe3->Draw("nostack");

    auto leg3=new TLegend(0.7,0.6,0.85,0.85);
    for(int ii=0;ii<11;ii++)
       leg3->AddEntry(hD2[ii],Form("%s kin%d","D2",kin[ii]),"L");
    leg3->Draw();



}
