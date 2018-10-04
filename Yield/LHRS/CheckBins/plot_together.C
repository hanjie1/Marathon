void plot_together()
{
     TFile *f1=new TFile("Xbj_new.root");
     TString target[4]={"H1","D2","He3","H3"};
     int kin[11]={0,1,2,3,4,5,7,9,11,13,15};

     TH1F *hH1[5];
     TH1F *hD2[11];
     TH1F *hHe3[11];
     TH1F *hH3[11];
     for(int ii=0;ii<11;ii++){
         if(ii<5){
            hH1[ii]=(TH1F *)f1->Get(Form("%s_kin%d","H1",kin[ii]));
         }
         hD2[ii]=(TH1F *)f1->Get(Form("%s_kin%d","D2",kin[ii]));
         hHe3[ii]=(TH1F *)f1->Get(Form("%s_kin%d","He3",kin[ii]));
         hH3[ii]=(TH1F *)f1->Get(Form("%s_kin%d","H3",kin[ii]));
     }

   //  int color[5]={1,2,3,4,6,7,8,9,11};

     



     THStack *hs[11];
     for(int ii=0;ii<11;ii++){
	TString hname = Form("kin%d",kin[ii]);
        hs[ii]=new THStack(hname.Data(),"xbj distribution");
        if(ii<5){
           hs[ii]->Add(hH1[ii]);
           hH1[ii]->SetLineColor(1);
	}
        hD2[ii]->SetLineColor(2);
	hs[ii]->Add(hD2[ii]);
        hH3[ii]->SetLineColor(4);
	hs[ii]->Add(hH3[ii]);
        hHe3[ii]->SetLineColor(3);
	hs[ii]->Add(hHe3[ii]);

         Double_t inte1=hD2[ii]->Integral();
         Double_t inte2=0;
         if(ii<5)inte2=hH1[ii]->Integral();
         Double_t inte3=hH3[ii]->Integral();
         Double_t inte4=hHe3[ii]->Integral();
         cout<<inte1<<"  "<<inte2<<"  "<<inte3<<"  "<<inte4<<endl;
         if(ii<5)hH1[ii]->Scale(inte1/inte2);
         hH3[ii]->Scale(inte1/inte3);
         hHe3[ii]->Scale(inte1/inte4);


     }

     TCanvas *c1[11];
     for(int ii=0;ii<11;ii++){     
         c1[ii]=new TCanvas(Form("c%d",ii+1));
         hs[ii]->Draw("nostack");

         auto leg=new TLegend(0.7,0.6,0.85,0.85);
         if(ii<5)leg->AddEntry(hH1[ii],Form("%s kin%d",target[0].Data(),kin[ii]),"L");
         leg->AddEntry(hD2[ii],Form("%s kin%d",target[1].Data(),kin[ii]),"L");
         leg->AddEntry(hHe3[ii],Form("%s kin%d",target[2].Data(),kin[ii]),"L");
         leg->AddEntry(hH3[ii],Form("%s kin%d",target[3].Data(),kin[ii]),"L");
         leg->Draw();
     }

     for(int ii=0;ii<11;ii++){
	 if(ii==0)c1[ii]->Print("All_xbj.pdf[");
	 c1[ii]->Print("All_xbj.pdf");
	 if(ii==10)c1[ii]->Print("All_xbj.pdf]");
     }


}
