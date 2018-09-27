void plot_all()
{
     TFile *f1=new TFile("Xbj_new.root");
     int kin[11]={0,1,2,3,4,5,7,9,11,13,15};

     TH1F *hH1[5];
     TH1F *hD2[11];
     TH1F *hHe3[11];
     TH1F *hH3[11];
     for(int ii=0;ii<11;ii++){
	 if(ii<5){
            hH1[ii]=(TH1F *)f1->Get(Form("%s_kin%d","H1",kin[ii]));
//	    hH1[ii]->Rebin(20);
         }
	 hD2[ii]=(TH1F *)f1->Get(Form("%s_kin%d","D2",kin[ii]));
//	 hD2[ii]->Rebin(20);
	 hHe3[ii]=(TH1F *)f1->Get(Form("%s_kin%d","He3",kin[ii]));
//	 hHe3[ii]->Rebin(20);
	 hH3[ii]=(TH1F *)f1->Get(Form("%s_kin%d","H3",kin[ii]));
//	 hH3[ii]->Rebin(20);
     }
     Double_t H1Fbin[5]={0.0},H1Lbin[5]={0.0};
     Double_t D2Fbin[11]={0.0},D2Lbin[11]={0.0};
     Double_t HeFbin[11]={0.0},HeLbin[11]={0.0};
     Double_t H3Fbin[11]={0.0},H3Lbin[11]={0.0};

     int color[11]={1,2,3,4,5,6,7,8,41,28,46};
     TCanvas *c1=new TCanvas("c1");
     for(int ii=0;ii<5;ii++){
	if(ii==0)hH1[ii]->Draw();
        else hH1[ii]->Draw("same");
	hH1[ii]->SetLineColor(color[ii]);
        Double_t tmp_max=hH1[ii]->GetBinContent(hH1[ii]->GetMaximumBin());
        Double_t tmp_th=tmp_max*0.25;
        H1Fbin[ii]=hH1[ii]->GetBinLowEdge(hH1[ii]->FindFirstBinAbove(tmp_th));
        H1Lbin[ii]=hH1[ii]->GetBinLowEdge(hH1[ii]->FindLastBinAbove(tmp_th));
     }

     TCanvas *c2=new TCanvas("c2");
     for(int ii=0;ii<11;ii++){
        if(ii==0)hD2[ii]->Draw();
        else hD2[ii]->Draw("same");
	hD2[ii]->SetLineColor(color[ii]);
        Double_t tmp_max=hD2[ii]->GetBinContent(hD2[ii]->GetMaximumBin());
        Double_t tmp_th=tmp_max*0.25;
        D2Fbin[ii]=hD2[ii]->GetBinLowEdge(hD2[ii]->FindFirstBinAbove(tmp_th));
        D2Lbin[ii]=hD2[ii]->GetBinLowEdge(hD2[ii]->FindLastBinAbove(tmp_th));
     }

     TCanvas *c3=new TCanvas("c3");
     for(int ii=0;ii<11;ii++){
        if(ii==0)hHe3[ii]->Draw();
        else hHe3[ii]->Draw("same");
	hHe3[ii]->SetLineColor(color[ii]);
        Double_t tmp_max=hHe3[ii]->GetBinContent(hHe3[ii]->GetMaximumBin());
        Double_t tmp_th=tmp_max*0.25;
        HeFbin[ii]=hHe3[ii]->GetBinLowEdge(hHe3[ii]->FindFirstBinAbove(tmp_th));
        HeLbin[ii]=hHe3[ii]->GetBinLowEdge(hHe3[ii]->FindLastBinAbove(tmp_th));
     }
     TCanvas *c4=new TCanvas("c4");
     for(int ii=0;ii<11;ii++){
        if(ii==0)hH3[ii]->Draw();
        else hH3[ii]->Draw("same");
	hH3[ii]->SetLineColor(color[ii]);
        Double_t tmp_max=hH3[ii]->GetBinContent(hH3[ii]->GetMaximumBin());
        Double_t tmp_th=tmp_max*0.25;
        H3Fbin[ii]=hH3[ii]->GetBinLowEdge(hH3[ii]->FindFirstBinAbove(tmp_th));
        H3Lbin[ii]=hH3[ii]->GetBinLowEdge(hH3[ii]->FindLastBinAbove(tmp_th));
     }
  
    ofstream myfile;
    myfile.open("Xrange_new_25per.txt");
    myfile<<"---------- H1 ----------"<<endl;
    for(int ii=0;ii<5;ii++){
        //double tmp_f=(H1Fbin[ii]-1)*0.02;
        myfile<<fixed<<setprecision(3);
        myfile<<H1Fbin[ii]<<",";
    }
    myfile<<endl;
    for(int ii=0;ii<5;ii++){
        //double tmp_l=(H1Lbin[ii]-1)*0.02;
        myfile<<H1Lbin[ii]<<",";
    }
    myfile<<endl;

    myfile<<"---------- D2 ----------"<<endl;
    for(int ii=0;ii<11;ii++){
       // double tmp_f=(D2Fbin[ii]-1)*0.02;
        myfile<<D2Fbin[ii]<<",";
    }
    myfile<<endl;
    for(int ii=0;ii<11;ii++){
        //double tmp_l=(D2Lbin[ii]-1)*0.02;
        myfile<<D2Lbin[ii]<<",";
    }
    myfile<<endl;

    myfile<<"---------- He ----------"<<endl;
    for(int ii=0;ii<11;ii++){
        //double tmp_f=(HeFbin[ii]-1)*0.02;
        myfile<<HeFbin[ii]<<",";
    }
    myfile<<endl;
    for(int ii=0;ii<11;ii++){
        //double tmp_l=(HeLbin[ii]-1)*0.02;
        myfile<<HeLbin[ii]<<",";
    }
    myfile<<endl;
    myfile<<"---------- H3 ----------"<<endl;
    for(int ii=0;ii<11;ii++){
        //double tmp_f=(H3Fbin[ii]-1)*0.02;
        myfile<<H3Fbin[ii]<<",";
    }
    myfile<<endl;
    for(int ii=0;ii<11;ii++){
        //double tmp_l=(H3Lbin[ii]-1)*0.02;
        myfile<<H3Lbin[ii]<<",";
    }
    myfile<<endl;
    myfile.close();

}
