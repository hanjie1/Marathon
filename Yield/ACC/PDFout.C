void PDFout()
{
        TFile *f1=new TFile("Nu_Theta_all.root");

	TH2F *hH1[5];
	TH2F *hD2[11];
	TH2F *hD2_2nd[4];
	TH2F *hHe[11];
	TH2F *hHe_2nd[4];
	TH2F *hH3[11];
	TH2F *hH3_2nd[4];

	int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
	TString target[4]={"H1","D2","He3","H3"};

	int maxkin=11;
        for(int ii=0;ii<4;ii++){
	    if(ii==0)maxkin=5;
	    else maxkin=11;
	    for(int jj=0;jj<maxkin;jj++){
		if(ii==0) hH1[jj]=(TH2F *)f1->Get(Form("H1_kin%d",kin[jj]));
		if(ii==1){
		   if(kin[jj]<7 || kin[jj]==15)hD2[jj]=(TH2F *)f1->Get(Form("D2_kin%d",kin[jj]));
		   else{
			hD2[jj]=(TH2F *)f1->Get(Form("D2_kin%d_1st",kin[jj]));
			hD2_2nd[jj-6]=(TH2F *)f1->Get(Form("D2_kin%d_2nd",kin[jj]));
		   }
		}

		if(ii==2){
		   if(kin[jj]<7 || kin[jj]==15)hHe[jj]=(TH2F *)f1->Get(Form("He3_kin%d",kin[jj]));
		   else{
			hHe[jj]=(TH2F *)f1->Get(Form("He3_kin%d_1st",kin[jj]));
			hHe_2nd[jj-6]=(TH2F *)f1->Get(Form("He3_kin%d_2nd",kin[jj]));
		   }
		}

		if(ii==3){
		   if(kin[jj]<7 || kin[jj]==15)hH3[jj]=(TH2F *)f1->Get(Form("H3_kin%d",kin[jj]));
		   else{
			hH3[jj]=(TH2F *)f1->Get(Form("H3_kin%d_1st",kin[jj]));
			hH3_2nd[jj-6]=(TH2F *)f1->Get(Form("H3_kin%d_2nd",kin[jj]));
		   }
		}
	    }
	}

        gStyle->SetOptStat(1111111);
        TCanvas *c1=new TCanvas("c1","c1",1500,1500);
	c1->Divide(3,2);
        for(int ii=0;ii<5;ii++){
	    c1->cd(ii+1);
	    hH1[ii]->Draw("COLZ");
	}
 
        TCanvas *c2=new TCanvas("c2","c2",1500,1500);
	c2->Divide(3,2);
        for(int ii=0;ii<6;ii++){
	    c2->cd(ii+1);
	    hD2[ii]->Draw("COLZ");
	}

        TCanvas *c3=new TCanvas("c3","c3",1500,1500);
	c3->Divide(3,2);
	int kk=0;
        for(int ii=6;ii<9;ii++){
	    c3->cd(kk+1);
	    hD2[ii]->Draw("COLZ");
	    kk++;
	    c3->cd(kk+1);
	    hD2_2nd[ii-6]->Draw("COLZ");
	    kk++;
	}

	TCanvas *c4=new TCanvas("c4","c4",1500,1500);	
	c4->Divide(2,2);
	c4->cd(1);
	hD2[9]->Draw("COLZ");
	c4->cd(2);
	hD2_2nd[3]->Draw("COLZ");
	c4->cd(3);
	hD2[10]->Draw("COLZ");

        TCanvas *c5=new TCanvas("c5","c5",1500,1500);
	c5->Divide(3,2);
        for(int ii=0;ii<6;ii++){
	    c5->cd(ii+1);
	    hHe[ii]->Draw("COLZ");
	}

        TCanvas *c6=new TCanvas("c6","c6",1500,1500);
	c6->Divide(3,2);
	kk=0;
        for(int ii=6;ii<9;ii++){
	    c6->cd(kk+1);
	    hHe[ii]->Draw("COLZ");
	    kk++;
	    c6->cd(kk+1);
	    hHe_2nd[ii-6]->Draw("COLZ");
	    kk++;
	}

	TCanvas *c7=new TCanvas("c7","c7",1500,1500);	
	c7->Divide(2,2);
	c7->cd(1);
	hHe[9]->Draw("COLZ");
	c7->cd(2);
	hHe_2nd[3]->Draw("COLZ");
	c7->cd(3);
	hHe[10]->Draw("COLZ");

        TCanvas *c8=new TCanvas("c8","c8",1500,1500);
	c8->Divide(3,2);
        for(int ii=0;ii<6;ii++){
	    c8->cd(ii+1);
	    hH3[ii]->Draw("COLZ");
	}

        TCanvas *c9=new TCanvas("c9","c9",1500,1500);
	c9->Divide(3,2);
	kk=0;
        for(int ii=6;ii<9;ii++){
	    c9->cd(kk+1);
	    hH3[ii]->Draw("COLZ");
	    kk++;
	    c9->cd(kk+1);
	    hH3_2nd[ii-6]->Draw("COLZ");
	    kk++;
	}

	TCanvas *c10=new TCanvas("c10","c10",1500,1500);	
	c10->Divide(2,2);
	c10->cd(1);
	hH3[9]->Draw("COLZ");
	c10->cd(2);
	hH3_2nd[3]->Draw("COLZ");
	c10->cd(3);
	hH3[10]->Draw("COLZ");

	c1->Print("Nu_Theta.pdf[");	
	c1->Print("Nu_Theta.pdf");	
	c2->Print("Nu_Theta.pdf");	
	c3->Print("Nu_Theta.pdf");	
	c4->Print("Nu_Theta.pdf");	
	c5->Print("Nu_Theta.pdf");	
	c6->Print("Nu_Theta.pdf");	
	c7->Print("Nu_Theta.pdf");	
	c8->Print("Nu_Theta.pdf");	
	c9->Print("Nu_Theta.pdf");	
	c10->Print("Nu_Theta.pdf");	
	c10->Print("Nu_Theta.pdf]");	

}
