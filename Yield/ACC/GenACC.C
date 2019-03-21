void GenACC()
{
        TFile *f1=new TFile("Nu_Theta_all.root");

	TH2F *hH1[5];
	TH2F *hD2[11];
	TH2F *hD2_2nd[4];
	TH2F *hHe[11];
	TH2F *hHe_2nd[4];
	TH2F *hH3[11];
	TH2F *hH3_2nd[4];
/*
	TH2F *hH1_new[5];
	TH2F *hD2_new[11];
	TH2F *hD2_new_2nd[4];
	TH2F *hHe_new[11];
	TH2F *hHe_new_2nd[4];
	TH2F *hH3_new[11];
	TH2F *hH3_new_2nd[4];
*/
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

	int nxbin=800/8+2;
	int nybin=400/10+2;
        for(int ii=0;ii<5;ii++){
	    hH1[ii]->Rebin2D(8,8);
            Int_t maxbin = hH1[ii]->GetMaximumBin();
            Double_t maxcontent = hH1[ii]->GetBinContent(maxbin);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hH1[ii]->GetBin(mm,nn);
                   Double_t content=hH1[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                      cout<<"There is underflow or overflow!"<<endl;
                      continue;
                   }

                   Double_t new_con=content/maxcontent;
                   new_con=(int)(new_con/0.001+0.5)*0.001;
                   new_con=new_con*100;
                   hH1[ii]->SetBinContent(nbin,new_con);
               }
            }
	}

        for(int ii=0;ii<11;ii++){
	    hD2[ii]->Rebin2D(8,8);
            Int_t maxbin = hD2[ii]->GetMaximumBin();
            Double_t maxcontent = hD2[ii]->GetBinContent(maxbin);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hD2[ii]->GetBin(mm,nn);
                   Double_t content=hD2[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                      cout<<"There is underflow or overflow!"<<endl;
                      continue;
                   }

                   Double_t new_con=content/maxcontent;
                   new_con=(int)(new_con/0.001+0.5)*0.001;
                   new_con=new_con*100;
                   hD2[ii]->SetBinContent(nbin,new_con);
               }
            }
	}

        for(int ii=0;ii<4;ii++){
	    hD2_2nd[ii]->Rebin2D(8,8);
            Int_t maxbin = hD2_2nd[ii]->GetMaximumBin();
            Double_t maxcontent = hD2_2nd[ii]->GetBinContent(maxbin);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hD2_2nd[ii]->GetBin(mm,nn);
                   Double_t content=hD2_2nd[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                      cout<<"There is underflow or overflow!"<<endl;
                      continue;
                   }

                   Double_t new_con=content/maxcontent;
                   new_con=(int)(new_con/0.001+0.5)*0.001;
                   new_con=new_con*100;
                   hD2_2nd[ii]->SetBinContent(nbin,new_con);
               }
            }
	}


        for(int ii=0;ii<11;ii++){
	    hHe[ii]->Rebin2D(8,8);
            Int_t maxbin = hHe[ii]->GetMaximumBin();
            Double_t maxcontent = hHe[ii]->GetBinContent(maxbin);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hHe[ii]->GetBin(mm,nn);
                   Double_t content=hHe[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                      cout<<"There is underflow or overflow!"<<endl;
                      continue;
                   }

                   Double_t new_con=content/maxcontent;
                   new_con=(int)(new_con/0.001+0.5)*0.001;
                   new_con=new_con*100;
                   hHe[ii]->SetBinContent(nbin,new_con);
               }
            }
	}

        for(int ii=0;ii<4;ii++){
	    hHe_2nd[ii]->Rebin2D(8,8);
            Int_t maxbin = hHe_2nd[ii]->GetMaximumBin();
            Double_t maxcontent = hHe_2nd[ii]->GetBinContent(maxbin);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hHe_2nd[ii]->GetBin(mm,nn);
                   Double_t content=hHe_2nd[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                      cout<<"There is underflow or overflow!"<<endl;
                      continue;
                   }

                   Double_t new_con=content/maxcontent;
                   new_con=(int)(new_con/0.001+0.5)*0.001;
                   new_con=new_con*100;
                   hHe_2nd[ii]->SetBinContent(nbin,new_con);
               }
            }
	}

        for(int ii=0;ii<11;ii++){
	    hH3[ii]->Rebin2D(8,8);
            Int_t maxbin = hH3[ii]->GetMaximumBin();
            Double_t maxcontent = hH3[ii]->GetBinContent(maxbin);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hH3[ii]->GetBin(mm,nn);
                   Double_t content=hH3[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                      cout<<"There is underflow or overflow!"<<endl;
                      continue;
                   }

                   Double_t new_con=content/maxcontent;
                   new_con=(int)(new_con/0.001+0.5)*0.001;
                   new_con=new_con*100;
                   hH3[ii]->SetBinContent(nbin,new_con);
               }
            }
	}

        for(int ii=0;ii<4;ii++){
	    hH3_2nd[ii]->Rebin2D(8,8);
            Int_t maxbin = hH3_2nd[ii]->GetMaximumBin();
            Double_t maxcontent = hH3_2nd[ii]->GetBinContent(maxbin);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hH3_2nd[ii]->GetBin(mm,nn);
                   Double_t content=hH3_2nd[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                      cout<<"There is underflow or overflow!"<<endl;
                      continue;
                   }

                   Double_t new_con=content/maxcontent;
                   new_con=(int)(new_con/0.001+0.5)*0.001;
                   new_con=new_con*100;
                   hH3_2nd[ii]->SetBinContent(nbin,new_con);
               }
            }
	}


        gStyle->SetOptStat(1111111);
        TCanvas *c1=new TCanvas("c1","c1",1500,1500);
	c1->Divide(3,2);
        for(int ii=0;ii<5;ii++){
	    c1->cd(ii+1);
	    hH1[ii]->Draw("lego");
	}
 
        TCanvas *c2=new TCanvas("c2","c2",1500,1500);
	c2->Divide(3,2);
        for(int ii=0;ii<6;ii++){
	    c2->cd(ii+1);
	    hD2[ii]->Draw("lego");
	}

        TCanvas *c3=new TCanvas("c3","c3",1500,1500);
	c3->Divide(3,2);
	int kk=0;
        for(int ii=6;ii<9;ii++){
	    c3->cd(kk+1);
	    hD2[ii]->Draw("lego");
	    kk++;
	    c3->cd(kk+1);
	    hD2_2nd[ii-6]->Draw("lego");
	    kk++;
	}

	TCanvas *c4=new TCanvas("c4","c4",1500,1500);	
	c4->Divide(2,2);
	c4->cd(1);
	hD2[9]->Draw("lego");
	c4->cd(2);
	hD2_2nd[3]->Draw("lego");
	c4->cd(3);
	hD2[10]->Draw("lego");

        TCanvas *c5=new TCanvas("c5","c5",1500,1500);
	c5->Divide(3,2);
        for(int ii=0;ii<6;ii++){
	    c5->cd(ii+1);
	    hHe[ii]->Draw("lego");
	}

        TCanvas *c6=new TCanvas("c6","c6",1500,1500);
	c6->Divide(3,2);
	kk=0;
        for(int ii=6;ii<9;ii++){
	    c6->cd(kk+1);
	    hHe[ii]->Draw("lego");
	    kk++;
	    c6->cd(kk+1);
	    hHe_2nd[ii-6]->Draw("lego");
	    kk++;
	}

	TCanvas *c7=new TCanvas("c7","c7",1500,1500);	
	c7->Divide(2,2);
	c7->cd(1);
	hHe[9]->Draw("lego");
	c7->cd(2);
	hHe_2nd[3]->Draw("lego");
	c7->cd(3);
	hHe[10]->Draw("lego");

        TCanvas *c8=new TCanvas("c8","c8",1500,1500);
	c8->Divide(3,2);
        for(int ii=0;ii<6;ii++){
	    c8->cd(ii+1);
	    hH3[ii]->Draw("lego");
	}

        TCanvas *c9=new TCanvas("c9","c9",1500,1500);
	c9->Divide(3,2);
	kk=0;
        for(int ii=6;ii<9;ii++){
	    c9->cd(kk+1);
	    hH3[ii]->Draw("lego");
	    kk++;
	    c9->cd(kk+1);
	    hH3_2nd[ii-6]->Draw("lego");
	    kk++;
	}

	TCanvas *c10=new TCanvas("c10","c10",1500,1500);	
	c10->Divide(2,2);
	c10->cd(1);
	hH3[9]->Draw("lego");
	c10->cd(2);
	hH3_2nd[3]->Draw("lego");
	c10->cd(3);
	hH3[10]->Draw("lego");

	c1->Print("ACC.pdf[");	
	c1->Print("ACC.pdf");	
	c2->Print("ACC.pdf");	
	c3->Print("ACC.pdf");	
	c4->Print("ACC.pdf");	
	c5->Print("ACC.pdf");	
	c6->Print("ACC.pdf");	
	c7->Print("ACC.pdf");	
	c8->Print("ACC.pdf");	
	c9->Print("ACC.pdf");	
	c10->Print("ACC.pdf");	
	c10->Print("ACC.pdf]");	

}
