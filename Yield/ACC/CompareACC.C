void CompareACC()
{
        TFile *f1=new TFile("ACC_rmXS.root");
        TFile *f2=new TFile("ACC_ratio.root","RECREATE");

        TH2F *hH1[5];
        TH2F *hD2[11];
        TH2F *hD2_2nd[4];
        TH2F *hHe[11];
        TH2F *hHe_2nd[4];
        TH2F *hH3[11];
        TH2F *hH3_2nd[4];

	TH2F *hDp[5];
	TH2F *hHeD[11];
	TH2F *hHeD_2nd[4];
	TH2F *hH3D[11];
	TH2F *hH3D_2nd[4];
	TH2F *hH3He[11];
	TH2F *hH3He_2nd[4];

        int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
        TString target[4]={"H1","D2","He3","H3"};

        int maxkin=11;
        for(int ii=0;ii<4;ii++){
            if(ii==0)maxkin=5;
            else maxkin=11;
            for(int jj=0;jj<maxkin;jj++){
                if(ii==0) hH1[jj]=(TH2F *)f1->Get(Form("H1_kin%d_new",kin[jj]));
                if(ii==1){
                   if(kin[jj]<7 || kin[jj]==15)hD2[jj]=(TH2F *)f1->Get(Form("D2_kin%d_new",kin[jj]));
                   else{
                        hD2[jj]=(TH2F *)f1->Get(Form("D2_kin%d_1st_new",kin[jj]));
                        hD2_2nd[jj-6]=(TH2F *)f1->Get(Form("D2_kin%d_2nd_new",kin[jj]));
                   }
                }

                if(ii==2){
                   if(kin[jj]<7 || kin[jj]==15)hHe[jj]=(TH2F *)f1->Get(Form("He3_kin%d_new",kin[jj]));
                   else{
                        hHe[jj]=(TH2F *)f1->Get(Form("He3_kin%d_1st_new",kin[jj]));
                        hHe_2nd[jj-6]=(TH2F *)f1->Get(Form("He3_kin%d_2nd_new",kin[jj]));
                   }
                }

                if(ii==3){
                   if(kin[jj]<7 || kin[jj]==15)hH3[jj]=(TH2F *)f1->Get(Form("H3_kin%d_new",kin[jj]));
                   else{
                        hH3[jj]=(TH2F *)f1->Get(Form("H3_kin%d_1st_new",kin[jj]));
                        hH3_2nd[jj-6]=(TH2F *)f1->Get(Form("H3_kin%d_2nd_new",kin[jj]));
                   }
                }
            }
        }

	int rbin=8;
        int nxbin=800/rbin+2;
        int nybin=400/rbin+2;

        for(int ii=0;ii<5;ii++){
	    hH1[ii]->Rebin2D(rbin,rbin);
	    hD2[ii]->Rebin2D(rbin,rbin);

	    hDp[ii]=new TH2F(Form("Dp_kin%d",kin[ii]),"D/p ACC ratio",800/rbin,-4.0,4.0,400/rbin,-0.2,0.2);
            Int_t maxbin1 = hH1[ii]->GetMaximumBin();
            Double_t maxcontent1 = hH1[ii]->GetBinContent(maxbin1);
            Int_t maxbin2 = hD2[ii]->GetMaximumBin();
            Double_t maxcontent2 = hD2[ii]->GetBinContent(maxbin2);

            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                  Double_t nbin1=hH1[ii]->GetBin(mm,nn);
                  Double_t content1=hH1[ii]->GetBinContent(nbin1);
                  Double_t nbin2=hD2[ii]->GetBin(mm,nn);
                  Double_t content2=hD2[ii]->GetBinContent(nbin2);

                  if(nbin1 != nbin2){cout<<"!!! Wrong bin !!!"<<endl; continue;}
                  if(content1==0||content2==0){
                     hDp[ii]->SetBinContent(nbin1,0);
                     continue;
                  }

                  Double_t new_con1=content1/maxcontent1;
                  Double_t new_con2=content2/maxcontent2;
                  Double_t ratio=new_con2/new_con1;
		 
		  Double_t H1_err=hH1[ii]->GetBinError(nbin1);
		  Double_t D2_err=hD2[ii]->GetBinError(nbin2);
 		  Double_t ratio_err=ratio*sqrt((H1_err/content1)*(H1_err/content1)+(D2_err/content2)*(D2_err/content2));

                  ratio=(int)(ratio/0.01+0.5)*0.01;
                  hDp[ii]->SetBinContent(nbin1,ratio);
                  hDp[ii]->SetBinError(nbin1,ratio_err);
             }
            }
	}

        for(int ii=0;ii<11;ii++){
	    hHe[ii]->Rebin2D(rbin,rbin);
            if(ii>4)hD2[ii]->Rebin2D(rbin,rbin);

	    if(kin[ii]<7||kin[ii]==15)	    	
		hHeD[ii]=new TH2F(Form("HeD_kin%d",kin[ii]),"He/D ACC ratio",800/rbin,-4.0,4.0,400/rbin,-0.2,0.2);
	    else
		hHeD[ii]=new TH2F(Form("HeD_kin%d_1st",kin[ii]),"He/D ACC ratio",800/rbin,-4.0,4.0,400/rbin,-0.2,0.2);

            Int_t maxbin1 = hD2[ii]->GetMaximumBin();
            Double_t maxcontent1 = hD2[ii]->GetBinContent(maxbin1);
            Int_t maxbin2 = hHe[ii]->GetMaximumBin();
            Double_t maxcontent2 = hHe[ii]->GetBinContent(maxbin2);

            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                  Double_t nbin1=hD2[ii]->GetBin(mm,nn);
                  Double_t content1=hD2[ii]->GetBinContent(nbin1);
                  Double_t nbin2=hHe[ii]->GetBin(mm,nn);
                  Double_t content2=hHe[ii]->GetBinContent(nbin2);

                  if(nbin1 != nbin2){cout<<"!!! Wrong bin !!!"<<endl; continue;}
                  if(content1==0||content2==0){
                     hHeD[ii]->SetBinContent(nbin1,0);
                     continue;
                  }

                  Double_t new_con1=content1/maxcontent1;
                  Double_t new_con2=content2/maxcontent2;
                  Double_t ratio=new_con2/new_con1;

		  Double_t D2_err=hD2[ii]->GetBinError(nbin1);
		  Double_t He_err=hHe[ii]->GetBinError(nbin2);
 		  Double_t ratio_err=ratio*sqrt((D2_err/content1)*(D2_err/content1)+(He_err/content2)*(He_err/content2));

                  ratio=(int)(ratio/0.01+0.5)*0.01;
                  hHeD[ii]->SetBinContent(nbin1,ratio);
                  hHeD[ii]->SetBinError(nbin1,ratio_err);
             }
            }
	}

        for(int ii=0;ii<4;ii++){
	    hD2_2nd[ii]->Rebin2D(rbin,rbin);
	    hHe_2nd[ii]->Rebin2D(rbin,rbin);
	    hHeD_2nd[ii]=new TH2F(Form("HeD_kin%d_2nd",kin[ii+6]),"He/D ACC ratio",800/rbin,-4.0,4.0,400/rbin,-0.2,0.2);

            Int_t maxbin1 = hD2_2nd[ii]->GetMaximumBin();
            Double_t maxcontent1 = hD2_2nd[ii]->GetBinContent(maxbin1);
            Int_t maxbin2 = hHe_2nd[ii]->GetMaximumBin();
            Double_t maxcontent2 = hHe_2nd[ii]->GetBinContent(maxbin2);

            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                  Double_t nbin1=hD2_2nd[ii]->GetBin(mm,nn);
                  Double_t content1=hD2_2nd[ii]->GetBinContent(nbin1);
                  Double_t nbin2=hHe_2nd[ii]->GetBin(mm,nn);
                  Double_t content2=hHe_2nd[ii]->GetBinContent(nbin2);

                  if(nbin1 != nbin2){cout<<"!!! Wrong bin !!!"<<endl; continue;}
                  if(content1==0||content2==0){
                     hHeD_2nd[ii]->SetBinContent(nbin1,0);
                     continue;
                  }

                  Double_t new_con1=content1/maxcontent1;
                  Double_t new_con2=content2/maxcontent2;
                  Double_t ratio=new_con2/new_con1;

		  Double_t D2_err=hD2_2nd[ii]->GetBinError(nbin1);
		  Double_t He_err=hHe_2nd[ii]->GetBinError(nbin2);
 		  Double_t ratio_err=ratio*sqrt((D2_err/content1)*(D2_err/content1)+(He_err/content2)*(He_err/content2));

                  ratio=(int)(ratio/0.01+0.5)*0.01;
                  hHeD_2nd[ii]->SetBinContent(nbin1,ratio);
                  hHeD_2nd[ii]->SetBinError(nbin1,ratio_err);
             }
            }
	}


        for(int ii=0;ii<11;ii++){
	    hH3[ii]->Rebin2D(rbin,rbin);

	    if(kin[ii]<7||kin[ii]==15)	    	
		hH3D[ii]=new TH2F(Form("H3D_kin%d",kin[ii]),"H3/D ACC ratio",800/rbin,-4.0,4.0,400/rbin,-0.2,0.2);
	    else
		hH3D[ii]=new TH2F(Form("H3D_kin%d_1st",kin[ii]),"H3/D ACC ratio",800/rbin,-4.0,4.0,400/rbin,-0.2,0.2);

            Int_t maxbin1 = hD2[ii]->GetMaximumBin();
            Double_t maxcontent1 = hD2[ii]->GetBinContent(maxbin1);
            Int_t maxbin2 = hH3[ii]->GetMaximumBin();
            Double_t maxcontent2 = hH3[ii]->GetBinContent(maxbin2);

            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                  Double_t nbin1=hD2[ii]->GetBin(mm,nn);
                  Double_t content1=hD2[ii]->GetBinContent(nbin1);
                  Double_t nbin2=hH3[ii]->GetBin(mm,nn);
                  Double_t content2=hH3[ii]->GetBinContent(nbin2);

                  if(nbin1 != nbin2){cout<<"!!! Wrong bin !!!"<<endl; continue;}
                  if(content1==0||content2==0){
                     hH3D[ii]->SetBinContent(nbin1,0);
                     continue;
                  }

                  Double_t new_con1=content1/maxcontent1;
                  Double_t new_con2=content2/maxcontent2;
                  Double_t ratio=new_con2/new_con1;

		  Double_t D2_err=hD2[ii]->GetBinError(nbin1);
		  Double_t H3_err=hH3[ii]->GetBinError(nbin2);
 		  Double_t ratio_err=ratio*sqrt((D2_err/content1)*(D2_err/content1)+(H3_err/content2)*(H3_err/content2));

                  ratio=(int)(ratio/0.01+0.5)*0.01;
                  hH3D[ii]->SetBinContent(nbin1,ratio);
                  hH3D[ii]->SetBinError(nbin1,ratio_err);
             }
            }
	}

        for(int ii=0;ii<4;ii++){
	    hH3_2nd[ii]->Rebin2D(rbin,rbin);
	    hH3D_2nd[ii]=new TH2F(Form("H3D_kin%d_2nd",kin[ii+6]),"H3/D ACC ratio",800/rbin,-4.0,4.0,400/rbin,-0.2,0.2);

            Int_t maxbin1 = hD2_2nd[ii]->GetMaximumBin();
            Double_t maxcontent1 = hD2_2nd[ii]->GetBinContent(maxbin1);
            Int_t maxbin2 = hH3_2nd[ii]->GetMaximumBin();
            Double_t maxcontent2 = hH3_2nd[ii]->GetBinContent(maxbin2);

            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                  Double_t nbin1=hD2_2nd[ii]->GetBin(mm,nn);
                  Double_t content1=hD2_2nd[ii]->GetBinContent(nbin1);
                  Double_t nbin2=hH3_2nd[ii]->GetBin(mm,nn);
                  Double_t content2=hH3_2nd[ii]->GetBinContent(nbin2);

                  if(nbin1 != nbin2){cout<<"!!! Wrong bin !!!"<<endl; continue;}
                  if(content1==0||content2==0){
                     hH3D_2nd[ii]->SetBinContent(nbin1,0);
                     continue;
                  }

                  Double_t new_con1=content1/maxcontent1;
                  Double_t new_con2=content2/maxcontent2;
                  Double_t ratio=new_con2/new_con1;

		  Double_t D2_err=hD2_2nd[ii]->GetBinError(nbin1);
		  Double_t H3_err=hH3_2nd[ii]->GetBinError(nbin2);
 		  Double_t ratio_err=ratio*sqrt((D2_err/content1)*(D2_err/content1)+(H3_err/content2)*(H3_err/content2));

                  ratio=(int)(ratio/0.01+0.5)*0.01;
                  hH3D_2nd[ii]->SetBinContent(nbin1,ratio);
                  hH3D_2nd[ii]->SetBinError(nbin1,ratio_err);
             }
            }
	}

        for(int ii=0;ii<11;ii++){

	    if(kin[ii]<7||kin[ii]==15)	    	
		hH3He[ii]=new TH2F(Form("H3He_kin%d",kin[ii]),"H3/He ACC ratio",800/rbin,-4.0,4.0,400/rbin,-0.2,0.2);
	    else
		hH3He[ii]=new TH2F(Form("H3He_kin%d_1st",kin[ii]),"H3/He ACC ratio",800/rbin,-4.0,4.0,400/rbin,-0.2,0.2);

            Int_t maxbin1 = hHe[ii]->GetMaximumBin();
            Double_t maxcontent1 = hHe[ii]->GetBinContent(maxbin1);
            Int_t maxbin2 = hH3[ii]->GetMaximumBin();
            Double_t maxcontent2 = hH3[ii]->GetBinContent(maxbin2);

            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                  Double_t nbin1=hHe[ii]->GetBin(mm,nn);
                  Double_t content1=hHe[ii]->GetBinContent(nbin1);
                  Double_t nbin2=hH3[ii]->GetBin(mm,nn);
                  Double_t content2=hH3[ii]->GetBinContent(nbin2);

                  if(nbin1 != nbin2){cout<<"!!! Wrong bin !!!"<<endl; continue;}
                  if(content1==0||content2==0){
                     hH3He[ii]->SetBinContent(nbin1,0);
                     continue;
                  }

                  Double_t new_con1=content1/maxcontent1;
                  Double_t new_con2=content2/maxcontent2;
                  Double_t ratio=new_con2/new_con1;

		  Double_t He_err=hHe[ii]->GetBinError(nbin1);
		  Double_t H3_err=hH3[ii]->GetBinError(nbin2);
 		  Double_t ratio_err=ratio*sqrt((He_err/content1)*(He_err/content1)+(H3_err/content2)*(H3_err/content2));

                  ratio=(int)(ratio/0.01+0.5)*0.01;
                  hH3He[ii]->SetBinContent(nbin1,ratio);
                  hH3He[ii]->SetBinError(nbin1,ratio_err);
             }
            }
	}

        for(int ii=0;ii<4;ii++){
	    hH3He_2nd[ii]=new TH2F(Form("H3He_kin%d_2nd",kin[ii+6]),"H3/He ACC ratio",800/rbin,-4.0,4.0,400/rbin,-0.2,0.2);

            Int_t maxbin1 = hHe_2nd[ii]->GetMaximumBin();
            Double_t maxcontent1 = hHe_2nd[ii]->GetBinContent(maxbin1);
            Int_t maxbin2 = hH3_2nd[ii]->GetMaximumBin();
            Double_t maxcontent2 = hH3_2nd[ii]->GetBinContent(maxbin2);

            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                  Double_t nbin1=hHe_2nd[ii]->GetBin(mm,nn);
                  Double_t content1=hHe_2nd[ii]->GetBinContent(nbin1);
                  Double_t nbin2=hH3_2nd[ii]->GetBin(mm,nn);
                  Double_t content2=hH3_2nd[ii]->GetBinContent(nbin2);

                  if(nbin1 != nbin2){cout<<"!!! Wrong bin !!!"<<endl; continue;}
                  if(content1==0||content2==0){
                     hH3He_2nd[ii]->SetBinContent(nbin1,0);
                     continue;
                  }

                  Double_t new_con1=content1/maxcontent1;
                  Double_t new_con2=content2/maxcontent2;
                  Double_t ratio=new_con2/new_con1;

		  Double_t He_err=hHe_2nd[ii]->GetBinError(nbin1);
		  Double_t H3_err=hH3_2nd[ii]->GetBinError(nbin2);
 		  Double_t ratio_err=ratio*sqrt((He_err/content1)*(He_err/content1)+(H3_err/content2)*(H3_err/content2));

                  ratio=(int)(ratio/0.01+0.5)*0.01;
                  hH3He_2nd[ii]->SetBinContent(nbin1,ratio);
                  hH3He_2nd[ii]->SetBinError(nbin1,ratio_err);
             }
            }
	}

        f2->Write();



}
