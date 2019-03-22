void GenXSTable()
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

        Double_t Nu_c=7.49;
        Double_t Theta_c[11]={16.8075,17.5717,19.1125,20.575,21.9401,23.2065,25.5858,27.7642,29.8087,31.7274,33.5552};
        Double_t Theta_c1[4]={25.5909,27.7744,29.8159,31.7325};//2nd run of kin 7,9,11,13; 2nd and 3rd kin15 are almost the same as 1st run, so use same center theta value;

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

	int rbin=8;
        Int_t nxbin = 800/rbin+2;
        Int_t nybin = 400/rbin+2;

        Double_t E0=10.589;
        Double_t dXBin=0.01*rbin;
        Double_t dYBin=0.001*rbin;
	for(int ii=0;ii<5;ii++){
	    ofstream outfile;
	    outfile.open(Form("Kin_table/H1_kin%d.inp",kin[ii]));

            Double_t th_min=100.0,th_max=0.0;
            Double_t Ep_min=100.0,Ep_max=0.0;

	    hH1[ii]->Rebin2D(rbin,rbin);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hH1[ii]->GetBin(mm,nn);
                   Double_t content=hH1[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                     cout<<"There is underflow or overflow!"<<endl;
                     continue;
                   }
                   Double_t theta=-4.0+(mm-1)*dXBin+Theta_c[ii];
                   Double_t nu=-0.2+(nn-1)*dYBin+Nu_c;
                   Double_t Ep=E0-nu;
                   if(theta<th_min)th_min=theta;
                   if(theta>th_max)th_max=theta;
                   if(Ep<Ep_min)Ep_min=Ep;
                   if(Ep>Ep_max)Ep_max=Ep;
               }
            }
cout<<"H1 kin"<<kin[ii]<<": "<<th_min<<" "<<th_max<<" "<<Ep_min<<" "<<Ep_max<<endl;
            outfile<<"Marathon"<<endl;
            outfile<<Form("H1_kin%d",kin[ii])<<endl;
            outfile<<endl;
            outfile<<endl;
            outfile<<endl;
            outfile<<"E    Ep    theta"<<endl;

            th_max=th_max+dXBin;
            Ep_min=Ep_min-dYBin;

            Double_t XXbin=(th_max-th_min)/dXBin+1;
            Double_t YYbin=(Ep_max-Ep_min)/dYBin+1;

            for(int ii=0;ii<XXbin;ii++){
               Double_t theta = th_min+dXBin*ii;
               for(int jj=0;jj<YYbin;jj++){
                  Double_t Ep=Ep_max-dYBin*jj;
                  outfile<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep<<" "<<theta<<endl;
               }
            }
            outfile.close();
	}	


	for(int ii=0;ii<11;ii++){
	    ofstream outfile;
	    if(kin[ii]<7||kin[ii]==15)outfile.open(Form("Kin_table/D2_kin%d.inp",kin[ii]));
	    else outfile.open(Form("Kin_table/D2_kin%d_1st.inp",kin[ii]));

            Double_t th_min=100.0,th_max=0.0;
            Double_t Ep_min=100.0,Ep_max=0.0;

	    hD2[ii]->Rebin2D(rbin,rbin);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hD2[ii]->GetBin(mm,nn);
                   Double_t content=hD2[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                     cout<<"There is underflow or overflow!"<<endl;
                     continue;
                   }
                   Double_t theta=-4.0+(mm-1)*dXBin+Theta_c[ii];
                   Double_t nu=-0.2+(nn-1)*dYBin+Nu_c;
                   Double_t Ep=E0-nu;
                   if(theta<th_min)th_min=theta;
                   if(theta>th_max)th_max=theta;
                   if(Ep<Ep_min)Ep_min=Ep;
                   if(Ep>Ep_max)Ep_max=Ep;
               }
            }
cout<<"D2 kin"<<kin[ii]<<": "<<th_min<<" "<<th_max<<" "<<Ep_min<<" "<<Ep_max<<endl;
            outfile<<"Marathon"<<endl;
            if(kin[ii]<7||kin[ii]==15)outfile<<Form("D2_kin%d",kin[ii])<<endl;
            else outfile<<Form("D2_kin%d_1st",kin[ii])<<endl;
            outfile<<endl;
            outfile<<endl;
            outfile<<endl;
            outfile<<"E    Ep    theta"<<endl;

            th_max=th_max+dXBin;
            Ep_min=Ep_min-dYBin;

            Double_t XXbin=(th_max-th_min)/dXBin+1;
            Double_t YYbin=(Ep_max-Ep_min)/dYBin+1;

            for(int ii=0;ii<XXbin;ii++){
               Double_t theta = th_min+dXBin*ii;
               for(int jj=0;jj<YYbin;jj++){
                  Double_t Ep=Ep_max-dYBin*jj;
                  outfile<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep<<" "<<theta<<endl;
               }
            }
            outfile.close();
	}	

	for(int ii=0;ii<11;ii++){
	    ofstream outfile;
	    if(kin[ii]<7 || kin[ii]==15)outfile.open(Form("Kin_table/He3_kin%d.inp",kin[ii]));
	    else outfile.open(Form("Kin_table/He3_kin%d_1st.inp",kin[ii]));

            Double_t th_min=100.0,th_max=0.0;
            Double_t Ep_min=100.0,Ep_max=0.0;

	    hHe[ii]->Rebin2D(rbin,rbin);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hHe[ii]->GetBin(mm,nn);
                   Double_t content=hHe[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                     cout<<"There is underflow or overflow!"<<endl;
                     continue;
                   }
                   Double_t theta=-4.0+(mm-1)*dXBin+Theta_c[ii];
                   Double_t nu=-0.2+(nn-1)*dYBin+Nu_c;
                   Double_t Ep=E0-nu;
                   if(theta<th_min)th_min=theta;
                   if(theta>th_max)th_max=theta;
                   if(Ep<Ep_min)Ep_min=Ep;
                   if(Ep>Ep_max)Ep_max=Ep;
               }
            }
            outfile<<"Marathon"<<endl;
            if(kin[ii]<7||kin[ii]==15)outfile<<Form("He3_kin%d",kin[ii])<<endl;
            else outfile<<Form("He3_kin%d_1st",kin[ii])<<endl;
            outfile<<endl;
            outfile<<endl;
            outfile<<endl;
            outfile<<"E    Ep    theta"<<endl;

            th_max=th_max+dXBin;
            Ep_min=Ep_min-dYBin;

            Double_t XXbin=(th_max-th_min)/dXBin+1;
            Double_t YYbin=(Ep_max-Ep_min)/dYBin+1;

            for(int ii=0;ii<XXbin;ii++){
               Double_t theta = th_min+dXBin*ii;
               for(int jj=0;jj<YYbin;jj++){
                  Double_t Ep=Ep_max-dYBin*jj;
                  outfile<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep<<" "<<theta<<endl;
               }
            }
            outfile.close();
	}	

	for(int ii=0;ii<11;ii++){
	    ofstream outfile;
	    if(kin[ii]<7||kin[ii]==15)outfile.open(Form("Kin_table/H3_kin%d.inp",kin[ii]));
	    else outfile.open(Form("Kin_table/H3_kin%d_1st.inp",kin[ii]));

            Double_t th_min=100.0,th_max=0.0;
            Double_t Ep_min=100.0,Ep_max=0.0;

	    hH3[ii]->Rebin2D(rbin,rbin);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hH3[ii]->GetBin(mm,nn);
                   Double_t content=hH3[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                     cout<<"There is underflow or overflow!"<<endl;
                     continue;
                   }
                   Double_t theta=-4.0+(mm-1)*dXBin+Theta_c[ii];
                   Double_t nu=-0.2+(nn-1)*dYBin+Nu_c;
                   Double_t Ep=E0-nu;
                   if(theta<th_min)th_min=theta;
                   if(theta>th_max)th_max=theta;
                   if(Ep<Ep_min)Ep_min=Ep;
                   if(Ep>Ep_max)Ep_max=Ep;
               }
            }
            outfile<<"Marathon"<<endl;
            if(kin[ii]<7||kin[ii]==15)outfile<<Form("H3_kin%d",kin[ii])<<endl;
            else outfile<<Form("H3_kin%d_1st",kin[ii])<<endl;
            outfile<<endl;
            outfile<<endl;
            outfile<<endl;
            outfile<<"E    Ep    theta"<<endl;

            th_max=th_max+dXBin;
            Ep_min=Ep_min-dYBin;

            Double_t XXbin=(th_max-th_min)/dXBin+1;
            Double_t YYbin=(Ep_max-Ep_min)/dYBin+1;

            for(int ii=0;ii<XXbin;ii++){
               Double_t theta = th_min+dXBin*ii;
               for(int jj=0;jj<YYbin;jj++){
                  Double_t Ep=Ep_max-dYBin*jj;
                  outfile<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep<<" "<<theta<<endl;
               }
            }
            outfile.close();
	}	


	for(int ii=0;ii<4;ii++){
	    ofstream outfile;
	    outfile.open(Form("Kin_table/D2_kin%d_2nd.inp",kin[ii+6]));

            Double_t th_min=100.0,th_max=0.0;
            Double_t Ep_min=100.0,Ep_max=0.0;

	    hD2_2nd[ii]->Rebin2D(rbin,rbin);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hD2_2nd[ii]->GetBin(mm,nn);
                   Double_t content=hD2_2nd[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                     cout<<"There is underflow or overflow!"<<endl;
                     continue;
                   }
                   Double_t theta=-4.0+(mm-1)*dXBin+Theta_c1[ii];
                   Double_t nu=-0.2+(nn-1)*dYBin+Nu_c;
                   Double_t Ep=E0-nu;
                   if(theta<th_min)th_min=theta;
                   if(theta>th_max)th_max=theta;
                   if(Ep<Ep_min)Ep_min=Ep;
                   if(Ep>Ep_max)Ep_max=Ep;
               }
            }
            outfile<<"Marathon"<<endl;
            outfile<<Form("D2_kin%d_2nd",kin[ii+6])<<endl;
            outfile<<endl;
            outfile<<endl;
            outfile<<endl;
            outfile<<"E    Ep    theta"<<endl;

            th_max=th_max+dXBin;
            Ep_min=Ep_min-dYBin;

            Double_t XXbin=(th_max-th_min)/dXBin+1;
            Double_t YYbin=(Ep_max-Ep_min)/dYBin+1;

            for(int ii=0;ii<XXbin;ii++){
               Double_t theta = th_min+dXBin*ii;
               for(int jj=0;jj<YYbin;jj++){
                  Double_t Ep=Ep_max-dYBin*jj;
                  outfile<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep<<" "<<theta<<endl;
               }
            }
            outfile.close();
	}	

	for(int ii=0;ii<4;ii++){
	    ofstream outfile;
	    outfile.open(Form("Kin_table/He3_kin%d_2nd.inp",kin[ii+6]));

            Double_t th_min=100.0,th_max=0.0;
            Double_t Ep_min=100.0,Ep_max=0.0;

	    hHe_2nd[ii]->Rebin2D(rbin,rbin);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hHe_2nd[ii]->GetBin(mm,nn);
                   Double_t content=hHe_2nd[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                     cout<<"There is underflow or overflow!"<<endl;
                     continue;
                   }
                   Double_t theta=-4.0+(mm-1)*dXBin+Theta_c1[ii];
                   Double_t nu=-0.2+(nn-1)*dYBin+Nu_c;
                   Double_t Ep=E0-nu;
                   if(theta<th_min)th_min=theta;
                   if(theta>th_max)th_max=theta;
                   if(Ep<Ep_min)Ep_min=Ep;
                   if(Ep>Ep_max)Ep_max=Ep;
               }
            }
            outfile<<"Marathon"<<endl;
            outfile<<Form("He3_kin%d_2nd",kin[ii+6])<<endl;
            outfile<<endl;
            outfile<<endl;
            outfile<<endl;
            outfile<<"E    Ep    theta"<<endl;

            th_max=th_max+dXBin;
            Ep_min=Ep_min-dYBin;

            Double_t XXbin=(th_max-th_min)/dXBin+1;
            Double_t YYbin=(Ep_max-Ep_min)/dYBin+1;

            for(int ii=0;ii<XXbin;ii++){
               Double_t theta = th_min+dXBin*ii;
               for(int jj=0;jj<YYbin;jj++){
                  Double_t Ep=Ep_max-dYBin*jj;
                  outfile<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep<<" "<<theta<<endl;
               }
            }
            outfile.close();
	}	

	for(int ii=0;ii<4;ii++){
	    ofstream outfile;
	    outfile.open(Form("Kin_table/H3_kin%d_2nd.inp",kin[ii+6]));

            Double_t th_min=100.0,th_max=0.0;
            Double_t Ep_min=100.0,Ep_max=0.0;

	    hH3_2nd[ii]->Rebin2D(rbin,rbin);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hH3_2nd[ii]->GetBin(mm,nn);
                   Double_t content=hH3_2nd[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                     cout<<"There is underflow or overflow!"<<endl;
                     continue;
                   }
                   Double_t theta=-4.0+(mm-1)*dXBin+Theta_c1[ii];
                   Double_t nu=-0.2+(nn-1)*dYBin+Nu_c;
                   Double_t Ep=E0-nu;
                   if(theta<th_min)th_min=theta;
                   if(theta>th_max)th_max=theta;
                   if(Ep<Ep_min)Ep_min=Ep;
                   if(Ep>Ep_max)Ep_max=Ep;
               }
            }
            outfile<<"Marathon"<<endl;
            outfile<<Form("H3_kin%d_2nd",kin[ii+6])<<endl;
            outfile<<endl;
            outfile<<endl;
            outfile<<endl;
            outfile<<"E    Ep    theta"<<endl;

            th_max=th_max+dXBin;
            Ep_min=Ep_min-dYBin;

            Double_t XXbin=(th_max-th_min)/dXBin+1;
            Double_t YYbin=(Ep_max-Ep_min)/dYBin+1;

            for(int ii=0;ii<XXbin;ii++){
               Double_t theta = th_min+dXBin*ii;
               for(int jj=0;jj<YYbin;jj++){
                  Double_t Ep=Ep_max-dYBin*jj;
                  outfile<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep<<" "<<theta<<endl;
               }
            }
            outfile.close();
	}	





}
