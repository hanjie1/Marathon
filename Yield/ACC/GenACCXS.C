#include "SearchXS.h"
int ReadXS(TString filename, Double_t Theta[nTh], Double_t Ep[nEp],Double_t XS_rad[nTh][nEp]){
	   for(int ii=0;ii<nTh;ii++){
		Theta[ii]=0.0;
		for(int jj=0;jj<nEp;jj++){
		    Ep[jj]=0.0;
		    XS_rad[ii][jj]=0.0;
		}
	    }

	   ifstream infile;
	   infile.open(Form("XS_table/%s",filename.Data()));
	   if(infile.is_open())cout<<filename.Data()<<" is open"<<endl;


           Ssiz_t from=0;
           TString content,tmp;
           Double_t Theta0=0.0,Ep0=0.0;
           int xx=0,yy=0;
           tmp.ReadLine(infile);
           while(tmp.ReadLine(infile)){
                tmp.Tokenize(content,from," ");
                tmp.Tokenize(content,from," ");
                tmp.Tokenize(content,from," ");
                Double_t aTheta=atof(content.Data());
                tmp.Tokenize(content,from," ");
                Double_t aEp=atof(content.Data());
                tmp.Tokenize(content,from," ");
                tmp.Tokenize(content,from," ");
                Double_t aXS=atof(content.Data());
                from=0;
		if(xx==0 && yy==0){
                   Theta0=aTheta; Ep0=aEp;
		   Theta[xx]=aTheta;  Ep[yy]=aEp;
		}
		if(aTheta!=Theta0){
		   xx++;
		   Theta0=aTheta;
		   Theta[xx]=aTheta;
		}

		if(aEp==Ep0)yy=0;
		else Ep[yy]=aEp;

		XS_rad[xx][yy]=aXS;
                yy++;
           }
           infile.close();

	   return 1;
}


void GenACCXS()
{
        TFile *f1=new TFile("Nu_Theta_all.root");
        TFile *f2=new TFile("ACC_rmXS.root","RECREATE");

	TH2F *hH1[5];
	TH2F *hD2[11];
	TH2F *hD2_2nd[4];
	TH2F *hHe[11];
	TH2F *hHe_2nd[4];
	TH2F *hH3[11];
	TH2F *hH3_2nd[4];

	TH2F *hH1new[5];
	TH2F *hD2new[11];
	TH2F *hD2new_2nd[4];
	TH2F *hHenew[11];
	TH2F *hHenew_2nd[4];
	TH2F *hH3new[11];
	TH2F *hH3new_2nd[4];

        Double_t Nu_c=7.49;
        Double_t Theta_c[11]={16.8075,17.5717,19.1125,20.575,21.9401,23.2065,25.5858,27.7642,29.8087,31.7274,33.5552};
        Double_t Theta_c1[4]={25.5909,27.7744,29.8159,31.7325};//2nd run of kin 7,9,11,13; 2nd and 3rd kin15 are almost the sam

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

	int nxbin=800+2;
	int nybin=400+2;
        Double_t E0=10.589;
        Double_t dXBin=0.01;
        Double_t dYBin=0.001;

        for(int ii=0;ii<5;ii++){
	    Double_t Ep[nEp],Theta[nTh],XS_rad[nTh][nEp];
	    TString filename=Form("H1_kin%d_xs.out",kin[ii]);	    
	    ReadXS(filename,Theta,Ep,XS_rad);

            hH1new[ii]=new TH2F(Form("H1_kin%d_new",kin[ii]),"Nu vs. Theta with XS removed",800,-4.0,4.0,400,-0.2,0.2);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hH1[ii]->GetBin(mm,nn);
                   Double_t content=hH1[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                      cout<<"There is underflow or overflow!"<<endl;
                      continue;
                   }
                   Double_t aTheta=-4.0+(mm-1)*dXBin+Theta_c[ii];
                   Double_t aNu=-0.2+(nn-1)*dYBin+Nu_c;
                   Double_t aEp=E0-aNu;

                   Double_t aXS=SearchXS(aEp,aTheta,Theta,Ep,XS_rad);

                   Double_t new_con=content/aXS;
		   Double_t new_err=sqrt(content)/aXS;
                   hH1new[ii]->SetBinContent(nbin,new_con);
                   hH1new[ii]->SetBinError(nbin,new_err);
               }
            }
	}


        for(int ii=0;ii<11;ii++){
	    Double_t Ep[nEp],Theta[nTh],XS_rad[nTh][nEp];
	    TString filename;
            if(kin[ii]<7||kin[ii]==15) filename=Form("D2_kin%d_xs.out",kin[ii]);	    
            else filename=Form("D2_kin%d_1st_xs.out",kin[ii]);	    
	    ReadXS(filename,Theta,Ep,XS_rad);

            if(kin[ii]<7||kin[ii]==15)
             hD2new[ii]=new TH2F(Form("D2_kin%d_new",kin[ii]),"Nu vs. Theta with XS removed",800,-4.0,4.0,400,-0.2,0.2);
	    else
	     hD2new[ii]=new TH2F(Form("D2_kin%d_1st_new",kin[ii]),"Nu vs. Theta with XS removed",800,-4.0,4.0,400,-0.2,0.2);

            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hD2[ii]->GetBin(mm,nn);
                   Double_t content=hD2[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                      cout<<"There is underflow or overflow!"<<endl;
                      continue;
                   }

                   Double_t aTheta=-4.0+(mm-1)*dXBin+Theta_c[ii];
                   Double_t aNu=-0.2+(nn-1)*dYBin+Nu_c;
                   Double_t aEp=E0-aNu;

                   Double_t aXS=SearchXS(aEp,aTheta,Theta,Ep,XS_rad);

                   Double_t new_con=content/aXS;
                   hD2new[ii]->SetBinContent(nbin,new_con);
		   Double_t new_err=sqrt(content)/aXS;
                   hD2new[ii]->SetBinError(nbin,new_err);
               }
            }
	}

        for(int ii=0;ii<4;ii++){
	    Double_t Ep[nEp],Theta[nTh],XS_rad[nTh][nEp];
	    TString filename;
            filename=Form("D2_kin%d_2nd_xs.out",kin[ii+6]);	    
	    ReadXS(filename,Theta,Ep,XS_rad);

	    hD2new_2nd[ii]=new TH2F(Form("D2_kin%d_2nd_new",kin[ii+6]),"Nu vs. Theta with XS removed",800,-4.0,4.0,400,-0.2,0.2);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hD2_2nd[ii]->GetBin(mm,nn);
                   Double_t content=hD2_2nd[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                      cout<<"There is underflow or overflow!"<<endl;
                      continue;
                   }

                   Double_t aTheta=-4.0+(mm-1)*dXBin+Theta_c1[ii];
                   Double_t aNu=-0.2+(nn-1)*dYBin+Nu_c;
                   Double_t aEp=E0-aNu;

                   Double_t aXS=SearchXS(aEp,aTheta,Theta,Ep,XS_rad);

                   Double_t new_con=content/aXS;
                   hD2new_2nd[ii]->SetBinContent(nbin,new_con);
		   Double_t new_err=sqrt(content)/aXS;
                   hD2new_2nd[ii]->SetBinError(nbin,new_err);
               }
            }
	}

        for(int ii=0;ii<11;ii++){
	    Double_t Ep[nEp],Theta[nTh],XS_rad[nTh][nEp];
	    TString filename;
            if(kin[ii]<7||kin[ii]==15) filename=Form("He3_kin%d_xs.out",kin[ii]);	    
            else filename=Form("He3_kin%d_1st_xs.out",kin[ii]);	    
	    ReadXS(filename,Theta,Ep,XS_rad);

            if(kin[ii]<7||kin[ii]==15)
             hHenew[ii]=new TH2F(Form("He3_kin%d_new",kin[ii]),"Nu vs. Theta with XS removed",800,-4.0,4.0,400,-0.2,0.2);
	    else
	     hHenew[ii]=new TH2F(Form("He3_kin%d_1st_new",kin[ii]),"Nu vs. Theta with XS removed",800,-4.0,4.0,400,-0.2,0.2);

            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hHe[ii]->GetBin(mm,nn);
                   Double_t content=hHe[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                      cout<<"There is underflow or overflow!"<<endl;
                      continue;
                   }

                   Double_t aTheta=-4.0+(mm-1)*dXBin+Theta_c[ii];
                   Double_t aNu=-0.2+(nn-1)*dYBin+Nu_c;
                   Double_t aEp=E0-aNu;

                   Double_t aXS=SearchXS(aEp,aTheta,Theta,Ep,XS_rad);

                   Double_t new_con=content/aXS;
                   hHenew[ii]->SetBinContent(nbin,new_con);
		   Double_t new_err=sqrt(content)/aXS;
                   hHenew[ii]->SetBinError(nbin,new_err);
               }
            }
	}

        for(int ii=0;ii<4;ii++){
	    Double_t Ep[nEp],Theta[nTh],XS_rad[nTh][nEp];
	    TString filename;
            filename=Form("He3_kin%d_2nd_xs.out",kin[ii+6]);	    
	    ReadXS(filename,Theta,Ep,XS_rad);

	    hHenew_2nd[ii]=new TH2F(Form("He3_kin%d_2nd_new",kin[ii+6]),"Nu vs. Theta with XS removed",800,-4.0,4.0,400,-0.2,0.2);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hHe_2nd[ii]->GetBin(mm,nn);
                   Double_t content=hHe_2nd[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                      cout<<"There is underflow or overflow!"<<endl;
                      continue;
                   }

                   Double_t aTheta=-4.0+(mm-1)*dXBin+Theta_c1[ii];
                   Double_t aNu=-0.2+(nn-1)*dYBin+Nu_c;
                   Double_t aEp=E0-aNu;

                   Double_t aXS=SearchXS(aEp,aTheta,Theta,Ep,XS_rad);

                   Double_t new_con=content/aXS;
                   hHenew_2nd[ii]->SetBinContent(nbin,new_con);
		   Double_t new_err=sqrt(content)/aXS;
                   hHenew_2nd[ii]->SetBinError(nbin,new_err);
               }
            }
	}

        for(int ii=0;ii<11;ii++){
	    Double_t Ep[nEp],Theta[nTh],XS_rad[nTh][nEp];
	    TString filename;
            if(kin[ii]<7||kin[ii]==15) filename=Form("H3_kin%d_xs.out",kin[ii]);	    
            else filename=Form("H3_kin%d_1st_xs.out",kin[ii]);	    
	    ReadXS(filename,Theta,Ep,XS_rad);

            if(kin[ii]<7||kin[ii]==15)
             hH3new[ii]=new TH2F(Form("H3_kin%d_new",kin[ii]),"Nu vs. Theta with XS removed",800,-4.0,4.0,400,-0.2,0.2);
	    else
	     hH3new[ii]=new TH2F(Form("H3_kin%d_1st_new",kin[ii]),"Nu vs. Theta with XS removed",800,-4.0,4.0,400,-0.2,0.2);

            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hH3[ii]->GetBin(mm,nn);
                   Double_t content=hH3[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                      cout<<"There is underflow or overflow!"<<endl;
                      continue;
                   }

                   Double_t aTheta=-4.0+(mm-1)*dXBin+Theta_c[ii];
                   Double_t aNu=-0.2+(nn-1)*dYBin+Nu_c;
                   Double_t aEp=E0-aNu;

                   Double_t aXS=SearchXS(aEp,aTheta,Theta,Ep,XS_rad);

                   Double_t new_con=content/aXS;
                   hH3new[ii]->SetBinContent(nbin,new_con);
		   Double_t new_err=sqrt(content)/aXS;
                   hH3new[ii]->SetBinError(nbin,new_err);
               }
            }
	}

        for(int ii=0;ii<4;ii++){
	    Double_t Ep[nEp],Theta[nTh],XS_rad[nTh][nEp];
	    TString filename;
            filename=Form("H3_kin%d_2nd_xs.out",kin[ii+6]);	    
	    ReadXS(filename,Theta,Ep,XS_rad);

	    hH3new_2nd[ii]=new TH2F(Form("H3_kin%d_2nd_new",kin[ii+6]),"Nu vs. Theta with XS removed",800,-4.0,4.0,400,-0.2,0.2);
            for(int mm=0;mm<nxbin;mm++){
               for(int nn=0;nn<nybin;nn++){
                   Double_t nbin=hH3_2nd[ii]->GetBin(mm,nn);
                   Double_t content=hH3_2nd[ii]->GetBinContent(nbin);
                   if(content==0)continue;

                   if(mm==0||nn==0||mm==nxbin-1||nn==nybin-1){
                      cout<<"There is underflow or overflow!"<<endl;
                      continue;
                   }

                   Double_t aTheta=-4.0+(mm-1)*dXBin+Theta_c1[ii];
                   Double_t aNu=-0.2+(nn-1)*dYBin+Nu_c;
                   Double_t aEp=E0-aNu;

                   Double_t aXS=SearchXS(aEp,aTheta,Theta,Ep,XS_rad);

                   Double_t new_con=content/aXS;
                   hH3new_2nd[ii]->SetBinContent(nbin,new_con);
		   Double_t new_err=sqrt(content)/aXS;
                   hH3new_2nd[ii]->SetBinError(nbin,new_err);
               }
            }
	}

        f2->Write();


        gStyle->SetOptStat(1111111);
        TCanvas *c1=new TCanvas("c1","c1",1500,1500);
	c1->Divide(3,2);
        for(int ii=0;ii<5;ii++){
	    c1->cd(ii+1);
	    hH1new[ii]->Draw("COLZ");
	}
 
        TCanvas *c2=new TCanvas("c2","c2",1500,1500);
	c2->Divide(3,2);
        for(int ii=0;ii<6;ii++){
	    c2->cd(ii+1);
	    hD2new[ii]->Draw("COLZ");
	}

        TCanvas *c3=new TCanvas("c3","c3",1500,1500);
	c3->Divide(3,2);
	int kk=0;
        for(int ii=6;ii<9;ii++){
	    c3->cd(kk+1);
	    hD2new[ii]->Draw("COLZ");
	    kk++;
	    c3->cd(kk+1);
	    hD2new_2nd[ii-6]->Draw("COLZ");
	    kk++;
	}

	TCanvas *c4=new TCanvas("c4","c4",1500,1500);	
	c4->Divide(2,2);
	c4->cd(1);
	hD2new[9]->Draw("COLZ");
	c4->cd(2);
	hD2new_2nd[3]->Draw("COLZ");
	c4->cd(3);
	hD2new[10]->Draw("COLZ");

        TCanvas *c5=new TCanvas("c5","c5",1500,1500);
	c5->Divide(3,2);
        for(int ii=0;ii<6;ii++){
	    c5->cd(ii+1);
	    hHenew[ii]->Draw("COLZ");
	}

        TCanvas *c6=new TCanvas("c6","c6",1500,1500);
	c6->Divide(3,2);
	kk=0;
        for(int ii=6;ii<9;ii++){
	    c6->cd(kk+1);
	    hHenew[ii]->Draw("COLZ");
	    kk++;
	    c6->cd(kk+1);
	    hHenew_2nd[ii-6]->Draw("COLZ");
	    kk++;
	}

	TCanvas *c7=new TCanvas("c7","c7",1500,1500);	
	c7->Divide(2,2);
	c7->cd(1);
	hHenew[9]->Draw("COLZ");
	c7->cd(2);
	hHenew_2nd[3]->Draw("COLZ");
	c7->cd(3);
	hHenew[10]->Draw("COLZ");

        TCanvas *c8=new TCanvas("c8","c8",1500,1500);
	c8->Divide(3,2);
        for(int ii=0;ii<6;ii++){
	    c8->cd(ii+1);
	    hH3new[ii]->Draw("COLZ");
	}

        TCanvas *c9=new TCanvas("c9","c9",1500,1500);
	c9->Divide(3,2);
	kk=0;
        for(int ii=6;ii<9;ii++){
	    c9->cd(kk+1);
	    hH3new[ii]->Draw("COLZ");
	    kk++;
	    c9->cd(kk+1);
	    hH3new_2nd[ii-6]->Draw("COLZ");
	    kk++;
	}

	TCanvas *c10=new TCanvas("c10","c10",1500,1500);	
	c10->Divide(2,2);
	c10->cd(1);
	hH3new[9]->Draw("COLZ");
	c10->cd(2);
	hH3new_2nd[3]->Draw("COLZ");
	c10->cd(3);
	hH3new[10]->Draw("COLZ");

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
