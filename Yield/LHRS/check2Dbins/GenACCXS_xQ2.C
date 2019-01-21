#define MAXBIN 100
#include "SearchXS.h"

using namespace TMath;
void GenACCXS_xQ2()
{
        TFile *f1=new TFile("X_Q2_all.root","UPDATE");
//        TFile *f2=new TFile("nu_th_removeXS.root","RECREATE");
        TH2F *hD2_kin0;
        TH2F *hD2_kin4;
        TH2F *hD2_kin15;

        TH2F *hD2_kin0_new=new TH2F("hD2_kin0_new","update histogram with XS removed",830,0.12,0.95,2400,1.5,13.5);
        TH2F *hD2_kin4_new=new TH2F("hD2_kin4_new","update histogram with XS removed",830,0.12,0.95,2400,1.5,13.5);
        TH2F *hD2_kin15_new=new TH2F("hD2_kin15_new","update histogram with XS removed",830,0.12,0.95,2400,1.5,13.5);

        hD2_kin0 = (TH2F *)f1->Get("D2_kin0");
        hD2_kin4 = (TH2F *)f1->Get("D2_kin4");
        hD2_kin15 = (TH2F *)f1->Get("D2_kin15");

        Int_t nxbin = 830+2;
        Int_t nybin = 2400+2;
        Double_t E0=10.589;
        Double_t Mp=0.93827;
        Double_t dXBin=0.001;
        Double_t dYBin=0.005;


	ifstream infile1;
	infile1.open("Table/D2_kin0_xs.out");

        Ssiz_t from=0;
        TString content,tmp;
	Double_t Theta[nTh]={0.0},Ep[nEp]={0.0};
	Double_t XS_rad[nTh][nEp];
        int xx=0,yy=0;
	tmp.ReadLine(infile1);
        while(tmp.ReadLine(infile1)){
              tmp.Tokenize(content,from," ");
              tmp.Tokenize(content,from," ");
              tmp.Tokenize(content,from," ");
              Theta[xx]=atof(content.Data());
              tmp.Tokenize(content,from," ");
              Ep[yy]=atof(content.Data());
              tmp.Tokenize(content,from," ");
              tmp.Tokenize(content,from," ");
              XS_rad[xx][yy]=atof(content.Data());
              from=0;
	      yy++;
	      if(yy%nEp==0){xx++;yy=0;}
        }
        infile1.close();



//	Int_t maxbin = hD2_kin0->GetMaximumBin();
//	Double_t maxcontent = hD2_kin0->GetBinContent(maxbin);
	for(int ii=0;ii<nxbin;ii++){
	    for(int jj=0;jj<nybin;jj++){
                Double_t nbin=hD2_kin0->GetBin(ii,jj);
                Double_t content=hD2_kin0->GetBinContent(nbin);
                if(content==0){
		   hD2_kin0_new->SetBinContent(nbin,0);
                   continue;
		}

                if(ii==0||jj==0||ii==nxbin-1||jj==nybin-1){
                   cout<<"There is underflow or overflow!"<<endl;
                   continue;
                }
	
                Double_t aXbj=0.12+(ii-1)*dXBin+dXBin/2.0;
                Double_t Q2=1.5+(jj-1)*dYBin+dYBin/2.0;
		Double_t nu=Q2/(2*Mp*aXbj);
                Double_t aEp=E0-nu;
		Double_t aTheta=2*ASin(sqrt(Q2/(4*E0*aEp)));
		aTheta=aTheta*180/pi;
		//cout<<aTheta<<"  "<<aEp<<endl;
		Double_t aXS=SearchXS(aEp,aTheta,Theta,Ep,XS_rad);	

		if(aXS==0)aXS=1;	
		Double_t new_con=content/aXS;
		//new_con=(int)(new_con/0.001+0.5)*0.001;
		//new_con=new_con*100;
		hD2_kin0_new->SetBinContent(nbin,new_con);
	    }
	}


	ifstream infile2;
	infile2.open("Table/D2_kin4_xs.out");

	for(int ii=0;ii<nTh;ii++){
	    Theta[ii]=0.0;
	    for(int jj=0;jj<nEp;jj++){
		Ep[jj]=0.0;
		XS_rad[ii][jj]=0.0;	
	    }
	}

        from=0;
        xx=0,yy=0;
	tmp.ReadLine(infile2);
        while(tmp.ReadLine(infile2)){
              tmp.Tokenize(content,from," ");
              tmp.Tokenize(content,from," ");
              tmp.Tokenize(content,from," ");
              Theta[xx]=atof(content.Data());
              tmp.Tokenize(content,from," ");
              Ep[yy]=atof(content.Data());
              tmp.Tokenize(content,from," ");
              tmp.Tokenize(content,from," ");
              XS_rad[xx][yy]=atof(content.Data());
              from=0;
	      yy++;
	      if(yy%nEp==0){xx++;yy=0;}
        }
        infile2.close();

	for(int ii=0;ii<nxbin;ii++){
	    for(int jj=0;jj<nybin;jj++){
                Double_t nbin=hD2_kin4->GetBin(ii,jj);
                Double_t content=hD2_kin4->GetBinContent(nbin);
                if(content==0){
		   hD2_kin4_new->SetBinContent(nbin,0);
                   continue;
		}

                if(ii==0||jj==0||ii==nxbin-1||jj==nybin-1){
                   cout<<"There is underflow or overflow!"<<endl;
                   continue;
                }
	
                Double_t aXbj=0.12+(ii-1)*dXBin+dXBin/2.0;
                Double_t Q2=1.5+(jj-1)*dYBin+dYBin/2.0;
		Double_t nu=Q2/(2*Mp*aXbj);
                Double_t aEp=E0-nu;
		Double_t aTheta=2*ASin(sqrt(Q2/(4*E0*aEp)));
		aTheta=aTheta*180/pi;
		//cout<<aTheta<<"  "<<aEp<<endl;
		Double_t aXS=SearchXS(aEp,aTheta,Theta,Ep,XS_rad);	
	
		if(aXS==0)aXS=1;	
		Double_t new_con=content/aXS;
		//new_con=(int)(new_con/0.001+0.5)*0.001;
		//new_con=new_con*100;
		hD2_kin4_new->SetBinContent(nbin,new_con);
	    }
	}

	for(int ii=0;ii<nTh;ii++){
	    Theta[ii]=0.0;
	    for(int jj=0;jj<nEp;jj++){
		Ep[jj]=0.0;
		XS_rad[ii][jj]=0.0;	
	    }
	}

	ifstream infile3;
	infile3.open("Table/D2_kin15_xs.out");

        from=0;
        xx=0,yy=0;
	tmp.ReadLine(infile3);
        while(tmp.ReadLine(infile3)){
              tmp.Tokenize(content,from," ");
              tmp.Tokenize(content,from," ");
              tmp.Tokenize(content,from," ");
              Theta[xx]=atof(content.Data());
              tmp.Tokenize(content,from," ");
              Ep[yy]=atof(content.Data());
              tmp.Tokenize(content,from," ");
              tmp.Tokenize(content,from," ");
              XS_rad[xx][yy]=atof(content.Data());
              from=0;
	      yy++;
	      if(yy%nEp==0){xx++;yy=0;}
        }
        infile3.close();

	for(int ii=0;ii<nxbin;ii++){
	    for(int jj=0;jj<nybin;jj++){
                Double_t nbin=hD2_kin15->GetBin(ii,jj);
                Double_t content=hD2_kin15->GetBinContent(nbin);
                if(content==0){
		   hD2_kin15_new->SetBinContent(nbin,0);
                   continue;
		}

                if(ii==0||jj==0||ii==nxbin-1||jj==nybin-1){
                   cout<<"There is underflow or overflow!"<<endl;
                   continue;
                }
	
                Double_t aXbj=0.12+(ii-1)*dXBin+dXBin/2.0;
                Double_t Q2=1.5+(jj-1)*dYBin+dYBin/2.0;
		Double_t nu=Q2/(2*Mp*aXbj);
                Double_t aEp=E0-nu;
		Double_t aTheta=2*ASin(sqrt(Q2/(4*E0*aEp)));
		aTheta=aTheta*180/pi;
		//cout<<aTheta<<"  "<<aEp<<endl;
		Double_t aXS=SearchXS(aEp,aTheta,Theta,Ep,XS_rad);	
	
		if(aXS==0)aXS=1;	
		Double_t new_con=content/aXS;
		//new_con=(int)(new_con/0.001+0.5)*0.001;
		//new_con=new_con*100;
		hD2_kin15_new->SetBinContent(nbin,new_con);
	    }
	}

	f1->Write();

}
