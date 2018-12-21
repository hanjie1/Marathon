#define MAXBIN 100

void GenACCXS()
{
        TFile *f1=new TFile("nu_th_all.root");
        TH2F *hD2_kin0;
        TH2F *hD2_kin4;
        TH2F *hD2_kin15;

        hD2_kin0 = (TH2F *)f1->Get("D2_kin0");
        hD2_kin4 = (TH2F *)f1->Get("D2_kin4");
        hD2_kin15 = (TH2F *)f1->Get("D2_kin15");

        Int_t nxbin = 2200+2;
        Int_t nybin = 1000+2;
        Double_t E0=10.589;

	ifstream infile1;
	infile1.open("Table/D2_kin0_xs.out");

	int nEp=28,nTh=40;

        Ssiz_t from=0;
        TString content,tmp;
	Double_t Theta[MAXBIN]={0.0},Ep[MAXBIN]={0.0};
	Double_t XS_rad[MAXBIN][MAXBIN]=0.0;
        int nn=0,xx=0,yy=0;
        while(tmp.ReadLine(infile1)){
              if(nn==0){nn++;continue;}
              tmp.Tokenize(content,from,", ");
              tmp.Tokenize(content,from,", ");
              tmp.Tokenize(content,from,", ");
              Theta[xx]=atof(content.Data());
              tmp.Tokenize(content,from,", ");
              Ep[yy]=atof(content.Data());
              tmp.Tokenize(content,from," ");
              tmp.Tokenize(content,from," ");
              XS_rad[xx][yy]=atof(content.Data());
              from=0;
              nn++;
	      yy++;
	      if(yy%nEp==0){xx++;yy=0;}
        }
        infile1.close();


//	Int_t maxbin = hD2_kin0->GetMaximumBin();
//	Double_t maxcontent = hD2_kin0->GetBinContent(maxbin);
	nn=0;
	for(int ii=0;ii<nxbin;ii++){
	    for(int jj=0;jj<nybin;jj++){
                Double_t nbin=hD2_kin0->GetBin(ii,jj);
                Double_t content=hD2_kin0->GetBinContent(nbin);
                if(content==0)continue;

                if(ii==0||jj==0||ii==nxbin-1||jj==nybin-1){
                   cout<<"There is underflow or overflow!"<<endl;
                   continue;
                }
	
                Double_t aTheta=14+(ii-1)*dXBin+dXBin/2.0;
                Double_t nu=7+(jj-1)*dYBin+dYBin/2.0;
                Double_t aEp=E0-nu;

		Double_t new_con=content/maxcontent;
		new_con=(int)(new_con/0.001+0.5)*0.001;
		new_con=new_con*100;
		hD2_kin0->SetBinContent(nbin,new_con);
	    }
	}
        
	maxbin = hD2_kin4->GetMaximumBin();
	maxcontent = hD2_kin4->GetBinContent(maxbin);
	for(int ii=0;ii<nxbin;ii++){
	    for(int jj=0;jj<nybin;jj++){
                Double_t nbin=hD2_kin4->GetBin(ii,jj);
                Double_t content=hD2_kin4->GetBinContent(nbin);
                if(content==0)continue;

                if(ii==0||jj==0||ii==nxbin-1||jj==nybin-1){
                   cout<<"There is underflow or overflow!"<<endl;
                   continue;
                }
		
		Double_t new_con=content/maxcontent;
		new_con=(int)(new_con/0.001+0.5)*0.001;
		new_con=new_con*100;
		hD2_kin4->SetBinContent(nbin,new_con);
	    }
	}
  
	maxbin = hD2_kin15->GetMaximumBin();
	maxcontent = hD2_kin15->GetBinContent(maxbin);
	for(int ii=0;ii<nxbin;ii++){
	    for(int jj=0;jj<nybin;jj++){
                Double_t nbin=hD2_kin15->GetBin(ii,jj);
                Double_t content=hD2_kin15->GetBinContent(nbin);
                if(content==0)continue;

                if(ii==0||jj==0||ii==nxbin-1||jj==nybin-1){
                   cout<<"There is underflow or overflow!"<<endl;
                   continue;
                }
		
		Double_t new_con=content/maxcontent;
		new_con=(int)(new_con/0.001+0.5)*0.001;
		new_con=new_con*100;
		hD2_kin15->SetBinContent(nbin,new_con);
	    }
	}

	TCanvas *c1=new TCanvas("c1");
        hD2_kin0->Draw("TEXT");

	TCanvas *c2=new TCanvas("c2");
        hD2_kin4->Draw("TEXT");

	TCanvas *c3=new TCanvas("c3");
        hD2_kin15->Draw("TEXT");
}
