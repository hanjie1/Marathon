void GenACC()
{
        TFile *f1=new TFile("nu_th_all.root");
        TH2F *hD2_kin0;
        TH2F *hD2_kin4;
        TH2F *hD2_kin15;

	Double_t xcen[3]={16.8075,21.9401,33.555};
	Double_t ycen=7.49;

        hD2_kin0 = (TH2F *)f1->Get("hD2_kin0_new");
        hD2_kin4 = (TH2F *)f1->Get("hD2_kin4_new");
        hD2_kin15 = (TH2F *)f1->Get("hD2_kin15_new");

	hD2_kin0->Rebin2D(20,20);
	hD2_kin4->Rebin2D(20,20);
	hD2_kin15->Rebin2D(20,20);

        Int_t nxbin = 2200/20+2;
        Int_t nybin = 1000/20+2;

	Int_t maxbin = hD2_kin0->GetMaximumBin();
	Double_t maxcontent = hD2_kin0->GetBinContent(maxbin);
	for(int ii=0;ii<nxbin;ii++){
	    for(int jj=0;jj<nybin;jj++){
                Double_t nbin=hD2_kin0->GetBin(ii,jj);
                Double_t content=hD2_kin0->GetBinContent(nbin);
                if(content==0)continue;

                if(ii==0||jj==0||ii==nxbin-1||jj==nybin-1){
                   cout<<"There is underflow or overflow!"<<endl;
                   continue;
                }
		
		//Double_t bin_x=hD2_kin0->GetXaxis()->GetBinLowEdge(ii);
		//Double_t bin_y=hD2_kin0->GetXaxis()->GetBinLowEdge(jj);

		//Double_t newx=bin_x-xcen[0];
		//Double_t newy=bin_y-ycen;

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
