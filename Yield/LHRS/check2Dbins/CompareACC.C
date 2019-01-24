void CompareACC()
{
        TFile *f1=new TFile("nu_th_all.root");
        TH2F *hH1_kin0;
        TH2F *hH1_kin4;
        TH2F *hD2_kin0;
        TH2F *hD2_kin4;


        hH1_kin0 = (TH2F *)f1->Get("hH1_kin0_new");
        hH1_kin4 = (TH2F *)f1->Get("hH1_kin4_new");

	hH1_kin0->Rebin2D(20,20);
	hH1_kin4->Rebin2D(20,20);

        hD2_kin0 = (TH2F *)f1->Get("hD2_kin0_new");
        hD2_kin4 = (TH2F *)f1->Get("hD2_kin4_new");

	hD2_kin0->Rebin2D(20,20);
	hD2_kin4->Rebin2D(20,20);

        Int_t nxbin = 2200/20+2;
        Int_t nybin = 1000/20+2;

	TH2F *hratio_kin0=new TH2F("hratio_kin0","D2 ACC/H1 ACC",110,14,36,50,7,8);
	TH2F *hratio_kin4=new TH2F("hratio_kin4","D2 ACC/H1 ACC",110,14,36,50,7,8);

        Int_t maxbin1 = hH1_kin0->GetMaximumBin();
        Double_t maxcontent1 = hH1_kin0->GetBinContent(maxbin1);
        Int_t maxbin2 = hD2_kin0->GetMaximumBin();
        Double_t maxcontent2 = hD2_kin0->GetBinContent(maxbin2);

	for(int ii=0;ii<nxbin;ii++){
	    for(int jj=0;jj<nybin;jj++){
                Double_t nbin1=hH1_kin0->GetBin(ii,jj);
                Double_t content1=hH1_kin0->GetBinContent(nbin1);
                Double_t nbin2=hD2_kin0->GetBin(ii,jj);
                Double_t content2=hD2_kin0->GetBinContent(nbin2);

		if(nbin1 != nbin2){cout<<"!!! Wrong bin !!!"<<endl; continue;}
                if(content1==0||content2==0){
		   hratio_kin0->SetBinContent(nbin1,0);
		   continue;
		}

                Double_t new_con1=content1/maxcontent1;
                Double_t new_con2=content2/maxcontent2;

		Double_t ratio=new_con2/new_con1;

                ratio=(int)(ratio/0.01+0.5)*0.01;
		hratio_kin0->SetBinContent(nbin1,ratio);
	    }
	}
        
        maxbin1 = hH1_kin4->GetMaximumBin();
        maxcontent1 = hH1_kin4->GetBinContent(maxbin1);
        maxbin2 = hD2_kin4->GetMaximumBin();
        maxcontent2 = hD2_kin4->GetBinContent(maxbin2);

	for(int ii=0;ii<nxbin;ii++){
	    for(int jj=0;jj<nybin;jj++){
                Double_t nbin1=hH1_kin4->GetBin(ii,jj);
                Double_t content1=hH1_kin4->GetBinContent(nbin1);
                Double_t nbin2=hD2_kin4->GetBin(ii,jj);
                Double_t content2=hD2_kin4->GetBinContent(nbin2);

		if(nbin1 != nbin2){cout<<"!!! Wrong bin !!!"<<endl; continue;}
                if(content1==0||content2==0){
		   hratio_kin4->SetBinContent(nbin1,0);
		   continue;
		}

                Double_t new_con1=content1/maxcontent1;
                Double_t new_con2=content2/maxcontent2;

		Double_t ratio=new_con2/new_con1;
                ratio=(int)(ratio/0.01+0.5)*0.01;
		hratio_kin4->SetBinContent(nbin1,ratio);
	    }
	}

	TCanvas *c1=new TCanvas("c1");
        hratio_kin0->Draw("TEXT");

	TCanvas *c2=new TCanvas("c2");
        hratio_kin4->Draw("TEXT");

}
