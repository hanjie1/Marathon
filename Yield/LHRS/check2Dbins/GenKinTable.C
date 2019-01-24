void GenKinTable(){
	TFile *f1=new TFile("nu_th_all.root");
	TH2F *hH1_kin0;
	TH2F *hH1_kin4;
//	TH2F *hH1_kin15;

	hH1_kin0 = (TH2F *)f1->Get("H1_kin0");
	hH1_kin4 = (TH2F *)f1->Get("H1_kin4");
//	hH1_kin15 = (TH2F *)f1->Get("H1_kin15");

	int rbin=10;
        hH1_kin0->Rebin2D(rbin,rbin);
        hH1_kin4->Rebin2D(rbin,rbin);
//        hH1_kin15->Rebin2D(rbin,rbin);

        Int_t nxbin = 2200/rbin+2;
        Int_t nybin = 1000/rbin+2;

	Double_t E0=10.589;
	Double_t dXBin=0.01*rbin;
	Double_t dYBin=0.001*rbin;

	Double_t th_min=100.0,th_max=0.0;
	Double_t Ep_min=100.0,Ep_max=0.0;
	for(int ii=0;ii<nxbin;ii++){
	    for(int jj=0;jj<nybin;jj++){
	        Double_t nbin=hH1_kin0->GetBin(ii,jj);
		Double_t content=hH1_kin0->GetBinContent(nbin);
		if(content==0)continue;
	   
		if(ii==0||jj==0||ii==nxbin-1||jj==nybin-1){
		   cout<<"There is underflow or overflow!"<<endl;
		   continue;
		}
		Double_t theta=14+(ii-1)*dXBin;
		Double_t nu=7+(jj-1)*dYBin;
		Double_t Ep=E0-nu;
	        if(theta<th_min)th_min=theta;
	        if(theta>th_max)th_max=theta;
	        if(Ep<Ep_min)Ep_min=Ep;
	        if(Ep>Ep_max)Ep_max=Ep;
	    }
	}
	ofstream outfile1;
	outfile1.open("Table/H1_kin0.inp");
        outfile1<<"Marathon"<<endl;
        outfile1<<"H1_kin0"<<endl;
        outfile1<<endl;
        outfile1<<endl;
        outfile1<<endl;
        outfile1<<"E    Ep    theta"<<endl;
        
	th_max=th_max+dXBin;
	Ep_min=Ep_min-dYBin;

	Double_t XXbin=(th_max-th_min)/dXBin+1;
	Double_t YYbin=(Ep_max-Ep_min)/dYBin+1;
    
        for(int ii=0;ii<XXbin;ii++){
	    Double_t theta = th_min+dXBin*ii;
	    for(int jj=0;jj<YYbin;jj++){
		Double_t Ep=Ep_max-dYBin*jj;
		outfile1<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep<<" "<<theta<<endl;
	    }
	}
	outfile1.close();

	th_min=100.0,th_max=0.0;
	Ep_min=100.0,Ep_max=0.0;
	for(int ii=0;ii<nxbin;ii++){
	    for(int jj=0;jj<nybin;jj++){
	        Double_t nbin=hH1_kin4->GetBin(ii,jj);
		Double_t content=hH1_kin4->GetBinContent(nbin);
		if(content==0)continue;
	   
		if(ii==0||jj==0||ii==nxbin-1||jj==nybin-1){
		   cout<<"There is underflow or overflow!"<<endl;
		   continue;
		}
		Double_t theta=14+(ii-1)*dXBin;
		Double_t nu=7+(jj-1)*dYBin;
		Double_t Ep=E0-nu;
	        if(theta<th_min)th_min=theta;
	        if(theta>th_max)th_max=theta;
	        if(Ep<Ep_min)Ep_min=Ep;
	        if(Ep>Ep_max)Ep_max=Ep;
	    }
	}

	ofstream outfile2;
	outfile2.open("Table/H1_kin4.inp");
        outfile2<<"Marathon"<<endl;
        outfile2<<"H1_kin4"<<endl;
        outfile2<<endl;
        outfile2<<endl;
        outfile2<<endl;
        outfile2<<"E    Ep    theta"<<endl;
        
	th_max=th_max+dXBin;
	Ep_min=Ep_min-dYBin;

	XXbin=(th_max-th_min)/dXBin+1;
	YYbin=(Ep_max-Ep_min)/dYBin+1;

        for(int ii=0;ii<XXbin;ii++){
	    Double_t theta = th_min+dXBin*ii;
	    for(int jj=0;jj<YYbin;jj++){
		Double_t Ep=Ep_max-dYBin*jj;
		outfile2<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep<<" "<<theta<<endl;
	    }
	}
	outfile2.close();
/*
	th_min=100.0,th_max=0.0;
	Ep_min=100.0,Ep_max=0.0;
	for(int ii=0;ii<nxbin;ii++){
	    for(int jj=0;jj<nybin;jj++){
	        Double_t nbin=hH1_kin15->GetBin(ii,jj);
		Double_t content=hH1_kin15->GetBinContent(nbin);
		if(content==0)continue;
	   
		if(ii==0||jj==0||ii==nxbin-1||jj==nybin-1){
		   cout<<"There is underflow or overflow!"<<endl;
		   continue;
		}
		Double_t theta=14+(ii-1)*dXBin;
		Double_t nu=7+(jj-1)*dYBin;
		Double_t Ep=E0-nu;
	        if(theta<th_min)th_min=theta;
	        if(theta>th_max)th_max=theta;
	        if(Ep<Ep_min)Ep_min=Ep;
	        if(Ep>Ep_max)Ep_max=Ep;
	    }
	}

	ofstream outfile3;
	outfile3.open("Table/H1_kin15.inp");
        outfile3<<"Marathon"<<endl;
        outfile3<<"H1_kin15"<<endl;
        outfile3<<endl;
        outfile3<<endl;
        outfile3<<endl;
        outfile3<<"E    Ep    theta"<<endl;
        
	th_max=th_max+dXBin;
	Ep_min=Ep_min-dYBin;

	XXbin=(th_max-th_min)/dXBin+1;
	YYbin=(Ep_max-Ep_min)/dYBin+1;

        for(int ii=0;ii<XXbin;ii++){
	    Double_t theta = th_min+dXBin*ii;
	    for(int jj=0;jj<YYbin;jj++){
		Double_t Ep=Ep_max-dYBin*jj;
		outfile3<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep<<" "<<theta<<endl;
	    }
	}
	outfile3.close();
*/
}
