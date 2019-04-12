void DefineCut()
{
     TFile *f1=new TFile("ACC_ratio.root");
     
     TH2F *hH3He;
     hH3He=(TH2F *)f1->Get("H3He_kin1");
     
     int rbin=8;
     int nxbin=800/rbin+2;
     int nybin=400/rbin+2;

     ofstream outfile;
     outfile.open("ACC_matrix4.dat");
     TH2F *hACC=new TH2F("hACC","ACC matrix",800/rbin,-4.0,4.0,400/rbin,-0.2,0.2); 
     for(int ii=0;ii<nxbin;ii++){
	for(int jj=0;jj<nybin;jj++){
           Double_t nbin=hH3He->GetBin(ii,jj);
           Double_t content=hH3He->GetBinContent(nbin);
	   Double_t err=hH3He->GetBinError(nbin);

	   if(content==0) continue;

	   Double_t xbin=hH3He->GetXaxis()->GetBinLowEdge(ii);
	   Double_t ybin=hH3He->GetYaxis()->GetBinLowEdge(jj);
           if(err/content<0.20 && content>=0.75 && content<=1.25){
              hACC->SetBinContent(nbin,1);
              outfile<<xbin<<"  "<<ybin<<"  "<<1<<endl;
	   }
           else{
	      hACC->SetBinContent(nbin,0);
              outfile<<xbin<<"  "<<ybin<<"  "<<0<<endl;
	   }
	}
     }

     outfile.close();
     hACC->Draw("TEXT");

}
