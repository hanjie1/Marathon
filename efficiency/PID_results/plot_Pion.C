void plot_Pion(){
    Double_t H1_pecErr[12]={0.0},D2_pecErr[12]={0.0},He_pecErr[12]={0.0},H3_pecErr[12]={0.0};
    Double_t H1_pec[12]={0.0},D2_pec[12]={0.0},He_pec[12]={0.0},H3_pec[12]={0.0};
    Double_t kin[12]={0,1,2,3,4,5,7,9,11,13,15,16};

    ifstream infile;
    infile.open("pitoe_Javier.txt");
    if(!infile.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;

    while(tmp.ReadLine(infile)){
          if(nn==0){nn++;continue;}
          if(tmp[0]=='#')continue;
          tmp.Tokenize(content,from,"  ");
          kin[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          H1_pec[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          H1_pecErr[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          D2_pec[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          D2_pecErr[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          He_pec[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          He_pecErr[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          H3_pec[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          H3_pecErr[nn-1]=atof(content.Data());
          from=0;
          nn++;
     }
    infile.close();

    TGraphErrors *gH1_pec=new TGraphErrors(5,kin,H1_pec,0,H1_pecErr);
    TGraphErrors *gD2_pec=new TGraphErrors(12,kin,D2_pec,0,D2_pecErr);
    TGraphErrors *gHe_pec=new TGraphErrors(12,kin,He_pec,0,He_pecErr);
    TGraphErrors *gH3_pec=new TGraphErrors(12,kin,H3_pec,0,H3_pecErr);

    TCanvas *c1=new TCanvas("c1","c1",1500,800);
    TMultiGraph *mg=new TMultiGraph();
    gH1_pec->SetMarkerStyle(8);
    gH1_pec->SetMarkerColor(1);
    gD2_pec->SetMarkerStyle(8);
    gD2_pec->SetMarkerColor(8);
    gHe_pec->SetMarkerStyle(8);
    gHe_pec->SetMarkerColor(4);
    gH3_pec->SetMarkerStyle(8);
    gH3_pec->SetMarkerColor(2);

    mg->Add(gH1_pec);
    mg->Add(gD2_pec);
    mg->Add(gHe_pec);
    mg->Add(gH3_pec);
    mg->Draw("AP");
    mg->SetTitle(";kinematic;#pi/e");

   auto leg1=new TLegend(0.65,0.7,0.75,0.85);
   leg1->AddEntry(gH1_pec,"{}^{1}H","P");
   leg1->AddEntry(gD2_pec,"{}^{2}H","P");
   leg1->AddEntry(gHe_pec,"{}^{3}He","P");
   leg1->AddEntry(gH3_pec,"{}^{3}H","P");
   leg1->Draw();

   mg->GetYaxis()->SetRangeUser(0,0.003);
   c1->Print("pitoe.pdf");
}
