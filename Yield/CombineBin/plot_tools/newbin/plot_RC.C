void plot_RC()
{
    ifstream file1;
    TString myfile="Ratio_Dp.dat";
    file1.open(myfile);
    if(!file1.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;

    Double_t x_Dp[19]={0.0},RC_Dp[19]={0.0},RC_HeD[46]={0.0};
    Double_t x_HeD[46]={0.0},RC_H3D[46]={0.0},RC_H3He[46]={0.0};

    while(tmp.ReadLine(file1)){
          tmp.Tokenize(content,from,"  ");
          x_Dp[nn]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          RC_Dp[nn]=atof(content.Data());
          from=0;
          nn++;
     }
    file1.close();

    ifstream file2;
    file2.open("Ratio_HeD.dat");
    if(!file2.is_open())return 0;

    from=0;
    nn=0;
    while(tmp.ReadLine(file2)){
          tmp.Tokenize(content,from,"  ");
          x_HeD[nn]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          RC_HeD[nn]=atof(content.Data());
          from=0;
          nn++;
     }
    file2.close();

    ifstream file3;
    file3.open("Ratio_H3D.dat");
    if(!file3.is_open())return 0;

    from=0;
    nn=0;
    while(tmp.ReadLine(file3)){
          tmp.Tokenize(content,from,"  ");
          x_HeD[nn]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          RC_H3D[nn]=atof(content.Data());
          from=0;
          nn++;
     }
    file3.close();

    ifstream file4;
    file4.open("Ratio_H3He.dat");
    if(!file4.is_open())return 0;

    from=0;
    nn=0;
    while(tmp.ReadLine(file4)){
          tmp.Tokenize(content,from,"  ");
          x_HeD[nn]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          RC_H3He[nn]=atof(content.Data());
          from=0;
          nn++;
     }
    file4.close();

    TGraph *g1=new TGraph(19,x_Dp,RC_Dp);
    TGraph *g2=new TGraph(46,x_HeD,RC_HeD);
    TGraph *g3=new TGraph(46,x_HeD,RC_H3D);
    TGraph *g4=new TGraph(46,x_HeD,RC_H3He);
    
    TMultiGraph *mg=new TMultiGraph();
    g1->SetMarkerStyle(20);
    g2->SetMarkerStyle(21);
    g3->SetMarkerStyle(22);
    g4->SetMarkerStyle(23);
    g1->SetMarkerColor(1);
    g2->SetMarkerColor(2);
    g3->SetMarkerColor(4);
    g4->SetMarkerColor(8);
    g1->SetMarkerSize(1.3);
    g2->SetMarkerSize(1.3);
    g3->SetMarkerSize(1.3);
    g4->SetMarkerSize(1.3);

    mg->Add(g1); 
    mg->Add(g2); 
    mg->Add(g3); 
    mg->Add(g4); 
    mg->Draw("AP");
    mg->SetTitle(";x;RC factor;");

   auto leg1=new TLegend(0.7,0.6,0.9,0.9);
   leg1->AddEntry(g1,"D/p","P");
   leg1->AddEntry(g2,"^{3}He/D","P");
   leg1->AddEntry(g3,"^{3}H/D","P");
   leg1->AddEntry(g4,"^{3}H/^{3}He","P");
   leg1->Draw();

}
