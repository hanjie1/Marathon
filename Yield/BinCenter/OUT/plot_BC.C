void plot_BC(){
    Double_t x[46]={0.0},BC[46]={0.0};
    TString target;

    cout<<"Which ratio?"<<endl;
    cin>>target;

    ifstream file;
    TString myfile=Form("%s_BCfac.dat",target.Data());
    file.open(myfile);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;

    while(tmp.ReadLine(file)){
          tmp.Tokenize(content,from," ");
          x[nn]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          BC[nn]=atof(content.Data());
          from=0;
cout<<x[nn]<<"  "<<BC[nn]<<endl;
          nn++;
     }
    file.close();

    TGraph *gBC=new TGraph(46,x,BC);
    gBC->SetMarkerStyle(8);
    gBC->SetMarkerColor(4);
    gBC->Draw("AP");
}
