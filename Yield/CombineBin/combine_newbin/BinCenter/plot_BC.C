#define MAXNUM 60
void plot_BC()
{
    ifstream infile;
    infile.open("H3He_BCfac.dat");
    Ssiz_t from=0;
    TString content,tmp;
    int mm=0;
    Double_t Xbc[MAXNUM]={0.0},BCfactor[MAXNUM]={0.0};
    while(tmp.ReadLine(infile)){
          tmp.Tokenize(content,from," ");
          Xbc[mm]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          BCfactor[mm]=atof(content.Data());
          from=0;
          mm++;
     }
    infile.close();

    TGraph *gBC=new TGraph();
    for(int ii=0;ii<mm;ii++){
	gBC->SetPoint(ii,Xbc[ii],BCfactor[ii]);
    }
    gBC->SetMarkerStyle(8);
    gBC->Draw("AP");
    gBC->SetTitle("D/p BC factor;xbj");

}
