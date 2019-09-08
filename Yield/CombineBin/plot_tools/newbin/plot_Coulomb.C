#define MAXNUM 60
void plot_Coulomb()
{
    ifstream infile;
    infile.open("Ratio_H3He.dat");
    Ssiz_t from=0;
    TString content,tmp;
    int mm=0;
    Double_t Xbc[MAXNUM]={0.0},CouFac[MAXNUM]={0.0};
    while(tmp.ReadLine(infile)){
          tmp.Tokenize(content,from," ");
          Xbc[mm]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          CouFac[mm]=atof(content.Data());
          from=0;
          mm++;
     }
    infile.close();

    TGraph *gCoul=new TGraph();
    for(int ii=0;ii<mm;ii++){
	gCoul->SetPoint(ii,Xbc[ii],CouFac[ii]);
    }
    gCoul->SetMarkerStyle(8);
    gCoul->Draw("AP");
    gCoul->SetTitle("D/p Coulomb correction factor;xbj");

}
