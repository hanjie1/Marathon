#define MAXNUM 60
TString Yieldpath="/w/halla-scifs17exp/triton/hanjie/MARATHON/analysis/Yield/CombineBin/plot_tools/";
int ReadFile(TString filename,Double_t xavg[MAXNUM],Double_t Yield[MAXNUM],Double_t Y_err[MAXNUM],Double_t RadCor[MAXNUM]){
    ifstream file;
    TString myfile=Yieldpath+filename;
    file.open(myfile);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    
    while(tmp.ReadLine(file)){
          tmp.Tokenize(content,from,"  ");
          xavg[nn]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          Yield[nn]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          Y_err[nn]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          RadCor[nn]=atof(content.Data());
          from=0;
          nn++;
     }
    file.close();
    return nn;

}

