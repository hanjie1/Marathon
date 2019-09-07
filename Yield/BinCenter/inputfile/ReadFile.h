#define MAXNUM 60
TString Yieldpath="/w/halla-scifs17exp/triton/hanjie/MARATHON/analysis/Yield/CombineBin/";
int ReadFile(TString filename,Double_t xavg[MAXNUM],Double_t Q2[MAXNUM],int kin[MAXNUM]){
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
          Q2[nn]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          kin[nn]=atoi(content.Data());

          from=0;
          nn++;
     }
    file.close();
    return nn;

}

