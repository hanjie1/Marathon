#define MAXNUM 60
int ReadFile(TString filename,Double_t xavg[MAXNUM],Double_t BCfac[MAXNUM]){
    ifstream file;
    file.open(filename);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    
    while(tmp.ReadLine(file)){
          tmp.Tokenize(content,from,"  ");
          xavg[nn]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          BCfac[nn]=atof(content.Data());

          from=0;
          nn++;
     }
    file.close();
    return nn;

}

