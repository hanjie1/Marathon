#define MAXBIN 26
TString Yieldpath="/w/halla-scifs17exp/triton/hanjie/MARATHON/analysis/Yield/CombineBin/Results/";
int ReadYield(TString filename,Double_t x[],Double_t Yield[],Double_t Y_err[]){
    ifstream file;
    TString myfile=Yieldpath+filename;
    file.open(myfile);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    
    while(tmp.ReadLine(file)){
          if(nn==0){nn++;continue;}
          if(tmp[0]=='#')continue; 
          tmp.Tokenize(content,from," ");
          x[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          Yield[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          Y_err[nn-1]=atof(content.Data());

          from=0;
          nn++;
     }
    file.close();
    return nn-1;

}

int ReadModel(TString filename,Double_t x[],Double_t Q2[],Double_t Ratio[],Double_t Rerr[]){
    ifstream file;
    TString myfile=Yieldpath+filename;
    file.open(myfile);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    
    while(tmp.ReadLine(file)){
          if(nn==0){nn++;continue;}
          if(tmp[0]=='#')continue; 
          tmp.Tokenize(content,from," ");
          x[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Q2[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Ratio[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Rerr[nn-1]=atof(content.Data());

          from=0;
          nn++;
     }
    file.close();
    return nn-1;

}

int ReadKP(TString filename,Double_t x[],Double_t F2P[],Double_t F2N[]){
    ifstream file;
    TString myfile=Yieldpath+filename;
    file.open(myfile);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    
    while(tmp.ReadLine(file)){
          if(nn==0){nn++;continue;}
          if(tmp[0]=='#')continue; 
          tmp.Tokenize(content,from," ");
          x[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          F2P[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          F2N[nn-1]=atof(content.Data());

          from=0;
          nn++;
     }
    file.close();
    return nn-1;

}

int ReadCJ(TString filename,Double_t x[],Double_t RD[]){
    ifstream file;
    TString myfile=Yieldpath+filename;
    file.open(myfile);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    
    while(tmp.ReadLine(file)){
          if(nn==0){nn++;continue;}
          if(tmp[0]=='#')continue; 
          tmp.Tokenize(content,from," ");
          x[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          RD[nn-1]=atof(content.Data());

          from=0;
          nn++;
     }
    file.close();
    return nn-1;

}
