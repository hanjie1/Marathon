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

int ReadModel(TString filename,Double_t x[],Double_t Yield[],Double_t Y_err[]){
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
          Y_err[nn-1]=atof(content.Data());

          from=0;
          nn++;
     }
    file.close();
    return nn-1;

}

int ReadNMC(TString filename,Double_t x[],Double_t Yield[],Double_t Y_loerr[],Double_t Y_hierr[]){
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
          Y_loerr[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Y_hierr[nn-1]=atof(content.Data());

          from=0;
          nn++;
     }
    file.close();
    return nn-1;

}

int ReadGlobal(TString filename,Double_t nexp[],Double_t x[],Double_t Q2[],Double_t Yield[],Double_t Yerr[],Double_t Ynorm[]){
    ifstream file;
    TString myfile=Yieldpath+filename;
    file.open(myfile);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    
    while(tmp.ReadLine(file)){
          if(tmp[0]=='#')continue; 
          tmp.Tokenize(content,from," ");
          nexp[nn]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          x[nn]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Q2[nn]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Yield[nn]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Yerr[nn]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Ynorm[nn]=atof(content.Data());

          from=0;
          nn++;
     }
    file.close();
    return nn;

}

int ReadKP(TString filename,Double_t x[],Double_t F2P[],Double_t F2D[]){
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
          tmp.Tokenize(content,from," ");
          F2D[nn-1]=atof(content.Data());

          from=0;
          nn++;
     }
    file.close();
    return nn-1;

}
                    
