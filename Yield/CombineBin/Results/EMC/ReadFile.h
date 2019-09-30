#define MAXBIN 26
TString Yieldpath="/w/halla-scifs17exp/triton/hanjie/MARATHON/analysis/Yield/CombineBin/Results/";
int ReadYield(TString filename,Double_t x[], Double_t Q2[],Double_t Yield[],Double_t Yerr[],Double_t ST[],Double_t SY[]){
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
          Yield[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          ST[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          SY[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
	  Yerr[nn-1]=atof(content.Data());

          from=0;
          nn++;
     }
    file.close();
    return nn-1;

}

int ReadNP(TString filename,Double_t x[],Double_t Yield[],Double_t Yerr[],Double_t ST[],Double_t SY[]){
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
          Yerr[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          ST[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          SY[nn-1]=atof(content.Data());

          from=0;
          nn++;
     }
    file.close();
    return nn;

}

int ReadHallC(TString filename,Double_t x[],Double_t R[],Double_t ST[],Double_t SY[]){
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
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          R[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          ST[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          SY[nn-1]=atof(content.Data());

          from=0;
          nn++;
     }
    file.close();
    return nn-1;

}

int ReadHERMES(TString filename,Double_t x[],Double_t R[],Double_t ST[],Double_t SY[]){
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
          R[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          ST[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          SY[nn-1]=atof(content.Data());
          from=0;
          nn++;
     }
    file.close();
    return nn-1;

}

int ReadISO(TString filename,Double_t x[],Double_t Yield[],Double_t ST[],Double_t SY[]){
    ifstream file;
    TString myfile=filename;
    file.open(myfile);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    
    while(tmp.ReadLine(file)){
          if(tmp[0]=='\\')continue; 
          tmp.Tokenize(content,from," & ");
          x[nn]=atof(content.Data());
          tmp.Tokenize(content,from," & ");
          tmp.Tokenize(content,from," & ");
          Yield[nn]=atof(content.Data());
          tmp.Tokenize(content,from," & ");
          ST[nn]=atof(content.Data());
          tmp.Tokenize(content,from," & ");
          SY[nn]=atof(content.Data());
          from=0;
          nn++;
     }
    file.close();
    return nn;

}
