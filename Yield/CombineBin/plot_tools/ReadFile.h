#define MAXBIN 26
TString Yieldpath="/w/halla-scifs17exp/triton/hanjie/MARATHON/analysis/Yield/CombineBin/Kin_Yield/";
int ReadYield(TString filename,int kin,Double_t x[][MAXBIN],Double_t xavg[][MAXBIN],Double_t Q2[][MAXBIN],Double_t Yield[][MAXBIN],Double_t Y_err[][MAXBIN]){
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
          tmp.Tokenize(content,from,", ");
          x[kin][nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,", ");
          xavg[kin][nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,", ");
          Q2[kin][nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,", ");
          Yield[kin][nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Y_err[kin][nn-1]=atof(content.Data());
          //cout<<x[kin][nn-1]<<"  "<<xavg[kin][nn-1]<<"  "<<Q2[kin][nn-1]<<"  "<<Yield[kin][nn-1]<<"  "<<Y_err[kin][nn-1]<<endl;
          from=0;
          nn++;
     }
    file.close();
    return 1;

}

int ReadRadCor(TString filename,int kin,Double_t x[][MAXBIN],Double_t Q2[][MAXBIN],Double_t RadCor[][MAXBIN]){
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
          x[kin][nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Q2[kin][nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          RadCor[kin][nn-1]=atof(content.Data());
          cout<<x[kin][nn-1]<<"  "<<Q2[kin][nn-1]<<"  "<<RadCor[kin][nn-1]<<endl;
          from=0;
          nn++;
     }
    file.close();
    return 1;
}
