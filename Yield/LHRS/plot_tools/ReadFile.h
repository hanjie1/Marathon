#define MAXBIN 35
TString Yieldpath="/home/hanjie/work/MARATHON/RadCor/T2_externals/RUNPLAN/datafile/";
TString RCpath="/home/hanjie/work/MARATHON/RadCor/T2_externals/OUT/";
int ReadYield(TString filename,int kin,Double_t x[][MAXBIN],Double_t Q2[][MAXBIN],Double_t Yield[][MAXBIN],Double_t Y_err[][MAXBIN]){
    ifstream file;
    TString myfile=Yieldpath+filename;
    file.open(myfile);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    int KKin=0;
    if(kin<=5)KKin=kin-1;
    if(kin==7)KKin=kin-2;
    if(kin==9)KKin=kin-3;
    if(kin==11)KKin=kin-4;
    if(kin==13)KKin=kin-5;
    if(kin==15)KKin=kin-6;
    
    while(tmp.ReadLine(file)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from,", ");
          int bin=atoi(content.Data());
          tmp.Tokenize(content,from,", ");
          x[KKin][bin]=atof(content.Data());
          tmp.Tokenize(content,from,", ");
          Q2[KKin][bin]=atof(content.Data());
          tmp.Tokenize(content,from,", ");
          Yield[KKin][bin]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Y_err[KKin][bin]=atof(content.Data());
          //cout<<x[KKin][bin]<<"  "<<Q2[KKin][bin]<<"  "<<Yield[KKin][bin]<<"  "<<Y_err[KKin][bin]<<endl;
          from=0;
          nn++;
     }
    file.close();
    return 1;

}

int ReadRC(TString filename1,int kin,Double_t x[][MAXBIN],Double_t Q2[][MAXBIN],Double_t BornXS[][MAXBIN],Double_t RadXS[][MAXBIN]){
    ifstream file;
    TString myfile=RCpath+filename1;
    file.open(myfile);
    if(!file.is_open())return 0;
    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    int KKin=0;
    if(kin<=5)KKin=kin-1;
    if(kin==7)KKin=kin-2;
    if(kin==9)KKin=kin-3;
    if(kin==11)KKin=kin-4;
    if(kin==13)KKin=kin-5;
    if(kin==15)KKin=kin-6;

    while(tmp.ReadLine(file)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from," ");
          x[KKin][nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Q2[KKin][nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          BornXS[KKin][nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          RadXS[KKin][nn-1]=atof(content.Data());
          from=0;
          //if(kin==7)cout<<kin<<"  "<<x[KKin][nn-1]<<"  "<<Q2[KKin][nn-1]<<"  "<<BornXS[KKin][nn-1]<<"  "<<RadXS[KKin][nn-1]<<endl;
          nn++;
     }
    file.close();
    return 1;

}
