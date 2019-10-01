#define MAXBIN 70
int ReadYield(TString filename,Double_t x[],Double_t Corr[]){
    ifstream file;
    TString myfile=filename;
    file.open(myfile);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    
    while(tmp.ReadLine(file)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from,",");
          x[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,",");
          Corr[nn-1]=atof(content.Data());

          from=0;
          nn++;
     }
    file.close();
    return nn-1;

}
