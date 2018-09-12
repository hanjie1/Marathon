TString rootpath="/lustre/expphy/cache/halla/triton/prod/marathon/pass1/";
//TString rootpath="/lustre/expphy/volatile/halla/triton/Marathon_Rootfiles/pass1_calibration/";
TString listpath="/w/halla-scifs17exp/triton/Runlist/";

TChain *GetTree(int run_number,int kin,TString TreeName){

        TChain *T=new TChain(TreeName.Data());
        //TString File=rootpath+Form("kin%d/tritium_%d.root",kin,run_number);
        TString File=rootpath+Form("tritium_%d.root",run_number);
        TString lastfile;
        if(!gSystem->AccessPathName(File)){
           if(TreeName!="EndLeft")T->Add(File);
	   lastfile=File;     

           int index=1;
	   //File=rootpath+Form("kin%d/tritium_%d_%d.root",kin,run_number,index);
	   File=rootpath+Form("tritium_%d_%d.root",run_number,index);
           while(!gSystem->AccessPathName(File)){
		  if(TreeName!="EndLeft")T->Add(File);
                  lastfile=File;
                  index++;
		  //File=rootpath+Form("kin%d/tritium_%d_%d.root",kin,run_number,index);
		  File=rootpath+Form("tritium_%d_%d.root",run_number,index);
           }
        }
        else {cout<<run_number<<" rootfile can't be found"<<endl; return 0;}

        if(TreeName=="EndLeft")
            T->Add(lastfile);

        return T;
}

inline TChain *GetFiles(TString filename,TString TreeName){
       
       ifstream infile;
       TString File = listpath+filename+".dat";
       infile.open(File);
       if(!infile.is_open()){cout<<"!!! run list file not found "<<endl;return 0;} 

     TString tmp;
     TString target;
     int kin=0;
     if(tmp.ReadToken(infile))target=tmp;
     else{
          cout<<"No target type!!!"<<endl;
          exit(0);
     }

     if(tmp.ReadToken(infile))kin=atoi(tmp);
     else{
          cout<<"No kinematic!!!"<<endl;
          exit(0);
     }

     vector< Int_t> runList;
     int run_number,success=0;
     Ssiz_t from=0;
     TString content;
     if(tmp.ReadLine(infile)){
        while(tmp.Tokenize(content,from,","))
         {
              run_number = atoi(content);
              if(run_number>0)runList.push_back(run_number);
         }
    }
    infile.close();

    TChain *T=new TChain("T");
    for(Int_t ii=0;ii<runList.size();ii++)
      {
        TString File=rootpath+Form("kin%d/tritium_%d.root",kin,runList[ii]);
        //TString File=rootpath+Form("tritium_%d.root",runList[ii]);
        TString lastfile;
        if(!gSystem->AccessPathName(File)){
           if(TreeName!="EndLeft")T->Add(File);
           lastfile=File;

           int index=1;
           File=rootpath+Form("kin%d/tritium_%d_%d.root",kin,runList[ii],index);
           //File=rootpath+Form("tritium_%d_%d.root",runList[ii],index);
           while(!gSystem->AccessPathName(File)){
                  if(TreeName!="EndLeft")T->Add(File);
                  lastfile=File;
                  index++;
                  File=rootpath+Form("kin%d/tritium_%d_%d.root",kin,runList[ii],index);
                  //File=rootpath+Form("tritium_%d_%d.root",runList[ii],index);
           }
        }
        else {cout<<run_number<<" rootfile can't be found"<<endl; exit(0);}

        if(TreeName=="EndLeft")
            T->Add(lastfile);

      }

    return T;

}

