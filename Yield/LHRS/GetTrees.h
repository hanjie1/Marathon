TString rootpath="/lustre/expphy/cache/halla/triton/prod/marathon/pass1/";
//TString rootpath="/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/";
//TString rootpath="/lustre/expphy/volatile/halla/triton/hanjie/Marathon/";
TString listpath="/home/hanjie/work/MARATHON/analysis/Yield/Runlist/";

TChain *GetTree(int run_number,int kin,TString TreeName){

        TChain *T=new TChain(TreeName.Data());
        TString File=rootpath+Form("kin%d/tritium_%d.root",kin,run_number);
        //TString File=rootpath+Form("tritium_%d.root",run_number);
        TString lastfile;
        if(!gSystem->AccessPathName(File)){
           if(TreeName!="EndLeft")T->Add(File);
	   lastfile=File;     

           int index=1;
	   File=rootpath+Form("kin%d/tritium_%d_%d.root",kin,run_number,index);
	   //File=rootpath+Form("tritium_%d_%d.root",run_number,index);
           while(!gSystem->AccessPathName(File)){
		  if(TreeName!="EndLeft")T->Add(File);
                  lastfile=File;
                  index++;
		  File=rootpath+Form("kin%d/tritium_%d_%d.root",kin,run_number,index);
		  //File=rootpath+Form("tritium_%d_%d.root",run_number,index);
           }
        }
        else {cout<<run_number<<" rootfile can't be found"<<endl; return 0;}

        if(TreeName=="EndLeft"){
            T->Add(lastfile);
	    if(T->GetFile()==NULL)return 0;
        }

        return T;
}

inline TChain *GetFiles(TString filename,int kin,TString TreeName){
       
       ifstream infile;
       TString File = listpath+filename+".dat";
       infile.open(File);
       if(!infile.is_open()){cout<<"!!! run list file not found "<<endl;return 0;} 

       TString tmp;
       int run_number;
       TChain *T=new TChain(TreeName.Data());
       while(tmp.ReadToken(infile)){
	     run_number=tmp.Atoi();
	     T=GetTree(run_number,kin,TreeName);             
       }
       
       infile.close();
       return T;

}

