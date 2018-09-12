#include <fstream>
static const char* ROOTFILE_FORMAT="%s/tritium_%d.root";
static const char* PATHS[] = {
  "/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/kin0",
  "/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/kin1",
  "/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/kin2",
  "/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/kin3",
  "/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/kin4",
  "/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/kin5",
  "/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/kin7",
  "/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/kin9",
  "/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/kin11",
  "/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/kin13",
  "/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/kin15",
  "/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/kin16",
  "/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/optics",
  "/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/positron",
  "/lustre/expphy/volatile/halla/triton/hanjie/Rootfiles",
   0
};

Bool_t IsFileExist(const Char_t * fname)
{
  fstream testfile;

  testfile.open(fname,ios_base::in);
  Bool_t isopen=testfile.is_open();
  testfile.close();

  return isopen;
}

TChain *GetTree(int run_number,TString TreeName){

        TChain *T=new TChain(TreeName.Data());
        TString File;
        const char** path=0;
        char filename[300];
        path = PATHS;
        Bool_t found=0;
        while ( path && *path ) {
           sprintf(filename,ROOTFILE_FORMAT,*path,run_number);
           if (IsFileExist(filename)) {
               found = 1;
               File=*path;
               cout <<filename<<endl;
               cout <<File<<endl;
               break;
           }
           path++;
        }
        if(found==0){
           cout<<"Root file "<<run_number<<" not found"<<endl;
           return 0;
        } 

        TString rootfile=File+Form("/tritium_%d.root",run_number);
        TString lastfile;
        if(!gSystem->AccessPathName(rootfile)){
           if(TreeName!="EndLeft")T->Add(rootfile);
	   lastfile=rootfile;     

           int index=1;
	   rootfile=File+Form("/tritium_%d_%d.root",run_number,index);
           while(!gSystem->AccessPathName(rootfile)){
		  if(TreeName!="EndLeft")T->Add(rootfile);
                  lastfile=rootfile;
                  index++;
	          rootfile=File+Form("/tritium_%d_%d.root",run_number,index);
           }
        }
        else {cout<<run_number<<" rootfile can't be found"<<endl; return 0;}

        if(TreeName=="EndLeft")
            T->Add(lastfile);

        return T;
}

