TString rootpath="/lustre/expphy/cache/halla/triton/prod/marathon/pass2/";
TString rootpath1="/lustre/expphy/volatile/halla/triton/hanjie/Rootfiles";
//TString rootpath="/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/";
//TString rootpath="/lustre/expphy/volatile/halla/triton/hanjie/Marathon/";

TChain *GetTree(int run_number,int kin,TString TreeName){

        TChain *T=new TChain(TreeName.Data());
        TString File[5];
        File[0]=rootpath+Form("kin%d/tritium_%d.root",kin,run_number);
        File[1]=rootpath+Form("kin%d_1st/tritium_%d.root",kin,run_number);
        File[2]=rootpath+Form("kin%d_2nd/tritium_%d.root",kin,run_number);
        File[3]=rootpath+Form("kin%d_3rd/tritium_%d.root",kin,run_number);
        File[4]=rootpath1+Form("/tritium_%d.root",run_number);

	TString kinfile;
        bool found=false;
        for(int ii=0;ii<5;ii++){
	    if(!gSystem->AccessPathName(File[ii])){
		found=true;
		if(ii==0)kinfile=Form("kin%d",kin);
		if(ii==1)kinfile=Form("kin%d_1st",kin);
		if(ii==2)kinfile=Form("kin%d_2nd",kin);
		if(ii==3)kinfile=Form("kin%d_3rd",kin);
		if(ii==4){kinfile="";rootpath=rootpath1;};
		break;
	    }
	}
        TString lastfile;
	TString filename=rootpath+Form("%s/tritium_%d.root",kinfile.Data(),run_number);
        if(found){
           if((TreeName!="EndLeft")&&(TreeName!="EndRight"))T->Add(filename);
	   lastfile=filename;     

           int index=1;
	   filename=rootpath+Form("%s/tritium_%d_%d.root",kinfile.Data(),run_number,index);
	   //filename=rootpath+Form("tritium_%d_%d.root",run_number,index);
           while(!gSystem->AccessPathName(filename)){
		  if((TreeName!="EndLeft")&&(TreeName!="EndRight"))T->Add(filename);
                  lastfile=filename;
                  index++;
		  filename=rootpath+Form("%s/tritium_%d_%d.root",kinfile.Data(),run_number,index);
		  //filename=rootpath+Form("tritium_%d_%d.root",run_number,index);
           }
        }
        else {cout<<run_number<<" rootfile can't be found"<<endl; return 0;}

        if((TreeName=="EndLeft")||(TreeName=="EndRight")){
            T->Add(lastfile);
	    if(T->GetFile()==NULL)return 0;
        }

        return T;
}
