TString rootpath="/lustre/expphy/cache/halla/triton/prod/marathon/pass2/";
//TString rootpath="/lustre/expphy/cache/halla/triton/prod/marathon/pass1_calibration/";
//TString rootpath="/lustre/expphy/volatile/halla/triton/hanjie/Marathon/";

TChain *GetChain(int run_number[],int size,int kin,TString TreeName){

    TChain *T=new TChain(TreeName.Data());

    for(int kk=0;kk<size;kk++){
        TString File[4];
        File[0]=rootpath+Form("kin%d/tritium_%d.root",kin,run_number[kk]);
        File[1]=rootpath+Form("kin%d_1st/tritium_%d.root",kin,run_number[kk]);
        File[2]=rootpath+Form("kin%d_2nd/tritium_%d.root",kin,run_number[kk]);
        File[3]=rootpath+Form("kin%d_3rd/tritium_%d.root",kin,run_number[kk]);

	TString kinfile;
        bool found=false;
        for(int ii=0;ii<4;ii++){
	    if(!gSystem->AccessPathName(File[ii])){
		found=true;
		if(ii==0)kinfile=Form("kin%d",kin);
		if(ii==1)kinfile=Form("kin%d_1st",kin);
		if(ii==2)kinfile=Form("kin%d_2nd",kin);
		if(ii==3)kinfile=Form("kin%d_3rd",kin);
		break;
	    }
	}
        TString lastfile;
	TString filename=rootpath+Form("%s/tritium_%d.root",kinfile.Data(),run_number[kk]);
        if(found){
           if((TreeName!="EndLeft")&&(TreeName!="EndRight"))T->Add(filename);
	   lastfile=filename;     

           int index=1;
	   filename=rootpath+Form("%s/tritium_%d_%d.root",kinfile.Data(),run_number[kk],index);
	   //filename=rootpath+Form("tritium_%d_%d.root",run_number,index);
           while(!gSystem->AccessPathName(filename)){
                  if((TreeName!="EndLeft")&&(TreeName!="EndRight"))T->Add(filename);
                  lastfile=filename;
                  index++;
		  filename=rootpath+Form("%s/tritium_%d_%d.root",kinfile.Data(),run_number[kk],index);
		  //filename=rootpath+Form("tritium_%d_%d.root",run_number,index);
           }
        }
        else {cout<<run_number[kk]<<" rootfile can't be found"<<endl; continue;}

        if(TreeName=="EndLeft"||TreeName=="EndRight"){
            T->Add(lastfile);
	    if(T->GetFile()==NULL)continue;
        }

     }
     return T;
}
