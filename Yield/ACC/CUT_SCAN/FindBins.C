// This program is to find the mini/max xbj value for each kinematics 
// which will be used to decide the bin x for Yield ratio

#include "SetCut.h"
#include "GetRunList.h"
#include "GetTrees.h"

void FindBins()
{
   TString target[4]={"H1","D2","He3","H3"};
   int kin[4]={1,4,9,15};
   Double_t x_min[4]={0.0},x_max[4]={0.0};


//first loop should be for kinematics and second layer should be for target;
//for each kinematicas, find the min/max xbj among all targets; 
   for(int ii=0;ii<4;ii++){
      Double_t xxmin=1000.0,xxmax=0.0;
      for(int jj=0;jj<4;jj++){
        if(jj==0 && kin[ii]>4)continue;
	int KKin=0;
        if(kin[ii]<6)KKin=kin[ii];
        if(kin[ii]==7)KKin=6;
        if(kin[ii]==9)KKin=7;
        if(kin[ii]==11)KKin=8;
        if(kin[ii]==13)KKin=9;
        if(kin[ii]==15)KKin=10;
        TCut VZ=Form("L.tr.vz<%f && L.tr.vz>%f",vz_max[KKin],vz_min[KKin]);

        vector<vector<Int_t> > runList;
        int run_number=0,nrun=0;
        nrun=GetRunList(runList,kin[ii],target[jj]);
        cout<<nrun<<" runs are added "<<endl;
        if(nrun==0)exit(0);

        Double_t run_xmin=1000,run_xmax=0.0; 
	TChain *T=new TChain("T");
	for(int kk=0;kk<nrun;kk++){
	    run_number=runList[kk][0];
            T=GetTree(run_number,kin[ii],"T");
            T->Draw(">>electron",trigger2+CK+Ep+beta+ACC+VZ+TRK+W2);
            TEventList *electron;
            gDirectory->GetObject("electron",electron);
            T->SetEventList(electron);

            Double_t axbj=0.0;
            T->SetBranchStatus("*",0);
            T->SetBranchStatus("EKLxe.x_bj",1);
            T->SetBranchAddress("EKLxe.x_bj",&axbj);

            Int_t nentries=electron->GetN();
	    Double_t tmpx_min=1000.0,tmpx_max=0.0;
            for(int mm=0;mm<nentries;mm++){
               T->GetEntry(electron->GetEntry(mm));
	       if(axbj<tmpx_min)tmpx_min=axbj;
	       if(axbj>tmpx_max)tmpx_max=axbj;
	    }

	    if(tmpx_min<run_xmin) run_xmin=tmpx_min;
	    if(tmpx_max>run_xmax) run_xmax=tmpx_max;
	}
            if(xxmin>run_xmin) xxmin=run_xmin;
            if(xxmax<run_xmax) xxmax=run_xmax;
      }
      x_min[ii]=(int)(xxmin/0.001+0.5)*0.001;
      x_max[ii]=(int)(xxmax/0.001+0.5)*0.001;
   }

   int nBin[4]={0};
   for(int ii=0;ii<4;ii++){
	nBin[ii]=(int)((x_max[ii]-x_min[ii])/0.03+0.5);
   }

   ofstream outfile;
   outfile.open("xbin.dat");
   outfile<<"xmin[4]={";
   outfile<<fixed<<setprecision(3);
   for(int ii=0;ii<4;ii++){
       if(ii<3)outfile<<x_min[ii]<<",";
       if(ii==3)outfile<<x_min[ii]<<"};"<<endl;
   }

   outfile<<"xmax[4]={";
   outfile<<fixed<<setprecision(3);
   for(int ii=0;ii<4;ii++){
       if(ii<3)outfile<<x_max[ii]<<",";
       if(ii==3)outfile<<x_max[ii]<<"};"<<endl;
   }

   outfile<<"nBin[4]={";
   for(int ii=0;ii<4;ii++){
       if(ii<3)outfile<<nBin[ii]<<",";
       if(ii==3)outfile<<nBin[ii]<<"};"<<endl;
   }
   outfile.close();

}
