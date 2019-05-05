#include "ReadFile.h"
void SortX(){
    Double_t x[MAXNUM]={0.0},Ratio[MAXNUM]={0.0},Rerr[MAXNUM]={0.0},RadCor[MAXNUM]={0.0};
    Double_t Q2[MAXNUM]={0.0};
    int kin[MAXNUM]={0};

    TString filename;
    filename="newbin/Ratio_Dp.dat";
    int totalN = ReadFile(filename,x,Q2,Ratio,Rerr,RadCor,kin);
    cout<<totalN<<endl;
    if(totalN==0){cout<<"No ratio !!"<<endl;exit(0);}

    int nmove=0;
    
    for(int ii=0;ii<totalN;ii++){
     int NN=totalN-ii;
     for(int jj=0;jj<NN-1;jj++){
	if(x[jj]>x[jj+1]){
	   Double_t tmp_x,tmp_Ratio,tmp_Rerr,tmp_RC,tmp_kin,tmp_Q2;
	   tmp_x=x[jj];
	   x[jj]=x[jj+1];
	   x[jj+1]=tmp_x;

	   tmp_Q2=Q2[jj];
	   Q2[jj]=Q2[jj+1];
	   Q2[jj+1]=tmp_Q2;

	   tmp_Ratio=Ratio[jj];
	   Ratio[jj]=Ratio[jj+1];
	   Ratio[jj+1]=tmp_Ratio;

	   tmp_Rerr=Rerr[jj];
	   Rerr[jj]=Rerr[jj+1];
	   Rerr[jj+1]=tmp_Rerr;

	   tmp_RC=RadCor[jj];
	   RadCor[jj]=RadCor[jj+1];
	   RadCor[jj+1]=tmp_RC;

	   tmp_kin=kin[jj];
	   kin[jj]=kin[jj+1];
	   kin[jj+1]=tmp_kin;

	   nmove++;
	}
     }
     if(nmove==0)break;
    }

    ofstream outfile;
    outfile.open("Xbj_sort_Dp.dat");
    for(int ii=0;ii<totalN;ii++){
	outfile<<x[ii]<<"  "<<Q2[ii]<<"  "<<Ratio[ii]<<"  "<<Rerr[ii]<<"  "<<RadCor[ii]<<"  "<<kin[ii]<<endl;
    }
    outfile.close();

}
