#include "SetCons.h"
#include "ReadFile.h"
using namespace std;

void genBCfile(){
    Double_t xavg[MAXNUM]={0.0},Q2[MAXNUM]={0.0};
    int kin[MAXNUM]={0};

    TString ratio="Dp";

    TString filename;
    filename=Form("combine_newbin/Xbj_sort_%s.dat",ratio.Data());
    int totalN = ReadFile(filename,xavg,Q2,kin);
    cout<<totalN<<endl;
    if(totalN==0){cout<<"No ratio !!"<<endl;exit(0);}

    ofstream outfile;
    outfile.open(Form("%s_x.dat",ratio.Data()));

    int totalBin=0;   
    if(ratio=="Dp")totalBin=8;
    else totalBin=19;
 
    int nn=0;
    int nkin[12]={0};

    for(int ii=0;ii<totalBin;ii++){
	int tmpN=0;
	Double_t xxBC=0.0;
        if(ratio=="Dp"){
	   xxBC=X_center_Dp[ii];
           tmpN=nn+BCnBin_Dp[ii];
	}
	else{
	   xxBC=X_center[ii];
           tmpN=nn+BCnBin[ii];
	}

	for(int jj=nn;jj<tmpN;jj++){
            int tmpKin=-1;
            if(kin[jj]<6)tmpKin=kin[jj];
            if(kin[jj]==7)tmpKin=6;
            if(kin[jj]==9)tmpKin=7;
            if(kin[jj]==11)tmpKin=8;
            if(kin[jj]==13)tmpKin=9;
            if(kin[jj]==15)tmpKin=10;
            if(kin[jj]==16)tmpKin=11;

	    Double_t tmpXmin=xmin[tmpKin];
	    Double_t tmpXmax=tmpXmin+dBin[tmpKin][0];
	    for(int kk=0;kk<nkin[tmpKin];kk++){
	        tmpXmin+=dBin[tmpKin][kk];
	        tmpXmax+=dBin[tmpKin][kk+1];
	    }
	   
	    outfile<<right<<setw(7)<<tmpXmin<<"  "<<right<<setw(7)<<tmpXmax<<"  "<<right<<setw(10)<<Q2[jj]<<"  ";
	    outfile<<right<<setw(7)<<xxBC<<"  "<<right<<setw(10)<<xavg[jj]<<"  "<<right<<setw(3)<<kin[jj]<<endl;
	    nkin[tmpKin]++;
	    nn++;

	}        
    }

    outfile.close();


}
