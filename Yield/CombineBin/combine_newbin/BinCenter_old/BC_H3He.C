#include "ReadFile.h"

void BC_H3He(){
     Double_t xbj[18]={0.19,0.22,0.25,0.29,0.33,0.36,0.385,0.43,0.48,0.51,0.55,0.59,0.63,0.67,0.7,0.74,0.78,0.82};
     int nBin[18]={2,3,4,4,4,2,2,3,2,2,2,2,2,2,2,3,2,2};

     Double_t Xbc_H3[MAXNUM]={0.0},Q2_bc_H3[MAXNUM]={0.0},Ybc_H3[MAXNUM]={0.0};
     Double_t Xbc_He3[MAXNUM]={0.0},Q2_bc_He3[MAXNUM]={0.0},Ybc_He3[MAXNUM]={0.0};
     Double_t Xi_H3[MAXNUM]={0.0},Q2_i_H3[MAXNUM]={0.0},Yi_H3[MAXNUM]={0.0};
     Double_t Xi_He3[MAXNUM]={0.0},Q2_i_He3[MAXNUM]={0.0},Yi_He3[MAXNUM]={0.0};
     Double_t BC_Corr[MAXNUM]={1.0};

     TString filename;
     filename="H3_Bincenter_xs.out";
     int nbc1=ReadFile(filename,Xbc_H3,Q2_bc_H3,Ybc_H3);
     filename="He3_Bincenter_xs.out";
     int nbc2=ReadFile(filename,Xbc_He3,Q2_bc_He3,Ybc_He3);
     if(nbc1!=nbc2){cout<<"Something wrong with BC file!"<<endl;exit(0);}

     filename="H3_H3He_xs.out";
     int n1=ReadFile(filename,Xi_H3,Q2_i_H3,Yi_H3);
     filename="He3_H3He_xs.out";
     int n2=ReadFile(filename,Xi_He3,Q2_i_He3,Yi_He3);
     if(n1!=n2){cout<<"Something wrong with all points file!"<<endl;exit(0);}

     nbc1=nbc1-1;
     n1=n1-1;

     int nn=1;
     for(int ii=0;ii<nbc1;ii++){
	 if(abs(Xbc_H3[ii]-Xbc_He3[ii])>0.001){cout<<"Xbj are different!!"<<endl;continue;}
	 Double_t H3He_BC=Ybc_H3[ii]/Ybc_He3[ii];

         int nbin=nn+nBin[ii];
	 for(int jj=nn;jj<nbin;jj++){
	    if(abs(Xi_H3[jj]-Xi_He3[jj])>0.001){cout<<"Xbj i are different!!"<<endl;continue;}
	    if(Yi_He3[jj]==0)continue;
	    Double_t H3He_i=Yi_H3[jj]/Yi_He3[jj];
	    BC_Corr[jj]=H3He_i/H3He_BC;
	    nn++;
	 }
     }

     ofstream outfile;
     outfile.open("BCfactor_H3He.dat");
     for(int ii=0;ii<n1;ii++){
	outfile<<Xi_H3[ii]<<"  "<<BC_Corr[ii]<<endl;
     }
     outfile.close();

}
