#include "ReadFile.h"

void BC_Dp(){
     Double_t xbj[6]={0.19,0.22,0.25,0.29,0.32,0.34};
     int nBin[6]={2,3,4,4,2,2};

     Double_t Xbc_D2[MAXNUM]={0.0},Q2_bc_D2[MAXNUM]={0.0},Ybc_D2[MAXNUM]={0.0};
     Double_t Xbc_H1[MAXNUM]={0.0},Q2_bc_H1[MAXNUM]={0.0},Ybc_H1[MAXNUM]={0.0};
     Double_t Xi_D2[MAXNUM]={0.0},Q2_i_D2[MAXNUM]={0.0},Yi_D2[MAXNUM]={0.0};
     Double_t Xi_H1[MAXNUM]={0.0},Q2_i_H1[MAXNUM]={0.0},Yi_H1[MAXNUM]={0.0};
     Double_t BC_Corr[MAXNUM]={1.0};

     TString filename;
     filename="D2_BC_Dp_xs.out";
     int nbc1=ReadFile(filename,Xbc_D2,Q2_bc_D2,Ybc_D2);
     filename="H1_BC_Dp_xs.out";
     int nbc2=ReadFile(filename,Xbc_H1,Q2_bc_H1,Ybc_H1);
     if(nbc1!=nbc2){cout<<"Something wrong with BC file!"<<endl;exit(0);}

     filename="D2_Dp_xs.out";
     int n1=ReadFile(filename,Xi_D2,Q2_i_D2,Yi_D2);
     filename="H1_Dp_xs.out";
     int n2=ReadFile(filename,Xi_H1,Q2_i_H1,Yi_H1);
     if(n1!=n2){cout<<"Something wrong with all points file!"<<endl;exit(0);}

     nbc1=nbc1-1;
     n1=n1-1;

     int nn=1;
     for(int ii=0;ii<nbc1;ii++){
	 if(abs(Xbc_D2[ii]-Xbc_H1[ii])>0.001){cout<<"Xbj are different!!"<<endl;continue;}
	 Double_t Dp_BC=Ybc_D2[ii]/Ybc_H1[ii];

         int nbin=nn+nBin[ii];
	 for(int jj=nn;jj<nbin;jj++){
	    if(abs(Xi_D2[jj]-Xi_H1[jj])>0.001){cout<<"Xbj i are different!!"<<endl;continue;}
	    if(Yi_H1[jj]==0)continue;
	    Double_t Dp_i=Yi_D2[jj]/Yi_H1[jj];
	    BC_Corr[jj]=Dp_i/Dp_BC;
cout<<Yi_D2[jj]/Ybc_D2[ii]<<"  "<<Yi_H1[jj]/Ybc_H1[ii]<<endl;
	    nn++;
	 }
     }

     ofstream outfile;
     outfile.open("BCfactor_Dp.dat");
     for(int ii=0;ii<n1-1;ii++){
	outfile<<Xi_D2[ii]<<"  "<<BC_Corr[ii]<<endl;
     }
     outfile<<Xi_D2[n1-1]<<"  "<<1.0<<endl;
     outfile.close();

}
