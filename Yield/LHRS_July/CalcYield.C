#include "GetTrees.h"
#include "SetCut.h"
#include "SetCons.h"
#include "CalcLum.h"
#include "CalcLT.h"
#include "SearchXS.h"
#include <TMath.h>

TRI_VAR GetLT(int run_number)
{
     TString table="LHRStest";
     TSQLServer* Server = TSQLServer::Connect("mysql://halladb/triton-work","triton-user","3He3Hdata");

     TString  query=Form("select * from %s where run_number=%d;",table.Data(),run_number);
     TSQLResult* result=Server->Query(query.Data());
     TSQLRow *row;
     int nrows = result->GetRowCount();
     TRI_VAR LT;
     if (nrows==0) {
        cout<<"run "<<run_number<<" is NOT in "<<table<<endl; 
        Server->Close();
        LT.value=0;
        LT.err=0;   
        return LT;
     }

     query=Form("select livetime,LT_err from %s where run_number=%d;",table.Data(),run_number);
     result=Server->Query(query.Data());
     row=result->Next();
     Double_t livetime=atof(row->GetField(0));
     Double_t livetime_err=atof(row->GetField(1));

     LT.value=livetime;
     LT.err=livetime_err;

     Server->Close(); 
     return LT;      
}

void CalcYield(){
     TString filename;
     cout<<"Input file name: "<<endl;
     cin>>filename;
     filename="/w/halla-scifs17exp/triton/Runlist/"+filename;
//     filename=""+filename;

     Double_t LUM=CalcLum(filename); //total luminosity get for this kinematics;
     cout<<"Get total Luminosity:  "<<LUM<<endl;

     ifstream infile;
     infile.open(filename);
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

     ofstream myfile;
     myfile.open(Form("RawYield/%s_kin%d.txt",target.Data(),kin));
     myfile<<"xbj  Yield"<<endl;

     vector<Int_t> runList;
     int run_number,nrun=0;
     Ssiz_t from=0;
     TString content;
     if(tmp.ReadLine(infile)){
        while(tmp.Tokenize(content,from,","))
         {
              run_number = atoi(content);
              if(run_number>0){runList.push_back(run_number);nrun++;}
         }
     }
 
     infile.close();

     Double_t xbj[17];
     for(int ii=0;ii<17;ii++){
         xbj[ii]=0.16+ii*0.02;
     }

     // read XS table
     Double_t xs_rad[nTh][nEp];
     Double_t xs_born[nTh][nEp];
     Double_t theta[nTh],Eprime[nEp];
     ifstream infile1;
     infile1.open(Form("RadCor/event/%s_kin%d_xs.out",target.Data(),kin));
     if(!infile1.is_open()){cout<<"!!! radiative correction file not found "<<endl;return 0;}
     tmp.ReadLine(infile1);
     from=0;
     int xx=0,yy=0,nline=0; 
     while(tmp.ReadLine(infile1)){
        for(int ii=0;ii<4;ii++){
           tmp.Tokenize(content,from,"  ");
           if(ii==0)theta[xx]=atof(content); 
           if(ii==1)Eprime[yy]=atof(content); 
           if(ii==2)xs_born[xx][yy]=atof(content); 
           if(ii==3)xs_rad[xx][yy]=atof(content); 
	}
        from=0;
       // cout<<theta[xx]<<"  "<<Eprime[yy]<<"  "<<xs_born[xx][yy]/xs_rad[xx][yy]<<endl;
        yy++;
        nline++;
        if(nline%nEp==0){xx++;yy=0;}
     }
     infile1.close();


     TString TreeName="T";
     TChain* T;
     Double_t totalNe[17]={0.0};
     Double_t RawNe[17]={0.0};
     Double_t totalQ2[17]={0.0};
     Double_t totalXbj[17]={0.0};
     Double_t totalNe_err[17]={0.0};
     for(int ii=0;ii<nrun;ii++){
         run_number=runList[ii];
         TRI_VAR LT=CalcLT(run_number,kin,1);
         Double_t livetime=LT.value; 
         Double_t livetime_err=LT.err; 
         cout<<"Get LT:  "<<livetime<<"  "<<livetime_err<<endl;
         T=GetTree(run_number,kin,TreeName);
         T->Draw(">>electron",trigger2+CK+Ep+beta+ACC+VZ+TRK);
         TEventList *electron;
         gDirectory->GetObject("electron",electron);
         T->SetEventList(electron);

         T->SetBranchStatus("*",0);
         T->SetBranchStatus("L.gold.p",1);
         T->SetBranchStatus("EKLx.angle",1);
         T->SetBranchStatus("EKLx.x_bj",1);
         T->SetBranchStatus("EKLx.Q2",1);
         T->SetBranchStatus("LeftBCMev.isrenewed",1);

         Double_t aEprime=0.0,aTheta=0.0,axbj=0.0,aQ2=0.0;
         Double_t isrenewed=0;
         T->SetBranchAddress("L.gold.p",&aEprime);
         T->SetBranchAddress("EKLx.angle",&aTheta);
         T->SetBranchAddress("EKLx.x_bj",&axbj);
         T->SetBranchAddress("EKLx.Q2",&aQ2);
         T->SetBranchAddress("LeftBCMev.isrenewed",&isrenewed);

         Double_t Radcor=1.0;
	 Int_t NNe[17]={0};
         Int_t nentries=electron->GetN();
         for(int jj=0;jj<nentries;jj++){
	     T->GetEntry(electron->GetEntry(jj));
         //    cout<<aTheta<<endl;
//	     Radcor=SearchXS(kin,aEprime,aTheta,theta,Eprime,xs_born,xs_rad);
             if(Radcor==0.0){
                cout<<"This event can't find Radcor"<<endl;
                Radcor=1.0;
             }
             for(int kk=0;kk<17;kk++){
		 Double_t dxbj=axbj-xbj[kk];
                 if(dxbj<0.02 && dxbj>=0){
		    totalNe[kk]+=1.0/livetime*Radcor;
		    NNe[kk]++;
                    RawNe[kk]++;
		    totalQ2[kk]+=aQ2;
		    totalXbj[kk]+=axbj;
		    break;
                 }
	     }
         } 

         for(int jj=0;jj<17;jj++){
           Double_t tmp_err=0.0;
           if(NNe[jj]!=0){
              tmp_err=sqrt(1.0/NNe[jj]+(livetime_err/livetime)*(livetime_err/livetime))*NNe[jj]/livetime;
           }
	   totalNe_err[jj]+=tmp_err*tmp_err;
         }
         delete T;
         cout<<"Get electron coutns"<<endl;
     }


    Double_t rawYield[17]={0.0};
    Double_t rawYield_err[17]={0.0};
    Double_t avgQ2[17]={0.0};
    Double_t avgXbj[17]={0.0};
    for(int ii=0;ii<17;ii++){
        totalNe_err[ii]=sqrt(totalNe_err[ii]);
	rawYield[ii]=totalNe[ii]/LUM;
	rawYield_err[ii]=totalNe_err[ii]/LUM;
        double tmp_x=0.17+ii*0.02;
        if(RawNe[ii]!=0){
           avgQ2[ii]=totalQ2[ii]/RawNe[ii];
	   avgXbj[ii]=totalXbj[ii]/RawNe[ii];
        }
        myfile<<avgXbj[ii]<<"   "<<rawYield[ii]<<"   "<<rawYield_err[ii]<<"  "<<avgQ2[ii]<<endl;
    }

   myfile.close();

}
