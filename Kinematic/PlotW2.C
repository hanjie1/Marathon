#include "GetTrees.h"
#include "GetRunList.h"
#include "SetCut.h"

void PlotW2(){
     TString target = "H3";
     int kin[12]={0,1,2,3,4,5,7,9,11,13,15,16};

     TGraph *hW2[12];
     TGraph *hQ2[12];

     for(int ii=0;ii<12;ii++){
         vector<Int_t> runList;
         int run_number=0,nrun=0;
         nrun=GetRunList(runList,kin[ii],target);
         cout<<nrun<<" runs are added "<<endl;
         if(nrun==0)exit(0);
	 TCut VZ;
	 if(ii<11)VZ=Form("L.tr.vz<%f && L.tr.vz>%f",vz_max[ii],vz_min[ii]);
	
	 hW2[ii]=new TGraph();
	 hQ2[ii]=new TGraph();
	 int npoints=0;
	 for(int jj=0;jj<nrun;jj++){
             run_number=runList[jj];
	     TChain* T;
             T=GetTree(run_number,kin[ii],"T");

	     if(ii<11){
                T->Draw(">>electron",trigger2+CK+Ep+beta+ACC+VZ+TRK);
                TEventList *electron;
                gDirectory->GetObject("electron",electron);
                T->SetEventList(electron);
   
                T->SetBranchStatus("*",0);
                T->SetBranchStatus("EKLxe.x_bj",1);
                T->SetBranchStatus("EKLxe.Q2",1);
                T->SetBranchStatus("EKLxe.W2",1);
   
                Double_t axbj=0.0,aQ2=0.0,aW2=0.0;
                T->SetBranchAddress("EKLxe.x_bj",&axbj);
                T->SetBranchAddress("EKLxe.Q2",&aQ2);
                T->SetBranchAddress("EKLxe.W2",&aW2);
   	    
                Int_t nentries=electron->GetN();
                for(int nn=0;nn<nentries;nn++){
                   T->GetEntry(electron->GetEntry(nn));
     	     	   hW2[ii]->SetPoint(npoints,axbj,aW2);
     	     	   hQ2[ii]->SetPoint(npoints,axbj,aQ2);
   	 	   npoints++;
       	        }
	     }

	     if(ii==11){
                T->Draw(">>electron",trigger2_R+CK_R+Ep_R+beta_R+ACC_R+VZ_R+TRK_R);
                TEventList *electron;
                gDirectory->GetObject("electron",electron);
                T->SetEventList(electron);
   
                T->SetBranchStatus("*",0);
                T->SetBranchStatus("EKRxe.x_bj",1);
                T->SetBranchStatus("EKRxe.Q2",1);
                T->SetBranchStatus("EKRxe.W2",1);
   
                Double_t axbj=0.0,aQ2=0.0,aW2=0.0;
                T->SetBranchAddress("EKRxe.x_bj",&axbj);
                T->SetBranchAddress("EKRxe.Q2",&aQ2);
                T->SetBranchAddress("EKRxe.W2",&aW2);
   	    
                Int_t nentries=electron->GetN();
                for(int nn=0;nn<nentries;nn++){
                   T->GetEntry(electron->GetEntry(nn));
     	     	   hW2[ii]->SetPoint(npoints,axbj,aW2);
     	     	   hQ2[ii]->SetPoint(npoints,axbj,aQ2);
   	 	   npoints++;
       	        }
	     }
	     delete T;
	 }

     }

     TCanvas *c1=new TCanvas("c1","c1",1500,1200);
     int color[12]={1,2,3,4,5,6,7,8,9,46,30,38};
     TMultiGraph *mg=new TMultiGraph();
     for(int ii=0;ii<12;ii++){
	hW2[ii]->SetMarkerColor(color[ii]);
	hW2[ii]->SetMarkerStyle(7);
	mg->Add(hW2[ii]);
     }
     mg->Draw("AP");
     mg->GetXaxis()->SetLimits(0,1);
     mg->GetYaxis()->SetRangeUser(0,14);
     mg->SetTitle(";x;{W}_{2};");

     auto leg1=new TLegend(0.15,0.15,0.25,0.65);
     for(int ii=0;ii<12;ii++){
      leg1->AddEntry(hW2[ii],Form("kin%d",kin[ii]),"P");
     }
     leg1->Draw();
     leg1->SetBorderSize(0);

     //c1->Print("W2_X.pdf");

     TCanvas *c2=new TCanvas("c2","c2",1500,1200);
     TMultiGraph *mg1=new TMultiGraph();
     for(int ii=0;ii<12;ii++){
	hQ2[ii]->SetMarkerColor(color[ii]);
	hQ2[ii]->SetMarkerStyle(7);
	mg1->Add(hQ2[ii]);
     }
     mg1->Draw("AP");
     mg1->GetXaxis()->SetLimits(0,1);
     mg1->GetYaxis()->SetRangeUser(0,14);
     mg1->SetTitle(";x;{Q}_{2};");

     auto leg2=new TLegend(0.15,0.35,0.25,0.85);
     for(int ii=0;ii<12;ii++){
      leg2->AddEntry(hQ2[ii],Form("kin%d",kin[ii]),"P");
     }
     leg2->Draw();
     leg2->SetBorderSize(0);

     //c2->Print("Q2_X.pdf");

}
