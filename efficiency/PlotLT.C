#include "GetTrees.h"
#include "SetCut.h"
#include "SetCons.h"
Double_t CalcLT_L(int run_number,int kin,int beamcut=0)
{
     TString TreeName="T";
     TChain *T=GetTree(run_number,kin,TreeName);
     TChain *EndT=GetTree(run_number,kin,"EndLeft");
     int noend=0;
     if(EndT==NULL)noend=1;

     Double_t Nscaler=0;
     Double_t Nmeas=0;
     Double_t Livetime=0;
     Double_t Livetime_err=0;
     if(beamcut==0){
	if(noend)Nscaler = T->GetMaximum("evLeftT2");
	else Nscaler = EndT->GetMaximum("EndLeftT2");
        Nmeas = T->GetEntries(trigger2);
        if(Nscaler!=0){
           Livetime=Nmeas/Nscaler;
	   Livetime_err=sqrt(1.0/Nmeas+1.0/Nscaler)*Livetime;
	}
        cout<<Livetime<<endl;
     }

     if(beamcut==1){
        T->SetBranchStatus("*",0);
        T->SetBranchStatus("evLeftT2",1);     
        T->SetBranchStatus("DL.evtypebits",1); 
        T->SetBranchStatus("evLeftdnew_r",1);
        T->SetBranchStatus("LeftBCMev.isrenewed",1);

	Double_t evT2,evtypebits,current_dnew,isrenewed;
        Double_t beamUp[5];
        T->SetBranchAddress("evLeftT2",&evT2);
        T->SetBranchAddress("DL.evtypebits",&evtypebits);
        T->SetBranchAddress("evLeftdnew_r",&current_dnew);
        T->SetBranchAddress("LeftBCMev.isrenewed",&isrenewed);

        Double_t maxT2;
        if(!noend)maxT2=EndT->GetMaximum("EndLeftT2");
        Double_t nentries=T->GetEntries();

	Double_t lastcount=0;
        Double_t entryT2=0;
        Nscaler=0.0;Nmeas=0.0;
        for(int ii=0;ii<nentries;ii++){
	    T->GetEntry(ii);
            Int_t trigger=(Int_t)evtypebits;                
	    if((trigger>>2)&1)entryT2=entryT2+1;

	    if(isrenewed){
               current_dnew=current_dnew*dnew_gain;
               if(current_dnew>0){
                  Nscaler=Nscaler+evT2-lastcount;
                  Nmeas=Nmeas+entryT2;
               }
               entryT2=0;
               lastcount=evT2;   
            }
        }
        
        if(maxT2!=lastcount && noend==0)Nscaler=Nscaler+maxT2-lastcount;
        Livetime=Nmeas/Nscaler;
        Livetime_err=sqrt(1.0/Nmeas+1.0/Nscaler)*Livetime;
        cout<<Livetime<<endl;

     }

     Livetime=1.0/Livetime;
     delete T;
     return Livetime;
}

Double_t CalcLT_R(int run_number,int kin,int beamcut=0)
{
     TString TreeName="T";
     TChain *T=GetTree(run_number,kin,TreeName);
     TChain *EndT=GetTree(run_number,kin,"EndRight");
     int noend=0;
     if(EndT==NULL)noend=1;

     Double_t Nscaler=0;
     Double_t Nmeas=0;
     Double_t Livetime=0;
     Double_t Livetime_err=0;
     if(beamcut==0){
	if(noend)Nscaler = T->GetMaximum("evRightT5");
	else Nscaler = EndT->GetMaximum("EndRightT5");
	TCut trigger5="(DR.evtypebits>>5)&1";
        Nmeas = T->GetEntries(trigger5);
        if(Nscaler!=0){
           Livetime=Nmeas/Nscaler;
	   Livetime_err=sqrt(1.0/Nmeas+1.0/Nscaler)*Livetime;
	}
        cout<<Livetime<<endl;
     }

     if(beamcut==1){
        T->SetBranchStatus("*",0);
        T->SetBranchStatus("evRightT5",1);     
        T->SetBranchStatus("DR.evtypebits",1); 
        T->SetBranchStatus("evRightdnew_r",1);
        T->SetBranchStatus("RightBCMev.isrenewed",1);

	Double_t evT2,evtypebits,current_dnew,isrenewed;
        Double_t beamUp[5];
        T->SetBranchAddress("evRightT5",&evT2);
        T->SetBranchAddress("DR.evtypebits",&evtypebits);
        T->SetBranchAddress("evRightdnew_r",&current_dnew);
        T->SetBranchAddress("RightBCMev.isrenewed",&isrenewed);

        Double_t maxT2;
        if(!noend)maxT2=EndT->GetMaximum("EndRightT5");
        Double_t nentries=T->GetEntries();

	Double_t lastcount=0;
        Double_t entryT2=0;
        Nscaler=0.0;Nmeas=0.0;
        for(int ii=0;ii<nentries;ii++){
	    T->GetEntry(ii);
            Int_t trigger=(Int_t)evtypebits;                
	    if((trigger>>5)&1)entryT2=entryT2+1;

	    if(isrenewed){
               current_dnew=current_dnew*dnew_gain;
               if(current_dnew>0){
                  Nscaler=Nscaler+evT2-lastcount;
                  Nmeas=Nmeas+entryT2;
               }
               entryT2=0;
               lastcount=evT2;   
            }
        }
        
        if(maxT2!=lastcount && noend==0)Nscaler=Nscaler+maxT2-lastcount;
        Livetime=Nmeas/Nscaler;
        Livetime_err=sqrt(1.0/Nmeas+1.0/Nscaler)*Livetime;
        cout<<kin<<"  "<<Livetime<<endl;

     }

     delete T;
     Livetime=1.0/Livetime;
     return Livetime;
}



void PlotLT(){
     Double_t kin[12]={0,1,2,3,4,5,7,9,11,13,15,16};

     Double_t H1_DT[5]={0.0},D2_DT[12]={0.0},He_DT[12]={0.0},H3_DT[12]={0.0};
     int H1_nrun[5]={2479,1213,1240,1286,2582};
     int D2_nrun[12]={2493,1214,1242,1288,2579,1350,1371,1405,1493,1639,1862,91280};
     int He_nrun[12]={2531,1219,1236,1283,2576,1353,1373,1408,1485,1632,1868,91365};
     int H3_nrun[12]={2503,1224,1245,1290,2577,1355,1377,1413,1489,1635,1880,91369};

     TGraph *gH1=new TGraph(5);
     TGraph *gD2=new TGraph(12);
     TGraph *gHe=new TGraph(12);
     TGraph *gH3=new TGraph(12);

     for(int ii=0;ii<11;ii++){
         if(ii<5){
	    H1_DT[ii]=CalcLT_L(H1_nrun[ii],kin[ii],1);
	    gH1->SetPoint(ii,kin[ii],H1_DT[ii]);
	 }
         D2_DT[ii]=CalcLT_L(D2_nrun[ii],kin[ii],1);
         He_DT[ii]=CalcLT_L(He_nrun[ii],kin[ii],1);
         H3_DT[ii]=CalcLT_L(H3_nrun[ii],kin[ii],1);

	 gD2->SetPoint(ii,kin[ii],D2_DT[ii]);
	 gHe->SetPoint(ii,kin[ii],He_DT[ii]);
	 gH3->SetPoint(ii,kin[ii],H3_DT[ii]);
     }

         D2_DT[11]=CalcLT_R(D2_nrun[11],kin[11],1);
         He_DT[11]=CalcLT_R(He_nrun[11],kin[11],1);
         H3_DT[11]=CalcLT_R(H3_nrun[11],kin[11],1);

	 gD2->SetPoint(11,kin[11],D2_DT[11]);
	 gHe->SetPoint(11,kin[11],He_DT[11]);
	 gH3->SetPoint(11,kin[11],H3_DT[11]);
    
	gH1->SetMarkerStyle(8); 
	gH1->SetMarkerColor(1); 
	gD2->SetMarkerStyle(8); 
	gD2->SetMarkerColor(3); 
	gHe->SetMarkerStyle(8); 
	gHe->SetMarkerColor(4); 
	gH3->SetMarkerStyle(8); 
	gH3->SetMarkerColor(2); 

	TCanvas *c1=new TCanvas("c1");
	TMultiGraph *mg=new TMultiGraph();
	mg->Add(gH1);
	mg->Add(gD2);
	mg->Add(gHe);
	mg->Add(gH3);
	mg->Draw("AP");
	mg->GetXaxis()->SetTitle("kinematic");
	mg->GetYaxis()->SetTitle("C_{DT}");

        auto leg1=new TLegend(0.7,0.6,0.85,0.85);
 	leg1->AddEntry(gH1,"H1","P");
        leg1->AddEntry(gD2,"D2","P");
        leg1->AddEntry(gHe,"He3","P");
        leg1->AddEntry(gH3,"H3","P");
        leg1->Draw();

	c1->Print("DT.pdf");
	
}



