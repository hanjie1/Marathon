#include "ReadFile.h"
#include "libEMCISO.h"
void He3_EMC(){
     Double_t x[19]={0.0},Q2[19]={0.0},F2np[19]={0.0},F2np_err[19]={0.0},F2np_ST[19]={0.0},F2np_SY[19]={0.0};
     Double_t HeD[19]={0.0},HeD_err[19]={0.0},HeD_ST[19]={0.0},HeD_SY[19]={0.0};
     Double_t HeD_iso[19]={0.0},HeD_err_iso[19]={0.0},HeD_iso_ST[19]={0.0},HeD_iso_SY[19]={0.0};
     Double_t pHeD[19]={0.0},pHeD_err[19]={0.0};
     Double_t pHeD_norm[19]={0.0},pHeD_err_norm[19]={0.0};
     Double_t xC[26]={0.0},HeD_C[26]={0.0},HeD_Cerr[26]={0.0},HeDC_ST[26]={0.0},HeDC_SY[26]={0.0};
     Double_t Hermes_x[12]={0.0},Hermes_HeD[12]={0.0},Hermes_HeDerr[12]={0.0},Hermes_HeDC_ST[12]={0.0},Hermes_HeDC_SY[12]={0.0};

     Double_t Z=2.0,A=3.0;
     Double_t Nc=0.024;

     TString Rfile;
     Rfile="F2n_F2p/results/F2np_final.dat";
     int nbin1=ReadNP(Rfile,x,F2np,F2np_err,F2np_ST,F2np_SY);
     Rfile="newbin/HeD_final.dat";
     int nbin2=ReadYield(Rfile,x,Q2,HeD,HeD_err,HeD_ST,HeD_SY);
     Rfile="Other_Data/HeD_HallC.dat";
     int nbin3=ReadHallC(Rfile,xC,HeD_C,HeDC_ST,HeDC_SY);
     Rfile="Other_Data/hermes_He3D.dat";
     int nbin4=ReadHERMES(Rfile,Hermes_x,Hermes_HeD,Hermes_HeDC_ST,Hermes_HeDC_SY);

     ofstream outfile;
     outfile.open("He_EMC_iso.dat");
     for(int ii=0;ii<19;ii++){
	 Double_t tmpIso=(1.0+F2np[ii])/(Z+(A-Z)*F2np[ii]);
  	 pHeD[ii]=HeD[ii]*2.0/A;	
  	 pHeD_err[ii]=HeD_err[ii]*2.0/A;	
  	 pHeD_norm[ii]=pHeD[ii]*(1.+Nc);	
  	 pHeD_err_norm[ii]=pHeD_err[ii]*(1.+Nc);	

	 HeD_iso[ii]=HeD[ii]*tmpIso*(1.0+Nc);
	 HeD_err_iso[ii]=HeD_iso[ii]*sqrt(pow(HeD_err[ii]/HeD[ii],2)
			+pow((2.0*Z-A)/((1+F2np[ii])*(Z+(A-Z)*F2np[ii]))*F2np_err[ii],2));
	 HeD_iso_ST[ii]=HeD_iso[ii]*sqrt(pow(HeD_ST[ii]/HeD[ii],2)
                        +pow((2.0*Z-A)/((1+F2np[ii])*(Z+(A-Z)*F2np[ii]))*F2np_ST[ii],2));
	 HeD_iso_SY[ii]=HeD_iso[ii]*sqrt(pow(HeD_SY[ii]/HeD[ii],2)
                        +pow((2.0*Z-A)/((1+F2np[ii])*(Z+(A-Z)*F2np[ii]))*F2np_SY[ii],2));
	 outfile<<fixed<<setprecision(2)<<x[ii]<<" & "<<Q2[ii]<<" & "<<setprecision(3)<<HeD_iso[ii]<<" & "<<HeD_iso_ST[ii]<<" & "<<HeD_iso_SY[ii]<<"  \\\\"<<endl;
	 outfile<<"\\hline"<<endl; 
     } 
     outfile.close();

     TGraphErrors *gHeD_C=new TGraphErrors();
     for(int ii=0;ii<nbin3;ii++){
	 if(xC[ii]>0.85)continue;
	 HeD_Cerr[ii]=sqrt(HeDC_ST[ii]*HeDC_ST[ii]+HeDC_SY[ii]*HeDC_SY[ii]);
	 gHeD_C->SetPoint(ii,xC[ii],HeD_C[ii]);
	 gHeD_C->SetPointError(ii,0.0,HeD_Cerr[ii]);
     }
   TGraphErrors *HallC_norm=new TGraphErrors(1);
   HallC_norm->SetPoint(0,xC[0],0.9);
   HallC_norm->SetPointError(0,0,HeD_C[0]*0.0184);


     TGraphErrors *gHeD_herme=new TGraphErrors();
     for(int ii=0;ii<nbin4;ii++){
	 Hermes_HeDerr[ii]=sqrt(pow(Hermes_HeDC_ST[ii],2)+pow(Hermes_HeDC_SY[ii],2));
	 gHeD_herme->SetPoint(ii,Hermes_x[ii],Hermes_HeD[ii]);
	 gHeD_herme->SetPointError(ii,0.0,Hermes_HeDerr[ii]);
     }

     TGraphErrors *gHeD=new TGraphErrors(19,x,pHeD,0,pHeD_err);
     TGraphErrors *gHeD_iso=new TGraphErrors(19,x,HeD_iso,0,HeD_err_iso);
     TGraphErrors *gHeD_norm=new TGraphErrors(19,x,pHeD_norm,0,pHeD_err_norm);

     auto f1_KP=new TF1("f1_KP","He_ISO(x)",0.16,0.85);
     auto f1_SLAC=new TF1("f1_SLAC","SLAC_EMC(x)",0.16,0.85);
     auto f1_SLAC_den=new TF1("f1_SLAC_den","SLAC_EMC_Den(x,3.0,2.0)",0.16,0.85);

     TLine *l1=new TLine(0,1,0.9,1);
     l1->SetLineColor(1);
     l1->SetLineStyle(7);

   gStyle->SetEndErrorSize(4);
     TCanvas *c1=new TCanvas("c1","c1",1500,1200);
     TMultiGraph *mg=new TMultiGraph();
     gHeD->SetMarkerStyle(8);
     gHeD->SetMarkerColor(4);
     gHeD->SetMarkerSize(1.5);
     gHeD->SetLineColor(4);
     gHeD_iso->SetMarkerStyle(8);
     gHeD_iso->SetMarkerColor(2);
     gHeD_iso->SetMarkerSize(1.5);
     gHeD_iso->SetLineColor(2);
     gHeD_norm->SetMarkerStyle(8);
     gHeD_norm->SetMarkerColor(4);
     gHeD_norm->SetMarkerSize(1.5);
     gHeD_norm->SetLineColor(4);
     gHeD_C->SetMarkerStyle(21);
     gHeD_C->SetMarkerColor(8);
     gHeD_C->SetMarkerSize(1.5);
     gHeD_C->SetLineColor(8);
     gHeD_herme->SetMarkerStyle(25);
     gHeD_herme->SetMarkerColor(1);
     gHeD_herme->SetMarkerSize(1.5);
     gHeD_herme->SetLineColor(1);
     HallC_norm->SetLineColor(8);
     HallC_norm->SetLineWidth(2);
//     mg->Add(gHeD,"P");
     mg->Add(gHeD_iso,"P");
     mg->Add(gHeD_C,"P");
     mg->Add(gHeD_norm,"P");
     mg->Add(gHeD_herme,"P");
     mg->Add(HallC_norm,"L");
     mg->Draw("A");
     mg->SetTitle(";Bjorken x;(F_{2}^{^{3}He}/F_{2}^{^{2}H})_{iso}");
     mg->GetYaxis()->SetRangeUser(0.85,1.25);
     mg->GetXaxis()->SetRangeUser(0,0.9);

     f1_KP->SetLineColor(12);
     f1_KP->SetLineStyle(1);
     f1_KP->Draw("same");

     f1_SLAC->SetLineColor(12);
     f1_SLAC->SetLineStyle(3);
     f1_SLAC->Draw("same");

     f1_SLAC_den->SetLineColor(12);
     f1_SLAC_den->SetLineStyle(7);
     f1_SLAC_den->Draw("same");

     l1->Draw("same");

     
   auto leg1=new TLegend(0.15,0.75,0.7,0.9);
   leg1->SetNColumns(3);
   leg1->AddEntry(gHeD_iso,"#scale[2]{MARATHON}","P");
   leg1->AddEntry(gHeD_norm,"#scale[2]{MARATHON no Iso. Cor.}","P");
   leg1->AddEntry(gHeD_C,"#scale[2]{E03-103}","P");
   leg1->AddEntry(gHeD_herme,"#scale[2]{HERMES}","P");
   leg1->AddEntry(f1_KP,"#scale[2]{KP model}","L");
   leg1->AddEntry(f1_SLAC,"#scale[2]{SLAC A fit}","L");
   leg1->AddEntry(f1_SLAC_den,"#scale[2]{SLAC density model}","L");
   leg1->SetMargin(0.3);
   leg1->Draw();

   TLatex latex;
   latex.SetTextSize(0.025);
   latex.DrawLatex(xC[0]+0.1,0.9,"E03-103 Norm. (1.84%)");

     c1->Print("He_EMC.pdf");
}
