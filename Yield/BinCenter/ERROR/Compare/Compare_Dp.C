#include "ReadFile.h"
void Compare_Dp(){
     Double_t x[MAXNUM]={0.0},BC111[MAXNUM]={0.0},BC211[MAXNUM]={0.0};
     Double_t BC121[MAXNUM]={0.0},BC122[MAXNUM]={0.0},BC123[MAXNUM]={0.0};

     TString filename;
     filename="../OUT/model111/Dp_BCfac.dat";
     ReadFile(filename,x,BC111);
     filename="../OUT/model211/Dp_BCfac.dat";
     ReadFile(filename,x,BC211);
     filename="../OUT/model121/Dp_BCfac.dat";
     ReadFile(filename,x,BC121);
     filename="../OUT/model122/Dp_BCfac.dat";
     ReadFile(filename,x,BC122);
     filename="../OUT/model123/Dp_BCfac.dat";
     ReadFile(filename,x,BC123);

     TGraph *gDp[4];
     int color[4]={1,2,4,8};
     for(int ii=0;ii<4;ii++){
	 gDp[ii]=new TGraph();
         gDp[ii]->SetMarkerStyle(8);
         gDp[ii]->SetMarkerColor(color[ii]);
     }
     for(int ii=0;ii<MAXNUM;ii++){
	 if(x[ii]==0)continue;
	 Double_t ratio=0.0;
	 ratio=BC211[ii]/BC111[ii];
	 gDp[0]->SetPoint(ii,x[ii],ratio);
	 ratio=BC121[ii]/BC111[ii];
	 gDp[1]->SetPoint(ii,x[ii],ratio);
	 ratio=BC122[ii]/BC111[ii];
	 gDp[2]->SetPoint(ii,x[ii],ratio);
	 ratio=BC123[ii]/BC111[ii];
	 gDp[3]->SetPoint(ii,x[ii],ratio);
     }

     TMultiGraph *mg=new TMultiGraph();
     for(int ii=0;ii<4;ii++) mg->Add(gDp[ii]);
     mg->Draw("AP");
     mg->SetTitle("D/p BC factor ratio;x");

     auto leg1=new TLegend(0.7,0.6,0.85,0.85);
     leg1->AddEntry(gDp[0],"model211","P");
     leg1->AddEntry(gDp[1],"model121","P");
     leg1->AddEntry(gDp[2],"model122","P");
     leg1->AddEntry(gDp[3],"model123","P");

     leg1->Draw();





}
