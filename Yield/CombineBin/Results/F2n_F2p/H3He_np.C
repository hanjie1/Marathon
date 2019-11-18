#include "ReadFile.h"
#include "libEMCR.h"

int H3He_np(){
   Double_t x[19]={0.0},Q2[19]={0.0};
   Double_t H3He[19]={0.0},H3He_err[19]={0.0},H3He_ST[19]={0.0},H3He_SY[19]={0.0};
   Double_t x_KP[110]={0.0},F2n_KP[110]={0.0},F2p_KP[110]={0.0};
   Double_t x_Wit[7]={0.0},Q2_Wit[7]={0.0},F2np_Wit_P[7]={0.0},F2np_Wit_R[7]={0.0},F2np_Wit_B[7]={0.0},F2np_Wit_E[7]={0.0};
   Double_t x_Wally[6]={0.0},F2np_Wally[6]={0.0};

   ifstream file_Wit;
   file_Wit.open("../Other_Data/Whitlow_NP.dat");
   if(!file_Wit.is_open())return 0;

   TString content,tmp;
   Ssiz_t from=0;
   int nn=0;
   while(tmp.ReadLine(file_Wit)){
          if(nn==0){nn++;continue;}
          if(tmp[0]=='#')continue;
          tmp.Tokenize(content,from," ");
          x_Wit[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Q2_Wit[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          F2np_Wit_P[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          F2np_Wit_R[nn-1]=F2np_Wit_P[nn-1]*atof(content.Data());
          tmp.Tokenize(content,from," ");
          F2np_Wit_B[nn-1]=F2np_Wit_P[nn-1]*atof(content.Data());
          tmp.Tokenize(content,from," ");
          F2np_Wit_E[nn-1]=atof(content.Data());

          from=0;
          nn++;
     }
    file_Wit.close();

   ifstream file_Wally;
   file_Wally.open("../Other_Data/Wally_NP.dat");
   if(!file_Wally.is_open())return 0;

   from=0;
   nn=0;
   while(tmp.ReadLine(file_Wally)){
          if(nn==0){nn++;continue;}
          if(tmp[0]=='#')continue;
          tmp.Tokenize(content,from,"  ");
          x_Wally[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          F2np_Wally[nn-1]=atof(content.Data());

          from=0;
          nn++;
     }
    file_Wally.close();


   TString Rfile="newbin/H3He_final.dat";
   int nbin=ReadYieldFinal(Rfile,x,Q2,H3He,H3He_ST,H3He_SY);
   Rfile="Model/F2dis_os1tm1ht1mec1_Dav18_He3Salme";
   int nbin_KP=ReadKP(Rfile,x_KP,F2p_KP,F2n_KP);

   Double_t Nc=0.976;

    ofstream outfile;
    outfile.open("results/F2np_final.dat");
    outfile<<"x   Q2   n/p    total_err     stat_err     sys_err"<<endl;;

    ofstream outfile1;
    outfile1.open("ForThesis.dat");

    Double_t H3He_np[19]={0.0},H3He_npErr[19]={0.0};
    Double_t np_ST[19]={0.0},np_SY[19]={0.0};
    for(int ii=0;ii<19;ii++){
	Double_t He_R=He_EMC(x[ii]);
	Double_t H3_R=H3_EMC(x[ii]);
	Double_t SR=H3_R/He_R;

	Double_t tmpH3He=H3He[ii]*Nc;
	Double_t tmpST=H3He_ST[ii]*Nc;
	Double_t tmpSY=H3He_SY[ii]*Nc;
	Double_t tmpH3He_err=sqrt(tmpST*tmpST+tmpSY*tmpSY);
 	
	H3He_np[ii]=(2.0*tmpH3He-SR)/(2.0*SR-tmpH3He);
	H3He_npErr[ii]=abs(3.0*SR/((2.0*SR-tmpH3He)*(2.0*SR-tmpH3He)))*tmpH3He_err;
	np_ST[ii]=abs(3.0*SR/((2.0*SR-tmpH3He)*(2.0*SR-tmpH3He)))*tmpST;
	np_SY[ii]=abs(3.0*SR/((2.0*SR-tmpH3He)*(2.0*SR-tmpH3He)))*tmpSY;
	outfile<<x[ii]<<"  "<<Q2[ii]<<"  "<<H3He_np[ii]<<"  "<<H3He_npErr[ii]<<"  "<<np_ST[ii]<<"  "<<np_SY[ii]<<endl;
	outfile1<<fixed<<setprecision(2)<<x[ii]<<" & "<<Q2[ii]<<" & "<<setprecision(3)<<H3He_np[ii]<<" & "<<np_ST[ii]<<" &  "<<np_SY[ii]<<"  \\\\"<<endl;
	outfile1<<"\\hline"<<endl;
    }
    outfile.close();
    outfile1.close();
    TGraphErrors *gH3He=new TGraphErrors(19,x,H3He_np,0,H3He_npErr);

    TGraphErrors *gWally=new TGraphErrors(6,x_Wally,F2np_Wally,0,0);

    TGraphErrors *gWit_P=new TGraphErrors(7);
    TGraphErrors *gWit_R=new TGraphErrors(7);
    TGraphErrors *gWit_B=new TGraphErrors(7);
    TGraphErrors *gWit_E=new TGraphErrors(7);

    for(int ii=0;ii<7;ii++){
	Double_t np_err;
	if(x_Wit[ii]<=0.75)np_err=0.01;
	else np_err=0.03;
        gWit_P->SetPoint(ii,x_Wit[ii],F2np_Wit_P[ii]);
        gWit_P->SetPointError(ii,0,np_err);
        gWit_R->SetPoint(ii,x_Wit[ii],F2np_Wit_R[ii]);
        gWit_R->SetPointError(ii,0,np_err);
        gWit_B->SetPoint(ii,x_Wit[ii],F2np_Wit_B[ii]);
        gWit_B->SetPointError(ii,0,np_err);
        gWit_E->SetPoint(ii,x_Wit[ii],F2np_Wit_E[ii]);
        gWit_E->SetPointError(ii,0,np_err);
    }

    TGraphErrors *gKP=new TGraphErrors();
    for(int ii=0;ii<nbin_KP;ii++){
	if(x_KP[ii]>0.85)continue;
	gKP->SetPoint(ii,x_KP[ii],F2n_KP[ii]/F2p_KP[ii]);
    }

   gStyle->SetEndErrorSize(4);

    TCanvas *c1=new TCanvas("c1","c1",1500,1200);
    TMultiGraph *mg1=new TMultiGraph();
    gH3He->SetMarkerColor(2);
    gH3He->SetMarkerStyle(8);
    gH3He->SetMarkerSize(2);
    gH3He->SetLineColor(2);
    gKP->SetLineColor(4);
    gKP->SetLineStyle(1);
    gKP->SetLineWidth(2);

    gWit_P->SetMarkerColor(1);
    gWit_P->SetMarkerStyle(24);
    gWit_P->SetMarkerSize(2);
    gWit_P->SetLineColor(1);
    gWit_R->SetMarkerColor(4);
    gWit_R->SetLineColor(4);
    gWit_R->SetMarkerStyle(25);
    gWit_R->SetMarkerSize(2);
    gWit_B->SetMarkerColor(8);
    gWit_B->SetLineColor(8);
    gWit_B->SetMarkerStyle(26);
    gWit_B->SetMarkerSize(2);
    gWit_E->SetMarkerColor(6);
    gWit_E->SetLineColor(6);
    gWit_E->SetMarkerStyle(27);
    gWit_E->SetMarkerSize(2);

    gWally->SetMarkerColor(46);
    gWally->SetLineColor(46);
    gWally->SetMarkerStyle(46);
    gWally->SetMarkerSize(2);

    mg1->Add(gH3He,"P");
    mg1->Add(gWit_P,"P");
    mg1->Add(gWit_R,"P");
    mg1->Add(gWit_B,"P");
//    mg1->Add(gWit_E,"P");
    mg1->Add(gWally,"P");
    //mg1->Add(gKP,"L");
    mg1->Draw("A"); 
    mg1->SetTitle(";Bjorken x;F_{2}^{n} / F_{2}^{p}");
    mg1->GetYaxis()->SetLabelOffset(0.0005);
//    mg1->GetYaxis()->SetRangeUser(0.3,0.85);
 
   auto leg1=new TLegend(0.48,0.7,0.9,0.9);
   leg1->SetNColumns(2);
   leg1->AddEntry(gH3He,"#scale[1.8]{MARATHON}","P");
   //leg1->AddEntry(gKP,"#scale[1]{KP model}","L");
   leg1->AddEntry(gWally,"#scale[1.8]{Melnitchouk and Thomas}","P");
   leg1->AddEntry(gWit_P,"#scale[1.8]{Whitlow et al. (Paris)}","P");
   leg1->AddEntry(gWit_R,"#scale[1.8]{Whitlow et al. (Reid)}","P");
   leg1->AddEntry(gWit_B,"#scale[1.8]{Whitlow et al. (Bonn)}","P");
//   leg1->AddEntry(gWit_E,"#scale[1.8]{Whitlow et al. (Frankfurt and Strikman)}","P");
   leg1->SetMargin(0.2);
   leg1->Draw();

   c1->Print("F2np_final_compare.pdf");

    return 0;
}
