#include "SetCut.h"
#include "GetTrees.h"
void Make_PID_plot(){
	TCanvas *c1=new TCanvas("c1","c1",1500,1000);
	gStyle->SetOptStat(0);
	gStyle->SetTitleY(0.05);
	gStyle->SetTitleX(0.5);

	int run_number=2503;
        TChain *T;
        T=GetTree(run_number,0,"T");

        TCut VZ=Form("L.tr.vz<%f && L.tr.vz>%f",vz_max[0],vz_min[0]);
        TCut CK_e="L.cer.asum_c>1500";
        TCut CK_pi="L.cer.asum_c<200";


        TH1F *hEp=new TH1F("hEp","electron E/p distribution",100,0,1.4);
        TH1F *hCer=new TH1F("hCer","electron Cherenkov distribution",200,0,10000);

	c1->Divide(2,1);
	c1->cd(1);
        T->Draw("L.cer.asum_c>>hCer",TRK+ACC+trigger2+VZ+beta+W2);
        hCer->SetLineWidth(1);
	hCer->SetTitle("(a);Cherenkov sum;counts;");

	TArrow *ar1=new TArrow(1500,100,1500,70,0.01);//,"|>");
	ar1->SetLineWidth(2);	
	ar1->Draw();

	c1->cd(2);
        T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>hEp",TRK+ACC+trigger2+VZ+beta+W2);
        hEp->SetLineWidth(1);
	hEp->SetTitle("(b);E/P;counts;");
	
	TArrow *ar2=new TArrow(0.7,600,0.7,300,0.01);//,"|>");	
	ar2->SetLineWidth(2);	
	ar2->Draw();

	c1->Print("PID_results/PIDplot_1.pdf");
}
