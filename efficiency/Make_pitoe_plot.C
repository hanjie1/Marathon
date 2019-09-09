#include "SetCut.h"
#include "GetTrees.h"
void Make_pitoe_plot(){
	TCanvas *c1=new TCanvas("c1","c1",1500,1200);
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
        TH1F *hEp_e=new TH1F("hEp_e","electron E/p distribution",100,0,1.4);
        TH1F *hEp_pi=new TH1F("hEp_pi","pion E/p distribution",100,0,1.4);

        Double_t ep_cut1=0.1;
        Double_t ep_cut2=0.7;

	c1->Divide(3,1);
	c1->cd(1);
        T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>hEp",TRK+ACC+trigger2+VZ+beta+W2);
        hEp->SetLineWidth(1);
	hEp->SetTitle("(a);E/P;counts;");
        gPad->SetLogy();

	c1->cd(2);
        T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>hEp_e",TRK+ACC+trigger2+VZ+beta+CK_e+W2);
        hEp_e->SetLineWidth(1);
	hEp_e->SetTitle("(b);E/P;counts;");
	

        int nbin1=hEp_e->FindBin(ep_cut1);
        Double_t ymax1 = hEp_e->GetBinContent(nbin1);	
	TLine *l1=new TLine(ep_cut1,0,ep_cut1,ymax1);
	l1->SetLineStyle(9);
	l1->SetLineWidth(1);
	l1->SetLineColor(9);
	l1->Draw();

        int nbin2=hEp_e->FindBin(ep_cut2);
        Double_t ymax2 = hEp_e->GetBinContent(nbin2);
	TLine *l2=new TLine(ep_cut2,0,ep_cut2,ymax2);
	l2->SetLineStyle(7);
	l2->SetLineWidth(1);
	l2->SetLineColor(9);
	l2->Draw();

	TLatex latex1;
	latex1.SetTextSize(0.07);
	//latex1.SetTextAlign(12);
	latex1.DrawLatex(0.4,4,"1");
	latex1.DrawLatex(0.9,4,"2");

        gPad->SetLogy();

	c1->cd(3);
        T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>hEp_pi",TRK+ACC+trigger1+VZ+beta+CK_pi+W2);

        hEp_pi->SetLineWidth(1);
	hEp_pi->SetTitle("(c);E/P;counts;");

        int nbin3=hEp_pi->FindBin(ep_cut1);
        Double_t ymax3 = hEp_pi->GetBinContent(nbin3);	
	TLine *l3=new TLine(ep_cut1,0,ep_cut1,ymax3);
	l3->SetLineStyle(9);
	l3->SetLineWidth(1);
	l3->SetLineColor(9);
	l3->Draw();

        int nbin4=hEp_pi->FindBin(ep_cut2);
        Double_t ymax4 = hEp_pi->GetBinContent(nbin4);
	TLine *l4=new TLine(ep_cut2,0,ep_cut2,ymax4);
	l4->SetLineStyle(7);
	l4->SetLineWidth(1);
	l4->SetLineColor(9);
	l4->Draw();
        gPad->SetLogy();

	TLatex latex2;
	latex2.SetTextSize(0.07);
//	latex2.SetTextAlign(12);
	latex2.DrawLatex(0.4,2,"3");
	latex2.DrawLatex(0.9,2,"4");

	c1->Print("PID_results/PIDplot.pdf");
}
