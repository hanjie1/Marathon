void PlotfH(){
     Double_t tau=4500;
     Double_t H3_t0=0.085099;
     Double_t He3_t0=0.0000226;

     TCanvas *c1=new TCanvas("c1","c1",1500,800);
     auto f1=new TF1("f1",Form("(%f+%f*(1.0-TMath::Exp(-x/%f)))/(%f+%f)",He3_t0,H3_t0,tau,He3_t0,H3_t0),83,175);  
     f1->GetHistogram()->GetXaxis()->SetTitle("Days");
     f1->GetHistogram()->GetYaxis()->SetTitle("f_{He}");
     f1->GetHistogram()->SetTitle("");

     f1->SetLineColor(8);

     f1->Draw();
     c1->Print("Thesis_plots/fH.pdf");

}
