void PlotfH(){
     Double_t tau=4500;
     Double_t H3_t0=0.085099;
     Double_t He3_t0=0.0000226;

     auto f1=new TF1("f1",Form("%f+%f*(1.0-TMath::Exp(-x/%f))",He3_t0,H3_t0,tau),83,175);  

     f1->SetLineColor(8);

     f1->Draw();

}
