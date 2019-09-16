void PlotPos(){
     auto f1=new TF1("f1","1-TMath::Exp(-10.1463193*x-2.38884408)",0.15,0.85);  
     auto f2=new TF1("f2","1-TMath::Exp(-9.15721674*x-2.46454072)",0.15,0.85);  
     auto f3=new TF1("f3","1-TMath::Exp(-8.42406685*x-2.63480349)",0.15,0.85);  
     auto f4=new TF1("f4","1-TMath::Exp(-8.42875968*x-2.61143195)",0.15,0.85);  

     f1->GetYaxis()->SetRangeUser(0.975,1.005);

     f1->SetLineColor(1);
     f2->SetLineColor(4);
     f3->SetLineColor(8);
     f4->SetLineColor(2);

     TCanvas *c1=new TCanvas("c1","c1",1500,1000);
     f1->Draw();
     f2->Draw("same");
     f3->Draw("same");
     f4->Draw("same");
     f1->GetHistogram()->GetXaxis()->SetTitle("Bjorken x");
     f1->GetHistogram()->GetYaxis()->SetTitle("C_{e^{+}}");
     f1->GetHistogram()->SetTitle("");

     auto legend = new TLegend(0.7,0.3,0.85,0.55);
     legend->AddEntry("f1","{}^{1}H","l");
     legend->AddEntry("f2","{}^{2}H","l");
     legend->AddEntry("f3","{}^{3}He","l");
     legend->AddEntry("f4","{}^{3}H","l");
     legend->Draw();

     c1->Print("Thesis_plots/positron.pdf");
}
