void PlotECC(){
     auto f1=new TF1("f1","1-TMath::Exp(-2.54412287*x-3.69437788)",0.15,0.85);  
     auto f2=new TF1("f2","1+TMath::Exp(-5.03250363*x-3.97644112)",0.15,0.85);  
     auto f3=new TF1("f3","1+TMath::Exp(-3.90944653*x-3.40802943)",0.15,0.85);  
     auto f4=new TF1("f4","(1-TMath::Exp(-2.54412287*x-3.69437788))*(1+TMath::Exp(-5.03250363*x-3.97644112))",0.15,0.85);  

     f1->GetYaxis()->SetRangeUser(0.95,1.05);

     f1->SetLineColor(8);
     f2->SetLineColor(2);
     f3->SetLineColor(1);
     f4->SetLineColor(4);

     f1->Draw();
     f2->Draw("same");
     f3->Draw("same");
     f4->Draw("same");

     auto legend = new TLegend(0.1,0.7,0.48,0.9);
     legend->AddEntry("f3","{}^{2}H/{}^{1}H","l");
     legend->AddEntry("f1","{}^{3}He/{}^{2}H","l");
     legend->AddEntry("f4","{}^{3}H/{}^{2}H","l");
     legend->AddEntry("f2","{}^{3}H/{}^{3}He","l");
     legend->Draw();

}
