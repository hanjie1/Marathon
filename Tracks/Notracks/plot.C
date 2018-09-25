void plot()
{
     ifstream filename;
     filename.open("H3_He_withcosmic.txt");

     Double_t kin[11],He_eff[11],H3_eff[11];
     Double_t He_eff_err[11],H3_eff_err[11];
     string line;
     
     int ii=0;
     Double_t tmp=0;
     while(getline(filename,line)){
           istringstream str(line);
           if(ii/11==0)str>>kin[ii]>>H3_eff[ii]>>H3_eff_err[ii];
           else str>>tmp>>He_eff[ii%11]>>He_eff_err[ii%11];
 	   ii++;
     }
   
    Double_t ratio[11],ratio_err[11];
    for(int jj=0;jj<11;jj++){
        He_eff[jj]=1-He_eff[jj];
        H3_eff[jj]=1-H3_eff[jj];
	ratio[jj]=He_eff[jj]/H3_eff[jj];
        ratio_err[jj]=sqrt(He_eff_err[jj]*He_eff_err[jj]/(H3_eff[jj]*H3_eff[jj])+(H3_eff_err[jj]*H3_eff_err[jj])/(He_eff[jj]*He_eff[jj]));
    }

 
   TGraphErrors *hHe=new TGraphErrors(11,kin,He_eff,0,He_eff_err);
   TGraphErrors *hH3=new TGraphErrors(11,kin,H3_eff,0,H3_eff_err);

   TCanvas *c1=new TCanvas();
   TMultiGraph *mg=new TMultiGraph();
   mg->Add(hHe);
   hHe->SetMarkerColor(1);
   hHe->SetMarkerStyle(8);
   mg->Add(hH3);
   hH3->SetMarkerColor(2);
   hH3->SetMarkerStyle(8);
   mg->Draw("AP");
   mg->SetTitle("track efficiency;kin;eff");

   auto leg = new TLegend(0.1,0.65,0.25,0.80);
   leg->AddEntry(hHe,"He3","P");
   leg->AddEntry(hH3,"Tritium","P");
   leg->Draw();

   TCanvas *c2=new TCanvas();
   TGraphErrors *hratio=new TGraphErrors(11,kin,ratio,0,ratio_err);
   hratio->SetMarkerColor(1);
   hratio->SetMarkerStyle(8);
   hratio->Draw("AP");
   hratio->SetTitle("He3/H3 track eff;kin;He3/H3");


}
