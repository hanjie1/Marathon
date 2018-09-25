void plot_cerped()
{
     ifstream file1;
     file1.open("../OUT/Rcer_ped.dat");

     Double_t nrun[1400]={0.0};
     Double_t cer1[1400]={0.0};
     Double_t cer2[1400]={0.0};
     Double_t cer3[1400]={0.0};
     Double_t cer4[1400]={0.0};
     Double_t cer5[1400]={0.0};
     Double_t cer6[1400]={0.0};
     Double_t cer7[1400]={0.0};
     Double_t cer8[1400]={0.0};
     Double_t cer9[1400]={0.0};
     Double_t cer10[1400]={0.0};
     Double_t flag[1400]={0.0};

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0,total=0;
     while(tmp.ReadLine(file1)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from,"  ");
          nrun[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          cer1[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          cer2[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          cer3[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          cer4[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          cer5[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          cer6[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          cer7[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          cer8[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          cer9[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          cer10[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          flag[nn-1]=atof(content.Data());
          nn++;
          total++;
          from=0;
     }
     file1.close();

     TGraph *hcer0=new TGraph();
     TGraph *hcer1=new TGraph();
     TGraph *hcer2=new TGraph();
     TGraph *hcer3=new TGraph();
     TGraph *hcer4=new TGraph();
     TGraph *hcer5=new TGraph();
     TGraph *hcer6=new TGraph();
     TGraph *hcer7=new TGraph();
     TGraph *hcer8=new TGraph();
     TGraph *hcer9=new TGraph();

     int mm=0;
     for(int ii=0;ii<total;ii++){
         if(flag[ii]==1||nrun[ii]==0)continue;
         hcer0->SetPoint(mm,nrun[ii],cer1[ii]);
         hcer1->SetPoint(mm,nrun[ii],cer2[ii]);
         hcer2->SetPoint(mm,nrun[ii],cer3[ii]);
         hcer3->SetPoint(mm,nrun[ii],cer4[ii]);
         hcer4->SetPoint(mm,nrun[ii],cer5[ii]);
         hcer5->SetPoint(mm,nrun[ii],cer6[ii]);
         hcer6->SetPoint(mm,nrun[ii],cer7[ii]);
         hcer7->SetPoint(mm,nrun[ii],cer8[ii]);
         hcer8->SetPoint(mm,nrun[ii],cer9[ii]);
         hcer9->SetPoint(mm,nrun[ii],cer10[ii]);
         mm++;
     }
     TCanvas *c1=new TCanvas("c1");
     hcer0->SetMarkerStyle(8);
     hcer0->SetMarkerColor(4);
     hcer0->Draw("AP");
     hcer0->SetTitle("RHRS cer1;run number;pedestal");

     TCanvas *c2=new TCanvas("c2");
     hcer1->SetMarkerStyle(8);
     hcer1->SetMarkerColor(4);
     hcer1->Draw("AP");
     hcer1->SetTitle("RHRS cer2;run number;pedestal");

     TCanvas *c3=new TCanvas("c3");
     hcer2->SetMarkerStyle(8);
     hcer2->SetMarkerColor(4);
     hcer2->Draw("AP");
     hcer2->SetTitle("RHRS cer3;run number;pedestal");

     TCanvas *c4=new TCanvas("c4");
     hcer3->SetMarkerStyle(8);
     hcer3->SetMarkerColor(4);
     hcer3->Draw("AP");
     hcer3->SetTitle("RHRS cer4;run number;pedestal");

     TCanvas *c5=new TCanvas("c5");
     hcer4->SetMarkerStyle(8);
     hcer4->SetMarkerColor(4);
     hcer4->Draw("AP");
     hcer4->SetTitle("RHRS cer5;run number;pedestal");

     TCanvas *c6=new TCanvas("c6");
     hcer5->SetMarkerStyle(8);
     hcer5->SetMarkerColor(4);
     hcer5->Draw("AP");
     hcer5->SetTitle("RHRS cer6;run number;pedestal");

     TCanvas *c7=new TCanvas("c7");
     hcer6->SetMarkerStyle(8);
     hcer6->SetMarkerColor(4);
     hcer6->Draw("AP");
     hcer6->SetTitle("RHRS cer7;run number;pedestal");

     TCanvas *c8=new TCanvas("c8");
     hcer7->SetMarkerStyle(8);
     hcer7->SetMarkerColor(4);
     hcer7->Draw("AP");
     hcer7->SetTitle("RHRS cer8;run number;pedestal");

     TCanvas *c9=new TCanvas("c9");
     hcer8->SetMarkerStyle(8);
     hcer8->SetMarkerColor(4);
     hcer8->Draw("AP");
     hcer8->SetTitle("RHRS cer9;run number;pedestal");

     TCanvas *c10=new TCanvas("c10");
     hcer9->SetMarkerStyle(8);
     hcer9->SetMarkerColor(4);
     hcer9->Draw("AP");
     hcer9->SetTitle("RHRS cer10;run number;pedestal");

     c1->Print("Rcer_ped.pdf[");
     c1->Print("Rcer_ped.pdf");
     c2->Print("Rcer_ped.pdf");
     c3->Print("Rcer_ped.pdf");
     c4->Print("Rcer_ped.pdf");
     c5->Print("Rcer_ped.pdf");
     c6->Print("Rcer_ped.pdf");
     c7->Print("Rcer_ped.pdf");
     c8->Print("Rcer_ped.pdf");
     c9->Print("Rcer_ped.pdf");
     c10->Print("Rcer_ped.pdf");
     c10->Print("Rcer_ped.pdf]");






}


