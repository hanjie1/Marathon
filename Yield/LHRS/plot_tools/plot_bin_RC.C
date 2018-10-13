void plot_bin_RC()
{
     ifstream H1_file,D2_file;
     ifstream H1_fileR,D2_fileR;

     Double_t H1_yield[17]={0.0},D2_yield[17]={0.0};
     Double_t H1_Q2[17]={0.0},D2_Q2[17]={0.0};
     Double_t H1_nQ2[17]={0.0},D2_nQ2[17]={0.0};
     Double_t H1_xbj[17]={0.0},D2_xbj[17]={0.0};
     Double_t H1_nxbj[17]={0.0},D2_nxbj[17]={0.0};
     Double_t H1_RC[17]={0.0},D2_RC[17]={0.0};
     Double_t H1_yield_err[17]={0.0},D2_yield_err[17]={0.0};
     Double_t H1_RC_err[17]={0.0},D2_RC_err[17]={0.0};
     Double_t xbj[17]={0.0};
     for(int ii=0;ii<17;ii++)
         xbj[ii]=0.17+ii*0.02;
     
     for(int ii=0;ii<4;ii++)
      {
         H1_file.open(Form("../RawYield/H1_kin%d.txt",ii+1));
         Ssiz_t from=0;
         TString content,tmp;
         int nn=0,nnx=0; 
         while(tmp.ReadLine(H1_file)){
	     if(nn==0){nn++;continue;}	
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             H1_yield[nn-1]+=atof(content.Data());
             tmp.Tokenize(content,from," ");
             H1_yield_err[nn-1]+=atof(content.Data())*atof(content.Data());
             nn++;
             from=0;
         }
	H1_file.close();

         H1_fileR.open(Form("../RC_Yield/H1_kin%d.txt",ii+1));
         from=0;
         nn=0,nnx=0;
         while(tmp.ReadLine(H1_fileR)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from," ");
             H1_xbj[nn-1]+=atof(content.Data());
             if(atof(content.Data())!=0)H1_nxbj[nn-1]+=1.0;
             tmp.Tokenize(content,from," ");
             H1_RC[nn-1]+=atof(content.Data());
             tmp.Tokenize(content,from," ");
             H1_RC_err[nn-1]+=atof(content.Data())*atof(content.Data());
             tmp.Tokenize(content,from," ");
             H1_Q2[nn-1]+=atof(content.Data());
             if(atof(content.Data())!=0)H1_nQ2[nn-1]+=1.0;
             nn++;
             from=0;
         }
        H1_fileR.close();


         D2_file.open(Form("../RawYield/D2_kin%d.txt",ii+1));
         from=0;
         nn=0,nnx=0;
         while(tmp.ReadLine(D2_file)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             D2_yield[nn-1]+=atof(content.Data());
             tmp.Tokenize(content,from," ");
             D2_yield_err[nn-1]+=atof(content.Data())*atof(content.Data());
             nn++;
             from=0;
         }
        D2_file.close();

         D2_fileR.open(Form("../RC_Yield/D2_kin%d.txt",ii+1));
         from=0;
         nn=0,nnx=0;
         while(tmp.ReadLine(D2_fileR)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from," ");
             D2_xbj[nn-1]+=atof(content.Data());
             if(atof(content.Data())!=0)D2_nxbj[nn-1]+=1.0;
             tmp.Tokenize(content,from," ");
             D2_RC[nn-1]+=atof(content.Data());
             tmp.Tokenize(content,from," ");
             D2_RC_err[nn-1]+=atof(content.Data())*atof(content.Data());
             tmp.Tokenize(content,from," ");
             D2_Q2[nn-1]+=atof(content.Data());
             if(atof(content.Data())!=0)D2_nQ2[nn-1]+=1.0;
             nn++;
             from=0;
         }
        D2_fileR.close();


      }

     ifstream H1_rcfile,D2_rcfile;
     Double_t H1_factor[17]={1.0},D2_factor[17]={1.0};

     H1_rcfile.open("../RadCor/H1_bin_xs.out");
         Ssiz_t from=0;
         TString content,tmp;
         int nn=0;
         Double_t sborn=0.0,srad=0.0;
         while(tmp.ReadLine(H1_rcfile)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             sborn=atof(content.Data());
             tmp.Tokenize(content,from," ");
             srad=atof(content.Data());
	     H1_factor[nn-1]=sborn/srad;
             //cout<<sborn<<"  "<<srad<<"  "<<H1_factor[nn-1]<<endl;
             nn++;
             from=0;
         }
        H1_rcfile.close();

     D2_rcfile.open("../RadCor/D2_bin_xs.out");
         from=0;
         content,tmp;
         nn=0;
         sborn=0.0,srad=0.0;
         while(tmp.ReadLine(D2_rcfile)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             sborn=atof(content.Data());
             tmp.Tokenize(content,from," ");
             srad=atof(content.Data());
             D2_factor[nn-1]=sborn/srad;
             nn++;
             from=0;
         }
        D2_rcfile.close();






     Double_t ratio[17]={0.0},RC_ratio[17]={0.0},Bin_ratio[17]={0.0};
     Double_t ratio_err[17]={0.0},RC_ratio_err[17]={0.0},Bin_ratio_err[17]={0.0};
     Double_t D2_RC_bin[17]={0.0},H1_RC_bin[17]={0.0};
     Double_t D2_RC_bin_err[17]={0.0},H1_RC_bin_err[17]={0.0};
     for(int ii=0;ii<17;ii++){
	 if(H1_yield[ii]!=0){
            ratio[ii]=D2_yield[ii]/H1_yield[ii];
            D2_yield_err[ii]=sqrt(D2_yield_err[ii])/D2_nQ2[ii];
            H1_yield_err[ii]=sqrt(H1_yield_err[ii])/H1_nQ2[ii];
            D2_RC_err[ii]=sqrt(D2_RC_err[ii])/D2_nQ2[ii];
            H1_RC_err[ii]=sqrt(H1_RC_err[ii])/H1_nQ2[ii];
	    ratio_err[ii]=ratio[ii]*sqrt((D2_yield_err[ii]/D2_yield[ii])*(D2_yield_err[ii]/D2_yield[ii])+(H1_yield_err[ii]/H1_yield[ii])*(H1_yield_err[ii]/H1_yield[ii]));

            Bin_ratio[ii]=(D2_yield[ii]*D2_factor[ii])/(H1_yield[ii]*H1_factor[ii]);
	    Bin_ratio_err[ii]=ratio_err[ii]*D2_factor[ii]/H1_factor[ii];

            D2_RC_bin[ii]=D2_yield[ii]*D2_factor[ii];
            D2_RC_bin_err[ii]=D2_yield_err[ii]*D2_factor[ii];
            H1_RC_bin[ii]=H1_yield[ii]*H1_factor[ii];
            H1_RC_bin_err[ii]=H1_yield_err[ii]*H1_factor[ii];
         }
         if(H1_RC[ii]!=0){
            RC_ratio[ii]=D2_RC[ii]/H1_RC[ii];
            RC_ratio_err[ii]=ratio[ii]*sqrt((D2_RC_err[ii]/D2_RC[ii])*(D2_RC_err[ii]/D2_RC[ii])+(H1_RC_err[ii]/H1_RC[ii])*(H1_RC_err[ii]/H1_RC[ii]));
         }
         if(H1_nQ2[ii]!=0)H1_Q2[ii]=H1_Q2[ii]/H1_nQ2[ii];
         if(D2_nQ2[ii]!=0)D2_Q2[ii]=D2_Q2[ii]/D2_nQ2[ii];
         if(H1_nxbj[ii]!=0)H1_xbj[ii]=H1_xbj[ii]/H1_nxbj[ii];
         if(D2_nxbj[ii]!=0)D2_xbj[ii]=D2_xbj[ii]/D2_nxbj[ii];
         //cout<<fixed<<setprecision(3);
         //cout<<D2_xbj[ii]<<" &  "<<ratio[ii]<<" & "<<RC_ratio[ii]<<" & "<<Bin_ratio[ii]<<endl;
         //cout<<D2_xbj[ii]<<",";
      }
     cout<<endl;

     Double_t xbj_err[17]={0.0};
     TCanvas *c1=new TCanvas();
     TGraphErrors *hH1=new TGraphErrors(17,D2_xbj,H1_yield,xbj_err,H1_yield_err);
     TGraphErrors *hH1_RC=new TGraphErrors(17,D2_xbj,H1_RC,xbj_err,H1_RC_err);
     TGraphErrors *hH1_RC_bin=new TGraphErrors(17,D2_xbj,H1_RC_bin,xbj_err,H1_RC_bin_err);
     TMultiGraph *mg_H1=new TMultiGraph();
     hH1->SetMarkerStyle(8);
     hH1->SetMarkerColor(1);
     hH1_RC->SetMarkerStyle(8);
     hH1_RC->SetMarkerColor(2);
     hH1_RC_bin->SetMarkerStyle(8);
     hH1_RC_bin->SetMarkerColor(4);
     mg_H1->Add(hH1);
     mg_H1->Add(hH1_RC);
     mg_H1->Add(hH1_RC_bin);
     mg_H1->Draw("AP");
     mg_H1->SetTitle("Hydrogen Yield;xbj;nb");

     auto leg1 = new TLegend(0.7,0.7,0.85,0.85);
     leg1->AddEntry(hH1,"before RC","P");
     leg1->AddEntry(hH1_RC,"after RC (event by event)","P");
     leg1->AddEntry(hH1_RC_bin,"after RC (bin by bin)","P");
     leg1->Draw(); 

     TCanvas *c3=new TCanvas();
     TGraphErrors *hD2=new TGraphErrors(17,D2_xbj,D2_yield,xbj_err,D2_yield_err);
     TGraphErrors *hD2_RC=new TGraphErrors(17,D2_xbj,D2_RC,xbj_err,D2_RC_err);
     TGraphErrors *hD2_RC_bin=new TGraphErrors(17,D2_xbj,D2_RC_bin,xbj_err,D2_RC_bin_err);
     TMultiGraph *mg_D2=new TMultiGraph();
     hD2->SetMarkerStyle(8);
     hD2->SetMarkerColor(1);
     hD2_RC->SetMarkerStyle(8);
     hD2_RC->SetMarkerColor(2);
     hD2_RC_bin->SetMarkerStyle(8);
     hD2_RC_bin->SetMarkerColor(4);
     mg_D2->Add(hD2);
     mg_D2->Add(hD2_RC);
     mg_D2->Add(hD2_RC_bin);
     mg_D2->Draw("AP");
     mg_D2->SetTitle("Deuterium Yield;xbj;nb");

     auto leg2 = new TLegend(0.7,0.7,0.85,0.85);
     leg2->AddEntry(hD2,"before RC","P");
     leg2->AddEntry(hD2_RC,"after RC (event by event)","P");
     leg2->AddEntry(hD2_RC_bin,"after RC (bin by bin)","P");
     leg2->Draw(); 

     Double_t SLAC_x[13]={0.1722,0.1898,0.2044,0.2206,0.2407,0.2585,0.276,0.293,0.312,0.3305,0.357,0.3852,0.4151};
     Double_t SLAC_Q2[13]={2.5482,2.7419,2.8948,3.0579,3.2515,3.4146,3.5675,3.7102,3.8631,4.0058,4.1995,4.3931,4.5868};
     Double_t SLAC_D2[13]={11700,11700,11300,11000,10300,9850,9490,8970,8590,8210,7390,6900,6480};
     Double_t SLAC_H1[13]={6520,6670,6670,6420,6330,5860,5550,5340,5080,5030,4770,4510,4020};

     Double_t SLAC_ratio[13];
     for(int ii=0;ii<13;ii++){
	 SLAC_ratio[ii]=SLAC_D2[ii]/SLAC_H1[ii];
     //    cout<<fixed<<setprecision(3);
      //   cout<<SLAC_x[ii]<<" & "<<SLAC_Q2[ii]<<" & "<<SLAC_ratio[ii]<<endl;
     }


     TCanvas *c2=new TCanvas();
     TGraphErrors *hratio=new TGraphErrors(17,D2_xbj,ratio,xbj_err,ratio_err);
     TGraphErrors *hratio_RC=new TGraphErrors(17,D2_xbj,RC_ratio,xbj_err,RC_ratio_err);
     TGraphErrors *hratio_RC_bin=new TGraphErrors(17,D2_xbj,Bin_ratio,xbj_err,Bin_ratio_err);
     TGraph *hSLAC=new TGraph(13,SLAC_x,SLAC_ratio);
     TMultiGraph *mg_ratio=new TMultiGraph();
     hratio->SetMarkerStyle(8);
     hratio->SetMarkerColor(1);
     hratio_RC->SetMarkerStyle(8);
     hratio_RC->SetMarkerColor(2);
     hSLAC->SetMarkerStyle(8);
     hSLAC->SetMarkerColor(4);
     hratio_RC_bin->SetMarkerStyle(8);
     hratio_RC_bin->SetMarkerColor(6);
     mg_ratio->Add(hratio);
     mg_ratio->Add(hratio_RC);
     mg_ratio->Add(hratio_RC_bin);
     mg_ratio->Add(hSLAC);
     mg_ratio->Draw("AP");
     mg_ratio->SetTitle("Deuterium/proton Data Yield;xbj");

     auto leg3 = new TLegend(0.7,0.7,0.85,0.85);
     leg3->AddEntry(hratio,"before RC","P");
     leg3->AddEntry(hratio_RC,"after RC (event by event)","P");
     leg3->AddEntry(hratio_RC_bin,"after RC (bin by bin)","P");
     leg3->AddEntry(hSLAC,"SLAC","P");
     leg3->Draw(); 

     TCanvas *c4=new TCanvas();
     TGraphErrors *hratio_data=new TGraphErrors(17,D2_xbj,ratio,xbj_err,ratio_err);
     hratio_data->SetMarkerStyle(8);
     hratio_data->SetMarkerColor(1);
     hratio_data->Draw("AP");
     hratio_data->SetTitle("Deuterium/proton Data Yield;xbj");


}
