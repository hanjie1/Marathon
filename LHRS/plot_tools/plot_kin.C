void plot_kin()
{
     ifstream H1_file,D2_file;
     ifstream H1_rcfile,D2_rcfile;
     ifstream H1_Rfile,D2_Rfile;

     Double_t H1_data[4][17],D2_data[4][17],H1_data_err[4][17],D2_data_err[4][17];
     Double_t H1_RC[4][17],D2_RC[4][17],H1_RC_err[4][17],D2_RC_err[4][17]; 
     Double_t H1_bin[4][17],D2_bin[4][17],H1_bin_err[4][17],D2_bin_err[4][17]; 
     Double_t H1_xbj[4][17],D2_xbj[4][17],H1_Q2[4][17],D2_Q2[4][17];

     Double_t H1_factor[4][17]={{0.0}},D2_factor[4][17]={{0.0}};

     for(int ii=0;ii<4;ii++){
	 for(int jj=0;jj<17;jj++){
	     H1_data[ii][jj]=0;D2_data[ii][jj]=0;H1_data_err[ii][jj]=0;D2_data_err[ii][jj]=0;
	     H1_RC[ii][jj]=0;  D2_RC[ii][jj]=0;  H1_RC_err[ii][jj]=0;  D2_RC_err[ii][jj]=0;
	     H1_bin[ii][jj]=0; D2_bin[ii][jj]=0; H1_bin_err[ii][jj]=0; D2_bin_err[ii][jj]=0;
	     H1_xbj[ii][jj]=0;    H1_Q2[ii][jj]=0;
	     D2_xbj[ii][jj]=0;    D2_Q2[ii][jj]=0;
	 }
     }

     for(int ii=0;ii<4;ii++){
         H1_file.open(Form("../RawYield/vz007/H1_kin%d.txt",ii+1));
         Ssiz_t from=0;
         TString content,tmp;
         int nn=0;
         while(tmp.ReadLine(H1_file)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from," ");
	     H1_xbj[ii][nn-1]=atof(content.Data());
             tmp.Tokenize(content,from," ");
             H1_data[ii][nn-1]=atof(content.Data());
             tmp.Tokenize(content,from," ");
             H1_data_err[ii][nn-1]=atof(content.Data());
             tmp.Tokenize(content,from," ");
             H1_Q2[ii][nn-1]=atof(content.Data());
             nn++;
             from=0;
         }
         H1_file.close();
 
         D2_file.open(Form("../RawYield/vz007/D2_kin%d.txt",ii+1));
         from=0;
         nn=0;
         while(tmp.ReadLine(D2_file)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from," ");
             D2_xbj[ii][nn-1]=atof(content.Data());
             tmp.Tokenize(content,from," ");
             D2_data[ii][nn-1]=atof(content.Data());
             tmp.Tokenize(content,from," ");
             D2_data_err[ii][nn-1]=atof(content.Data());
             tmp.Tokenize(content,from," ");
             D2_Q2[ii][nn-1]=atof(content.Data());
             nn++;
             from=0;
          }
         D2_file.close();

         H1_rcfile.open(Form("../RC_Yield/vz007/H1_kin%d.txt",ii+1));
         from=0;
         nn=0;
         while(tmp.ReadLine(H1_rcfile)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             H1_RC[ii][nn-1]=atof(content.Data());
             tmp.Tokenize(content,from," ");
             H1_RC_err[ii][nn-1]=atof(content.Data());
             tmp.Tokenize(content,from," ");
             nn++;
             from=0;
         }
         H1_rcfile.close();

         D2_rcfile.open(Form("../RC_Yield/vz007/D2_kin%d.txt",ii+1));
         from=0;
         nn=0;
         while(tmp.ReadLine(D2_rcfile)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             D2_RC[ii][nn-1]=atof(content.Data());
             tmp.Tokenize(content,from," ");
             D2_RC_err[ii][nn-1]=atof(content.Data());
             tmp.Tokenize(content,from," ");
             nn++;
             from=0;
         }
         D2_rcfile.close();

         H1_Rfile.open(Form("../RadCor/bin/vz007/H1_kin%d_xs.out",ii+1));
         from=0;
         nn=0;
         Double_t H1_sigb=0.0,H1_sigr=0.0;
         Double_t D2_sigb=0.0,D2_sigr=0.0;
         while(tmp.ReadLine(H1_Rfile)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             H1_sigb=atof(content.Data());
             tmp.Tokenize(content,from," ");
             H1_sigr=atof(content.Data());
             if(H1_sigr!=0)H1_factor[ii][nn-1]=H1_sigb/H1_sigr;
             nn++;
             from=0;
         }
         H1_Rfile.close();

         D2_Rfile.open(Form("../RadCor/bin/vz007/D2_kin%d_xs.out",ii+1));
         from=0;
         nn=0;
         while(tmp.ReadLine(D2_Rfile)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             D2_sigb=atof(content.Data());
             tmp.Tokenize(content,from," ");
             D2_sigr=atof(content.Data());
             if(D2_sigr!=0)D2_factor[ii][nn-1]=D2_sigb/D2_sigr;
             nn++;
             from=0;
         }
         D2_Rfile.close();

     }

    Double_t Raw_ratio[4][17]={{0.0}},Raw_ratio_err[4][17]={{0.0}};
    Double_t RC_ratio[4][17]={{0.0}},RC_ratio_err[4][17]={{0.0}};
    Double_t bin_ratio[4][17]={{0.0}},bin_ratio_err[4][17]={{0.0}};

    Double_t H1_posA=0.079031,H1_posB=8.587641;
    Double_t D2_posA=0.087557,D2_posB=8.446762;

    Double_t H1_end[4]={0.020565,0.01805,0.01798, 0.01635};
    Double_t D2_end[4]={0.008531,0.007625,0.00721, 0.00699};

    Double_t Raw_H1corr[4][17]={{0.0}};
    Double_t Raw_D2corr[4][17]={{0.0}};
    Double_t Raw_Rcorr1[4][17]={{0.0}},Raw_Rcorr_err1[4][17]={{0.0}};
    Double_t Raw_Rcorr2[4][17]={{0.0}},Raw_Rcorr_err2[4][17]={{0.0}};

    Double_t du[4][17]={{0.0}},du_err[4][17]={{0.0}};

    ofstream pfile;
    pfile.open("RadCor.dat");
    pfile<<"kin xbj    Q2    H1_factor     D2_factor"<<endl;    
    for(int ii=0;ii<4;ii++){
        for(int jj=0;jj<17;jj++){
	    if(H1_xbj[ii][jj]==0||D2_xbj[ii][jj]==0)continue;
            Raw_ratio[ii][jj]=D2_data[ii][jj]/H1_data[ii][jj];
            Raw_ratio_err[ii][jj]=Raw_ratio[ii][jj]*sqrt(H1_data_err[ii][jj]*H1_data_err[ii][jj]/(H1_data[ii][jj]*H1_data[ii][jj])+
                                                 D2_data_err[ii][jj]*D2_data_err[ii][jj]/(D2_data[ii][jj]*D2_data[ii][jj]));
          
            Double_t H1_corr=1-H1_posA*exp(-1.0*H1_posB*H1_xbj[ii][jj]);
            Double_t D2_corr=1-D2_posA*exp(-1.0*D2_posB*D2_xbj[ii][jj]);
            Raw_H1corr[ii][jj]=H1_data[ii][jj]*H1_corr;
            Raw_D2corr[ii][jj]=D2_data[ii][jj]*D2_corr;
            Raw_Rcorr1[ii][jj]=Raw_D2corr[ii][jj]/Raw_H1corr[ii][jj];
            Raw_Rcorr_err1[ii][jj]=Raw_ratio_err[ii][jj]*D2_corr/H1_corr;
 
            Raw_Rcorr2[ii][jj]=Raw_Rcorr1[ii][jj]*(1.0-D2_end[ii])/(1.0-H1_end[ii]);
            Raw_Rcorr_err2[ii][jj]=Raw_Rcorr_err1[ii][jj]*(1.0-D2_end[ii])/(1.0-H1_end[ii]);

            RC_ratio[ii][jj]=D2_RC[ii][jj]/H1_RC[ii][jj]*D2_corr*(1.0-D2_end[ii])/((1.0-H1_end[ii])*H1_corr);
            RC_ratio_err[ii][jj]=RC_ratio[ii][jj]*sqrt(H1_RC_err[ii][jj]*H1_RC_err[ii][jj]/(H1_RC[ii][jj]*H1_RC[ii][jj])+
                                                       D2_RC_err[ii][jj]*D2_RC_err[ii][jj]/(D2_RC[ii][jj]*D2_RC[ii][jj]))*D2_corr*(1.0-D2_end[ii])/((1.0-H1_end[ii])*H1_corr);
            bin_ratio[ii][jj]=(D2_data[ii][jj]*D2_factor[ii][jj])/(H1_data[ii][jj]*H1_factor[ii][jj]);
            bin_ratio_err[ii][jj]=Raw_ratio_err[ii][jj]*D2_factor[ii][jj]/H1_factor[ii][jj];

            Double_t RR=1.0095-0.0109*H1_xbj[ii][jj]-0.0821*H1_xbj[ii][jj]*H1_xbj[ii][jj]; 
            du[ii][jj]=RC_ratio[ii][jj]/RR-1.0;
            du_err[ii][jj]=RC_ratio_err[ii][jj]/RR;
            pfile<<fixed;
            pfile<<ii+1<<"    "<<setprecision(3)<<H1_xbj[ii][jj]<<"  "<<H1_Q2[ii][jj]<<"  "<<setprecision(3)<<H1_factor[ii][jj]<<"  "<<D2_factor[ii][jj]<< endl;    
        }
        
    }
    pfile.close();



    TGraphErrors *gRaw_kin[4];
    TGraphErrors *gRaw_corr1[4];
    TGraphErrors *gRaw_corr2[4];
    TGraphErrors *gRC_kin[4];
    TGraphErrors *gBin_kin[4];
    TGraphErrors *gDU[4];
         
    for(int ii=0;ii<4;ii++){
	gRaw_kin[ii]=new TGraphErrors();
	gRaw_corr1[ii]=new TGraphErrors();
	gRaw_corr2[ii]=new TGraphErrors();
	gRC_kin[ii]=new TGraphErrors();
	gBin_kin[ii]=new TGraphErrors();
	gDU[ii]=new TGraphErrors();
        int nn=0;
        for(int jj=0;jj<17;jj++){
           if(H1_xbj[ii][jj]==0||D2_xbj[ii][jj]==0)continue;
           if(ii==0 && 0<jj&&jj<5){
              gRaw_kin[ii]->SetPoint(nn,H1_xbj[ii][jj],Raw_ratio[ii][jj]);
              gRaw_kin[ii]->SetPointError(nn,0,Raw_ratio_err[ii][jj]);
              gRaw_corr1[ii]->SetPoint(nn,H1_xbj[ii][jj],Raw_Rcorr1[ii][jj]);
              gRaw_corr1[ii]->SetPointError(nn,0,Raw_Rcorr_err1[ii][jj]);
              gRaw_corr2[ii]->SetPoint(nn,H1_xbj[ii][jj],Raw_Rcorr2[ii][jj]);
              gRaw_corr2[ii]->SetPointError(nn,0,Raw_Rcorr_err2[ii][jj]);
              gRC_kin[ii]->SetPoint(nn,H1_xbj[ii][jj],RC_ratio[ii][jj]);
              gRC_kin[ii]->SetPointError(nn,0,RC_ratio_err[ii][jj]);
              gBin_kin[ii]->SetPoint(nn,H1_xbj[ii][jj],bin_ratio[ii][jj]);
              gBin_kin[ii]->SetPointError(nn,0,bin_ratio_err[ii][jj]);
              gDU[ii]->SetPoint(nn,H1_xbj[ii][jj],du[ii][jj]);
              gDU[ii]->SetPointError(nn,0,du_err[ii][jj]);
     	      nn++;
	   }
           if(ii==1 && 1<jj&&jj<8){
              gRaw_kin[ii]->SetPoint(nn,H1_xbj[ii][jj],Raw_ratio[ii][jj]);
              gRaw_kin[ii]->SetPointError(nn,0,Raw_ratio_err[ii][jj]);
              gRaw_corr1[ii]->SetPoint(nn,H1_xbj[ii][jj],Raw_Rcorr1[ii][jj]);
              gRaw_corr1[ii]->SetPointError(nn,0,Raw_Rcorr_err1[ii][jj]);
              gRaw_corr2[ii]->SetPoint(nn,H1_xbj[ii][jj],Raw_Rcorr2[ii][jj]);
              gRaw_corr2[ii]->SetPointError(nn,0,Raw_Rcorr_err2[ii][jj]);
              gRC_kin[ii]->SetPoint(nn,H1_xbj[ii][jj],RC_ratio[ii][jj]);
              gRC_kin[ii]->SetPointError(nn,0,RC_ratio_err[ii][jj]);
              gBin_kin[ii]->SetPoint(nn,H1_xbj[ii][jj],bin_ratio[ii][jj]);
              gBin_kin[ii]->SetPointError(nn,0,bin_ratio_err[ii][jj]);
              gDU[ii]->SetPoint(nn,H1_xbj[ii][jj],du[ii][jj]);
              gDU[ii]->SetPointError(nn,0,du_err[ii][jj]);
              nn++;
           }
           if(ii==2 && 3<jj&&jj<10){
              gRaw_kin[ii]->SetPoint(nn,H1_xbj[ii][jj],Raw_ratio[ii][jj]);
              gRaw_kin[ii]->SetPointError(nn,0,Raw_ratio_err[ii][jj]);
              gRaw_corr1[ii]->SetPoint(nn,H1_xbj[ii][jj],Raw_Rcorr1[ii][jj]);
              gRaw_corr1[ii]->SetPointError(nn,0,Raw_Rcorr_err1[ii][jj]);
              gRaw_corr2[ii]->SetPoint(nn,H1_xbj[ii][jj],Raw_Rcorr2[ii][jj]);
              gRaw_corr2[ii]->SetPointError(nn,0,Raw_Rcorr_err2[ii][jj]);
              gRC_kin[ii]->SetPoint(nn,H1_xbj[ii][jj],RC_ratio[ii][jj]);
              gRC_kin[ii]->SetPointError(nn,0,RC_ratio_err[ii][jj]);
              gBin_kin[ii]->SetPoint(nn,H1_xbj[ii][jj],bin_ratio[ii][jj]);
              gBin_kin[ii]->SetPointError(nn,0,bin_ratio_err[ii][jj]);
              gDU[ii]->SetPoint(nn,H1_xbj[ii][jj],du[ii][jj]);
              gDU[ii]->SetPointError(nn,0,du_err[ii][jj]);
              nn++;
           }
           if(ii==3 && 5<jj&&jj<12){
              gRaw_kin[ii]->SetPoint(nn,H1_xbj[ii][jj],Raw_ratio[ii][jj]);
              gRaw_kin[ii]->SetPointError(nn,0,Raw_ratio_err[ii][jj]);
              gRaw_corr1[ii]->SetPoint(nn,H1_xbj[ii][jj],Raw_Rcorr1[ii][jj]);
              gRaw_corr1[ii]->SetPointError(nn,0,Raw_Rcorr_err1[ii][jj]);
              gRaw_corr2[ii]->SetPoint(nn,H1_xbj[ii][jj],Raw_Rcorr2[ii][jj]);
              gRaw_corr2[ii]->SetPointError(nn,0,Raw_Rcorr_err2[ii][jj]);
              gRC_kin[ii]->SetPoint(nn,H1_xbj[ii][jj],RC_ratio[ii][jj]);
              gRC_kin[ii]->SetPointError(nn,0,RC_ratio_err[ii][jj]);
              gBin_kin[ii]->SetPoint(nn,H1_xbj[ii][jj],bin_ratio[ii][jj]);
              gBin_kin[ii]->SetPointError(nn,0,bin_ratio_err[ii][jj]);
              gDU[ii]->SetPoint(nn,H1_xbj[ii][jj],du[ii][jj]);
              gDU[ii]->SetPointError(nn,0,du_err[ii][jj]);
              nn++;
           }

        }
    }

     Double_t M_x[11]={0.0},M_ratio[11]={0.0},M_ratio_err[11]={0.0};
     M_x[0]=H1_xbj[0][1];
     M_ratio[0]=RC_ratio[0][1];
     M_ratio_err[0]=RC_ratio_err[0][1];
     
     M_x[1]=(H1_xbj[0][2]+H1_xbj[1][2])/2.0;
     M_ratio[1]=(RC_ratio[0][2]+RC_ratio[1][2])/2.0;
     M_ratio_err[1]=sqrt(RC_ratio_err[0][2]*RC_ratio_err[0][2]+RC_ratio_err[1][2]*RC_ratio_err[1][2])/2.0;

     M_x[2]=(H1_xbj[0][3]+H1_xbj[1][3])/2.0;
     M_ratio[2]=(RC_ratio[0][3]+RC_ratio[1][3])/2.0;
     M_ratio_err[2]=sqrt(RC_ratio_err[0][3]*RC_ratio_err[0][3]+RC_ratio_err[1][3]*RC_ratio_err[1][3])/2.0;

     M_x[3]=(H1_xbj[0][4]+H1_xbj[1][4]+H1_xbj[2][4])/3.0;
     M_ratio[3]=(RC_ratio[0][4]+RC_ratio[1][4]+RC_ratio[2][4])/3.0;
     M_ratio_err[3]=sqrt(RC_ratio_err[0][4]*RC_ratio_err[0][4]+RC_ratio_err[1][4]*RC_ratio_err[1][4]+RC_ratio_err[2][4]*RC_ratio_err[2][4])/3.0;

     M_x[4]=(H1_xbj[1][5]+H1_xbj[2][5])/2.0;
     M_ratio[4]=(RC_ratio[1][5]+RC_ratio[2][5])/2.0;
     M_ratio_err[4]=sqrt(RC_ratio_err[1][5]*RC_ratio_err[1][5]+RC_ratio_err[2][5]*RC_ratio_err[2][5])/2.0;

     M_x[5]=(H1_xbj[1][6]+H1_xbj[2][6]+H1_xbj[3][6])/3.0;
     M_ratio[5]=(RC_ratio[1][6]+RC_ratio[2][6]+RC_ratio[3][6])/3.0;
     M_ratio_err[5]=sqrt(RC_ratio_err[1][6]*RC_ratio_err[1][6]+RC_ratio_err[2][6]*RC_ratio_err[2][6]+RC_ratio_err[3][6]*RC_ratio_err[3][6])/3.0;

     M_x[6]=(H1_xbj[1][7]+H1_xbj[2][7]+H1_xbj[3][7])/3.0;
     M_ratio[6]=(RC_ratio[1][7]+RC_ratio[2][7]+RC_ratio[3][7])/3.0;
     M_ratio_err[6]=sqrt(RC_ratio_err[1][7]*RC_ratio_err[1][7]+RC_ratio_err[2][7]*RC_ratio_err[2][7]+RC_ratio_err[3][7]*RC_ratio_err[3][7])/3.0;

     M_x[7]=(H1_xbj[2][8]+H1_xbj[3][8])/2.0;
     M_ratio[7]=(RC_ratio[2][8]+RC_ratio[3][8])/2.0;
     M_ratio_err[7]=sqrt(RC_ratio_err[2][8]*RC_ratio_err[2][8]+RC_ratio_err[3][8]*RC_ratio_err[3][8])/2.0;

     M_x[8]=(H1_xbj[2][9]+H1_xbj[3][9])/2.0;
     M_ratio[8]=(RC_ratio[2][9]+RC_ratio[3][9])/2.0;
     M_ratio_err[8]=sqrt(RC_ratio_err[2][9]*RC_ratio_err[2][9]+RC_ratio_err[3][9]*RC_ratio_err[3][9])/2.0;

     M_x[9]=H1_xbj[3][10];
     M_ratio[9]=RC_ratio[3][10];
     M_ratio_err[9]=RC_ratio_err[3][10];

     M_x[10]=H1_xbj[3][11];
     M_ratio[10]=RC_ratio[3][11];
     M_ratio_err[10]=RC_ratio_err[3][11];





     Double_t SLAC_x[12]={0.1722,0.1898,0.2044,0.2206,0.2407,0.2585,0.276,0.293,0.312,0.3305,0.357,0.3852};
     Double_t SLAC_Q2[12]={2.5482,2.7419,2.8948,3.0579,3.2515,3.4146,3.5675,3.7102,3.8631,4.0058,4.1995,4.3931};
     Double_t SLAC_D2[12]={11700,11700,11300,11000,10300,9850,9490,8970,8590,8210,7390,6900};
     Double_t SLAC_D2_err[12]={297,230,214,212,202,200,174,153,145,142,133,127};
     Double_t SLAC_H1[12]={6520,6670,6670,6420,6330,5860,5550,5340,5080,5030,4770,4510};
     Double_t SLAC_H1_err[12]={171,162,152,163,165,144,141,99,108,112,113,101};

     Double_t Sdu_x[12]={0.175,0.205,0.235,0.265,0.295,0.325,0.355,0.385,0.305,0.335,0.365,0.395};
     Double_t Sdu[12]={0.853,0.774,0.726,0.751,0.687,0.671,0.646,0.633,0.647,0.663,0.62,0.643};
     Double_t Sdu_err[12]={0.027,0.025,0.025,0.026,0.022,0.023,0.023,0.023,0.058,0.033,0.031,0.032};

     Double_t Sdu1[12]={0.804,0.754,0.733,0.704,0.678,0.642,0.638,0.623,0.631,0.656,0.618,0.596};
     Double_t Sdu_err1[12]={0.025,0.025,0.025,0.022,0.022,0.022,0.023,0.025,0.038,0.031,0.031,0.028};

     Double_t SLAC_ratio[12],SLAC_ratio_err[12];
     for(int ii=0;ii<12;ii++){
         SLAC_ratio[ii]=SLAC_D2[ii]/SLAC_H1[ii];
         SLAC_ratio_err[ii]=SLAC_ratio[ii]*sqrt(SLAC_H1_err[ii]*SLAC_H1_err[ii]/(SLAC_H1[ii]*SLAC_H1[ii])+
                                                 SLAC_D2_err[ii]*SLAC_D2_err[ii]/(SLAC_D2[ii]*SLAC_D2[ii]));

     //    cout<<fixed<<setprecision(3);
      //   cout<<SLAC_x[ii]<<" & "<<SLAC_Q2[ii]<<" & "<<SLAC_ratio[ii]<<endl;
     }
    TGraphErrors *gSLAC=new TGraphErrors(12,SLAC_x,SLAC_ratio,0,SLAC_ratio_err);
    TGraphErrors *gSDU=new TGraphErrors(12,Sdu_x,Sdu,0,Sdu_err);
    TGraphErrors *gSDU1=new TGraphErrors(12,Sdu_x,Sdu1,0,Sdu_err1);



    TCanvas *c1=new TCanvas("c1","c1",1500,1500);
    TMultiGraph *mg1=new TMultiGraph();
    gRaw_kin[2]->SetMarkerStyle(8);
    gRaw_kin[2]->SetMarkerColor(1);
    gRaw_corr1[2]->SetMarkerStyle(8);
    gRaw_corr1[2]->SetMarkerColor(7);
    gRaw_corr2[2]->SetMarkerStyle(8);
    gRaw_corr2[2]->SetMarkerColor(4);
    gRC_kin[2]->SetMarkerStyle(8);
    gRC_kin[2]->SetMarkerColor(2);
    mg1->Add(gRaw_kin[2]);
    mg1->Add(gRaw_corr1[2]);
    mg1->Add(gRaw_corr2[2]);
    mg1->Add(gRC_kin[2]);
    mg1->Draw("AP");
    mg1->SetTitle("Kin 3 D2/H1;xbj");

    auto leg2 = new TLegend(0.7,0.7,0.85,0.85);
    leg2->AddEntry(gRaw_kin[2],"data","P");
    leg2->AddEntry(gRaw_corr1[2],"+positron correction","P");
    leg2->AddEntry(gRaw_corr2[2],"+end cup subtraction","P");
    leg2->AddEntry(gRC_kin[2],"+radiative correction","P");
    leg2->Draw();

/*
   TCanvas *c2=new TCanvas("c2","c2");
   TMultiGraph *mg=new TMultiGraph();

   for(int ii=0;ii<4;ii++){
       int mm=0;
       if(ii==0)mm=8;
       if(ii==1)mm=22;
       if(ii==2)mm=21;
       if(ii==3)mm=34;
       gRaw_kin[ii]->SetMarkerStyle(mm);
       gRaw_corr1[ii]->SetMarkerStyle(mm);
       gRaw_corr2[ii]->SetMarkerStyle(mm);
       gRC_kin[ii]->SetMarkerStyle(mm);
       gRaw_kin[ii]->SetMarkerColor(1);
       gRaw_corr1[ii]->SetMarkerColor(3);
       gRaw_corr2[ii]->SetMarkerColor(7);
       gRC_kin[ii]->SetMarkerColor(2);
//       mg->Add(gRaw_kin[ii]);
//       mg->Add(gRaw_corr1[ii]);
//       mg->Add(gRaw_corr2[ii]);
       mg->Add(gRC_kin[ii]);
   }

   gSLAC->SetMarkerStyle(8);
   gSLAC->SetMarkerColor(1);
   mg->Add(gSLAC);
   mg->Draw("AP");
   mg->GetYaxis()->SetTitle("#sigma_{d}/#sigma_{p}");
   mg->GetYaxis()->CenterTitle("true");
   mg->GetXaxis()->SetTitle("x");
   mg->GetXaxis()->CenterTitle("true");

    auto leg1 = new TLegend(0.7,0.7,0.85,0.85);
    leg1->AddEntry(gRC_kin[0],"MARATHON kin1","P");
    leg1->AddEntry(gRC_kin[1],"MARATHON kin2","P");
    leg1->AddEntry(gRC_kin[2],"MARATHON kin3","P");
    leg1->AddEntry(gRC_kin[3],"MARATHON kin4","P");
    leg1->AddEntry(gSLAC,"SLAC Bodek et al.","P");
    leg1->Draw();


   TCanvas *c3=new TCanvas("c3","c3");
    TMultiGraph *mg2=new TMultiGraph();
   Double_t SLAC_np[12]={0.785216,0.746311,0.687669,0.7081,0.62373,0.678879,0.709504,0.681022,0.694151,0.637233,0.556824,0.540533};
   Double_t SLAC_np_err[12]={0.06516,0.054566,0.0500066,0.0544469,0.0529663,0.0535172,0.0535599,0.0423485,0.0459898,0.0461614,0.0463167,0.0446567};
   TGraphErrors *gSLAC_np=new TGraphErrors(12,SLAC_x,SLAC_np,0,SLAC_np_err);

   Double_t NMC_x[3]={0.175,0.25,0.35};
   Double_t NMC_np[3]={0.774,0.700,0.588};
   Double_t NMC_np_err[3]={0.0};
   NMC_np_err[0]=sqrt(0.008*0.008+0.005*0.005);
   NMC_np_err[1]=sqrt(0.007*0.007+0.007*0.007);
   NMC_np_err[2]=sqrt(0.011*0.011+0.009*0.009);
   TGraphErrors *gNMC_np=new TGraphErrors(3,NMC_x,NMC_np,0,NMC_np_err);

   Double_t NMC_npcorr[3],NMC_npcorr_err[3];
   for(int ii=0;ii<3;ii++){
       Double_t RR=1.0095-0.0109*NMC_x[ii]-0.0821*NMC_x[ii]*NMC_x[ii];
       NMC_npcorr[ii]=(NMC_np[ii]+1.0)/RR-1;
       NMC_npcorr_err[ii]=NMC_np_err[ii]/RR;
   }
   TGraphErrors *gNMC_npcorr=new TGraphErrors(3,NMC_x,NMC_npcorr,0,NMC_npcorr_err);

   for(int ii=0;ii<4;ii++){
       int mm=0;
       if(ii==0)mm=8;
       if(ii==1)mm=22;
       if(ii==2)mm=21;
       if(ii==3)mm=34;
       gDU[ii]->SetMarkerStyle(mm);
       gDU[ii]->SetMarkerColor(2);
       mg2->Add(gDU[ii]);
   }

   gSLAC_np->SetMarkerStyle(8);
   gSLAC_np->SetMarkerColor(1);
   gNMC_np->SetMarkerStyle(8);
   gNMC_np->SetMarkerColor(4);
   gNMC_npcorr->SetMarkerStyle(8);
   gNMC_npcorr->SetMarkerColor(6);
   mg2->Add(gSLAC_np);
   mg2->Add(gNMC_np);
//   mg2->Add(gNMC_npcorr);
   mg2->Draw("AP");
   mg2->GetXaxis()->SetTitle("x");
   mg2->GetXaxis()->CenterTitle(true);
   mg2->GetYaxis()->SetTitle("F2n/F2p");
   mg2->GetYaxis()->CenterTitle(true);

    auto leg2 = new TLegend(0.7,0.7,0.85,0.85);
    leg2->AddEntry(gDU[0],"MARATHON kin1","P");
    leg2->AddEntry(gDU[1],"MARATHON kin2","P");
    leg2->AddEntry(gDU[2],"MARATHON kin3","P");
    leg2->AddEntry(gDU[3],"MARATHON kin4","P");
    leg2->AddEntry(gSLAC_np,"SLAC Bodek et al.","P");
    leg2->AddEntry(gNMC_np,"NMC","P");
//    leg2->AddEntry(gNMC_npcorr,"NMC with R*","P");
    leg2->Draw();


   TCanvas *c4=new TCanvas("c4","c4",1500,1500);
   TGraphErrors *M_dp=new TGraphErrors(11,M_x,M_ratio,0,M_ratio_err);
   TMultiGraph *mg4=new TMultiGraph();
   M_dp->SetMarkerStyle(8);
   M_dp->SetMarkerColor(2);
   gSLAC->SetMarkerStyle(8);
   gSLAC->SetMarkerColor(1);
   mg4->Add(gSLAC);
   mg4->Add(M_dp);
   mg4->SetTitle("sigma_d/sigma_p ;xbj");
   mg4->Draw("AP");

    auto leg4 = new TLegend(0.7,0.7,0.85,0.85);
    leg4->AddEntry(M_dp,"MARATHON","P");
    leg4->AddEntry(gSLAC,"SLAC Bodek et al.","P");
    leg4->Draw();

   TCanvas *c5=new TCanvas("c5","c5");
   Double_t factor_ratio[6]={0.0};
   Double_t factor_x[6]={0.0};
   for(int ii=2;ii<8;ii++){
	factor_ratio[ii-2]=D2_factor[1][ii]/H1_factor[1][ii];	           
        factor_x[ii-2]=H1_xbj[1][ii];
   }
   
   TGraph *gRatio=new TGraph(6,factor_x,factor_ratio);
   gRatio->SetMarkerStyle(8);
   gRatio->Draw("AP");
*/
}
