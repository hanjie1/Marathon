Double_t superR(Double_t xbj){
	 Double_t R=1.0;
	 R=1.02308-0.712485*xbj+8.53441*pow(xbj,2)-55.0103*pow(xbj,3)+211.738*pow(xbj,4)
	   -509.626*pow(xbj,5)+778.786*pow(xbj,6)-737.467*pow(xbj,7)+395.695*pow(xbj,8)-92.0484*pow(xbj,9);

	 return R;
}

Double_t D2_EMC(Double_t xbj){
	 Double_t R=1.0;
	 R=0.294856+5.50251*xbj-60.2084*pow(xbj,2)+364.398*pow(xbj,3)-1353.58*pow(xbj,4)+3208.92*pow(xbj,5)
	   -4876.38*pow(xbj,6)+4601.35*pow(xbj,7)-2456.78*pow(xbj,8)+567.901*pow(xbj,9);

	 return R;
}

int F2np(){
    ifstream infile1;
    infile1.open("../bin003/H3He_final.dat");
 
    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    Double_t H3He_xx[22]={0.0},H3He_Ratio[22]={0.0},H3He_Rerr[22]={0.0};
    while(tmp.ReadLine(infile1)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from,"  ");
          H3He_xx[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          H3He_Ratio[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          H3He_Rerr[nn-1]=atof(content.Data());
          from=0;
          nn++;
     }
    infile1.close();

    nn=0;
    ifstream infile2;
    infile2.open("../bin003/Dp_final.dat");
    Double_t Dp_xx[7]={0.0},Dp_Ratio[7]={0.0},Dp_Rerr[7]={0.0};
    while(tmp.ReadLine(infile2)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from,"  ");
          Dp_xx[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          Dp_Ratio[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from,"  ");
          Dp_Rerr[nn-1]=atof(content.Data());
          from=0;
          nn++;
     }
    infile2.close();

    Double_t H3He_np[22]={0.0},H3He_npErr[22]={0.0};
    for(int ii=0;ii<22;ii++){
	Double_t SR=superR(H3He_xx[ii]);
	H3He_np[ii]=(2.0*SR*H3He_Ratio[ii]-1.0)/(2.0-SR*H3He_Ratio[ii]);
	H3He_npErr[ii]=3.0*SR/((2.0-SR*H3He_Ratio[ii])*(2.0-SR*H3He_Ratio[ii]))*H3He_Rerr[ii];
    }

    Double_t Dp_np[7]={0.0},Dp_npErr[7]={0.0};
    for(int ii=0;ii<7;ii++){
	Double_t D2_R=D2_EMC(Dp_xx[ii])*2.0;
	Dp_np[ii]=Dp_Ratio[ii]/D2_R-1.0;
	Dp_npErr[ii]=1.0/D2_R*Dp_Rerr[ii];
    } 

    TGraphErrors *gH3He=new TGraphErrors(22,H3He_xx,H3He_np,0,H3He_npErr);
    TGraphErrors *gDp=new TGraphErrors(7,Dp_xx,Dp_np,0,Dp_npErr);

    TCanvas *c1=new TCanvas("c1","c1",1500,1500);
    gH3He->SetMarkerStyle(8);
    gH3He->SetMarkerColor(4);
    gH3He->Draw("AP");

    TCanvas *c2=new TCanvas("c2","c2",1500,1500);
    gDp->SetMarkerStyle(8);
    gDp->SetMarkerColor(4);
    gDp->Draw("AP");

    TCanvas *c3=new TCanvas("c3","c3",1500,1500);
    TMultiGraph *mg1=new TMultiGraph();
    gDp->SetMarkerColor(2);
    mg1->Add(gDp);
    mg1->Add(gH3He);
    mg1->Draw("AP"); 

    return 0;
}
