void CalcRatioBoilingError()
{
     Double_t H1_re[1]={0.0},D2_re[1]={0.0},H3_re[1]={0.0},He3_re[1]={0.0};

     ifstream infile1;
     infile1.open("Results/H1_boiling.dat");
     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(infile1)){
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from," ");
          H1_re[nn]=atof(content.Data());
          from=0;
          nn++;
     }
     infile1.close();

     ifstream infile2;
     infile2.open("Results/D2_boiling.dat");
     from=0;
     nn=0;
     while(tmp.ReadLine(infile2)){
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from," ");
          D2_re[nn]=atof(content.Data());
          from=0;
          nn++;
     }
     infile2.close();

     ifstream infile3;
     infile3.open("Results/He3_boiling.dat");
     from=0;
     nn=0;
     while(tmp.ReadLine(infile3)){
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from," ");
          He3_re[nn]=atof(content.Data());
          from=0;
          nn++;
     }
     infile3.close();

     ifstream infile4;
     infile4.open("Results/H3_boiling.dat");
     from=0;
     nn=0;
     while(tmp.ReadLine(infile4)){
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from,"  ");
          tmp.Tokenize(content,from," ");
          H3_re[nn]=atof(content.Data());
          from=0;
          nn++;
     }
     infile4.close();

     ofstream outfile1;
     outfile1.open("Results/Dp_boiling.dat");
     Double_t Dp_re[1]={0.0},HeD_re[1]={0.0},H3D_re[1]={0.0},H3He_re[1]={0.0};
     for(int ii=0;ii<1;ii++){
	Dp_re[ii]=sqrt(D2_re[ii]*D2_re[ii]+H1_re[ii]*H1_re[ii]);
	outfile1<<ii<<"  "<<Dp_re[ii]<<endl;
     }
     outfile1.close();

     int kin[1]={16};
     ofstream outfile2;
     outfile2.open("Results/HeD_boiling.dat");
     ofstream outfile3;
     outfile3.open("Results/H3D_boiling.dat");
     ofstream outfile4;
     outfile4.open("Results/H3He_boiling.dat");
     for(int ii=0;ii<1;ii++){
	HeD_re[ii]=sqrt(D2_re[ii]*D2_re[ii]+He3_re[ii]*He3_re[ii]);
	H3D_re[ii]=sqrt(D2_re[ii]*D2_re[ii]+H3_re[ii]*H3_re[ii]);
	H3He_re[ii]=sqrt(H3_re[ii]*H3_re[ii]+He3_re[ii]*He3_re[ii]);
	outfile2<<kin[ii]<<"  "<<HeD_re[ii]<<endl;
	outfile3<<kin[ii]<<"  "<<H3D_re[ii]<<endl;
	outfile4<<kin[ii]<<"  "<<H3He_re[ii]<<endl;
     }
     outfile2.close();
     outfile3.close();
     outfile4.close();

}
