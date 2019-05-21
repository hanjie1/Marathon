int GetRunList(vector<vector<Int_t> >& runList,const int kin,TString target)
{
     int run_number,nrun=0;
     Ssiz_t from=0;
     TString content;
     TString rundat=Form("/w/halla-scifs17exp/triton/Runlist/%s_kin%d.dat",target.Data(),kin);

     if(!gSystem->AccessPathName(rundat)){
        ifstream infile;
        infile.open(rundat);
        TString tmp,content;
        tmp.ReadLine(infile);
        tmp.ReadLine(infile);
        if(tmp.ReadLine(infile)){
           while(tmp.Tokenize(content,from,","))
           {
              run_number = atoi(content);
              if(run_number>0){
                runList.push_back(vector<Int_t>());
                runList[nrun].push_back(run_number);
                runList[nrun].push_back(0);
                nrun++;
	      }
           }
        }
        infile.close();
     }
     else{
	  rundat=Form("/w/halla-scifs17exp/triton/Runlist/%s_kin%d_1st.dat",target.Data(),kin);
          if(!gSystem->AccessPathName(rundat)){
              ifstream infile;
              infile.open(rundat);
              TString tmp,content;
              tmp.ReadLine(infile);
              tmp.ReadLine(infile);
              if(tmp.ReadLine(infile)){
                while(tmp.Tokenize(content,from,","))
                {
                   run_number = atoi(content);
                   if(run_number>0){
                      runList.push_back(vector<Int_t>());
                      runList[nrun].push_back(run_number);
                      runList[nrun].push_back(1);
                      nrun++;
	           }
                }
              }
              infile.close();
          } 
          else{
                cout<<"Can't find the runlist !!"<<endl;
                return 0;
          }
          rundat=Form("/w/halla-scifs17exp/triton/Runlist/%s_kin%d_2nd.dat",target.Data(),kin);
          if(!gSystem->AccessPathName(rundat)){
              ifstream infile;
              infile.open(rundat);
              TString tmp,content;
              tmp.ReadLine(infile);
              tmp.ReadLine(infile);
              from=0;
              if(tmp.ReadLine(infile)){
                while(tmp.Tokenize(content,from,","))
                {
                   run_number = atoi(content);
                   if(run_number>0){
                      runList.push_back(vector<Int_t>());
                      runList[nrun].push_back(run_number);
                      runList[nrun].push_back(2);
                      nrun++;
	           }
                }
              }
              infile.close();
          }  
          rundat=Form("/w/halla-scifs17exp/triton/Runlist/%s_kin%d_3rd.dat",target.Data(),kin);
          if(!gSystem->AccessPathName(rundat)){
              ifstream infile;
              infile.open(rundat);
              TString tmp,content;
              tmp.ReadLine(infile);
              tmp.ReadLine(infile);
              from=0;
              if(tmp.ReadLine(infile)){
                while(tmp.Tokenize(content,from,","))
                {
                   run_number = atoi(content);
                   if(run_number>0){
                      runList.push_back(vector<Int_t>());
                      runList[nrun].push_back(run_number);
                      runList[nrun].push_back(3);
                      nrun++;
	           }
                }
              }
              infile.close();
          }


     }
     return nrun;
}
