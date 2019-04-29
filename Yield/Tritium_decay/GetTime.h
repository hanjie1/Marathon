Int_t GetTime(TChain *T){

      THaRun* run = 0;
      T->LoadTree(0);
      TDirectory* fDir = T->GetDirectory();
      if (fDir) fDir->GetObject("Run_Data",run);
      if (!run) return 0;

      TDatime datetime = run->GetDate();
      Int_t time2 = datetime.GetDate();
      Int_t year2=time2/10000;
      Int_t month2=(time2-year2*10000)/100;
      Int_t day2=time2-year2*10000-month2*100;

      Int_t time1=20171023;
      Int_t year1=time1/10000;
      Int_t month1=(time1-year1*10000)/100;
      Int_t day1=time1-year1*10000-month1*100;

      if(year2<=year1){
         cout<<"This is 2017 run ?!"<<endl;
         return 0;
      }      
      if(month2>=month1){
         cout<<"This is not a Spring run ?!"<<endl;
         return 0;
      }

      Int_t m_day[12]={31,28,31,30,31,30,31,31,30,31,30,31};
      Int_t totaldays=0;

      Int_t inidays=m_day[month1-1]-day1+m_day[10]+m_day[11];
      
      switch(month2){
	case 1:
	   totaldays=inidays+day2;
	   break;
	case 2:
	   totaldays=inidays;
           for(int ii=0;ii<month2-1;ii++)
	       totaldays+=m_day[ii];
	   totaldays=totaldays+day2;
	   break;
	case 3:
	   totaldays=inidays;
           for(int ii=0;ii<month2-1;ii++)
	       totaldays+=m_day[ii];
	   totaldays=totaldays+day2;
	   break;
	case 4:
	   totaldays=inidays;
           for(int ii=0;ii<month2-1;ii++)
	       totaldays+=m_day[ii];
	   totaldays=totaldays+day2;
	   break;
	case 5:
	   totaldays=inidays;
           for(int ii=0;ii<month2-1;ii++)
	       totaldays+=m_day[ii];
	   totaldays=totaldays+day2;
	   break;
      } 

      return totaldays;

}
