Double_t D2_EMC(Double_t xbj){
         Double_t R=1.0;
         R=0.959252+1.70218*xbj-21.3065*pow(xbj,2)+141.459*pow(xbj,3)-568.77*pow(xbj,4)+1443.77*pow(xbj,5)
           -2338.81*pow(xbj,6)+2353.76*pow(xbj,7)-1344.9*pow(xbj,8)+334.502*pow(xbj,9);

         return R;
}

Double_t He_EMC(Double_t x){
        Double_t R=1.0;
        R=0.949891+2.38613*x-30.205*pow(x,2)+200.737*pow(x,3)-807.611*pow(x,4)+2050.04*pow(x,5)-3319.42*pow(x,6)
          +3338.26*pow(x,7)-1906.05*pow(x,8)+473.934*pow(x,9);
        return R;
}

Double_t H3_EMC(Double_t x){
        Double_t R=1.0;
        R=0.939068+2.78422*x-35.2978*pow(x,2)+234.697*pow(x,3)-940.068*pow(x,4)+2368.52*pow(x,5)-3802.3*pow(x,6)
          +3791.3*pow(x,7)-2147.29*pow(x,8)+529.779*pow(x,9);

        return R;
}

Double_t D2EMC_ERROR(Double_t x){
	 Double_t R=0.0;
         R=1.69123-45.5375*x+585.576*pow(x,2)-4174.74*pow(x,3)+17904.8*pow(x,4)-47620*pow(x,5)+78775.9*pow(x,6)
          -78722.1*pow(x,7)+43482.4*pow(x,8)-10191.3*pow(x,9);
	 R=R/100.0;

        return R;
}

Double_t H3EMC_ERROR(Double_t x){
	 Double_t R=0.0;
         R=6.74195-180.328*x+2182.81*pow(x,2)-14749.7*pow(x,3)+60718.2*pow(x,4)-157320*pow(x,5)+257010*pow(x,6)
          -256507*pow(x,7)+142755*pow(x,8)-33945.1*pow(x,9);
	 R=R/100.0;

        return R;
}

Double_t HeEMC_ERROR(Double_t x){
	 Double_t R=0.0;
         R=12.7752-290.548*x+2969.05*pow(x,2)-17442.9*pow(x,3)+64366.9*pow(x,4)-153392*pow(x,5)+234949*pow(x,6)
          -222825*pow(x,7)+118926*pow(x,8)-27288.7*pow(x,9);
	 R=R/100.0;

        return R;
}
