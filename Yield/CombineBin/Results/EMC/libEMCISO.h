Double_t H3_ISO(Double_t x){
        Double_t R=1.0;
        R=1.09527-1.73012*x+14.9157*pow(x,2)-72.0011*pow(x,3)+209.873*pow(x,4)-384.816*pow(x,5)+434.79*pow(x,6)
          -276.485*pow(x,7)+75.879*pow(x,8);
        return R;
}

Double_t He_ISO(Double_t x){
        Double_t R=1.0;
        R=1.02927-0.336485*x+2.5182*pow(x,2)-11.9978*pow(x,3)+35.3705*pow(x,4)-69.2632*pow(x,5)+87.672*pow(x,6)
          -64.3681*pow(x,7)+20.8184*pow(x,8);
        return R;
}
 
Double_t SLAC_EMC( Double_t x){
	 Double_t A=3.0;
         Double_t Cx=0.017+0.018*log(x)+0.005*log(x)*log(x);
	 Double_t alphax=-0.07+2.189*x-24.667*pow(x,2)+145.291*pow(x,3)-497.237*pow(x,4)+1013.129*pow(x,5)
			-1208.393*pow(x,6)+775.767*pow(x,7)-205.872*pow(x,8);
	 Double_t R=exp(Cx)*pow(A,alphax);
	 return R;
}
