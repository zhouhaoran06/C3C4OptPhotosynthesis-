#' Calculate C4 RuBP carboxylation and PEP carboxylation photosynthesis rate
#'
#' Calculate C4 RuBP carboxylation and PEP carboxylation photosynthesis rate using von Caemmerer photosynthesis model for C4 pathway (von Caemmerer S, 2000) with optimality to stomatal resistance and energy allocation to leaf. The use of the function is as follows: C4RuBPCarboxyPEPCarboxy(Vcmax25,J25,Ca,Tleaf,Phis,VPD,StartV). The package "rootSolve" should be installed before using the function.
#' @param Vcmax25 Maximal velocity of Rubisco carboxylation at 25°C (μmol m-2 s-1).
#' @param J25 Electron transport at a specific light intensity at 25°C. The value should be adjusted at different light intensities (μmol m-2 s-1).
#' @param Ca CO2 concentration in the atmosphere (ppm).
#' @param Tleaf Leaf temperature (oC).
#' @param Phis Soil water potential (Mpa).
#' @param VPD Vapour pressure deficit (kPa).
#' @param StartV A vector that gives the start points for the estimation (c(A,rs,fl)).A is photosynthesis rate (μmol m-2 s-1),rs is stomatal resistance (μmol-1 m2 s1, the reciprocal of stomatal conductance) and fl is the energy allocation ratio of leaf to the total plant.
#'
#' @return This package will return a dataframe of the numerical solution of the model using multiroot in the package "rootSolve".
#' @return Model$root: A vector or solution of c(A,rs,fl).
#'
#' @export

C4RuBPCarboxyPEPCarboxy <- function(Vcmax25,J25,Ca,Tleaf,Phis,VPD,StartV){

  #Temperature adjustment for Vcmax,K(Kc and Ko) and gammastar from 25°C to Tleaf
  Vcmax0<-Vcmax25*exp(26.37-65.33/0.008314/(273+Tleaf))
  k<-302*exp(32.06-79.43/0.008314/(273+Tleaf))*(1+210/(256*exp(14.68-36.38/0.008314/(273+Tleaf))))
  gammastar<-36*exp(15.27-37.83/0.008314/(273+Tleaf))
  rm<-(1+exp((1.4*(273+Tleaf)-437.4)/0.008314/(273+Tleaf)))/exp(18.81-49.6/0.008314/(273+Tleaf))
  Kp<-26.29*2.1^(Tleaf/10-1)
  Vpmax<-exp(33.66-70.37/0.008314/(273+Tleaf))/(1+exp((0.4*(273+Tleaf)-117.91)/0.008314/(273+Tleaf)))

  #Define other constants
  Delta<-VPD*0.01*1.6 #change the unit
  dv<-2
  bv<-3.8
  kc<-0.001044
  rho<-18.3
  gbs<-0.003

  #Define the model
  Model<-function(x){
    PhistoV <- -(-(x[3]*Delta/(1-x[3])/kc/x[2]/rho)+Phis)/dv
    F1<-x[1]^3*(-rm-x[2])+x[1]*(-Ca*exp(-PhistoV^bv)*Vcmax0-exp(-PhistoV^bv)*Kp*Vcmax0)+x[1]^2*(Ca+Kp+exp(-PhistoV^bv)*rm*Vcmax0+exp(-PhistoV^bv)*x[2]*Vcmax0)-x[1]^3*(gbs*rm^2+2*gbs*rm*x[2]+gbs*x[2]^2)+Ca^2*exp(-PhistoV^bv)*gbs*Vcmax0+Ca*exp(-PhistoV^bv)*gbs*Kp*Vcmax0+Ca*exp(-PhistoV^bv)*Vcmax0*Vpmax-x[1]^2*(-2*Ca*gbs*rm-gbs*k*rm-gbs*Kp*rm-2*Ca*gbs*x[2]-gbs*k*x[2]-gbs*Kp*x[2]-exp(-PhistoV^bv)*gbs*rm^2*Vcmax0-2*exp(-PhistoV^bv)*gbs*rm*x[2]*Vcmax0-exp(-PhistoV^bv)*gbs*x[2]^2*Vcmax0-rm*Vpmax-x[2]*Vpmax)-Ca*exp(-PhistoV^bv)*gbs*Vcmax0*gammastar-exp(-PhistoV^bv)*gbs*Kp*Vcmax0*gammastar-x[1]*(Ca^2*gbs+Ca*gbs*k+Ca*gbs*Kp+gbs*k*Kp+2*Ca*exp(-PhistoV^bv)*gbs*rm*Vcmax0+exp(-PhistoV^bv)*gbs*Kp*rm*Vcmax0+2*Ca*exp(-PhistoV^bv)*gbs*x[2]*Vcmax0+exp(-PhistoV^bv)*gbs*Kp*x[2]*Vcmax0+Ca*Vpmax+exp(-PhistoV^bv)*rm*Vcmax0*Vpmax+exp(-PhistoV^bv)*x[2]*Vcmax0*Vpmax-exp(-PhistoV^bv)*gbs*rm*Vcmax0*gammastar-exp(-PhistoV^bv)*gbs*x[2]*Vcmax0*gammastar)
    F2<-x[1]*(x[1]*(gbs*(2*Ca+k+Kp)-x[1]*(1+2*gbs*(rm+x[2]))+Vpmax)+exp(-PhistoV^bv)*Vcmax0*(x[1]*(1+2*gbs*(rm+x[2]))-Vpmax+gbs*(-2*Ca-Kp+gammastar)))+1/dv/(1-x[3])/kc/x[2]^2/rho*bv*exp(-PhistoV^bv)*x[3]*Vcmax0*(x[1]^2*(rm+x[2])*(1+gbs*(rm+x[2]))+Ca*(gbs*(Ca+Kp)+Vpmax)-x[1]*(Kp+Ca*(1+2*gbs*(rm+x[2]))+(rm+x[2])*(Vpmax+gbs*(Kp-gammastar)))-gbs*(Ca+Kp)*gammastar)*Delta*(PhistoV^(-1+bv))
    F3<-x[1]*(gbs*(Ca+k)*(Ca+Kp)+Ca*Vpmax+x[1]*(3*x[1]*(rm+x[2])*(1+gbs*(rm+x[2]))-2*(Kp+Ca*(1+2*gbs*(rm+x[2]))+(rm+x[2])*(gbs*(k+Kp)+Vpmax)))+exp(-PhistoV^bv)*Vcmax0*(Kp-2*x[1]*(rm+x[2])*(1+gbs*(rm+x[2]))+Ca*(1+2*gbs*(rm+x[2]))+(rm+x[2])*(Vpmax+gbs*(Kp-gammastar))))+1/dv*bv*exp(-PhistoV^bv)*x[3]*Vcmax0*(x[1]^2*(rm+x[2])*(1+gbs*(rm+x[2]))+Ca*(gbs*(Ca+Kp)+Vpmax)-x[1]*(Kp+Ca*(1+2*gbs*(rm+x[2]))+(rm+x[2])*(Vpmax+gbs*(Kp-gammastar)))-gbs*(Ca+Kp)*gammastar)*(-Delta/(1-x[3])/kc/x[2]/rho-x[3]*Delta/(1-x[3])^2/kc/x[2]/rho)*(PhistoV^(-1+bv))
    c(F1=F1,F2=F2,F3=F3)
  }
  ss<-multiroot(f=Model,start=StartV)
  return(ss)
}
