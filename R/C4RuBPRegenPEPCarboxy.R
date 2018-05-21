#' Calculate C4 RuBP regeneration and PEP carboxylation photosynthesis rate
#'
#' Calculate C4 RuBP regeneration and PEP carboxylation photosynthesis rate using von Caemmerer photosynthesis model for C4 pathway (von Caemmerer S, 2000) with optimality to stomatal resistance and energy allocation to leaf. The use of the function is as follows: C4RuBPRegenPEPCarboxy(Vcmax25,J25,Ca,Tleaf,Phis,VPD,StartV). The package "rootSolve" should be installed before using the function.
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

C4RuBPRegenPEPCarboxy <- function(Vcmax25,J25,Ca,Tleaf,Phis,VPD,StartV){

  #Temperature adjustment for Vcmax,K(Kc and Ko) and gammastar from 25°C to Tleaf
  c_J0<-log(J25*(1+exp((0.65*298-200)/0.008314/298)))+43.9/0.008314/298
  J0<-exp(c_J0-43.9/0.008314/(273+Tleaf))/(1+exp((0.65*(273+Tleaf)-200)/0.008314/(273+Tleaf)))
  gammastar<-36*exp(15.27-37.83/0.008314/(273+Tleaf))
  rm<-(1+exp((1.4*(273+Tleaf)-437.4)/0.008314/(273+Tleaf)))/exp(18.81-49.6/0.008314/(273+Tleaf))
  Kp<-26.29*2.1^(Tleaf/10-1)
  Vpmax<-exp(34.18-70.37/0.008314/(273+Tleaf))/(1+exp((0.4*(273+Tleaf)-117.91)/0.008314/(273+Tleaf)))

  #Define other constants
  Delta<-VPD*0.01*1.6 #change the unit
  dj<-2.5
  bj<-3
  kc<-0.001044
  rho<-18.3
  gbs<-0.003
  e1<-4.5
  e2<-10.5
  x0<-0.6

  #Define the model
  Model<-function(x){
    PhistoJ <- -(-(x[3]*Delta/(1-x[3])/kc/x[2]/rho)+Phis)/dj
    F1<-x[1]^3*(-e1*rm-e1*x[2])+x[1]*(-Ca*exp(-PhistoJ^bj)*J0*x0-exp(-PhistoJ^bj)*J0*Kp*x0)+x[1]^2*(Ca*e1+e1*Kp+exp(-PhistoJ^bj)*J0*rm*x0+exp(-PhistoJ^bj)*J0*x[2]*x0)-x[1]^3*(e1*gbs*rm^2+2*e1*gbs*rm*x[2]+e1*gbs*x[2]^2)+Ca^2*exp(-PhistoJ^bj)*gbs*J0*x0+Ca*exp(-PhistoJ^bj)*gbs*J0*Kp*x0+Ca*exp(-PhistoJ^bj)*J0*Vpmax*x0-Ca*exp(-PhistoJ^bj)*gbs*J0*x0*gammastar-exp(-PhistoJ^bj)*gbs*J0*Kp*x0*gammastar-x[1]^2*(-2*Ca*e1*gbs*rm-e1*gbs*Kp*rm-2*Ca*e1*gbs*x[2]-e1*gbs*Kp*x[2]-e1*rm*Vpmax-e1*x[2]*Vpmax-exp(-PhistoJ^bj)*gbs*J0*rm^2*x0-2*exp(-PhistoJ^bj)*gbs*J0*rm*x[2]*x0-exp(-PhistoJ^bj)*gbs*J0*x[2]^2*x0-e2*gbs*rm*gammastar-e2*gbs*x[2]*gammastar)-x[1]*(Ca^2*e1*gbs+Ca*e1*gbs*Kp+Ca*e1*Vpmax+2*Ca*exp(-PhistoJ^bj)*gbs*J0*rm*x0+exp(-PhistoJ^bj)*gbs*J0*Kp*rm*x0+2*Ca*exp(-PhistoJ^bj)*gbs*J0*x[2]*x0+exp(-PhistoJ^bj)*gbs*J0*Kp*x[2]*x0+exp(-PhistoJ^bj)*J0*rm*Vpmax*x0+exp(-PhistoJ^bj)*J0*x[2]*Vpmax*x0+Ca*e2*gbs*gammastar+e2*gbs*Kp*gammastar-exp(-PhistoJ^bj)*gbs*J0*rm*x0*gammastar-exp(-PhistoJ^bj)*gbs*J0*x[2]*x0*gammastar)
    F2<-x[1]*(x[1]*(-x[1]*e1*(1+2*gbs*(rm+x[2]))+e1*(2*Ca*gbs+gbs*Kp+Vpmax)+e2*gbs*gammastar)+exp(-PhistoJ^bj)*J0*x0*(x[1]*(1+2*gbs*(rm+x[2]))-Vpmax+gbs*(-2*Ca-Kp+gammastar)))+1/dj/(1-x[3])/kc/x[2]^2/rho*bj*exp(-PhistoJ^bj)*x[3]*J0*x0*(x[1]^2*(rm+x[2])*(1+gbs*(rm+x[2]))+Ca*(gbs*(Ca+Kp)+Vpmax)-x[1]*(Kp+Ca*(1+2*gbs*(rm+x[2]))+(rm+x[2])*(Vpmax+gbs*(Kp-gammastar)))-gbs*(Ca+Kp)*gammastar)*Delta*PhistoJ^(bj-1)
    F3<-x[1]*(Ca*e1*(gbs*(Ca+Kp)+Vpmax)+exp(-PhistoJ^bj)*J0*x0*(Kp-2*x[1]*(rm+x[2])*(1+gbs*(rm+x[2]))+Ca*(1+2*gbs*(rm+x[2]))+(rm+x[2])*(Vpmax+gbs*(Kp-gammastar)))+e2*gbs*(Ca+Kp)*gammastar+x[1]*(3*x[1]*e1*(rm+x[2])*(1+gbs*(rm+x[2]))-2*e1*(Kp+Ca*(1+2*gbs*(rm+x[2]))+(rm+x[2])*(gbs*Kp+Vpmax))-2*e2*gbs*(rm+x[2])*gammastar))+1/dj*bj*exp(-PhistoJ^bj)*x[3]*J0*x0*(x[1]^2*(rm+x[2])*(1+gbs*(rm+x[2]))+Ca*(gbs*(Ca+Kp)+Vpmax)-x[1]*(Kp+Ca*(1+2*gbs*(rm+x[2]))+(rm+x[2])*(Vpmax+gbs*(Kp-gammastar)))-gbs*(Ca+Kp)*gammastar)*(-Delta/(1-x[3])/kc/x[2]/rho-x[3]*Delta/(1-x[3])^2/kc/x[2]/rho)*PhistoJ^(bj-1)
    c(F1=F1,F2=F2,F3=F3)
    }
  ss<-multiroot(f=Model,start=StartV)
  return(ss)
}
