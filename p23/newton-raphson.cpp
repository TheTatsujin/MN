#include "mn_calculo_de_ceros.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

/// PARAMETROS DE LA DISTRIBUCION GAMMA
real alfa,beta_,d;

/// DISTRIBUCI�N GAMMA
real Gamma(real x){
  if(x<0) return 0;
  return d*pow(x,alfa-1.)*exp(-beta_*x);
}

/// DERIVADA DE LA DISTRIBUCI�N GAMMA
real Gammap(real x){
  if(x<0) return 0;
  return d*(alfa-1)*pow(x,alfa-2.)*exp(-beta_*x)-d*beta_*pow(x,alfa-1.)*exp(-beta_*x);
}

/// CALCULO DE LOS PAR�METROS DE LA DISTRIBUCI�N GAMMA A PARTIR DE LA
/// MEDIA Y VARIANZA MUESTRAL
void calculo_parametros_Gamma(real media, real varianza){
  /// calculo de alfa y beta_
  beta_=media/varianza;
  alfa=media*beta_;

  /// calculo de d a trav�s de la integral (se ver� en el tema 5)
  real h=0.001;
  real suma=0;
  for(real x=0;x<100;x+=h) suma+=h*pow(x,alfa-1.)*exp(-beta_*x);
  d=1./suma;

  cout << "\nalfa = " << alfa << " beta_ = " << beta_ << " d = " << d << "\n";

}

/// FUNCION QUE IMPLEMENTA EL METODO DE NEWTON-RAPHSON APROXIMANDO LA FUNCION DERIVADA
/// LA FUNCI�N DEVUELVE EL N�MERO DE ITERACIONES REALIZADAS SI TERMINA BIEN Y DEVUELVE -1
/// EN CASO CONTRARIO
int mn_newton_raphson (
real (*f)( real), /// funcion sobre la que se calcula el cero
real &x0, /// ra�z inicial que actualiza la funci�n
int NiterMax, /// n�mero de iteraciones m�ximo
real TOL) /// tolerancia para parar el algoritmo
{
    if (f(x0) == 0.) return 0;
    real df_x0 = mn_derivada1(f, x0);

    if (df_x0 == 0.) return -1;
    real x_nexti = x0 - f(x0) / df_x0;

    for(int i=0; i<NiterMax; i++){
        if (f(x_nexti) == 0.) return i;
        if (mn_distancia(x0, x_nexti) <= TOL) return i;
        x0 = x_nexti;
        df_x0 = mn_derivada1(f, x0);
        if (df_x0 == 0.) break;
        x_nexti = x0 - f(x0) / df_x0;
    }
    return -1;

}



/// APROXIMACI�N DERIVADA PRIMERA DE UNA FUNCI�N
real mn_derivada1(
real (*f)( real), /// funci�n que se deriva
real x) /// punto donde se eval�a la derivada primera
{
   /// CALCULO DE LA RAIZ CUADRADA DE LA UNIDAD DE REDONDEO u
   static real sqrt_u = sqrt( (double) mn_precision_aritmetica());

   /// CALCULO DESPLAZAMIENTO DE x PARA CALCULAR LA DERIVADA
   /// NOS ALEJAMOS DE x CON LA MITAD DE BITS QUE PERMITE LA ARITM�TICA
   real h=(mn_abs(x)+1.)*sqrt_u;

   return (f(x+h)-f(x-h))/(2.*h);

}
