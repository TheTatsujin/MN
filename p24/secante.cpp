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


/// FUNCION QUE IMPLEMENTA EL METODO DE LA SECANTE
/// LA FUNCI�N DEVUELVE EL N�MERO DE ITERACIONES REALIZADAS SI TERMINA BIEN Y DEVUELVE -1
/// EN CASO CONTRARIO
int mn_secante (
real (*f)( real), /// funcion sobre la que se calcula el cero
real &x0, /// primera aproximaci�n ra�z que actualiza la funci�n
real &x1, /// segunda aproximaci�n ra�z que actualiza la funci�n
int NiterMax, /// n�mero de iteraciones m�ximo
real TOL) /// tolerancia para para el algoritmo
{
  if (f(x0) == 0. || f(x1) == 0.) return 0;
  real x2 = x1 - f(x1) / mn_slope(f, x0, x1);
  for(int i=0; i<NiterMax; i++){
    if(f(x2) == 0. || mn_distancia(x2, x1) <= TOL) return i;
    x0 = x1;
    x1 = x2;
    x2 = x1 - f(x1) / mn_slope(f, x0, x1);
  }
  return -1;
}

real mn_slope(
real (*f)( real),
real x0,
real x1)
{ return (f(x1) - f(x0)) / (x1 - x0); }
