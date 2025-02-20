#include "mn_calculo_de_ceros.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

/// PARAMETROS DE LA DISTRIBUCION GAMMA
real alfa,beta_,d;

/// DISTRIBUCIÓN GAMMA
real Gamma(real x){
  if(x<0) return 0;
  return d*pow(x,alfa-1.)*exp(-beta_*x);
}

/// DERIVADA DE LA DISTRIBUCIÓN GAMMA
real Gammap(real x){
  if(x<0) return 0;
  return d*(alfa-1)*pow(x,alfa-2.)*exp(-beta_*x)-d*beta_*pow(x,alfa-1.)*exp(-beta_*x);
}

/// CALCULO DE LOS PARÁMETROS DE LA DISTRIBUCIÓN GAMMA A PARTIR DE LA
/// MEDIA Y VARIANZA MUESTRAL
void calculo_parametros_Gamma(real media, real varianza){
  /// calculo de alfa y beta_
  beta_=media/varianza;
  alfa=media*beta_;

  /// calculo de d a través de la integral (se verá en el tema 5)
  real h=0.001;
  real suma=0;
  for(real x=0;x<100;x+=h) suma+=h*pow(x,alfa-1.)*exp(-beta_*x);
  d=1./suma;

  cout << "\nalfa = " << alfa << " beta_ = " << beta_ << " d = " << d << "\n";

}


/// FUNCION QUE IMPLEMENTA EL METODO DE LA REGULA FALSI  Y DEVUELVE EL NÚMERO DE ITERACIONES
/// SI ALGO VA MAL DEVUELVE -1.  LA RAÍZ SE DEVUELVE COMO PARÁMETRO
int mn_regula_falsi (
real (*f)( real), /// función a la cual se calcula un cero
real &a, real &b, /// intervalo inicial para buscar la raíz
real &x, /// valor de salida de la raíz
real TOL,  /// tolerancia para parar las iteraciones del algoritmo
int NiterMax) /// número máximo de iteraciones permitidas
{
  if (f(a)*f(b) > 0.) return -1;
  real previous_x;
  x = a - ((b-a)/(f(b)-f(a)))*f(a);
  for (int i = 0; i < NiterMax; i++){
    if (f(x) == 0.) return i;
    (f(x)*f(a) < 0.) ? b = x : a = x;
    previous_x = x;
    x = a - ((b-a)/(f(b)-f(a)))*f(a);
    if (mn_distancia(previous_x, x) <= TOL) return i+1;
  }
  return -1;
}



