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

/// FUNCION QUE IMPLEMENTA EL CALCULO DE LAS RAICES DE UN POLINOMIO DE GRADO 2
/// DEVUELVE EL NÚMERO DE RAÍCES OBTENIDAS
/// LAS RAÍCES SALEN ORDENADAS POR VALOR ABSOLUTO. ES DECIR |x1|<=|x2|
int mn_ceros_pol_grado_2(
real a, real b, real c, // coeficientes polinomio de grado 2
real &x1, // primera raíz
real &x2) // segunda raíz
{
  if(a==0.) return(0);
  real dis=b*b-4*a*c;
  if(dis<0.) return(0);
  else if(dis==0.){
    x1=-b/(2*a);
    return(1);
  }
  dis=sqrt(dis);
  if(b>0){
    x1=(-b+dis)/(2*a);
    x2=(-b-dis)/(2*a);
  }
  else{
    x1=(-b-dis)/(2*a);
    x2=(-b+dis)/(2*a);
  }
  return(2);
}

/// FUNCION QUE IMPLEMENTA EL METODO DE MULLER USANDO UNA APROXIMACIÓN DE LAS DERIVADAS
/// LA FUNCIÓN DEVUELVE EL NÚMERO DE ITERACIONES REALIZADAS SI TERMINA BIEN Y DEVUELVE -1
/// EN CASO CONTRARIO
int mn_muller (
real (*f)( real), /// funcion sobre la que se calcula el cero
real &x0, /// raíz inicial que actualiza la función
int NiterMax, /// número de iteraciones máximo
real TOL) /// tolerancia para parar el algoritmo
{
  if(f(x0) == 0.) return 0;
  real quadratic_solution = x0;
  if (!mn_quadratic(f, x0, quadratic_solution)) return -1;

  real x1 = x0 + quadratic_solution;
  for(int i=0;i<NiterMax;i++){
    if(!mn_quadratic(f, x0, quadratic_solution)) break;
    if(f(x1) == 0. || mn_distancia(x1, x0) <= TOL) return i;
    x0 = x1;
    if (!mn_quadratic(f, x0, quadratic_solution)) break;
    x1 = x0 + quadratic_solution;
  }
  return -1;
}

bool mn_quadratic(
real (*f)(real),
real x0,
real &closest_solution)
{
    real solution_b = 0.;
    return (mn_ceros_pol_grado_2(
            mn_derivada2(f, x0)/2.,
            mn_derivada1(f, x0),
            f(x0),
            closest_solution, solution_b) > 0);
}

/// APROXIMACIÓN DERIVADA PRIMERA DE UNA FUNCIÓN
real mn_derivada1(
real (*f)( real), /// función que se deriva
real x) /// punto donde se evalúa la derivada primera
{
   /// CALCULO DE LA RAIZ CUADRADA DE LA UNIDAD DE REDONDEO u
   static real sqrt_u = sqrt( (double) mn_precision_aritmetica());

   /// CALCULO DESPLAZAMIENTO DE x PARA CALCULAR LA DERIVADA
   /// NOS ALEJAMOS DE x CON LA MITAD DE BITS QUE PERMITE LA ARITMÉTICA
   real h=(mn_abs(x)+1.)*sqrt_u;

   return (f(x+h)-f(x-h))/(2.*h);

}

/// APROXIMACIÓN DERIVADA SEGUNDA DE UNA FUNCIÓN
real mn_derivada2(
real (*f)( real), /// función que se deriva
real x) /// punto donde se evalúa la derivada primera
{
   /// CALCULO DE LA RAIZ CUADRADA DE LA UNIDAD DE REDONDEO u
   static real sqrt_u = sqrt( (double) mn_precision_aritmetica());

   /// CALCULO DESPLAZAMIENTO DE x PARA CALCULAR LA DERIVADA
   /// NOS ALEJAMOS DE x CON LA MITAD DE BITS QUE PERMITE LA ARITMÉTICA
   real h=(mn_abs(x)+1.)*sqrt_u;

   return (f(x+h)+f(x-h)-2*f(x))/(h*h);

}
