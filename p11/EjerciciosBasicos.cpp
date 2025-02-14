/// INCLUSION DE LIBRERIAS NECESARIAS
#include <stdio.h>
#include "EjerciciosBasicos.h"

/// P11.1 FUNCIÓN QUE CALCULA LA MEDIA DE UN VECTOR
real mn_media(Array1D< real > &u){
  /// Variable para acumular el total.
  real total = 0;
  /// Bucle para sumar todos los elementos.
  for(int i=0; i<u.dim(); i++) total += u[i];
  /// Dividimos entre el total para cacular la media.
  return total / u.dim();
}

/// P11.2 FUNCIÓN QUE CALCULA EL MAXIMO DE UN VECTOR
real mn_max(Array1D< real > &u){
  /// Asumimos el primer valor como el máximo.
  real current_max = u[0];
  /// Bucle para comparar todos los elementos.
  for (int i=0; i<u.dim(); i++) {
    /// Si econtramos otro elemento mayor entonces actualizamos el máximo.
    if (current_max < u[i]) current_max = u[i];
  }
  /// Al llegar al final tenemos el máximo valor.
  return current_max;
}

/// P11.3 FUNCIÓN QUE CALCULA EL MINIMO DE UN VECTOR
real mn_min(Array1D< real > &u){
  /// Asumimos el primer valor como el mínimo.
  real current_min = u[0];
  /// Bucle para comparar todos los elementos.
  for (int i=0; i<u.dim(); i++) {
    /// Si econtramos otro elemento menor entonces actualizamos el mínimo.
    if (current_min > u[i]) current_min = u[i];
  }
  /// Al llegar al final tenemos el menor valor.
  return current_min;
}

/// P11.4 FUNCIÓN QUE ORDENA UN VECTOR DE MENOR A MAYOR
void mn_ordenar(Array1D< real > &u){
  /// Se recorre todo el vector para garantizar orden.
  for (int i = u.dim() - 1; i >= 0; i--){
    /// Otro cursor para elegir la nueva posición del elemento.
    int j = i;
    /// Buscamos un número mayor y cambiamos la posición.
    while (--j >=0 && u[j] > u[i]) std::swap(u[j], u[i]);

  }
}

/// P11.5 FUNCIÓN PARA MULTIPLICAR UNA MATRIZ POR UN VECTOR
Array1D< real > mn_multiplicacion_matriz_vector(Array2D< real > &A,Array1D< real > &u){
  return A*u;
}

/// P11.6 FUNCIÓN QUE DETERMINA SI UN NÚMERO ENTERO ES PRIMO
bool mn_es_primo(int i){
    /// Primero se verifica si es 1 o 2
    if (i <= 2 && i > 0) return true;
    /// Verificamos si es par
    if (i % 2 == 0) return false;

    /// Hallamos la raiz de i para realizar el cálculo
    int sqroot = sqrt(i);

    ///Bucle donde recorremos todos lo posibles números que podrían
    /// ser factores, saltando los valores pares.

    for (int j = 3; j <= sqroot; j+=2) {
        if (i % j == 0) return false;
    }

    /// Si un numero no es divisible por su raiz, sabemos que es primo
    return true;
}

/// P11.7 FUNCIÓN QUE CALCULA EL FACTORIAL DE UN NÚMERO NATURAL
real mn_factorial(int n){
  if (n == 1) return 1;
  return n*mn_factorial(n - 1);
}


/// P11.8 FUNCIÓN QUE CALCULA UNA POTENCIA CON UN NÚMERO NATURAL
/// NO SE PUEDE USAR LA FUNCIÓN pow()
real mn_potencia(real x,int n){
    real result = 1.;
    for (int i = 0; i < n; i++) result *= x;
    return result;
}

/// P11.9 FUNCIÓN QUE CALCULA EL DESARROLLO DE TAYLOR DE e^x
/// e^x = 1 + x + x^2/2! + ...... +x^n/n!
real mn_exp(real x,int n){
  real factorial = 1;
  real exponent = 1;
  real result = 1;

  for (int i = 1; i <= n; i++){
    exponent *= x;
    factorial *= i;
    result += exponent / factorial;
  }

  return result;
}

/// P11.10 FUNCIÓN QUE CALCULA EL DESARROLLO DE TAYLOR DE cos(x)
///  cos(x) = 1 - x^2/2! + x^4/4! - x^6/6!+...... +- x^(2n)/(2n)!
real mn_cos(real x,int n){
  real factorial = 1;
  real exponent = 1;
  real result = 1;
  real sign = 1;

  for (int i = 2; i < 2*n; i+=2){
    exponent *= x*x;
    factorial *= i * (i-1);
    sign *= -1;
    result += sign*exponent/factorial;
  }

  return result;
}

/// P11.11 FUNCIÓN QUE CALCULA EL DESARROLLO DE TAYLOR DE sin(x)
///  sin(x) = x - x^3/3! + x^5/5! - x^7/7!+...... +- x^(2n+1)/(2n+1)!
real mn_sin(real x,int n){
  real result = x;
  real exponent = x;
  real factorial = 1;
  real sign = 1;

  for (int i = 3; i <= 2*n + 1; i+=2){
    exponent *= x*x;
    factorial *= i * (i-1);
    sign *= -1;

    result += sign * exponent / factorial;
  }

  return result;
}

/// P11.12 FUNCIÓN QUE CALCULA EL DESARROLLO DE TAYLOR DE ln(x)
/// ln(x) = (x-1) - ((x-1)^2)/2 + ((x-1)^3)/3 - ((x-1)^4)/4+...... +- ((x-1)^n)/n
real mn_ln(real x,int n){
  real sign = 1;
  real exponent = x - 1;
  real result = x - 1;

  for (int i = 2; i <= n; i++) {
    exponent *= (x - 1);
    sign *= -1;
    result += sign * exponent / i;
  }
  return result;
}

/// P11.13 FUNCIÓN QUE CALCULA y^x DONDE y,x SON NÚMERO REALES
/// USAR LAS FUNCIONES IMPLEMENTADAS mn_exp() y mn_ln() TENIENDO EN CUENTA y^x=e^(x*ln(y))
real mn_pow(real y,real x,int n){
  return mn_exp(x*mn_ln(y, n), n);
}

/// P11.13 FUNCIÓN QUE CALCULA EL LIMITE DE LA SECUENCIA  yn=(1.+1./n).^n CUANDO n TIENDE A
/// INFINITO EL ALGORITMO PARA CUANDO LA DIFERENCIA EN VALOR ABSOLUTO DE LA DIFERENCIA
/// ENTRE UN VALOR DE LA SECUENCIA Y EL ANTERIOR ES INFERIOR AL PARAMETRO tolerancia
/// EL LIMITE DE LA SECUENCIA ES EL NUMERO e=2.71828182846
/// IMPORTANTE : PARA QUE LAS CONSTANTES LAS TRATE COMO NÚMEROS REALES HAY QUE AÑADIR UN .,
/// ES DECIR, POR EJEMPLO  1. (EN LUGAR DE 1). SI HACEMOS 1/2 EL RESULTADO ES CERO PORQUE HACE
/// LA DIVISIÓN EN PRECISIÓN ENTERA. SIN EMBARGO  1./2.=1./2=1/2.=0.5
real mn_limite1(real tolerancia){
 real n = 2;
 real current = 2.;
 real previous = 1.;

 while(fabs(current - previous) >= tolerancia){
    previous = current;
    current = mn_potencia(1. + 1./n, n);
    n++;
 }
 return current;
}

/// P11.14 FUNCIÓN QUE CALCULA EL LIMITE DE LA FUNCIÓN f(x)=sin(x)/x CUANDO x TIENDE HACIA 0.
/// EL PARAMETRO tolerancia SE UTILIZA PARA PARAR EL ALGORITMO CUANDO ESTAMOS CERCA DEL LÍMITE
/// EL VALOR DEL LÍMITE ES 1.
real mn_limite2(real tolerancia){
  real limit = 1.;
  real previous = 0.;
  real n = 1;
  while (fabs(limit - previous) >= tolerancia) {
    previous = limit;
    limit = sin(1./n) / (1./n);
    n++;
  }
  return limit;
}

/// P11.15 FUNCIÓN QUE CALCULA EL LIMITE DE LA SECUENCIA  yn=X(n+1)/X(n) DONDE X(n) ES LA
/// SUCESIÓN DE FIBONACCI DEFINIDA COMO X(n+1)=X(n)+X(n-1) EMPEZANDO POR X(1)=X(2)=1
/// EL ALGORITMO PARA CUANDO LA DIFERENCIA EN VALOR ABSOLUTO
/// ENTRE UN VALOR DE LA SECUENCIA yn Y EL ANTERIOR ES INFERIOR AL PARAMETRO tolerancia
/// EL LIMITE DE LA SECUENCIA yn ES EL NÚMERO AÚREO IGUAL A (1+SQRT(5))/2 = 1.618033988....
real mn_limite3(real tolerancia){
  real current_fibonacci = 2.;
  real previous_fibonacci = 1.;
  real current_ratio = 2.;
  real previous_ratio = 1.;

  while (fabs(current_ratio - previous_ratio) >= tolerancia) {
    real temporal = current_fibonacci;
    current_fibonacci = current_fibonacci + previous_fibonacci;
    previous_fibonacci = temporal;
    previous_ratio = current_ratio;
    current_ratio = current_fibonacci / previous_fibonacci;
  }

  return current_ratio;
}

/// P11.16 CÁLCULO DEL NÚMERO PI POR EL MÉTODO DE MONTECARLO. EL ÁREA DEL CÍRCULO DE RADIO
/// 1 ES PI. Y EL AREA DEL CUADRADO DE LADO 2 DONDE SE INSCRIBE EL CÍRCULO ES 4. POR TANTO
/// SI SE ELIGE UN PUNTO AL AZAR EN EL CUADRADO, LA PROBABILIDAD DE QUE CAIGA EN
/// EL CÍRCULO ES PI/4. EL MÉTODO DE MONTECARLO APROXIMA PI COGIENDO PUNTOS AL AZAR EN
/// EL CUADRADO [-1,1]x[-1,1] Y VIENDO QUE PROPORCIÓN CAE EN EL CÍRCULO.
/// NOTA : LA FUNCIÓN rand() DEVUELVE UN VALOR ENTERO ALEATORIO ENTRE 0 Y RAND_MAX
bool is_inside_circle(real x, real y){
    return x*x + y*y <= 1.;
}

real rand_number() {
    real sign = 1.;
    if (rand()*1. / RAND_MAX <= 0.5) sign *= -1.;
    return sign*rand()*1. / RAND_MAX;
}

real calculo_pi_montecarlo(int Nintentos){
    if (Nintentos == 0) return 0.;
    real inside = 0.;

    for (int i = 0; i < Nintentos; i++){
        if (is_inside_circle(rand_number(), rand_number())) inside++;
    }

    return 4. * inside / Nintentos;
}
