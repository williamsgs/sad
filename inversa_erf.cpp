#include <iostream>
#include <cmath>
#include <math.h>
#include <iomanip>
double simpson(int n, double a, double b);
double gauss (double x);
double erf (double x);
double erf1 (double x, double xinit);
double d_erf1 (double x);
double erfrational (double x);
double erftaylor (double x);
double NR(double xinit, double x, double eps);
double erfsymp (double x);
double factorial (double a);
double potencia (double x,double y);
double factorial2 (double f2);
const double raizdepi = 1.772453850905516027298167483341145;
const double euler = 2.71828182845904523536028747135;
const double q = 67.0;
int main ()
{
  double xinit = 0.0;
  double eps = 1.0e-4;
  double x = -0.8;
  double inverse = NR(xinit, x, eps);
  double i=0.0;
  double j = 2.0*i;
  double h = (2.0*i)-1.0;
  double a = 0.0;
  double b = x;
  double y = b;
  double n = 50000.0;
  double w = 0.0;
  for (w= -0.99999999999999999999999999999999999999; w <= 0.99999999999999999999999999999999999999;w += 0.0001)
    {
      std::cout << w << "\t" << std::setprecision(17)<< NR(xinit, w, eps)  <<std::endl;
    } 
  return 0;
}

double factorial (double a)
{
  if (a>1.0)
    return (a * factorial (a-1.0));
  else
    return 1.0;
}
double factorial2 (double f2)
{

  if (f2>1.0)
    return (f2 * factorial2 (f2-2.0));
  else
    return 1.0;
}

double potencia(double x, double y)
{
  double p=1.0;
  for(double i=1.0; i<=y; i += 1.0)
    p=p*x;
  return p;
}

double erftaylor (double x)
{
  double w =0.0;
  double suma =0.0;
  double pi = 3.142502653589793;

  for (double i = 0.0; i <= q ;i += 1.0){
    double j= 2.0*i;
    suma = suma + (potencia (-1.0, i)*(potencia(x,j)/( factorial(i)*(2.0*i+1.0))));
  }
  w = (suma * (2.0*x))/raizdepi;
  return w;
}
double erfsymp (double x)
{
  double w =0.0;
  double suma =0.0;
  double z =0.0;
  double g = -1.0 * (x*x);

  for (double i = 0.0; i <= q ;i += 1.0)
    {
      double t = x;
      double h = (2.0*i)-1.0;
      double a = 2.0*(t*t);

      suma = suma + (potencia (-1.0, i)*((factorial2 (h))/(potencia (a ,i))));

    }
  w = ((pow (euler,g))/(x*raizdepi)) * (1.0 + suma);
  z = 1.0 - w ;
  return z;
}
double erfrational (double x)
{
  double ak0 =1.000000000000013;
  double ak1 =1.474885681937094;
  double ak2 =1.089127207353042;
  double ak3 =0.481934851516365;
  double ak4 =0.133025422885837;
  double ak5 =0.021627200301105;
  double ak6 =0.001630015433745;
  double ak70 =0.000000000566405;
  double ak7 = ak70 * -1.0 ;
  double bk0 =1.000000000000000;
  double bk1 =2.603264849035166;
  double bk2 =3.026597029346489;
  double bk3 =2.046071816911715;
  double bk4 =0.873486411474986;
  double bk5 =0.237214006125950;
  double bk6 =0.038334123870994;
  double bk7 =0.002889083295887;
  double w =0.0;
  double suma1 =0.0;
  double suma2 =0.0;
  double z =0.0;
  double g = -1.0 * (x*x);
  suma1 = (((pow (x,0.0))* ak0) + ((pow (x,1.0))* ak1)+ ((pow (x,2.0))* ak2) + ((pow (x,3.0))* ak3) + ((pow (x,4.0))* ak4) + ((pow (x,5.0))* ak5) + ((pow (x,6.0))* ak6) + ((pow (x,7.0))* ak7));
  suma2 =(((pow (x,0.0))* bk0)+((pow (x,1.0))* bk1) + ((pow (x,2.0))* bk2) + ((pow (x,3.0))* bk3) + ((pow (x,4.0))* bk4) + ((pow (x,5.0))* bk5)+ ((pow (x,6.0))* bk6)+ ((pow (x,7.0))* bk7));
  w = pow (euler,g) * ( suma1 / suma2);
  z = 1.0 - w; 
  return z;
}
 
double erf (double x)
{
  double i = std::abs (x);
  double a =(erfsymp (i))*-1.0;
  double b =(erfrational (i))*-1.0;
  double c =(erftaylor (i))*-1.0;
  double d = erftaylor (i);
  double w =erfrational (i);
  double g = erfsymp (i);
  double y =0.0;

  if (x <= -6.0)
    {
      y = a;
    }
  else if ((x <= -1.5) && ( x > -6.0))
    {
      y =  b;
    }
  else if ((x < 0.0) && ( x > -1.5))
    {
      y =  c;
    }
  else if ((x > 0.0) && ( x <= 1.5))
    {
      y =  d;
    }
  else if ((x > 1.5) && ( x <= 6.0))
    {
      y =  w;
    }
  else if (x > 6.0)
    {
      y =  g;
    }
  else if (x==0.0)
    {
      y =  0.0;
    }
  return y;
}
double NR(double xinit, double x, double eps)
{
  double temp;
  double f = 0.0;
  double fdash = 0.0;
  f = erf1 ( x,xinit);
  fdash = d_erf1 ( xinit);
  double xnplus1=xinit-(f/fdash);
  temp=xinit;
  xinit=xnplus1;
  
  if(std::abs(temp-xinit)>=eps){
    return NR(xinit,x,eps);
    
  }
  return xinit;
}
double d_erf1 (double xinit)
{
  double a=0.0;
  double g = -1.0 * (xinit*xinit);
  a = ((2.0*(pow (euler,g)))/(raizdepi));
  return a;
}
double erf1 (double x,double xinit)
{
  double a=0.0;
  a = (erf(xinit)) - x;
  return a;
}
double simpson(int n, double a, double b)
{
  double m = std::abs (b);
  int i;
  double dx, sum,w,y, z, k;

  dx = (m - a) / n;
  sum = gauss(a) + gauss(m);
  for (i = 1.0; i < n; i++) {
    w = a + dx * i;
    sum += 2.0 * (1.0 + i%2) * gauss(w);
  }
  sum *= dx/3.0;
  y = sum * (2.0/raizdepi);
  z = y * -1.0;
  if (b>0.0)
    {
      k = y ;
    }
  else if (b<0.0)
    {
      k = y*-1.0;
    }
  else if (b == 0.0)
    {
      k= 0.0;
    }

  return k;
}
double gauss(double x)
{
  double a=0.0;
  double g = -1.0 * (x*x);
  a=pow (euler,g);
  return a;
}
