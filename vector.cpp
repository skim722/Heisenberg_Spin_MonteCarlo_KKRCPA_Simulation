#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "magnet.h"
#define PI 3.14159265

//#include "vector.h"
// det, (vector)dot, cross, transpose, (matrix)dot, inverse
////////////////////////////////////////////////////////////////////////
double det(const double a[3][3])
{
 return (
          a[0][0] * (a[1][1]*a[2][2] - a[1][2]*a[2][1])
        - a[0][1] * (a[1][0]*a[2][2] - a[1][2]*a[2][0])
        + a[0][2] * (a[1][0]*a[2][1] - a[1][1]*a[2][0])
        );
}
////////////////////////////////////////////////////////////////////////
double dot(const double a[3], const double b[3])
{
 return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
////////////////////////////////////////////////////////////////////////
void cross(const double a[3], const double b[3], double c[3])
{
c[0] = a[1]*b[2] - a[2]*b[1];
c[1] =-a[0]*b[2] + a[2]*b[0];
c[2] = a[0]*b[1] - a[1]*b[0];
}
////////////////////////////////////////////////////////////////////////
void transpose(const double a[3][3], double at[3][3])
{
 at[0][0] = a[0][0];
 at[0][1] = a[1][0];
 at[0][2] = a[2][0];

 at[1][0] = a[0][1];
 at[1][1] = a[1][1];
 at[1][2] = a[2][1];

 at[2][0] = a[0][2];
 at[2][1] = a[1][2];
 at[2][2] = a[2][2];
}
////////////////////////////////////////////////////////////////////////
void dot(const double a[3][3], const double b[3], double ab[3])
{
 ab[0] = dot(a[0],b);
 ab[1] = dot(a[1],b);
 ab[2] = dot(a[2],b);
}

////////////////////////////////////////////////////////////////////////
void dot(const double a[3][3], const double b[3][3], double ab[3][3])
{
 double bt[3][3];

 transpose(b,bt);

 ab[0][0] = dot(a[0],bt[0]);
 ab[0][1] = dot(a[0],bt[1]);
 ab[0][2] = dot(a[0],bt[2]);

 ab[1][0] = dot(a[1],bt[0]);
 ab[1][1] = dot(a[1],bt[1]);
 ab[1][2] = dot(a[1],bt[2]);

 ab[2][0] = dot(a[2],bt[0]);
 ab[2][1] = dot(a[2],bt[1]);
 ab[2][2] = dot(a[2],bt[2]);
}

////////////////////////////////////////////////////////////////////////
void inverse(const double a[3][3], double ai[3][3])
{
 double ac[3][3];
 cross(a[1],a[2],ac[0]);
 cross(a[2],a[0],ac[1]); // minus sign
 cross(a[0],a[1],ac[2]);

 double di = 1./det(a);
 
 ac[0][0] *= di;
 ac[0][1] *= di;
 ac[0][2] *= di;
 ac[1][0] *= di;
 ac[1][1] *= di;
 ac[1][2] *= di;
 ac[2][0] *= di;
 ac[2][1] *= di;
 ac[2][2] *= di;

 transpose(ac,ai);
}
////////////////////////////////////////////////////////////////////////
int lattice2metric(const double a[3][3], double g[3][3])
{
 for (int i=0;i<3;i++)
     for (int j=0;j<3;j++)
         g[i][j] = dot(a[i],a[j]);
}
////////////////////////////////////////////////////////////////////////
int normalize(double a[3]) {double r=1./sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
                       a[0]*=r;a[1]*=r;a[2]*=r;}

////////////////////////////////////////////////////////////////////////
int vscopy(int n, double a, double x[], double y[])
          {for (int i=0;i<n;i++) y[i]=a*x[i];}

////////////////////////////////////////////////////////////////////////
int vaxpy(int n, double a, double x[], double y[])
          {for (int i=0;i<n;i++) y[i]+=a*x[i];}

////////////////////////////////////////////////////////////////////////
inline double mydot(double v[], double w[])
{
 double s = 0.;
 for (int i=0;i<DIMENSIONS;i++) s += v[i]*w[i];
 return s;
}

////////////////////////////////////////////////////////////////////////
inline double mynorm(double v[]) 
{
 double s=0.0;
 for (int i=0;i<DIMENSIONS;i++) s +=v[i]*v[i]; 
return sqrt(s); 
}
////////////////////////////////////////////////////////////////////////
void conv_spherical( double xyz[3],double spherical[3])
{
  spherical[0]=sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
                    spherical[1]=180/PI*atan(xyz[1]/xyz[0]);
                    spherical[2]=180/PI*atan(sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1])/xyz[2]);
}
///////////////////////////////////////////////////////////////////////
void conv_cartesian(double spherical[3],double xyz[3])
{
  xyz[0]=spherical[0]*cos(spherical[1]*PI/180)*sin(spherical[2]*PI/180);
  xyz[1]=spherical[0]*sin(spherical[1]*PI/180)*sin(spherical[2]*PI/180);
  xyz[2]=spherical[0]*cos(spherical[2]*PI/180);
}

////////////////////////////////////////////////////////////////////
double angle( double v[], double u[])
{
  double param=mydot(u,v)/(mynorm(u)*mynorm(v));
  return acos(param) * 180.0 / PI;
}
///////////////////////////////////////////////////////////////////////
void rotation(const double u[3], double theta, double s[3],double rot_s[3])
{
   rot_s[0]=(u[0]*u[0]*(1-cos(theta*PI/180))+cos(theta*PI/180))*s[0]+(u[0]*u[1]*(1-cos(theta*PI/180))-u[2]*sin(theta*PI/180))*s[1]+(u[0]*u[2]*(1-cos(theta*PI/180))+u[1]*sin(theta*PI/180))*s[2];
   rot_s[1]=(u[0]*u[1]*(1-cos(theta*PI/180))+u[2]*sin(theta*PI/180))*s[0]+(u[1]*u[1]*(1-cos(theta*PI/180))+cos(theta*PI/180))*s[1]+(u[1]*u[2]*(1-cos(theta*PI/180))-u[0]*sin(theta*PI/180))*s[2];
   rot_s[2]=(u[0]*u[2]*(1-cos(theta*PI/180))-u[1]*sin(theta*PI/180))*s[0]+(u[1]*u[2]*(1-cos(theta*PI/180))+u[0]*sin(theta*PI/180))*s[1]+(u[2]*u[2]*(1-cos(theta*PI/180))-cos(theta*PI/180))*s[2];
}
