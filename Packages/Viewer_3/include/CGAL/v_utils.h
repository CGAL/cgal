// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Francois Rebufat <Francois.Rebufat@sophia.inria.fr>

#include <CGAL/IO/Color.h>
#include <utility>
// #include <vector>
#ifndef DRAWABLE
#include <CGAL/Drawable_object.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#define V_UTILS

CGAL_BEGIN_NAMESPACE

// typedef int           Size;
// typedef unsigned char Precision;


void invert(double *, double *);

void invert(double *);
void add_mat(double *, double *);

void set_matrix(double *, double *);
std::vector<double> apply_mat(double *, double , double , double );

std::vector<double> compute_plan(double , double , double ,double ,
			    double , double ,double , double ,
			    double );

double intersect_plan(const std::vector<double> &, double , double );

std::vector<double> translate(const std::vector<double> &, double,double,double );

std::vector<double> normal(double , double , double , double , double ,
	    double , double , double , double);




void set_color(Color );



double rad2deg(double );

std::pair<double,double> get_angles(const double & ,const double &  ,const
			       double &,const double & ,const
			       double  &,const double & );


void draw_tube(double , double , double , double , double ,
	       double , Size, Precision );


void draw_sphere(double , double , double , Size , Precision );

void draw_triangle(double , double , double ,double , double
		   , double ,double , double , double );

void draw_triangle_nice(double , double , double ,double , double
		   , double ,double , double , double );

void draw_triangle_2(double, double, double, double
		   ,double , double);

void draw_shrink_triangle_2(double, double, double, double
		   ,double , double);



void add_mat(double m[16], double res[16])
{
  for (int i=0; i<16 ; i++)
    res[i]=res[i] + m[i];
}

void set_matrix(double* m, double* n)
{
  for (int i=0; i<16 ; i++)
    m[i]=n[i];
}

void invert(double m[16], double res[16])
{

  double x[4][4];
  int i, j, k;
  double out[4][4];
  for (i=0; i<=3; i++)
    for(j=0; j <=3; j++)
      if (i==j) 
	out[i][j]=1;
      else
	out[i][j]=0;

  for (i=0; i<=3; i++)
    x[0][i]=m[i];  
  for (i=4; i<=7; i++)
    x[1][i-4]=m[i];
  for (i=8; i<=11; i++)
    x[2][i-8]=m[i];
  for (i=12; i<=15; i++)
    x[3][i-12]=m[i];

  for (i = 0; i < 4; i++) {
	if (x[i][i] != 1.0) {
	    double divby = x[i][i];
	    for (j = 0; j < 4; j++) {
		out[i][j] /= divby;
		x[i][j] /= divby;
	    }
	}
	for (j = 0; j < 4; j++) {
	    if (j != i) {
		if (x[j][i] != 0.0) {
		    double mulby = x[j][i];
		      for (k = 0; k < 4; k++) {
			  x[j][k] -= mulby * x[i][k];
			  out[j][k] -= mulby * out[i][k];
		      }
		}
	    }
	}
  }
  for (i=0; i<=3; i++)
    res[i]=out[0][i];  
  for (i=4; i<=7; i++)
    res[i]=out[1][i-4]; 
  for (i=8; i<=11; i++)
    res[i]=out[2][i-8]; 
  for (i=12; i<=15; i++)
    res[i]=out[3][i-12]; 

}

void invert(double m[16])
{
  invert(m,m);
}
 


std::vector<double> apply_mat(double* m, double x, double y, double z)
{

  std::vector<double> v(3);
  v[0]=m[0]*x + m[4]*y + m[8]*z + m[12];
  v[1]=m[1]*x + m[5]*y + m[9]*z + m[13];
  v[2]=m[2]*x + m[6]*y + m[10]*z + m[14];
  return(v);
}

std::vector<double> compute_plan(double x1, double y1, double z1,double x2,
			    double y2, double z2,double x3, double y3,
			    double z3)
{
  std::vector<double> v(4);
  v[0] = y1*(z2-z3) + y2*(z3-z1) + y3*(z1-z2);
  v[1] = z1*(x2-x3) + z2*(x3-x1) + z3*(x1-x2);
  v[2] = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);
  v[3] = -(x1*(y2*z3-y3*z2) + x2*(y3*z1-y1*z3) + x3*(y1*z2 - y2*z1));
  return v;
}

double intersect_plan(const std::vector<double> &v, double x, double y)
{
  return((v[0]*x + v[1]*y + v[3])/v[2]);
}

  
std::vector<double> translate(const std::vector<double> &v, double x,double y,double z)
{
  std::vector<double> vr(3);
  vr[0]= v[0] + x;
  vr[1]= v[1] + y;
  vr[2]= v[2] + z;
  return vr;
}

std::vector<double> normal(double x1, double y1, double z1, double x2, double y2,
	    double z2, double x3, double y3, double z3)
{
 std::vector<double> v(3);
  v[0] = y1*(z2-z3) + y2*(z3-z1) + y3*(z1-z2);
  v[1] = z1*(x2-x3) + z2*(x3-x1) + z3*(x1-x2);
  v[2] = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);
  double den = sqrt(pow(v[0],2) + pow(v[1],2) + pow(v[2],2) );
  v[0] = v[0]/den; v[1] = v[1]/den; v[2] = v[2]/den; 
  return v;
}



const double PI=3.14159265358979323844;

void set_color(Color c)
{
  glColor4ub(c.red(),c.green(),c.blue(),c.alpha());
}

double rad2deg(double a)
{
return((a*180)/PI);
}

std::pair<double,double> get_angles(const double &x1,const double &y1,const
			       double &z1,const double &x2,const
			       double &y2,const double &z2)
{
  double X, Y, aux ;
  double l=sqrt(pow(x1-x2,2) +pow(y1-y2,2) + pow(z1-z2,2));
  
  Y=asin((x2-x1)/l);
  aux=-((y2-y1)/l)/cos(asin((x2-x1)/l));
  if (aux > 1.0)  aux = 1.0;
  if (aux < -1.0) aux = -1.0;
  X= asin(aux);
  if ((z2-z1) <0) X=PI-X;

  return std::pair<double,double>(rad2deg(X),rad2deg(Y));
}


void draw_tube(double x1, double y1, double z1, double x2, double y2,
	       double z2, Size s, Precision p)
  {
    

   GLUquadricObj *q= gluNewQuadric();
   std::pair<double, double> ag=get_angles(x1,y1,z1,x2,y2,z2);
   double l=sqrt(pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2));
   glPushMatrix();
   glTranslatef(x1, y1, z1);

   glRotatef(ag.first,1,0,0);
   glRotatef(ag.second,0,1,0);
   gluQuadricNormals(q, GL_SMOOTH);
   gluCylinder(q, s, s, l,p, 1);
   glPopMatrix();
   gluDeleteQuadric(q);
  }

void draw_sphere(double x, double y, double z, Size s, Precision p)
{
  glPushMatrix();
  glTranslatef(x,y,z);
  GLUquadricObj *q = gluNewQuadric();
  gluQuadricNormals(q, GL_SMOOTH);
  gluQuadricTexture(q, GL_FALSE);
  gluSphere(q,s,p,p);
  glPopMatrix();
  gluDeleteQuadric(q);
}


void draw_triangle(double x1, double y1, double z1,double x2, double
		   y2, double z2,double x3, double y3, double z3)

{
  glBegin(GL_TRIANGLES);

  std::vector<double> v1(3);
  v1 = normal(x1,y1,z1,x2,y2,z2,x3,y3,z3);
       
  glNormal3d(v1[0],v1[1],v1[2]);
  glVertex3f(x1,y1,z1);
 glNormal3d(v1[0],v1[1],v1[2]);
  glVertex3f(x2,y2,z2);
 glNormal3d(v1[0],v1[1],v1[2]);
  glVertex3f(x3,y3,z3);
  glEnd();
}


void draw_triangle_nice(double x1, double y1, double z1,double x2, double
		   y2, double z2,double x3, double y3, double z3)

{
  glBegin(GL_TRIANGLES);

  std::vector<double> v1(3);
  v1 = normal(x1,y1,z1,x2,y2,z2,x3,y3,z3);
       
  glNormal3d(v1[0],v1[1],v1[2]);
  glVertex3f(x1,y1,z1);
  glNormal3d(-v1[0],-v1[1],-v1[2]);
  glVertex3f(x2,y2,z2);
  glNormal3d(v1[0],v1[1],v1[2]);
  glVertex3f(x3,y3,z3);
  glEnd();
}



void draw_triangle_2(double x1, double y1, double x2, double
		   y2,double x3, double y3)

{
  glBegin(GL_TRIANGLES);
  glVertex2f(x1,y1);
  glVertex2f(x2,y2);
  glVertex2f(x3,y3);
  glEnd();
}

void draw_shrink_triangle_2(double x1, double y1, double x2, double
		   y2,double x3, double y3)
{
  double bx= (x1+x2+x3)/3;
  double by= (y1+y2+y3)/3;

  x1= x1+(bx - x1)/3;
  y1= y1+(by - y1)/3;
  x2= x2+(bx - x2)/3;
  y2= y2+(by - y2)/3;
  x3= x3+(bx - x3)/3;
  y3= y3+(by - y3)/3;
  glBegin(GL_TRIANGLES);
  glVertex2f(x1,y1);
  glVertex2f(x2,y2);
  glVertex2f(x3,y3);
  glEnd();
}

CGAL_END_NAMESPACE
