// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#include <QObject>
#include "Implicit_function_interface.h"


void puiss(double& x, double& y, int n);

class Knot3_implicit_function :
  public QObject,
  public Implicit_function_interface
{
  Q_OBJECT
  Q_INTERFACES(Implicit_function_interface);
  
public:
  virtual QString name() const { return "Knot3 function"; }
  
  virtual double operator()(double a, double b, double c) const
  {
    double e=0.025;
    
    double x, y, z, t, den;
    den=1+a*a+b*b+c*c;
    x=2*a/den;
    y=2*b/den;
    z=2*c/den;
    t=(1-a*a-b*b-c*c)/den;
    
    double x19=x, y19=y, z17=z, t17=t;
    puiss(x19,y19,19);
    puiss(z17,t17,17);
    
    double f1 = z17-x19;
    double f2 = t17-y19;
    
    f1 = f1*f1;
    f2 = f2*f2;
    e=e*e/(den-1);
    
    return f1+f2-e;
  }
  
  virtual Bbox bbox() const
  {
    const double radius = 4.;
    double r = radius * 1.2;
    return Bbox(-r,-r,-r,r,r,r);
  }
};



void puiss(double& x, double& y, int n) {
  
  double xx = 1, yy = 0;
  
  while(n>0) {
    if (n&1) {
      double xxx = xx, yyy = yy;
      xx = xxx*x - yyy*y;
      yy = xxx*y + yyy*x;
    }
    
    double xxx = x, yyy = y;
    x=xxx*xxx-yyy*yyy;
    y=2*xxx*yyy;
    
    n/=2;
  }
  
  x = xx;
  y = yy;
}



#include <QtPlugin>
Q_EXPORT_PLUGIN2(Knot3_implicit_function, Knot3_implicit_function);
#include "Knot_implicit_function.moc"
