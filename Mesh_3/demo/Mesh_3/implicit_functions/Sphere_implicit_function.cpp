// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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


const double radius = 1.;

class Sphere_implicit_function :
  public QObject,
  public Implicit_function_interface
{
  Q_OBJECT
  Q_INTERFACES(Implicit_function_interface)
  
public:
  virtual QString name() const { return "Sphere function"; }
  
  virtual double operator()(double x, double y, double z) const
  {
    return (x*x + y*y + z*z - radius);
  }
  
  virtual Bbox bbox() const
  {
    double r = radius * 1.2;
    return Bbox(-r,-r,-r,r,r,r);
  }
};



#include <QtPlugin>
Q_EXPORT_PLUGIN2(Sphere_implicit_function, Sphere_implicit_function)
#include "Sphere_implicit_function.moc"
