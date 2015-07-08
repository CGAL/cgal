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

const double radius = 4.;

class Tanglecube_implicit_function :
  public QObject,
  public Implicit_function_interface
{
  Q_OBJECT
  Q_INTERFACES(Implicit_function_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.Mesh3Demo.Implicit_function_interface/1.0")

public:
  virtual QString name() const { return "Tanglecube function"; }
  
  virtual double operator()(double x, double y, double z) const
  {
    double x2=x*x, y2=y*y, z2=z*z;
    double x4=x2*x2, y4=y2*y2, z4=z2*z2;
    return x4 - 5*x2 + y4 - 5*y2 + z4 - 5*z2 + 11.8;
  }
  
  virtual Bbox bbox() const
  {
    double r = radius * 1.2;
    return Bbox(-r,-r,-r,r,r,r);
  }
};

#include "Tanglecube_implicit_function.moc"
