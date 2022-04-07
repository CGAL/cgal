// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
  Q_PLUGIN_METADATA(IID "com.geometryfactory.Mesh3Demo.Implicit_function_interface/1.0" FILE "sphere_implicit_function.json")

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

#include "Sphere_implicit_function.moc"
