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

class Klein_implicit_function :
  public QObject,
  public Implicit_function_interface
{
  Q_OBJECT
  Q_INTERFACES(Implicit_function_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.Mesh3Demo.Implicit_function_interface/1.0" FILE "klein_implicit_function.json")

public:
  virtual QString name() const { return "Klein function"; }

  virtual double operator()(double x, double y, double z) const
  {
    return   (x*x+y*y+z*z+2*y-1)
           * ( (x*x+y*y+z*z-2*y-1) *(x*x+y*y+z*z-2*y-1)-8*z*z)
           + 16*x*z* (x*x+y*y+z*z-2*y-1);
  }

  virtual Bbox bbox() const
  {
    const double radius = 6.;
    double r = radius * 1.1;
    return Bbox(-r,-r,-r,r,r,r);
  }
};

#include "Klein_implicit_function.moc"
