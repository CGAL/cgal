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


const double radius = 4.;

class Tanglecube_implicit_function :
  public QObject,
  public Implicit_function_interface
{
  Q_OBJECT
  Q_INTERFACES(Implicit_function_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.Mesh3Demo.Implicit_function_interface/1.0" FILE "tanglecube_implicit_function.json")

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
