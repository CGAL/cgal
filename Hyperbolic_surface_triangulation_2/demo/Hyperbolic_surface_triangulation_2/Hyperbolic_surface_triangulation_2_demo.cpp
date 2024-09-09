// Copyright (c) 2024
// INRIA Nancy (France), and Université Gustave Eiffel Marne-la-Vallee (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Vincent Despré, Loïc Dubois, Monique Teillaud

#include "window.h"

#include <CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_surface_traits_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>
#include <CGAL/Hyperbolic_surface_triangulation_2.h>

using namespace CGAL;

typedef Simple_cartesian<Gmpq>                                                 Kernel;
typedef Hyperbolic_Delaunay_triangulation_traits_2<Kernel>              ParentTraits;
typedef Hyperbolic_surface_traits_2<ParentTraits>                      Traits;
typedef Hyperbolic_fundamental_domain_2<Traits>                         Domain;
typedef Hyperbolic_fundamental_domain_factory_2<Traits>                 Factory;
typedef Hyperbolic_surface_triangulation_2<Traits>                      Triangulation;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){
  // 1. Generate the triangulation
  Factory factory;
  Domain domain = factory.make_hyperbolic_fundamental_domain_g2(time(NULL));
  Triangulation triangulation = Triangulation(domain);
  triangulation.make_Delaunay();

  // 2. Draw the triangulation
  QApplication app(argc, argv);
  app.setApplicationName("Hyperbolic surfaces triangulation 2 Demo");

  DemoWindow window;
  window.item().draw_triangulation(triangulation);
  window.show();

  QStringList args = app.arguments();
  args.removeAt(0);
  return app.exec();
}
