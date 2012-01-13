// Copyright (c) 2002  Max-Planck-Institute Saarbruecken (Germany)
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
// Author(s)     : Peter Hachenberger

#include <CGAL/basic.h>
#ifdef CGAL_VSH_FILTERED
  #include <CGAL/Simple_cartesian.h>
  #include <CGAL/Lazy_kernel.h>
  #include <CGAL/Gmpq.h>
#else
  #ifdef CGAL_USE_LEDA
    #include <CGAL/leda_integer.h>
  #else
    #include <CGAL/Gmpz.h>
  #endif
  #include <CGAL/Cartesian.h>
  #include <CGAL/Homogeneous.h>
  #include <CGAL/Quotient.h>
  #include <CGAL/cartesian_homogeneous_conversion.h>
#endif

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/visual_hull_creator.h>
#include <CGAL/Nef_nary_intersection_3.h>
#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>
#include <fstream>
#include <list>

#ifdef CGAL_VSH_FILTERED
  typedef CGAL::Simple_cartesian<CGAL::Gmpq> EKernel;
  typedef CGAL::Lazy_kernel<EKernel> Kernel;
#else
  #ifdef CGAL_USE_LEDA
    typedef leda_integer NT;
  #else
    typedef CGAL::Gmpz NT;
  #endif
  typedef CGAL::Quotient<NT> CNT;
  typedef CGAL::Cartesian<CNT> CKernel;
  typedef CGAL::Homogeneous<NT> Kernel;
  typedef CKernel::FT FT;
  typedef CKernel::Point_3 CPoint;
#endif

typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 Plane_3;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef CGAL::Nef_nary_intersection_3<Nef_polyhedron> NaryInt;

#ifdef CGAL_VSH_FILTERED
Point_3 read_point(std::ifstream& in) {
  Point_3 p;
  in >> p;
  return p;
}
#else
Point_3 read_point(std::ifstream& in) {
  double x,y,z;
  in >> x;
  in >> y;
  in >> z;
  CPoint p(x,y,z);
  return quotient_cartesian_to_homogeneous(p);
}
#endif

int main(int argc, char* argv[]) {
  // We've put the typedefs here as VC7 gives us an ICE if they are global typedefs
  typedef Nef_polyhedron::SNC_structure SNC_structure;
  typedef CGAL::visual_hull_creator<SNC_structure> VHC;

  if(argc!=2) {
    std::cerr << "Usage: visual_hull file" << std::endl;
    std::cerr << "For more information read the README file" << std::endl;
    return 1;
  }

  std::ifstream in(argv[1]);
  
  NaryInt ni;

  CGAL::Timer t;

  Point_3 room_min = read_point(in);
  Point_3 room_max = read_point(in);

  int ncameras;
  in >> ncameras;
  for(int cam=0; cam<ncameras; ++cam) {

    Point_3 camera(read_point(in));

    int npolygons;
    in >> npolygons;

    std::list<std::list<Point_3> > polygon_list;
    for(int poly=0; poly<npolygons; ++poly) {

      int npoints;
      in >> npoints;

      std::list<Point_3> input_points;
      for(int pnt=0; pnt<npoints; ++pnt)
	input_points.push_back(read_point(in));
      polygon_list.push_back(input_points);
    }

    std::list<std::list<Point_3> >::iterator li;
    for(li=polygon_list.begin(); li!=polygon_list.end(); ++li) {
      std::list<Point_3>::iterator pi(li->begin()), pimin(pi), pi_next,pi_prev;
      for(; pi!=li->end(); ++pi) {
	if(CGAL::lexicographically_xyz_smaller(*pi,*pimin))
	  pimin=pi;
      }
      pi_next=pi_prev=pimin;
      ++pi_next;
      if(pi_next==li->end()) pi_next=li->begin();
      if(pi_prev==li->begin()) pi_prev=li->end();
      --pi_prev;
      if(CGAL::orientation(*pi_prev,*pimin,*pi_next,camera)
	 == CGAL::POSITIVE)
	li->reverse();
    }

    t.start();
    Nef_polyhedron N;
    VHC vhc(room_min, room_max, camera, polygon_list);
    N.delegate(vhc,true);
    CGAL_assertion(N.is_valid());
    t.stop();
    std::cerr << "create view " << t.time() << std::endl;
    t.reset();
    ni.add_polyhedron(N);
  }

  Nef_polyhedron result = ni.get_intersection();

  QApplication a(argc,argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron>* w =
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(result);
  a.setMainWidget(w);
  w->show();
  a.exec();
}
