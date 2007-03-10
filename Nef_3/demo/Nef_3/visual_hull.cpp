// Copyright (c) 2002  Max-Planck-Institute Saarbruecken (Germany)
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
// Author(s)     : Peter Hachenberger

#include <CGAL/basic.h>
#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*){
  std::cout << "Sorry, this demo needs QT..." << std::endl; return 0;}
#else
#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Quotient.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/cartesian_homogeneous_conversion.h>
#include <CGAL/Nef_3/visual_hull_creator.h>
#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>
#include <fstream>
#include <list>

typedef CGAL::Gmpz NT;
typedef CGAL::Quotient<NT> CNT;
typedef CGAL::Cartesian<CNT> CKernel;
typedef CGAL::Homogeneous<NT> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;

typedef CKernel::FT FT;
typedef CKernel::Point_3 CPoint;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 Plane_3;

Point_3 read_point(std::ifstream& in) {
  double x,y,z;
  in >> x;
  in >> y;
  in >> z;
  CPoint p(x,y,z);
  return quotient_cartesian_to_homogeneous(p);
}

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

  typedef std::multimap<Nef_polyhedron::Size_type,Nef_polyhedron> PQ;
  typedef PQ::iterator      PQ_iterator;
  PQ pq;

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
    pq.insert(std::make_pair(N.number_of_vertices(),N));
    t.stop();
  }

  t.start();
  PQ_iterator i1, i2;
  while(pq.size() > 1) {
    i1 = i2 = pq.begin();
    ++i2;
    Nef_polyhedron Ntmp(i1->second * i2->second);
    pq.erase(i1);
    pq.erase(i2);
    pq.insert(std::make_pair(Ntmp.number_of_vertices(),Ntmp));
  }
  t.stop();

  Nef_polyhedron result(pq.begin()->second);

  QApplication a(argc,argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron>* w =
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(result);
  a.setMainWidget(w);
  w->show();
  a.exec();
}
#endif
