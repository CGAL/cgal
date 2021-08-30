// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot, Sylvain Pion

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/CGAL_Ipelet_base.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <boost/utility.hpp>

namespace CGAL_pca{

//~ typedef CGAL::Exact_predicates_exact_constructions_kernel                 Kernel;
typedef CGAL::Simple_cartesian<double>                 Kernel;

const std::string Slab[] = {
  "PCA","Help"
};

const std::string Hmsg[] = {
  "(Principal Component Analysis) given a set of points, draw a segment that is on the line defined by the eigen vector associated to the highest eigen value of the covariance matrix of the input points"
};

class pcaIpelet
  : public CGAL::Ipelet_base<Kernel,2>{
public:
  pcaIpelet()
    :CGAL::Ipelet_base<Kernel,2>("PCA",Slab,Hmsg){}
  void protected_run(int);
};


void pcaIpelet::protected_run(int fn)
{
  if (fn==1) {
    show_help();
    return;
  }


  std::list<Point_2> pt_list;
  std::list<Circle_2> cir_list;
  std::list<Polygon_2> poly_list;
  std::list<Kernel::Triangle_2> tri_list;
  std::list<Segment_2> sg_list;

  Iso_rectangle_2 bbox=
    read_active_objects(
      CGAL::dispatch_or_drop_output<Point_2,Polygon_2,Circle_2,Segment_2>(
        std::back_inserter(pt_list),
        std::back_inserter(poly_list),
        std::back_inserter(cir_list),
        std::back_inserter(sg_list)
      )
    );



  for (std::list<Polygon_2>::iterator it=poly_list.begin();it!=poly_list.end();++it)
    if (it->size()==3){
      tri_list.push_back(Kernel::Triangle_2(*(it->vertices_begin()),
                                            *boost::next(it->vertices_begin()),
                                            *boost::next(it->vertices_begin(),2)
                                            ));
    }
    else{
      print_error_message("This implementation is limited to triangles");
      return;
    }


  int s=0;

  if (!pt_list.empty()) s=1;
  if (!cir_list.empty()) s+=2;
  if (!tri_list.empty()) s+=4;
  if (!sg_list.empty()) s+=8;


  if (s==0) {
    print_error_message("Nothing is selected");
    return;
  }


  Kernel::Line_2 line;
  Kernel::Point_2 centroid;

  switch (s){
    case 1://points
      linear_least_squares_fitting_2(pt_list.begin(),pt_list.end(),line,centroid,CGAL::Dimension_tag<0>());
      break;
    case 2://circles
      linear_least_squares_fitting_2(cir_list.begin(),cir_list.end(),line,centroid,CGAL::Dimension_tag<2>());
      break;
    case 4://triangles
      linear_least_squares_fitting_2(tri_list.begin(),tri_list.end(),line,centroid,CGAL::Dimension_tag<2>());
      break;
    case 8://segments
      linear_least_squares_fitting_2(sg_list.begin(),sg_list.end(),line,centroid,CGAL::Dimension_tag<1>());
      break;
    default:
      print_error_message("Please select a set of points or segments or triangles or circles");
      return;
  }



  CGAL::Object obj_cgal = CGAL::intersection(line,bbox);
  Segment_2 seg;
  if (CGAL::assign(seg, obj_cgal))
    draw_in_ipe(seg);



}

}

CGAL_IPELET(CGAL_pca::pcaIpelet)

