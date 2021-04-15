// Copyright (c) 2005-2009  INRIA Sophia-Antipolis (France).
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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/CGAL_Ipelet_base.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <boost/format.hpp>


namespace CGAL_alpha_shapes{

typedef CGAL::Exact_predicates_inexact_constructions_kernel       Kernel;
typedef CGAL::Regular_triangulation_vertex_base_2<Kernel>         Rvb;
typedef CGAL::Alpha_shape_vertex_base_2<Kernel,Rvb>               Vb;
typedef CGAL::Regular_triangulation_face_base_2<Kernel>           Rf;
typedef CGAL::Alpha_shape_face_base_2<Kernel,Rf>                  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>               Tds;
typedef CGAL::Regular_triangulation_2<Kernel,Tds>                 Regular;
typedef CGAL::Alpha_shape_2<Regular>                              Alpha_shape_2;

const std::string  sublabel[] = {
  "k-th Alpha-shape", "Help"
};

const std::string  helpmsg[] = {
  "Draw alpha-shape for the k-th critical alpha value"
};


class ASphapeIpelet
  : public CGAL::Ipelet_base<Kernel,2> {
public:
  ASphapeIpelet()
    : CGAL::Ipelet_base<Kernel,2>("Alpha-shapes",sublabel,helpmsg){}
  void protected_run(int);
};


void ASphapeIpelet::protected_run(int fn)
{
  if (fn==1) {
    show_help();
    return;
  }

  std::list<Weighted_point_2> LWP;

  std::list<Circle_2> cir_list;
  std::list<Point_2> pt_list;




  read_active_objects(
    CGAL::dispatch_or_drop_output<Point_2,Circle_2>(std::back_inserter(pt_list),
                                                    std::back_inserter(cir_list))
  );


  if (pt_list.empty() && cir_list.empty()) {
    print_error_message("No circle nor point selected");
    return;
  }


  for (std::list<Point_2>::iterator it=pt_list.begin();it!=pt_list.end();++it)
    LWP.push_back(Weighted_point_2(*it,0));
  for (std::list<Circle_2>::iterator it=cir_list.begin();it!=cir_list.end();++it)
    LWP.push_back(Weighted_point_2(it->center(),it->squared_radius()));

  Alpha_shape_2 A(LWP.begin(),LWP.end());
  int alpha=-1;
  int nb_ret;
  boost::tie(nb_ret,alpha)=request_value_from_user<int>((boost::format("# Spectral critical value (0-%d)") % A.number_of_alphas()).str() );
  if (nb_ret == -1) return;

  if(alpha<0 || (std::size_t) alpha>A.number_of_alphas()){
    print_error_message("Not a good value");
    return;
  }


  A.set_alpha(alpha==0?(std::max)(std::numeric_limits<double>::epsilon(),A.get_nth_alpha(0)/2.):
              (std::size_t) alpha==A.number_of_alphas()?A.get_nth_alpha(alpha-1)+1:A.get_nth_alpha(alpha-1)/2.+A.get_nth_alpha(alpha)/2.);
  for ( Alpha_shape_2::Alpha_shape_edges_iterator it=A.alpha_shape_edges_begin();it!=A.alpha_shape_edges_end();++it)
    draw_in_ipe(A.segment(*it));

  for (Alpha_shape_2::Finite_faces_iterator it=A.finite_faces_begin();it!=A.finite_faces_end();++it){
    if (A.classify(it)==Alpha_shape_2::INTERIOR){
      std::list<Point_2> LP;
      LP.push_back(Point_2(it->vertex(0)->point()));
      LP.push_back(Point_2(it->vertex(1)->point()));
      LP.push_back(Point_2(it->vertex(2)->point()));
      draw_polyline_in_ipe(LP.begin(),LP.end(),true,false,true);
    }
  }
  group_selected_objects_();
  return;
}

}


CGAL_IPELET(CGAL_alpha_shapes::ASphapeIpelet)


