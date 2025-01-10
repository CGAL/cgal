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
// Author(s)     : Vincent Despré, Loïc Dubois, Marc Pouget, Monique Teillaud

// This file contains the declaration and the implementation of the class Delaunay_triangulation_on_hyperbolic_surface_2

#ifndef CGAL_DELAUNAY_HYPERBOLIC_SURFACE_TRIANGULATION_2
#define CGAL_DELAUNAY_HYPERBOLIC_SURFACE_TRIANGULATION_2

#include <CGAL/Triangulation_on_hyperbolic_surface_2.h>
#include <CGAL/basic.h>
#include <CGAL/Combinatorial_map.h>

#include <map>
#include <vector>
#include <queue>

namespace CGAL {

/*
Represents a geodesic Delaunay triangulation of a closed orientable hyperbolic surface.
*/

template<class Traits, class Attributes = Combinatorial_map_with_cross_ratios_item<Traits>>
  class Delaunay_triangulation_on_hyperbolic_surface_2: public Triangulation_on_hyperbolic_surface_2<Traits, Attributes>{
  public:
  typedef Triangulation_on_hyperbolic_surface_2<Traits, Attributes>    Base;//or T_on_HS_2
  typedef typename Base::Combinatorial_map_with_cross_ratios Combinatorial_map_with_cross_ratios;
  typedef typename Base::Anchor Anchor;

  typedef typename Combinatorial_map_with_cross_ratios::Dart_descriptor                                 Dart_descriptor;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<0>             Vertex_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<1>             Edge_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<2>             Face_range;

  typedef typename Combinatorial_map_with_cross_ratios::Dart_const_descriptor                           Dart_const_descriptor;
  typedef typename Combinatorial_map_with_cross_ratios::Dart_const_range                                Dart_const_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_const_range<0>       Vertex_const_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_const_range<1>       Edge_const_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_const_range<2>       Face_const_range;

  typedef typename Traits::FT                                                                Number;
  typedef typename Traits::Complex                                                           Complex_number;
  typedef typename Traits::Hyperbolic_point_2                                                Point;
  typedef Hyperbolic_isometry_2<Traits>                                                      Isometry;
  typedef Hyperbolic_fundamental_domain_2<Traits>                                            Domain;

  Delaunay_triangulation_on_hyperbolic_surface_2() {};
  Delaunay_triangulation_on_hyperbolic_surface_2(const Hyperbolic_fundamental_domain_2<Traits>& domain);
  Delaunay_triangulation_on_hyperbolic_surface_2(Combinatorial_map_with_cross_ratios& cmap, Anchor& anchor);

  bool is_Delaunay_flippable(Dart_const_descriptor dart) const;
  void flip(Dart_descriptor dart);
  bool is_Delaunay() const;
  int make_Delaunay();

protected:
  Dart_descriptor pick_edge_to_flip();
  Dart_const_descriptor pick_edge_to_flip() const;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class Traits, class Attributes>
Delaunay_triangulation_on_hyperbolic_surface_2<Traits,Attributes>::Delaunay_triangulation_on_hyperbolic_surface_2(const Domain& domain)
  : Base(domain)
{
  make_Delaunay();
}

template<class Traits, class Attributes>
  Delaunay_triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Delaunay_triangulation_on_hyperbolic_surface_2(Combinatorial_map_with_cross_ratios& cmap, Anchor& anchor)
  : Base(cmap, anchor)
  {
    make_Delaunay();
  }

////////////////////////////////////////////////////////////////////////////////
template<class Traits, class Attributes>
bool Delaunay_triangulation_on_hyperbolic_surface_2<Traits, Attributes>::is_Delaunay_flippable(Dart_const_descriptor dart) const{
  CGAL_precondition(Base::is_valid());
  return ( Base::get_cross_ratio(dart).imag()>Number(0) );
}

template<class Traits, class Attributes>
void Delaunay_triangulation_on_hyperbolic_surface_2<Traits, Attributes>::flip(Dart_descriptor dart){
  CGAL_precondition(Base::is_valid());
  // First gather all the information needed

  Dart_descriptor a = Base::opposite(dart); // Get a fresh descriptor
  Dart_descriptor b = Base::ccw(a);
  Dart_descriptor c = Base::cw(a);

   Dart_descriptor d = Base::opposite(a);
   Dart_descriptor e = Base::ccw(d);
   Dart_descriptor f = Base::cw(d);

   Complex_number cross_ratio_AB = Base::get_cross_ratio(e);
   Complex_number cross_ratio_BC = Base::get_cross_ratio(f);
   Complex_number cross_ratio_CD = Base::get_cross_ratio(b);
   Complex_number cross_ratio_DA = Base::get_cross_ratio(c);
   Complex_number cross_ratio_AC = Base::get_cross_ratio(a);

   // Modify the anchor

   if (this->_anchor.dart == a){
     this->_anchor.dart = e;
     this->_anchor.vertices[1] = Point(Base::fourth_point_from_cross_ratio(this->_anchor.vertices[1], this->_anchor.vertices[2], this->_anchor.vertices[0], cross_ratio_AC));
   } else if (this->_anchor.dart == b){
     this->_anchor.vertices[2] = Point(Base::fourth_point_from_cross_ratio(this->_anchor.vertices[0], this->_anchor.vertices[1], this->_anchor.vertices[2], cross_ratio_AC));
   } else if (this->_anchor.dart == c){
     this->_anchor.vertices[2] = Point(Base::fourth_point_from_cross_ratio(this->_anchor.vertices[2], this->_anchor.vertices[0], this->_anchor.vertices[1], cross_ratio_AC));
   } else if (this->_anchor.dart == d){
     this->_anchor.dart = b;
     this->_anchor.vertices[1] = Point(Base::fourth_point_from_cross_ratio(this->_anchor.vertices[1], this->_anchor.vertices[2], this->_anchor.vertices[0], cross_ratio_AC));
   } else if (this->_anchor.dart == e){
     this->_anchor.vertices[2] = Point(Base::fourth_point_from_cross_ratio(this->_anchor.vertices[0], this->_anchor.vertices[1], this->_anchor.vertices[2], cross_ratio_AC));
   } else if (this->_anchor.dart == f){
     this->_anchor.vertices[2] = Point(Base::fourth_point_from_cross_ratio(this->_anchor.vertices[2], this->_anchor.vertices[0], this->_anchor.vertices[1], cross_ratio_AC));
   }

   // Compute the new cross ratios

   Complex_number one (Number(1), Number(0));
   Complex_number cross_ratio_BD = (cross_ratio_AC) / ((cross_ratio_AC) - one) ;
   Complex_number cross_ratio_AB_2 = one - (one - (cross_ratio_AB)) * (cross_ratio_AC) ;
   Complex_number cross_ratio_BC_2 = one - (one - (cross_ratio_BC)) / (cross_ratio_BD) ;
   Complex_number cross_ratio_CD_2 = one - (one - (cross_ratio_CD)) * (cross_ratio_AC) ;
   Complex_number cross_ratio_DA_2 = one - (one - (cross_ratio_DA)) / (cross_ratio_BD) ;

   // Make the topological flip

   this->_combinatorial_map.template unlink_beta<1>(a);
   this->_combinatorial_map.template unlink_beta<1>(b);
   this->_combinatorial_map.template unlink_beta<1>(c);

   this->_combinatorial_map.template unlink_beta<1>(d);
   this->_combinatorial_map.template unlink_beta<1>(e);
   this->_combinatorial_map.template unlink_beta<1>(f);


   this->_combinatorial_map.template link_beta<1>(b, a);
   this->_combinatorial_map.template link_beta<1>(a, f);
   this->_combinatorial_map.template link_beta<1>(f, b);

   this->_combinatorial_map.template link_beta<1>(e, d);
   this->_combinatorial_map.template link_beta<1>(d, c);
   this->_combinatorial_map.template link_beta<1>(c, e);

   // And give the new cross ratios to the edges

   this->_combinatorial_map.template info<1>(a) = cross_ratio_BD;
   this->_combinatorial_map.template info<1>(e) = cross_ratio_AB_2;
   this->_combinatorial_map.template info<1>(f) = cross_ratio_BC_2;
   this->_combinatorial_map.template info<1>(b) = cross_ratio_CD_2;
   this->_combinatorial_map.template info<1>(c) = cross_ratio_DA_2;

   // Take care of the particular cases where we need to "flip again"

   if (Base::opposite(e) == b){
     this->_combinatorial_map.template info<1>(e) = one - (one - cross_ratio_AB_2) * (cross_ratio_AC) ;
   }

   if (Base::opposite(f) == c){
     this->_combinatorial_map.template info<1>(f) = one - (one - cross_ratio_BC_2) / (cross_ratio_BD) ;
   }
}

template<class Traits, class Attributes>
bool Delaunay_triangulation_on_hyperbolic_surface_2<Traits, Attributes>::is_Delaunay() const{
  if (! Base::is_valid()){
    return false;
  }
  return (pick_edge_to_flip() == nullptr);
}

template<class Traits, class Attributes>
int Delaunay_triangulation_on_hyperbolic_surface_2<Traits, Attributes>::make_Delaunay(){
  CGAL_precondition(this->is_valid());
  int number_of_flips_done = 0;

  Dart_descriptor edge_to_flip = pick_edge_to_flip();
  while (edge_to_flip != nullptr){
    flip(edge_to_flip);
    edge_to_flip = pick_edge_to_flip();
    number_of_flips_done++;
  }

  return number_of_flips_done;
}

////////////////////////////////////////////////////////////////////////////////
template<class Traits, class Attributes>
  typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Dart_descriptor Delaunay_triangulation_on_hyperbolic_surface_2<Traits, Attributes>::pick_edge_to_flip(){
  auto &cm=this->_combinatorial_map.darts();
  for (auto it = cm.begin(); it != cm.end(); ++it){
    if ( is_Delaunay_flippable(it) ){
      return it;
    }
  }
  return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits, class Attributes>
  typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Dart_const_descriptor Delaunay_triangulation_on_hyperbolic_surface_2<Traits, Attributes>::pick_edge_to_flip() const{
  const auto &cm=this->_combinatorial_map.darts();
  for (auto it = cm.begin(); it != cm.end(); ++it){
    if ( is_Delaunay_flippable(it) ){
      return it;
    }
  }
  return nullptr;
}

} // namespace CGAL

#endif // CGAL_DELAUNAY_HYPERBOLIC_SURFACE_TRIANGULATION_2
