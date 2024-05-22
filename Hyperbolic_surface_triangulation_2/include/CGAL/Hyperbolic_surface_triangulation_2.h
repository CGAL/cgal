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

// This file contains the declaration and the implementation of the class Hyperbolic_surface_triangulation_2

#ifndef CGAL_HYPERBOLIC_SURFACE_TRIANGULATION_2
#define CGAL_HYPERBOLIC_SURFACE_TRIANGULATION_2

#include "Complex_without_sqrt.h"
#include "Hyperbolic_isometry_2.h"
#include "Hyperbolic_fundamental_domain_2.h"

#include <CGAL/basic.h>
#include <CGAL/Combinatorial_map.h>

#include <map>
#include <vector>
#include <queue>

namespace CGAL {

/*
Represents a geodesic triangulation of a closed orientable hyperbolic surface.
The triangulation is stored as combinatorial map decorated with one cross-ratio per edge.
It is also possible to specify an anchor for the triangulation. An anchor consists in 1) a dart of the combinatorial map, belonging by definition to a vertex V and a triangle T, together with
2) three points A,B,C in the hyperbolic plane. The points A,B,C are the three vertices in counter-clockwise order of a triangle. This triangle is a lift
of T, and A is a lift of V.
*/
template<class Traits>
class Hyperbolic_surface_triangulation_2{
public:
  struct Combinatorial_map_with_cross_ratios_item{
    template <class CMap>
    struct Dart_wrapper{
        typedef Cell_attribute<CMap, Complex_without_sqrt<typename Traits::FT>> Edge_attrib;
        typedef std::tuple<void,Edge_attrib,void>   Attributes;
    };
  };
  typedef Combinatorial_map<2,Combinatorial_map_with_cross_ratios_item>                                 Combinatorial_map_with_cross_ratios;

  struct Anchor{
    typename Combinatorial_map_with_cross_ratios::Dart_handle dart;
    typename Traits::Hyperbolic_point_2 vertices[3];
  };

  typedef typename Combinatorial_map_with_cross_ratios::Dart_handle                                     Dart_handle;
  typedef typename Combinatorial_map_with_cross_ratios::Dart_range                                      Dart_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<0>             Vertex_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<1>             Edge_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<2>             Face_range;

  typedef typename Combinatorial_map_with_cross_ratios::Dart_const_handle                               Dart_const_handle;
  typedef typename Combinatorial_map_with_cross_ratios::Dart_const_range                                Dart_const_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_const_range<1>       Edge_const_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_const_range<2>       Face_const_range;

  typedef typename Traits::FT                                                                Number;
  typedef typename Traits::Complex                                                           ComplexNumber;
  typedef typename Traits::Hyperbolic_point_2                                                Point;
  typedef Hyperbolic_isometry_2<Traits>                                                      Isometry;
  typedef Hyperbolic_fundamental_domain_2<Traits>                                            Domain;

  Hyperbolic_surface_triangulation_2() {};
  Hyperbolic_surface_triangulation_2(const Hyperbolic_fundamental_domain_2<Traits>& domain);
  Hyperbolic_surface_triangulation_2(const Combinatorial_map_with_cross_ratios& cmap);
  Hyperbolic_surface_triangulation_2(const Combinatorial_map_with_cross_ratios& cmap, const Anchor& anchor);

  //Hyperbolic_surface_triangulation_2& operator=(Hyperbolic_surface_triangulation_2&& other);
  Hyperbolic_surface_triangulation_2& operator=(Hyperbolic_surface_triangulation_2 other);

  Combinatorial_map_with_cross_ratios& get_combinatorial_map_ref();
  bool has_anchor() const;
  Anchor& get_anchor_ref();

  void to_stream(std::ostream& s) const;
  void from_stream(std::istream& s);

  bool is_delaunay_flippable(Dart_handle dart) const;
  void flip(Dart_handle dart);
  bool is_delaunay() const;
  int make_delaunay();
  std::vector<std::tuple<Dart_const_handle,Point,Point,Point>> lift(bool center=true) const;

  bool is_valid() const;

private:
  Combinatorial_map_with_cross_ratios _combinatorial_map;
  bool _has_anchor = false;
  Anchor _anchor;

  Dart_handle ccw(Dart_handle dart);
  Dart_handle cw(Dart_handle dart);
  Dart_handle opposite(Dart_handle dart);
  Dart_const_handle const_ccw(Dart_const_handle dart) const;
  Dart_const_handle const_cw(Dart_const_handle dart) const;
  Dart_const_handle const_opposite(Dart_const_handle dart) const;

  ComplexNumber get_cross_ratio(Dart_const_handle dart) const;

  Dart_handle pick_edge_to_flip();

  void copy_from(const Combinatorial_map_with_cross_ratios& cmap);
  void copy_from(const Combinatorial_map_with_cross_ratios& cmap, const Anchor& anchor);

  // Returns the cross ratio of the points a,b,c,d
  ComplexNumber cross_ratio(const Point& a, const Point& b, const Point& c, const Point& d) const;
  // Returns the point d such that the cross ratio of a,b,c,d is cratio
  Point fourth_point_from_cross_ratio(const Point& a, const Point& b, const Point& c, const ComplexNumber& cratio) const;
};

template<class Traits> std::ostream& operator<<(std::ostream& s, Hyperbolic_surface_triangulation_2<Traits>& Hyperbolic_surface_triangulation_2);
template<class Traits> void operator>>(std::istream& s, Hyperbolic_surface_triangulation_2<Traits>& triangulation);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class Traits>
Hyperbolic_surface_triangulation_2<Traits>::Hyperbolic_surface_triangulation_2(const Domain& domain){
  // (Triangulates by adding an internal edge between domain.vertex(size-1) and the other vertices)
  _combinatorial_map.clear();
  int size = domain.size();

  // Make the triangles
  Dart_handle dart_of_triangle[size-2];
  for (int k=0; k<size-2; k++){
    dart_of_triangle[k] = _combinatorial_map.make_combinatorial_polygon(3);
  }

  // Sew the internal edges and set their cross ratios
  Dart_handle dart_1, dart_2;
  Point p0,p1,p2,p3;

  for (int k=1; k<size-2; k++){
    dart_1 = dart_of_triangle[k];
    dart_2 = cw(dart_of_triangle[k-1]);

    p0 = domain.vertex(size-1);
    p1 = domain.vertex(k-1);
    p2 = domain.vertex(k);
    p3 = domain.vertex(k+1);

    _combinatorial_map.template sew<2>(dart_1, dart_2);
    _combinatorial_map.template set_attribute<1>(dart_1, _combinatorial_map.template create_attribute<1>(cross_ratio(p0,p1,p2,p3)));
  }

  // Sew the boundary edges and set their cross ratios
  for (int k1=0; k1<size; k1++){
    int k2 = domain.paired_side(k1);

    p0 = domain.vertex((k1+1)%size);
    p2 = domain.vertex(k1);

    if (k1 == size-1){
      dart_1 = dart_of_triangle[0];
      p1 = domain.vertex(1);
    } else if (k1 == size-2) {
      dart_1 = cw(dart_of_triangle[size-3]);
      p1 = domain.vertex(size-3);
    } else {
      dart_1 = ccw(dart_of_triangle[k1]);
      p1 = domain.vertex(size-1);
    }

    if (k2 == size-1){
      dart_2 = dart_of_triangle[0];
      p3 = domain.vertex(1);
    } else if (k2 == size-2) {
      dart_2 = cw(dart_of_triangle[size-3]);
      p3 = domain.vertex(size-3);
    } else {
      dart_2 = ccw(dart_of_triangle[k2]);
      p3 = domain.vertex(size-1);
    }

    p3 = domain.side_pairing(k1).evaluate(p3);

    if (_combinatorial_map.template is_sewable<2>(dart_1, dart_2)){
      _combinatorial_map.template sew<2>(dart_1, dart_2);
      _combinatorial_map.template set_attribute<1>(dart_1, _combinatorial_map.template create_attribute<1>(cross_ratio(p0,p1,p2,p3)));
    }
  }

  // Set the anchor
  _anchor.dart = dart_of_triangle[0];
  _anchor.vertices[0] = domain.vertex(size-1);
  _anchor.vertices[1] = domain.vertex(0);
  _anchor.vertices[2] = domain.vertex(1);
  _has_anchor = true;
}

template<class Traits>
Hyperbolic_surface_triangulation_2<Traits>::Hyperbolic_surface_triangulation_2(const Combinatorial_map_with_cross_ratios& cmap){
  copy_from(cmap);
}

template<class Traits>
Hyperbolic_surface_triangulation_2<Traits>::Hyperbolic_surface_triangulation_2(const Combinatorial_map_with_cross_ratios& cmap, const Anchor& anchor){
  copy_from(cmap, anchor);
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
//Hyperbolic_surface_triangulation_2<Traits>& Hyperbolic_surface_triangulation_2<Traits>::operator=(Hyperbolic_surface_triangulation_2<Traits>&& other){
Hyperbolic_surface_triangulation_2<Traits>& Hyperbolic_surface_triangulation_2<Traits>::operator=(Hyperbolic_surface_triangulation_2<Traits> other){
  if (other.has_anchor()){
    copy_from(other.get_combinatorial_map_ref(), other.get_anchor_ref());
  }
  else {
    copy_from(other.get_combinatorial_map_ref());
  }
  return *this;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
typename Hyperbolic_surface_triangulation_2<Traits>::Combinatorial_map_with_cross_ratios& Hyperbolic_surface_triangulation_2<Traits>::get_combinatorial_map_ref(){
  return _combinatorial_map;
}

template<class Traits>
bool Hyperbolic_surface_triangulation_2<Traits>::has_anchor() const {
  return _has_anchor;
}

template<class Traits>
typename Hyperbolic_surface_triangulation_2<Traits>::Anchor& Hyperbolic_surface_triangulation_2<Traits>::get_anchor_ref(){
  return _anchor;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
bool Hyperbolic_surface_triangulation_2<Traits>::is_delaunay_flippable(Dart_handle dart) const{
  return ( get_cross_ratio(dart).imaginary_part()>Number(0) );
}

template<class Traits>
void Hyperbolic_surface_triangulation_2<Traits>::flip(Dart_handle dart){

  // First gather all the information needed

   Dart_handle a = opposite(dart); // Get a fresh handle
   Dart_handle b = ccw(a);
   Dart_handle c = cw(a);

   Dart_handle d = opposite(a);
   Dart_handle e = ccw(d);
   Dart_handle f = cw(d);

   ComplexNumber cross_ratio_AB = get_cross_ratio(e);
   ComplexNumber cross_ratio_BC = get_cross_ratio(f);
   ComplexNumber cross_ratio_CD = get_cross_ratio(b);
   ComplexNumber cross_ratio_DA = get_cross_ratio(c);
   ComplexNumber cross_ratio_AC = get_cross_ratio(a);

   // Modify the anchor

   if (_anchor.dart == a){
     _anchor.dart = e;
     _anchor.vertices[1] = Point(fourth_point_from_cross_ratio(_anchor.vertices[1], _anchor.vertices[2], _anchor.vertices[0], cross_ratio_AC));
   } else if (_anchor.dart == b){
     _anchor.vertices[2] = Point(fourth_point_from_cross_ratio(_anchor.vertices[0], _anchor.vertices[1], _anchor.vertices[2], cross_ratio_AC));
   } else if (_anchor.dart == c){
     _anchor.vertices[2] = Point(fourth_point_from_cross_ratio(_anchor.vertices[2], _anchor.vertices[0], _anchor.vertices[1], cross_ratio_AC));
   } else if (_anchor.dart == d){
     _anchor.dart = b;
     _anchor.vertices[1] = Point(fourth_point_from_cross_ratio(_anchor.vertices[1], _anchor.vertices[2], _anchor.vertices[0], cross_ratio_AC));
   } else if (_anchor.dart == e){
     _anchor.vertices[2] = Point(fourth_point_from_cross_ratio(_anchor.vertices[0], _anchor.vertices[1], _anchor.vertices[2], cross_ratio_AC));
   } else if (_anchor.dart == f){
     _anchor.vertices[2] = Point(fourth_point_from_cross_ratio(_anchor.vertices[2], _anchor.vertices[0], _anchor.vertices[1], cross_ratio_AC));
   }

   // Compute the new cross ratios

   ComplexNumber one (Number(1), Number(0));
   ComplexNumber cross_ratio_BD = (cross_ratio_AC) / ((cross_ratio_AC) - one) ;
   ComplexNumber cross_ratio_AB_2 = one - (one - (cross_ratio_AB)) * (cross_ratio_AC) ;
   ComplexNumber cross_ratio_BC_2 = one - (one - (cross_ratio_BC)) / (cross_ratio_BD) ;
   ComplexNumber cross_ratio_CD_2 = one - (one - (cross_ratio_CD)) * (cross_ratio_AC) ;
   ComplexNumber cross_ratio_DA_2 = one - (one - (cross_ratio_DA)) / (cross_ratio_BD) ;

   // Make the topological flip

   _combinatorial_map.template unlink_beta<1>(a);
   _combinatorial_map.template unlink_beta<1>(b);
   _combinatorial_map.template unlink_beta<1>(c);

   _combinatorial_map.template unlink_beta<1>(d);
   _combinatorial_map.template unlink_beta<1>(e);
   _combinatorial_map.template unlink_beta<1>(f);


   _combinatorial_map.template link_beta<1>(b, a);
   _combinatorial_map.template link_beta<1>(a, f);
   _combinatorial_map.template link_beta<1>(f, b);

   _combinatorial_map.template link_beta<1>(e, d);
   _combinatorial_map.template link_beta<1>(d, c);
   _combinatorial_map.template link_beta<1>(c, e);

   // And give the new cross ratios to the edges

   _combinatorial_map.template info<1>(a) = cross_ratio_BD;
   _combinatorial_map.template info<1>(e) = cross_ratio_AB_2;
   _combinatorial_map.template info<1>(f) = cross_ratio_BC_2;
   _combinatorial_map.template info<1>(b) = cross_ratio_CD_2;
   _combinatorial_map.template info<1>(c) = cross_ratio_DA_2;

   // Take care of the particular cases where we need to "flip again"

   if (opposite(e) == b){
     _combinatorial_map.template info<1>(e) = one - (one - cross_ratio_AB_2) * (cross_ratio_AC) ;
   }

   if (opposite(f) == c){
     _combinatorial_map.template info<1>(f) = one - (one - cross_ratio_BC_2) / (cross_ratio_BD) ;
   }
}

template<class Traits>
bool Hyperbolic_surface_triangulation_2<Traits>::is_delaunay() const{
  if (! is_valid()){
    return false;
  }
  return (pick_edge_to_flip() == nullptr);
}

template<class Traits>
int Hyperbolic_surface_triangulation_2<Traits>::make_delaunay(){
  int number_of_flips_done = 0;

  Dart_handle edge_to_flip = pick_edge_to_flip();
  while (edge_to_flip != nullptr){
    flip(edge_to_flip);
    edge_to_flip = pick_edge_to_flip();
    number_of_flips_done++;
  }

  return number_of_flips_done;
}


template<class Traits>
std::vector<std::tuple<typename Hyperbolic_surface_triangulation_2<Traits>::Dart_const_handle, typename Hyperbolic_surface_triangulation_2<Traits>::Point, typename Hyperbolic_surface_triangulation_2<Traits>::Point, typename Hyperbolic_surface_triangulation_2<Traits>::Point>> Hyperbolic_surface_triangulation_2<Traits>::lift(bool center) const{
  std::vector<std::tuple<Dart_const_handle,Point,Point,Point>> realizations;

  size_t visited_darts_mark = _combinatorial_map.get_new_mark();
  _combinatorial_map.unmark_all(visited_darts_mark);

  struct Compare {
    bool operator()(std::pair<Dart_const_handle,double> const & x, std::pair<Dart_const_handle,double> const & y) {
        return x.second > y.second;
    }
  };
  std::priority_queue<std::pair<Dart_const_handle,double>, std::vector<std::pair<Dart_const_handle,double>>, Compare> queue;

  std::map<Dart_const_handle, Point> positions;

  Dart_const_range darts = _combinatorial_map.darts();

  _combinatorial_map.mark(_anchor.dart, visited_darts_mark);
  _combinatorial_map.mark(const_ccw(_anchor.dart), visited_darts_mark);
  _combinatorial_map.mark(const_cw(_anchor.dart), visited_darts_mark);

  if (center){
    Isometry center_the_drawing = hyperbolic_translation<Traits>(_anchor.vertices[0]);
    positions[_anchor.dart] = center_the_drawing.evaluate(_anchor.vertices[0]);
    positions[const_ccw(_anchor.dart)] = center_the_drawing.evaluate(_anchor.vertices[1]);
    positions[const_cw(_anchor.dart)] = center_the_drawing.evaluate(_anchor.vertices[2]);
  } else {
    positions[_anchor.dart] = _anchor.vertices[0];
    positions[const_ccw(_anchor.dart)] = _anchor.vertices[1];
    positions[const_cw(_anchor.dart)] = _anchor.vertices[2];
  }

  std::tuple<Dart_const_handle,Point,Point,Point> value = std::make_tuple(_anchor.dart, positions[_anchor.dart], positions[const_ccw(_anchor.dart)], positions[const_cw(_anchor.dart)]);
  realizations.push_back(value);

  ComplexNumber anchor_z0 (_anchor.vertices[0].x(), _anchor.vertices[0].y());
  ComplexNumber anchor_z1 (_anchor.vertices[1].x(), _anchor.vertices[1].y());
  ComplexNumber anchor_z2 (_anchor.vertices[2].x(), _anchor.vertices[2].y());

  double weight_of_anchor_dart = CGAL::to_double(anchor_z0.squared_modulus() + anchor_z1.squared_modulus());
  double weight_of_ccw_anchor_dart = CGAL::to_double(anchor_z1.squared_modulus() + anchor_z2.squared_modulus());
  double weight_of_cw_anchor_dart = CGAL::to_double(anchor_z2.squared_modulus() + anchor_z0.squared_modulus());

  queue.push(std::make_pair(_anchor.dart, weight_of_anchor_dart));
  queue.push(std::make_pair(const_ccw(_anchor.dart), weight_of_ccw_anchor_dart));
  queue.push(std::make_pair(const_cw(_anchor.dart), weight_of_cw_anchor_dart));



  while( ! queue.empty() ){
    Dart_const_handle invader = queue.top().first;
    queue.pop();

    Dart_const_handle invaded = const_opposite(invader);

    if (!_combinatorial_map.is_marked(invaded, visited_darts_mark)){
      _combinatorial_map.mark(invaded, visited_darts_mark);
      _combinatorial_map.mark(const_ccw(invaded), visited_darts_mark);
      _combinatorial_map.mark(const_cw(invaded), visited_darts_mark);

      Point a = positions[const_ccw(invader)];
      Point b = positions[const_cw(invader)];
      Point c = positions[invader];
      ComplexNumber cross_ratio = get_cross_ratio(invader);

      positions[invaded] = a;
      positions[const_ccw(invaded)] = c;
      Point d = fourth_point_from_cross_ratio(a, b, c, cross_ratio);
      positions[const_cw(invaded)] = d;

      ComplexNumber za (a.x(), a.y());
      ComplexNumber zc (c.x(), c.y());
      double invaded_distance_to_zero = CGAL::to_double(za.squared_modulus());
      double invaded_ccw_distance_to_zero = CGAL::to_double(zc.squared_modulus());
      ComplexNumber znew (positions[const_cw(invaded)].x(), positions[const_cw(invaded)].y());
      double invaded_cw_distance_to_zero = CGAL::to_double(znew.squared_modulus());

      double invaded_ccw_weight = invaded_ccw_distance_to_zero + invaded_cw_distance_to_zero;
      double invaded_cw_weight = invaded_cw_distance_to_zero + invaded_distance_to_zero;

      queue.push( std::make_pair(const_ccw(invaded), invaded_ccw_weight) );
      queue.push( std::make_pair(const_cw(invaded), invaded_cw_weight) );

      value = std::make_tuple(invaded, Point(a), Point(c), Point(d));
      realizations.push_back(value);
    }
  }

  _combinatorial_map.free_mark(visited_darts_mark);

  return realizations;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
bool Hyperbolic_surface_triangulation_2<Traits>::is_valid() const{
  // 1. Check the combinatorial map

  // Check that the combinatorial map is valid
  if ( !_combinatorial_map.is_valid() ){
    return false;
  }

  // Check that the combinatorial map has no 1,2-boundary
  for (int k=1; k<3; k++){
    if ( !_combinatorial_map.is_without_boundary(k) ){
      return false;
    }
  }

  // 2. Check the anchor, if any

  if  (_has_anchor){
    // Check that the dart handle of the anchor points to a dart of the combinatorial map
    if ( !_combinatorial_map.is_dart_used(_anchor.dart) ){
      return false;
    }

    // Check that the three vertices of the anchor lie within the open unit disk
    for (int k=0; k<3; k++){
      if (_anchor.vertices[k].get_z().squared_modulus() >= Number(1)){
        return false;
      }
    }
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
void Hyperbolic_surface_triangulation_2<Traits>::to_stream(std::ostream& s) const{
    // Give indices to the darts
    std::map<Dart_const_handle, int> darts_indices;
    int current_dart_index = 0;
    for (typename Dart_const_range::const_iterator it=_combinatorial_map.darts().begin(); it!=_combinatorial_map.darts().end(); ++it){
      darts_indices[it] = current_dart_index;
      current_dart_index++;
    }

    // Store the number of darts
    s << current_dart_index << std::endl;

    // Store the anchor, if any
    if (_has_anchor){
      s << "yes" << std::endl;
      s << darts_indices[_anchor.dart] << std::endl;
      s << _anchor.vertices[0] << std::endl;
      s << _anchor.vertices[1] << std::endl;
      s << _anchor.vertices[2] << std::endl;
    } else {
      s << "no" << std::endl;
    }

    // Store the triangles
    for (typename Face_const_range::const_iterator it = _combinatorial_map.template one_dart_per_cell<2>().begin(); it != _combinatorial_map.template one_dart_per_cell<2>().end(); ++it){
      s << darts_indices[it] << std::endl;
      s << darts_indices[const_cw(it)] << std::endl;
      s << darts_indices[const_ccw(it)] << std::endl;
    }

    // Store the edges
    for (typename Edge_range::const_iterator it = _combinatorial_map.template one_dart_per_cell<1>().begin(); it != _combinatorial_map.template one_dart_per_cell<1>().end(); ++it){
      s << darts_indices[it] << std::endl;
      s << darts_indices[const_opposite(it)] << std::endl;
      s << get_cross_ratio(it);
    }
}

template<class Traits>
void Hyperbolic_surface_triangulation_2<Traits>::from_stream(std::istream& s){
  _combinatorial_map.clear();

  // Load the number of darts
  std::string line;
  s >> line;
  int nb_darts = std::stoi(line);

  // Load the anchor
  int anchor_dart_id;
  s >> line;
  if (!line.compare("yes")){
    _has_anchor = true;

    s >> line;
    anchor_dart_id = std::stoi(line); // (*) _anchor.dart_id is set at the end of the function

    s >> _anchor.vertices[0];
    s >> _anchor.vertices[1];
    s >> _anchor.vertices[2];
  } else {
    _has_anchor = false;
  }

  // Load the triangles
  std::vector<Dart_handle> darts_by_id (nb_darts);
  int index1, index2, index3;
  for (int k=0; k<nb_darts/3; k++){
    Dart_handle triangle_dart = _combinatorial_map.make_combinatorial_polygon(3);

    s >> line;
    index1 = std::stoi(line);
    s >> line;
    index2 = std::stoi(line);
    s >> line;
    index3 = std::stoi(line);

    darts_by_id[index1] = triangle_dart;
    darts_by_id[index2] = cw(triangle_dart);
    darts_by_id[index3] = ccw(triangle_dart);
  }

  // Load the edges
  Dart_handle dart_1, dart_2;
  ComplexNumber cross_ratio;
  for (int k=0; k<nb_darts/2; k++){
    s >> line;
    index1 = std::stoi(line);
    s >> line;
    index2 = std::stoi(line);
    dart_1 = darts_by_id[index1];
    dart_2 = darts_by_id[index2];
    _combinatorial_map.template sew<2>(dart_1, dart_2);
    s >> cross_ratio;
    _combinatorial_map.template set_attribute<1>(dart_1, _combinatorial_map.template create_attribute<1>(cross_ratio));
  }

  // (*) here
  if (_has_anchor){
    _anchor.dart = darts_by_id[anchor_dart_id];
  }
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
std::ostream& operator<<(std::ostream& s, Hyperbolic_surface_triangulation_2<Traits>& triangulation){
  triangulation.to_stream(s);
  return s;
}

template<class Traits>
void operator>>(std::istream& s, Hyperbolic_surface_triangulation_2<Traits>& triangulation){
  triangulation.from_stream(s);
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
typename Hyperbolic_surface_triangulation_2<Traits>::Dart_handle Hyperbolic_surface_triangulation_2<Traits>::ccw(Dart_handle dart){
  return _combinatorial_map.beta(dart, 1);
}

template<class Traits>
typename Hyperbolic_surface_triangulation_2<Traits>::Dart_handle Hyperbolic_surface_triangulation_2<Traits>::cw(Dart_handle dart){
  return _combinatorial_map.beta(dart, 0);
}

template<class Traits>
typename Hyperbolic_surface_triangulation_2<Traits>::Dart_handle Hyperbolic_surface_triangulation_2<Traits>::opposite(Dart_handle dart){
  return _combinatorial_map.opposite(dart);
}

template<class Traits>
typename Hyperbolic_surface_triangulation_2<Traits>::Dart_const_handle Hyperbolic_surface_triangulation_2<Traits>::const_ccw(Dart_const_handle dart) const{
  return _combinatorial_map.beta(dart, 1);
}

template<class Traits>
typename Hyperbolic_surface_triangulation_2<Traits>::Dart_const_handle Hyperbolic_surface_triangulation_2<Traits>::const_cw(Dart_const_handle dart)  const{
  return _combinatorial_map.beta(dart, 0);
}

template<class Traits>
typename Hyperbolic_surface_triangulation_2<Traits>::Dart_const_handle Hyperbolic_surface_triangulation_2<Traits>::const_opposite(Dart_const_handle dart) const{
  return _combinatorial_map.opposite(dart);
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
typename Hyperbolic_surface_triangulation_2<Traits>::ComplexNumber Hyperbolic_surface_triangulation_2<Traits>::get_cross_ratio(Dart_const_handle dart) const{
  return _combinatorial_map.template info_of_attribute<1>(_combinatorial_map.template attribute<1>(dart));
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
typename Hyperbolic_surface_triangulation_2<Traits>::Dart_handle Hyperbolic_surface_triangulation_2<Traits>::pick_edge_to_flip(){
  for (typename Dart_range::iterator it = _combinatorial_map.darts().begin(); it != _combinatorial_map.darts().end(); ++it){
    if ( is_delaunay_flippable(it) ){
      return it;
    }
  }
  return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
void Hyperbolic_surface_triangulation_2<Traits>::copy_from(const Combinatorial_map_with_cross_ratios& cmap){
  //_combinatorial_map.copy_from_const(cmap);
  _combinatorial_map.copy(cmap);
  _has_anchor = false;
}

template<class Traits>
void Hyperbolic_surface_triangulation_2<Traits>::copy_from(const Combinatorial_map_with_cross_ratios& cmap, const Anchor& anchor){
  // Because of the anchor, we must operate the copy ourself
  _combinatorial_map.clear();

  // Copy the triangles and fill the darts conversion table
  std::map<Dart_const_handle, Dart_handle> darts_table;
  for (typename Face_const_range::const_iterator it=cmap.template one_dart_per_cell<2>().begin(); it!=cmap.template one_dart_per_cell<2>().end(); ++it){
    Dart_handle new_dart = _combinatorial_map.make_combinatorial_polygon(3);
    darts_table[it] = new_dart;
    darts_table[cmap.beta(it,0)] = _combinatorial_map.beta(new_dart,0);
    darts_table[cmap.beta(it,1)] = _combinatorial_map.beta(new_dart,1);
  }

  // Sew the edges and set their cross-ratios
  for (typename Edge_const_range::const_iterator it=cmap.template one_dart_per_cell<1>().begin(); it!=cmap.template one_dart_per_cell<1>().end(); ++it){
    Dart_handle dart_1 = darts_table[it];
    Dart_handle dart_2 = darts_table[cmap.opposite(it)];
    ComplexNumber cratio = cmap.template info_of_attribute<1>(cmap.template attribute<1>(it));

    _combinatorial_map.template sew<2>(dart_1, dart_2);
    _combinatorial_map.template set_attribute<1>(dart_1, _combinatorial_map.template create_attribute<1>(cratio));
  }

  cmap.opposite(anchor.dart);

  // Set the anchor
  _anchor.dart = darts_table[anchor.dart];
  for (int k=0; k<3; k++){
    _anchor.vertices[k] = anchor.vertices[k];
  }
  _has_anchor = true;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
typename Traits::Complex Hyperbolic_surface_triangulation_2<Traits>::cross_ratio(const Point& a, const Point& b, const Point& c, const Point& d) const{
  ComplexNumber za (a.x(), a.y());
  ComplexNumber zb (b.x(), b.y());
  ComplexNumber zc (c.x(), c.y());
  ComplexNumber zd (d.x(), d.y());
  return (zd-zb)*(zc-za) / ((zd-za)*(zc-zb));
}

template<class Traits>
typename Hyperbolic_surface_triangulation_2<Traits>::Point Hyperbolic_surface_triangulation_2<Traits>::fourth_point_from_cross_ratio(const Point& a, const Point& b, const Point& c, const ComplexNumber& cratio) const{
  ComplexNumber za (a.x(), a.y());
  ComplexNumber zb (b.x(), b.y());
  ComplexNumber zc (c.x(), c.y());
  ComplexNumber result = ( cratio*za*(zc-zb) + zb*(za-zc) ) / ( cratio*(zc-zb) + (za-zc) );
  return Point(result.real_part(), result.imaginary_part());
}

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_SURFACE_TRIANGULATION_2
