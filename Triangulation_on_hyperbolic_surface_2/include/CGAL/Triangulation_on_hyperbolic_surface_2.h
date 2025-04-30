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

#ifndef CGAL_TRIANGULATION_ON_HYPERBOLIC_SURFACE_2_H
#define CGAL_TRIANGULATION_ON_HYPERBOLIC_SURFACE_2_H

#include <CGAL/license/Triangulation_on_hyperbolic_surface_2.h>

#include <CGAL/Combinatorial_map.h>
#include <CGAL/Hyperbolic_fundamental_domain_2.h>

#include <CGAL/assertions.h>

#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <utility>

#include <optional>

namespace CGAL {

/*
Represents a geodesic triangulation of a closed orientable hyperbolic surface.
The triangulation is stored as combinatorial map decorated with one cross-ratio per edge.
It is also possible to specify an anchor for the triangulation. An anchor consists in 1) a dart of the combinatorial map, belonging by definition to a vertex V and a triangle T, together with
2) three points A,B,C in the hyperbolic plane. The points A,B,C are the three vertices in counter-clockwise order of a triangle. This triangle is a lift
of T, and A is a lift of V.
*/
template<class Traits>
struct Combinatorial_map_with_cross_ratios_item
{
  template <class CMap>
  struct Dart_wrapper
  {
    typedef Cell_attribute<CMap, Complex_number<typename Traits::FT> >    Edge_attrib;
    typedef std::tuple<void, Edge_attrib, void>                           Attributes;
  };
};

template<class Traits, class Attributes = Combinatorial_map_with_cross_ratios_item<Traits> >
class Triangulation_on_hyperbolic_surface_2
{
public:
  typedef Combinatorial_map<2, Attributes> Combinatorial_map_with_cross_ratios;

  struct Anchor
  {
    typename Combinatorial_map_with_cross_ratios::Dart_descriptor dart;
    typename Traits::Hyperbolic_point_2 vertices[3];

    Anchor(){};

    Anchor(typename Combinatorial_map_with_cross_ratios::Dart_descriptor dart,
           typename Traits::Hyperbolic_point_2 a, typename Traits::Hyperbolic_point_2 b, typename Traits::Hyperbolic_point_2 c)
    {
      this->dart = dart;
      vertices[0] = a;
      vertices[1] = b;
      vertices[2] = c;
    }
  };

  typedef typename Combinatorial_map_with_cross_ratios::Dart_descriptor                           Dart_descriptor;
  // typedef typename Combinatorial_map_with_cross_ratios::Dart_range                               Dart_range;
  // typedef typename Combinatorial_map_with_cross_ratios::Dart_range::iterator                     Dart_iterator;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<0>       Vertex_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<1>       Edge_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_range<2>       Face_range;

  typedef typename Combinatorial_map_with_cross_ratios::Dart_const_descriptor                     Dart_const_descriptor;
  typedef typename Combinatorial_map_with_cross_ratios::Dart_const_range                          Dart_const_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_const_range<0> Vertex_const_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_const_range<1> Edge_const_range;
  typedef typename Combinatorial_map_with_cross_ratios::template One_dart_per_cell_const_range<2> Face_const_range;

  typedef typename Traits::FT                                                                     Number;
  typedef typename Traits::Complex                                                                Complex_number;
  typedef typename Traits::Hyperbolic_point_2                                                     Point;
  typedef Hyperbolic_isometry_2<Traits>                                                           Isometry;
  typedef Hyperbolic_fundamental_domain_2<Traits>                                                 Domain;

  Triangulation_on_hyperbolic_surface_2() {};
  Triangulation_on_hyperbolic_surface_2(const Hyperbolic_fundamental_domain_2<Traits>& domain);
// Triangulation_on_hyperbolic_surface_2(Combinatorial_map_with_cross_ratios& cmap);
  Triangulation_on_hyperbolic_surface_2(Combinatorial_map_with_cross_ratios& cmap, Anchor& anchor);

  Combinatorial_map_with_cross_ratios& combinatorial_map();
  bool has_anchor() const;
  Anchor& anchor();
  const Anchor& anchor() const;

  void to_stream(std::ostream& s) const;
  void from_stream(std::istream& s);

  bool is_Delaunay_flippable(Dart_const_descriptor dart) const;
  void flip(Dart_descriptor dart);
  bool is_Delaunay() const;
  int make_Delaunay();
  std::vector<std::tuple<Dart_const_descriptor, Point, Point, Point> > lift(bool center=true) const;

  bool is_valid() const;

  // The following methods are not documented but they are non private for internal future use.
  Dart_descriptor ccw(Dart_descriptor dart);
  Dart_descriptor cw(Dart_descriptor dart);
  Dart_descriptor opposite(Dart_descriptor dart);
  Dart_const_descriptor const_ccw(Dart_const_descriptor dart) const;
  Dart_const_descriptor const_cw(Dart_const_descriptor dart) const;
  Dart_const_descriptor const_opposite(Dart_const_descriptor dart) const;

  Complex_number get_cross_ratio(Dart_const_descriptor dart) const;

  // returns the cross ratio of the points a,b,c,d
  Complex_number cross_ratio(const Point& a, const Point& b, const Point& c, const Point& d) const;
  // returns the point d such that the cross ratio of a,b,c,d is cratio
  Point fourth_point_from_cross_ratio(const Point& a, const Point& b, const Point& c, const Complex_number& cratio) const;

// Wrapper around the Cmap for iterating over vertices, edges or faces.
Vertex_range vertices_range() {
  return combinatorial_map_.template one_dart_per_cell<0>();
}
Edge_range edges_range() {
  return combinatorial_map_.template one_dart_per_cell<1>();
}
Face_range faces_range() {  
  return combinatorial_map_.template one_dart_per_cell<2>();
}
Vertex_const_range vertices_const_range() const {
  return combinatorial_map_.template one_dart_per_cell<0>();
}
Edge_const_range edges_const_range() const {
  return combinatorial_map_.template one_dart_per_cell<1>();
}
Face_const_range faces_const_range() const {
  return combinatorial_map_.template one_dart_per_cell<2>();
}

protected:
  Combinatorial_map_with_cross_ratios combinatorial_map_;
  // bool has_anchor_ = false;
  std::optional<Anchor> anchor_ = std::nullopt;

  Dart_descriptor pick_edge_to_flip();
  Dart_const_descriptor pick_edge_to_flip() const;

  void copy_from(Combinatorial_map_with_cross_ratios& cmap);
  void copy_from(Combinatorial_map_with_cross_ratios& cmap, const Anchor& anchor);
};

// template<class Traits, class Attributes> std::ostream& operator<<(std::ostream& s, const Triangulation_on_hyperbolic_surface_2<Traits, Attributes>& triangulation);
// template<class Traits, class Attributes> void operator>>(std::istream& s, Triangulation_on_hyperbolic_surface_2<Traits, Attributes>& triangulation);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class Traits, class Attributes>
Triangulation_on_hyperbolic_surface_2<Traits,Attributes>::
Triangulation_on_hyperbolic_surface_2(const Domain& domain)
{
  // (Triangulates by adding an internal edge between domain.vertex(size-1) and the other vertices)
  combinatorial_map_.clear();
  std::size_t size = domain.size();

  // Make the triangles
  std::vector<Dart_descriptor> dart_of_triangle(size-2);
  for (std::size_t  k=0; k<size-2; ++k) {
    dart_of_triangle[k] = combinatorial_map_.make_combinatorial_polygon(3);
  }

  // Sew the internal edges and set their cross ratios
  Dart_descriptor dart_1, dart_2;
  Point p0, p1, p2, p3;

  for (std::size_t  k=1; k<size-2; ++k) {
    dart_1 = dart_of_triangle[k];
    dart_2 = cw(dart_of_triangle[k-1]);

    p0 = domain.vertex(size-1);
    p1 = domain.vertex(k-1);
    p2 = domain.vertex(k);
    p3 = domain.vertex(k+1);

    combinatorial_map_.template sew<2>(dart_1, dart_2);
    combinatorial_map_.template set_attribute<1>(dart_1, combinatorial_map_.template create_attribute<1>(cross_ratio(p0,p1,p2,p3)));
  }

  // Sew the boundary edges and set their cross ratios
  for (std::size_t  k1=0; k1<size; k1++) {
    std::size_t k2 = domain.paired_side(k1);

    p0 = domain.vertex((k1+1)%size);
    p2 = domain.vertex(k1);

    if (k1 == size-1) {
      dart_1 = dart_of_triangle[0];
      p1 = domain.vertex(1);
    } else if (k1 == size-2) {
      dart_1 = cw(dart_of_triangle[size-3]);
      p1 = domain.vertex(size-3);
    } else {
      dart_1 = ccw(dart_of_triangle[k1]);
      p1 = domain.vertex(size-1);
    }

    if (k2 == size-1) {
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

    if (combinatorial_map_.template is_sewable<2>(dart_1, dart_2)) {
      combinatorial_map_.template sew<2>(dart_1, dart_2);
      combinatorial_map_.template set_attribute<1>(dart_1, combinatorial_map_.template create_attribute<1>(cross_ratio(p0,p1,p2,p3)));
    }
  }

  // Set the anchor
  anchor_ = Anchor(dart_of_triangle[0], domain.vertex(size-1), domain.vertex(0), domain.vertex(1));
}

/* template<class Traits, class Attributes> */
/*   Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Triangulation_on_hyperbolic_surface_2(Combinatorial_map_with_cross_ratios& cmap) { */
/*   copy_from(cmap); */
/* } */

template<class Traits, class Attributes>
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
Triangulation_on_hyperbolic_surface_2(Combinatorial_map_with_cross_ratios& cmap,
                                      Anchor& anchor)
{
  copy_from(cmap, anchor);
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits, class Attributes>
typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Combinatorial_map_with_cross_ratios&
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
combinatorial_map()
{
  return combinatorial_map_;
}

template<class Traits, class Attributes>
bool
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
has_anchor() const
{
  CGAL_precondition(is_valid());
  // return has_anchor_;
  return anchor_.has_value();
}

template<class Traits, class Attributes>
typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Anchor&
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
anchor()
{
  CGAL_precondition(is_valid() && has_anchor());
  return anchor_.value();
}

template<class Traits, class Attributes>
const typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Anchor&
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
anchor() const
{
  CGAL_precondition(is_valid() && has_anchor());
  return anchor_.value();
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits, class Attributes>
bool
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
is_Delaunay_flippable(Dart_const_descriptor dart) const
{
  CGAL_precondition(is_valid());
  return (get_cross_ratio(dart).imag() > Number(0));
}

template<class Traits, class Attributes>
void
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
flip(Dart_descriptor dart)
{
  CGAL_precondition(is_valid());

  // First gather all the information needed

   Dart_descriptor a = opposite(dart); // Get a fresh descriptor
   Dart_descriptor b = ccw(a);
   Dart_descriptor c = cw(a);

   Dart_descriptor d = opposite(a);
   Dart_descriptor e = ccw(d);
   Dart_descriptor f = cw(d);

   Complex_number cross_ratio_AB = get_cross_ratio(e);
   Complex_number cross_ratio_BC = get_cross_ratio(f);
   Complex_number cross_ratio_CD = get_cross_ratio(b);
   Complex_number cross_ratio_DA = get_cross_ratio(c);
   Complex_number cross_ratio_AC = get_cross_ratio(a);

   // Modify the anchor
   if(has_anchor()){
     if (anchor_.value().dart == a) {
       anchor_.value().dart = e;
       anchor_.value().vertices[1] = Point(fourth_point_from_cross_ratio(anchor_.value().vertices[1], anchor_.value().vertices[2], anchor_.value().vertices[0], cross_ratio_AC));
     } else if (anchor_.value().dart == b) {
       anchor_.value().vertices[2] = Point(fourth_point_from_cross_ratio(anchor_.value().vertices[0], anchor_.value().vertices[1], anchor_.value().vertices[2], cross_ratio_AC));
     } else if (anchor_.value().dart == c) {
       anchor_.value().vertices[2] = Point(fourth_point_from_cross_ratio(anchor_.value().vertices[2], anchor_.value().vertices[0], anchor_.value().vertices[1], cross_ratio_AC));
     } else if (anchor_.value().dart == d) {
       anchor_.value().dart = b;
       anchor_.value().vertices[1] = Point(fourth_point_from_cross_ratio(anchor_.value().vertices[1], anchor_.value().vertices[2], anchor_.value().vertices[0], cross_ratio_AC));
     } else if (anchor_.value().dart == e) {
       anchor_.value().vertices[2] = Point(fourth_point_from_cross_ratio(anchor_.value().vertices[0], anchor_.value().vertices[1], anchor_.value().vertices[2], cross_ratio_AC));
     } else if (anchor_.value().dart == f) {
       anchor_.value().vertices[2] = Point(fourth_point_from_cross_ratio(anchor_.value().vertices[2], anchor_.value().vertices[0], anchor_.value().vertices[1], cross_ratio_AC));
     }
  }

   // Compute the new cross ratios

   Complex_number one (Number(1), Number(0));
   Complex_number cross_ratio_BD = (cross_ratio_AC) / ((cross_ratio_AC) - one);
   Complex_number cross_ratio_AB_2 = one - (one - (cross_ratio_AB)) * (cross_ratio_AC);
   Complex_number cross_ratio_BC_2 = one - (one - (cross_ratio_BC)) / (cross_ratio_BD);
   Complex_number cross_ratio_CD_2 = one - (one - (cross_ratio_CD)) * (cross_ratio_AC);
   Complex_number cross_ratio_DA_2 = one - (one - (cross_ratio_DA)) / (cross_ratio_BD);

   // Make the topological flip

   combinatorial_map_.template unlink_beta<1>(a);
   combinatorial_map_.template unlink_beta<1>(b);
   combinatorial_map_.template unlink_beta<1>(c);

   combinatorial_map_.template unlink_beta<1>(d);
   combinatorial_map_.template unlink_beta<1>(e);
   combinatorial_map_.template unlink_beta<1>(f);


   combinatorial_map_.template link_beta<1>(b, a);
   combinatorial_map_.template link_beta<1>(a, f);
   combinatorial_map_.template link_beta<1>(f, b);

   combinatorial_map_.template link_beta<1>(e, d);
   combinatorial_map_.template link_beta<1>(d, c);
   combinatorial_map_.template link_beta<1>(c, e);

   // And give the new cross ratios to the edges

   combinatorial_map_.template info<1>(a) = cross_ratio_BD;
   combinatorial_map_.template info<1>(e) = cross_ratio_AB_2;
   combinatorial_map_.template info<1>(f) = cross_ratio_BC_2;
   combinatorial_map_.template info<1>(b) = cross_ratio_CD_2;
   combinatorial_map_.template info<1>(c) = cross_ratio_DA_2;

   // Take care of the particular cases where we need to "flip again"

   if (opposite(e) == b) {
     combinatorial_map_.template info<1>(e) = one - (one - cross_ratio_AB_2) * (cross_ratio_AC);
   }

   if (opposite(f) == c) {
     combinatorial_map_.template info<1>(f) = one - (one - cross_ratio_BC_2) / (cross_ratio_BD);
   }
}

template<class Traits, class Attributes>
bool
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
is_Delaunay() const
{
  if (! is_valid()) {
    return false;
  }
  return (pick_edge_to_flip() == nullptr);
}

template<class Traits, class Attributes>
int
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
make_Delaunay()
{
  CGAL_precondition(is_valid());
  int number_of_flips_done = 0;

  Dart_descriptor edge_to_flip = pick_edge_to_flip();
  while (edge_to_flip != nullptr) {
    flip(edge_to_flip);
    edge_to_flip = pick_edge_to_flip();
    number_of_flips_done++;
  }

  return number_of_flips_done;
}


template<class Traits, class Attributes>
std::vector<std::tuple<typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Dart_const_descriptor,
                       typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Point,
                       typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Point,
                       typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Point> >
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
lift(bool center) const
{
  CGAL_precondition(is_valid() && has_anchor());

  std::vector<std::tuple<Dart_const_descriptor, Point, Point, Point> > realizations;

  size_t visited_darts_mark = combinatorial_map_.get_new_mark();
  combinatorial_map_.unmark_all(visited_darts_mark);

  struct Compare
  {
    bool operator()(const std::pair<Dart_const_descriptor,double>& x,
                    const std::pair<Dart_const_descriptor,double>& y)
    {
      return x.second > y.second;
    }
  };

  std::priority_queue<std::pair<Dart_const_descriptor, double>,
                      std::vector<std::pair<Dart_const_descriptor, double> >, Compare> queue;

  std::unordered_map<Dart_const_descriptor, Point> positions;

  Dart_const_range darts = combinatorial_map_.darts();

  combinatorial_map_.mark(anchor_.value().dart, visited_darts_mark);
  combinatorial_map_.mark(const_ccw(anchor_.value().dart), visited_darts_mark);
  combinatorial_map_.mark(const_cw(anchor_.value().dart), visited_darts_mark);

  if (center) {
    Isometry center_the_drawing = hyperbolic_translation<Traits>(anchor_.value().vertices[0]);
    positions[anchor_.value().dart] = center_the_drawing.evaluate(anchor_.value().vertices[0]);
    positions[const_ccw(anchor_.value().dart)] = center_the_drawing.evaluate(anchor_.value().vertices[1]);
    positions[const_cw(anchor_.value().dart)] = center_the_drawing.evaluate(anchor_.value().vertices[2]);
  } else {
    positions[anchor_.value().dart] = anchor_.value().vertices[0];
    positions[const_ccw(anchor_.value().dart)] = anchor_.value().vertices[1];
    positions[const_cw(anchor_.value().dart)] = anchor_.value().vertices[2];
  }

  std::tuple<Dart_const_descriptor, Point, Point, Point> value =
    std::make_tuple(anchor_.value().dart,
                    positions[anchor_.value().dart],
                    positions[const_ccw(anchor_.value().dart)],
                    positions[const_cw(anchor_.value().dart)]);
  realizations.push_back(value);

  Complex_number anchor_z0(anchor_.value().vertices[0].x(), anchor_.value().vertices[0].y());
  Complex_number anchor_z1(anchor_.value().vertices[1].x(), anchor_.value().vertices[1].y());
  Complex_number anchor_z2(anchor_.value().vertices[2].x(), anchor_.value().vertices[2].y());

  double weight_of_anchor_dart = CGAL::to_double(norm(anchor_z0) + norm(anchor_z1));
  double weight_of_ccw_anchor_dart = CGAL::to_double(norm(anchor_z1) + norm(anchor_z2));
  double weight_of_cw_anchor_dart = CGAL::to_double(norm(anchor_z2) + norm(anchor_z0));

  queue.push(std::make_pair(anchor_.value().dart, weight_of_anchor_dart));
  queue.push(std::make_pair(const_ccw(anchor_.value().dart), weight_of_ccw_anchor_dart));
  queue.push(std::make_pair(const_cw(anchor_.value().dart), weight_of_cw_anchor_dart));

  while (! queue.empty()) {
    Dart_const_descriptor invader = queue.top().first;
    queue.pop();

    Dart_const_descriptor invaded = const_opposite(invader);

    if (!combinatorial_map_.is_marked(invaded, visited_darts_mark)) {
      combinatorial_map_.mark(invaded, visited_darts_mark);
      combinatorial_map_.mark(const_ccw(invaded), visited_darts_mark);
      combinatorial_map_.mark(const_cw(invaded), visited_darts_mark);

      const Point& a = positions[const_ccw(invader)];
      const Point& b = positions[const_cw(invader)];
      const Point& c = positions[invader];
      Complex_number cross_ratio = get_cross_ratio(invader);

      positions[invaded] = a;
      positions[const_ccw(invaded)] = c;
      Point d = fourth_point_from_cross_ratio(a, b, c, cross_ratio);
      positions[const_cw(invaded)] = d;

      Complex_number za(a.x(), a.y());
      Complex_number zc(c.x(), c.y());
      double invaded_distance_to_zero = CGAL::to_double(norm(za));
      double invaded_ccw_distance_to_zero = CGAL::to_double(norm(zc));
      Complex_number znew(positions[const_cw(invaded)].x(), positions[const_cw(invaded)].y());
      double invaded_cw_distance_to_zero = CGAL::to_double(norm(znew));

      double invaded_ccw_weight = invaded_ccw_distance_to_zero + invaded_cw_distance_to_zero;
      double invaded_cw_weight = invaded_cw_distance_to_zero + invaded_distance_to_zero;

      queue.push(std::make_pair(const_ccw(invaded), invaded_ccw_weight));
      queue.push(std::make_pair(const_cw(invaded), invaded_cw_weight));

      value = std::make_tuple(invaded, Point(a), Point(c), Point(d));
      realizations.push_back(value);
    }
  }

  combinatorial_map_.free_mark(visited_darts_mark);

  return realizations;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits, class Attributes>
bool
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
is_valid() const
{
  // 1. Check the combinatorial map

  // Check that the combinatorial map is valid
  if (!combinatorial_map_.is_valid()) {
    return false;
  }

  // Check that the combinatorial map has no 1,2-boundary
  for (int k=1; k<3; ++k) {
    if (!combinatorial_map_.is_without_boundary(k)) {
      return false;
    }
  }

  // 2. Check the anchor, if any

  if (has_anchor()) {
    // Check that the dart descriptor of the anchor points to a dart of the combinatorial map
    if (!combinatorial_map_.is_dart_used(anchor_.value().dart)) {
      return false;
    }

    // Check that the three vertices of the anchor lie within the open unit disk
    for (int k=0; k<3; ++k) {
      // if (anchor_.vertices[k].get_z() >= Number(1)) {
      if (norm(Complex_number(anchor_.value().vertices[k].x(),anchor_.value().vertices[k].y())) >= Number(1)) {
        return false;
      }
    }
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits, class Attributes>
void
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
to_stream(std::ostream& s) const
{
  CGAL_precondition(is_valid() && has_anchor());

  // Give indices to the darts
  std::map<Dart_const_descriptor, int> darts_indices;
  int current_dart_index = 0;
  for (typename Dart_const_range::const_iterator it=combinatorial_map_.darts().begin(); it!=combinatorial_map_.darts().end(); ++it) {
    darts_indices[it] = current_dart_index;
    current_dart_index++;
  }

  // Store the number of darts
  s << current_dart_index << std::endl;

  // Store the anchor, if any
  if (has_anchor()) {
    s << "yes" << std::endl;
    s << darts_indices[anchor_.value().dart] << std::endl;
    s << anchor_.value().vertices[0] << std::endl;
    s << anchor_.value().vertices[1] << std::endl;
    s << anchor_.value().vertices[2] << std::endl;
  } else {
    s << "no" << std::endl;
  }

  // Store the triangles
  for (typename Face_const_range::const_iterator it = faces_const_range().begin(); it != faces_const_range().end(); ++it) {
    s << darts_indices[it] << std::endl;
    s << darts_indices[const_cw(it)] << std::endl;
    s << darts_indices[const_ccw(it)] << std::endl;
  }

  // Store the edges
  for (typename Edge_const_range::const_iterator it = edges_const_range().begin(); it != edges_const_range().end(); ++it) {
    s << darts_indices[it] << std::endl;
    s << darts_indices[const_opposite(it)] << std::endl;
    s << get_cross_ratio(it);
  }
}

template<class Traits, class Attributes>
void
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
from_stream(std::istream& s)
{
  combinatorial_map_.clear();

  // Load the number of darts
  std::string line;
  s >> line;
  int nb_darts = std::stoi(line);

  // Load the anchor
  int anchor_dart_id = 0;
  s >> line;
  if (!line.compare("yes")) {
    anchor_ = Anchor();

    s >> line;
    anchor_dart_id = std::stoi(line); // (*) anchor_.dart_id is set at the end of the function

    s >> anchor_.value().vertices[0];
    s >> anchor_.value().vertices[1];
    s >> anchor_.value().vertices[2];
  } else {
    // has_anchor_ = false;
  }

  // Load the triangles
  std::vector<Dart_descriptor> darts_by_id (nb_darts);
  int index1, index2, index3;
  for (int k=0; k<nb_darts/3; ++k) {
    Dart_descriptor triangle_dart = combinatorial_map_.make_combinatorial_polygon(3);

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
  Dart_descriptor dart_1, dart_2;
  Complex_number cross_ratio;
  for (int k=0; k<nb_darts/2; ++k) {
    s >> line;
    index1 = std::stoi(line);
    s >> line;
    index2 = std::stoi(line);
    dart_1 = darts_by_id[index1];
    dart_2 = darts_by_id[index2];
    combinatorial_map_.template sew<2>(dart_1, dart_2);
    s >> cross_ratio;
    combinatorial_map_.template set_attribute<1>(dart_1, combinatorial_map_.template create_attribute<1>(cross_ratio));
  }

  // (*) here
  if (has_anchor()) {
    anchor_.value().dart = darts_by_id[anchor_dart_id];
  }
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits, class Attributes>
typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Dart_descriptor Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::ccw(Dart_descriptor dart) {
  return combinatorial_map_.beta(dart, 1);
}

template<class Traits, class Attributes>
typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Dart_descriptor Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::cw(Dart_descriptor dart) {
  return combinatorial_map_.beta(dart, 0);
}

template<class Traits, class Attributes>
typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Dart_descriptor Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::opposite(Dart_descriptor dart) {
  return combinatorial_map_.opposite(dart);
}

template<class Traits, class Attributes>
typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Dart_const_descriptor Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::const_ccw(Dart_const_descriptor dart) const {
  return combinatorial_map_.beta(dart, 1);
}

template<class Traits, class Attributes>
typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Dart_const_descriptor Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::const_cw(Dart_const_descriptor dart) const {
  return combinatorial_map_.beta(dart, 0);
}

template<class Traits, class Attributes>
typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Dart_const_descriptor Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::const_opposite(Dart_const_descriptor dart) const {
  return combinatorial_map_.opposite(dart);
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits, class Attributes>
typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Complex_number Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::get_cross_ratio(Dart_const_descriptor dart) const {
  return combinatorial_map_.template info_of_attribute<1>(combinatorial_map_.template attribute<1>(dart));
}

////////////////////////////////////////////////////////////////////////////////
template<class Traits, class Attributes>
typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Dart_descriptor Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
pick_edge_to_flip()
{
  auto& cm = combinatorial_map_.darts();
  for (auto it=cm.begin(); it!=cm.end(); ++it) {
    if (is_Delaunay_flippable(it)) {
      return it;
    }
  }
  return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits, class Attributes>
typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Dart_const_descriptor
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
pick_edge_to_flip() const
{
  const auto& cm = combinatorial_map_.darts();
  for (auto it=cm.begin(); it!=cm.end(); ++it) {
    if (is_Delaunay_flippable(it) ) {
      return it;
    }
  }
  return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits, class Attributes>
void
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
copy_from(Combinatorial_map_with_cross_ratios& cmap)
{
  //combinatorial_map_.copy_from_const(cmap);
  combinatorial_map_.copy(cmap);
  // has_anchor_ = false;
}

template<class Traits, class Attributes>
void
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
copy_from(Combinatorial_map_with_cross_ratios& cmap,
          const Anchor& anchor)
{
  // Because of the anchor, we must operate the copy ourself
  combinatorial_map_.clear();

  // Copy the triangles and fill the darts conversion table
  std::map<Dart_const_descriptor, Dart_descriptor> darts_table;
  for (typename Face_const_range::const_iterator it=cmap.template one_dart_per_cell<2>().begin(); it!=cmap.template one_dart_per_cell<2>().end(); ++it) {
    Dart_descriptor new_dart = combinatorial_map_.make_combinatorial_polygon(3);
    darts_table[it] = new_dart;
    darts_table[cmap.beta(it,0)] = combinatorial_map_.beta(new_dart,0);
    darts_table[cmap.beta(it,1)] = combinatorial_map_.beta(new_dart,1);
  }

  // Sew the edges and set their cross-ratios
  for (typename Edge_const_range::const_iterator it=cmap.template one_dart_per_cell<1>().begin(); it!=cmap.template one_dart_per_cell<1>().end(); ++it) {
    Dart_descriptor dart_1 = darts_table[it];
    Dart_descriptor dart_2 = darts_table[cmap.opposite(it)];
    Complex_number cratio = cmap.template info_of_attribute<1>(cmap.template attribute<1>(it));

    combinatorial_map_.template sew<2>(dart_1, dart_2);
    combinatorial_map_.template set_attribute<1>(dart_1, combinatorial_map_.template create_attribute<1>(cratio));
  }

  cmap.opposite(anchor.dart);

  // Set the anchor
  anchor_ = Anchor();
  anchor_.value().dart = darts_table[anchor.dart];
  for (int k=0; k<3; ++k) {
    anchor_.value().vertices[k] = anchor.vertices[k];
  }
  // has_anchor_ = true;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits, class Attributes>
typename Traits::Complex
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
cross_ratio(const Point& a, const Point& b, const Point& c, const Point& d) const
{
  Complex_number za (a.x(), a.y());
  Complex_number zb (b.x(), b.y());
  Complex_number zc (c.x(), c.y());
  Complex_number zd (d.x(), d.y());
  return (zd-zb)*(zc-za) / ((zd-za)*(zc-zb));
}

template<class Traits, class Attributes>
typename Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::Point
Triangulation_on_hyperbolic_surface_2<Traits, Attributes>::
fourth_point_from_cross_ratio(const Point& a, const Point& b, const Point& c,
                              const Complex_number& cratio) const
{
  Complex_number za (a.x(), a.y());
  Complex_number zb (b.x(), b.y());
  Complex_number zc (c.x(), c.y());
  Complex_number result = ( cratio*za*(zc-zb) + zb*(za-zc) ) / ( cratio*(zc-zb) + (za-zc));
  return Point(result.real(), result.imag());
}

} // namespace CGAL

#endif // CGAL_TRIANGULATION_ON_HYPERBOLIC_SURFACE_2_H

