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

#ifndef CGAL_HYPERBOLIC_FUNDAMENTAL_DOMAIN_2_H
#define CGAL_HYPERBOLIC_FUNDAMENTAL_DOMAIN_2_H

#include <CGAL/license/Triangulation_on_hyperbolic_surface_2.h>

#include <CGAL/Hyperbolic_surface_traits_2.h>
#include <CGAL/Hyperbolic_isometry_2.h>
#include <CGAL/assertions.h>
#include <iostream>
#include <vector>

namespace CGAL {

/*
Represents a convex geodesic hyperbolic domain D of a closed orientable hyperbolic surface.
The domain D is given as a convex geodesic hyperbolic polygon P given by the list of its vertices in the hyperbolic plane,
together with a pairing of the sides of P, such that every two paired sides have the same length, and such that
identifying every two paired sides in a way that respects the orientation of P would result in a closed
orientable hyperbolic surface.
*/
template<class Traits>
class Hyperbolic_fundamental_domain_2
{
public:
  typedef typename Traits::Hyperbolic_point_2                    Point;
  typedef typename Traits::FT                          FT;
  typedef Complex_number<FT>                                          Complex_number;

  Hyperbolic_fundamental_domain_2() {};

  //TODO DOC added precondition
  template<class PointRange, class PairingRange>
  Hyperbolic_fundamental_domain_2(PointRange & vertices, PairingRange & pairings)
  {
    vertices_ = std::vector<Point>(vertices.begin(), vertices.end());
    pairings_ = std::vector<std::size_t>(pairings.begin(), pairings.end());
    CGAL_precondition(is_valid());
    CGAL_precondition(is_valid_length_pairing());
  }

  // returns the number of vertices (equivalently, the number of sides)
  std::size_t size() const;

  // returns the index-th vertex
  const Point& vertex(std::size_t index) const;

  // returns the index of the side paired to side A, where A is the index-th side
  std::size_t paired_side(std::size_t index) const;

  // returns the isometry that maps side A to side B, where B is the index-th side, and A is the side paired to B
  Hyperbolic_isometry_2<Traits> side_pairing(std::size_t index) const;

  std::istream& from_stream(std::istream& s);
  std::ostream& to_stream(std::ostream& s) const;

  std::ostream& to_json(std::ostream& s) const;

  bool is_valid() const;
  bool is_valid_length_pairing() const;

private:
  std::vector<Point> vertices_;
  std::vector<std::size_t> pairings_;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class Traits>
std::size_t
Hyperbolic_fundamental_domain_2<Traits>::
size() const
{
  CGAL_precondition(is_valid());
  return vertices_.size();
}

template<class Traits>
const typename Hyperbolic_fundamental_domain_2<Traits>::Point&
Hyperbolic_fundamental_domain_2<Traits>::
vertex(std::size_t index) const
{
  CGAL_precondition(is_valid());
  return vertices_[index];
}

template<class Traits>
std::size_t
Hyperbolic_fundamental_domain_2<Traits>::
paired_side(std::size_t index) const
{
  CGAL_precondition(is_valid());
  return pairings_[index];
}

template<class Traits>
Hyperbolic_isometry_2<Traits>
Hyperbolic_fundamental_domain_2<Traits>::
side_pairing(std::size_t index) const
{
  CGAL_precondition(is_valid());
  std::size_t n = size();
  std::size_t paired_index = paired_side(index);

  //const Point& p1,p2,q1,q2;
  const Point& q1 = vertex(index);
  const Point& q2 = vertex((index+1)%n);
  const Point& p2 = vertex(paired_index);
  const Point& p1 = vertex((paired_index+1)%n);

  Hyperbolic_isometry_2<Traits> isom = isometry_pairing_the_sides<Traits>(p1,p2,q1,q2);
  return isom;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
std::ostream&
Hyperbolic_fundamental_domain_2<Traits>::
to_stream(std::ostream& s) const
{
  std::size_t n = size();

  s << std::to_string(n) << std::endl;

  for (std::size_t k=0; k<n; ++k) {
    s << paired_side(k) << std::endl;
  }

  for (std::size_t k=0; k<n; ++k) {
    s << vertex(k) << std::endl;
  }
  return s;
}

template<class Traits>
std::istream&
Hyperbolic_fundamental_domain_2<Traits>::
from_stream(std::istream& s)
{
  vertices_.clear();
  pairings_.clear();

  std::string line;
  s >> line;
  std::size_t size = std::stoi(line);
  vertices_.reserve(size);
  pairings_.reserve(size);

  for (std::size_t k=0; k<size; ++k) {
    s >> line;
    pairings_.push_back(std::stoi(line));
  }

  for (std::size_t k=0; k<size; ++k) {
    Point p;
    s >> p;
    vertices_.push_back(p);
  }

  CGAL_precondition(is_valid());
  CGAL_precondition(is_valid_length_pairing());

  return s;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
bool
Hyperbolic_fundamental_domain_2<Traits>::
is_valid()const
{
  // Get the number of vertices
  std::size_t n = vertices_.size();

  // Check that the number of vertices is even
  if (n%2) {
    return false;
  }

  // Check that there are as many side pairings as vertices
  if (pairings_.size() != n) {
    return false;
  }

  // Check that the pairings_ vector encodes a perfect matching of the set {0,1,\dots,n-1}
  std::vector<bool> already_paired(n);
  for (std::size_t k=0; k<n; ++k) {
    already_paired[k] = false;
  }
  for (std::size_t k=0; k<n; ++k) {
    std::size_t paired_side = pairings_[k];
    if (paired_side>=n) {
      return false;
    }
    if (already_paired[paired_side]) {
      return false;
    }
    already_paired[paired_side] = true;
  }

  // Check that the vertices all lie within the open unit disk
  for (std::size_t k=0; k<n; ++k) {
    if (norm(Complex_number(vertices_[k].x(), vertices_[k].y())) >= typename Traits::FT(1)) {
      return false;
    }
  }

  return true;
}

//TODO DOC
template<class Traits>
bool
Hyperbolic_fundamental_domain_2<Traits>::
is_valid_length_pairing()const
{
  // Get the number of vertices (= nb of sides)
  std::size_t n = vertices_.size();
  Point v1,v2,v1p,v2p;

  for (std::size_t k=0; k<n; ++k) {
    v1 = Point(vertices_[k%n].x(), vertices_[k%n].y());
    v2 = Point(vertices_[(k+1)%n].x(), vertices_[(k+1)%n].y());
    std::size_t kp = pairings_[k];
    v1p = Point(vertices_[kp%n].x(), vertices_[kp%n].y());
    v2p = Point(vertices_[(kp+1)%n].x(), vertices_[(kp+1)%n].y());

    //Try object design? typename Traits::cosh_hd chd = gt.cosh_hd_object();
    //if (!(chd(v1,v2) == chd(v1p,v2p))) {
    if (!(Traits::cosh_hd(v1,v2) == Traits::cosh_hd(v1p,v2p))) {
      return false;
    }
  }

 return true;
}

//////////////////////////////////////////////////////
//       TO JSON OUTPUT
//////////////////////////////////////////////////////
template<class Traits>
std::ostream&
Hyperbolic_fundamental_domain_2<Traits>::
to_json(std::ostream& s) const
{
    const std::size_t n = size();

    s << "{\n";
    s << "  \"type\": " << "\"domain\"" << ",\n";
    s << "  \"size\": " << n << ",\n";

    s << "  \"paired_side\": [";
    for (std::size_t k = 0; k < n; ++k)
    {
        if (k > 0) s << ", ";
        s << paired_side(k);
    }
    s << "],\n";

    s << "  \"vertices\": [";
    for (std::size_t k = 0; k < n; ++k)
    {
        if (k > 0) s << "," << std::endl;
        s << "[\"" << vertex(k).x() << "\", \"" << vertex(k).y() << "\"]" ;
    }
    s << "]\n";

    s << "}";

    return s;
}

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_FUNDAMENTAL_DOMAIN_2_H
