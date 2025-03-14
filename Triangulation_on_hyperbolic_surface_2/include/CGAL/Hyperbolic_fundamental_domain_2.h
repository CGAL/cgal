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

  Hyperbolic_fundamental_domain_2() {};

  template<class PointRange, class PairingRange>
  Hyperbolic_fundamental_domain_2(PointRange & vertices, PairingRange & pairings)
  {
    _vertices = std::vector<Point>(vertices.begin(), vertices.end());
    _pairings = std::vector<int>(pairings.begin(), pairings.end());
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

  bool is_valid() const;

private:
  std::vector<Point> _vertices;
  std::vector<int> _pairings;
};

//template<class Traits> std::ostream& operator<<(std::ostream& s, const Hyperbolic_fundamental_domain_2<Traits>& domain);
//template<class Traits> std::istream& operator>>(std::istream& s, Hyperbolic_fundamental_domain_2<Traits>& domain);

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
  return _vertices.size();
}

template<class Traits>
const typename Hyperbolic_fundamental_domain_2<Traits>::Point&
Hyperbolic_fundamental_domain_2<Traits>::
vertex(std::size_t index) const
{
  CGAL_precondition(is_valid());
  return _vertices[index];
}

template<class Traits>
std::size_t
Hyperbolic_fundamental_domain_2<Traits>::
paired_side(std::size_t index) const
{
  CGAL_precondition(is_valid());
  return _pairings[index];
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
  _vertices.clear();
  _pairings.clear();

  std::string line;
  s >> line;
  std::size_t size = std::stoi(line);
  _vertices.reserve(size);
  _pairings.reserve(size);

  for (std::size_t k=0; k<size; ++k) {
    s >> line;
    _pairings.push_back(std::stoi(line));
  }

  for (std::size_t k=0; k<size; ++k) {
    Point p;
    s >> p;
    _vertices.push_back(p);
  }
  return s;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
bool
Hyperbolic_fundamental_domain_2<Traits>::
is_valid()const
{
  // Get the number of vertices
  std::size_t n = _vertices.size();

  // Check that the number of vertices is even
  if (n%2) {
    return false;
  }

  // Check that there are as many side pairings as vertices
  if (_pairings.size() != n) {
    return false;
  }

  // Check that the _pairings vector encodes a perfect matching of the set {0,1,\dots,n-1}
  std::vector<bool> already_paired(n);
  for (std::size_t k=0; k<n; ++k) {
    already_paired[k] = false;
  }
  for (std::size_t k=0; k<n; ++k) {
    std::size_t paired_side = _pairings[k];
    if ((paired_side<0) || (paired_side>=n)) {
      return false;
    }
    if (already_paired[paired_side]) {
      return false;
    }
    already_paired[paired_side] = true;
  }

  // Check that the vertices all lie within the open unit disk
  for (std::size_t k=0; k<n; ++k) {
    if (norm(Complex_number(_vertices[k].x(), _vertices[k].y())) >= typename Traits::FT(1)) {
      return false;
    }
  }

  return true;
}

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_FUNDAMENTAL_DOMAIN_2_H
