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

// This file contains the declaration and the implementation of the class Hyperbolic_fundamental_domain_2

#ifndef CGAL_HYPERBOLIC_FUNDAMENTAL_DOMAIN_2
#define CGAL_HYPERBOLIC_FUNDAMENTAL_DOMAIN_2

#include <CGAL/Complex_without_sqrt.h>
#include "Hyperbolic_isometry_2.h"

#include <vector>
#include <iostream>

namespace CGAL {

/*
Represents a convex geodesic hyperbolic domain D of a closed orientable hyperbolic surface.
The domain D is given as a convex geodesic hyperbolic polygon P given by the list of its vertices in the hyperbolic plane,
together with a pairing of the sides of P, such that every two paired sides have the same length, and such that
identifying every two paired sides in a way that respects the orientation of P would result in a closed
orientable hyperbolic surface.
*/
template<class Traits>
class Hyperbolic_fundamental_domain_2 {

public:
  typedef typename Traits::Hyperbolic_point_2                    Point;

  Hyperbolic_fundamental_domain_2();

  Hyperbolic_fundamental_domain_2(typename std::vector<Point>::iterator vfirst,
				  typename std::vector<Point>::iterator vlast,
				  typename std::vector<int>::iterator pfirst,
				  typename std::vector<int>::iterator plast);
  
  int size() const; // Returns the number of vertices (equivalently, the number of sides)
  Point vertex(int index) const; // Returns the index-th vertex
  int paired_side(int index) const; // Returns the index of the side paired to side A, where A is the index-th side
  Hyperbolic_isometry_2<Traits> side_pairing(int index) const;// Returns the isometry that maps side A to side B, where B is the index-th side, and A is the side paired to B

  void from_stream(std::istream& s);
  void to_stream(std::ostream& s) const;

  bool is_valid() const;

private:
  std::vector<Point>                                               _vertices;
  std::vector<int>                                                  _pairings;
};

template<class Traits> std::ostream& operator<<(std::ostream& s, const Hyperbolic_fundamental_domain_2<Traits>& domain);
template<class Traits> void operator>>(std::istream& s, Hyperbolic_fundamental_domain_2<Traits>& domain);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class Traits>
  Hyperbolic_fundamental_domain_2<Traits>::Hyperbolic_fundamental_domain_2(){}

template<class Traits>
  Hyperbolic_fundamental_domain_2<Traits>::Hyperbolic_fundamental_domain_2(typename std::vector<Point>::iterator vfirst,
					       typename std::vector<Point>::iterator vlast,
					       typename std::vector<int>::iterator pfirst,
					       typename std::vector<int>::iterator plast){
  _vertices = std::vector<Point>(vfirst, vlast);
  _pairings = std::vector<int>(pfirst, plast);
 }

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
int Hyperbolic_fundamental_domain_2<Traits>::size() const{
  return _vertices.size();
}

template<class Traits>
typename Hyperbolic_fundamental_domain_2<Traits>::Point Hyperbolic_fundamental_domain_2<Traits>::vertex(int index) const{
  return _vertices[index];
}

template<class Traits>
int Hyperbolic_fundamental_domain_2<Traits>::paired_side(int index) const{
  return _pairings[index];
}

template<class Traits>
Hyperbolic_isometry_2<Traits> Hyperbolic_fundamental_domain_2<Traits>::side_pairing(int index) const{
  int n = size();
  int paired_index = paired_side(index);

  Point p1,p2,q1,q2;
  q1 = vertex(index);
  q2 = vertex((index+1)%n);
  p2 = vertex(paired_index);
  p1 = vertex((paired_index+1)%n);

  Hyperbolic_isometry_2<Traits> isom = isometry_pairing_the_sides<Traits>(p1,p2,q1,q2);
  return isom;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
void Hyperbolic_fundamental_domain_2<Traits>::to_stream(std::ostream& s) const{
  int n = size();

  s << std::to_string(n) << std::endl;

  for (int k=0; k<n; k++){
    s << paired_side(k) << std::endl;
  }

  for (int k=0; k<n; k++){
    s << vertex(k) << std::endl;
  }
}

template<class Traits>
void Hyperbolic_fundamental_domain_2<Traits>::from_stream(std::istream& s){
  _vertices.clear();
  _pairings.clear();

  std::string line;
  s >> line;
  int size = std::stoi(line);
  for (int k=0; k<size; k++){
    s >> line;
    _pairings.push_back(std::stoi(line));
  }

  for (int k=0; k<size; k++){
    Point p;
    s >> p;
    _vertices.push_back(p);
  }
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
std::ostream& operator<<(std::ostream& s, const Hyperbolic_fundamental_domain_2<Traits>& domain){
  domain.to_stream(s);
  return s;
}

template<class Traits>
void operator>>(std::istream& s, Hyperbolic_fundamental_domain_2<Traits>& domain){
  domain.from_stream(s);
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
bool Hyperbolic_fundamental_domain_2<Traits>::is_valid()const{
  // Get the number of vertices
  int n = _vertices.size();

  // Check that the number of vertices is even
  if (n%2){
    return false;
  }

  // Check that there are as many side pairings as vertices
  if (_pairings.size() != n){
    return false;
  }

  // Check that the _pairings vector encodes a perfect matching of the set {0,1,\dots,n-1}
  std::vector<bool> already_paired(n);
  for (int k=0; k<n; k++){
    already_paired[k] = false;
  }
  for (int k=0; k<n; k++){
    int paired_side = _pairings[k];
    if ((paired_side<0) || (paired_side>=n)){
      return false;
    }
    if (already_paired[paired_side]){
      return false;
    }
    already_paired[paired_side] = true;
  }

  // Check that the vertices all lie within the open unit disk
  for (int k=0; k<n; k++){
    if (_vertices[k].get_z().squared_modulus() >= typename Traits::FT(1)){
      return false;
    }
  }

  return true;
}

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_FUNDAMENTAL_DOMAIN_2
