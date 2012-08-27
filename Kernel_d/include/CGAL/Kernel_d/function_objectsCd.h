// Copyright (c) 2000,2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Michael Seel, Kurt Mehlhorn

#ifndef CGAL_FUNCTION_OBJECTSCD_H
#define CGAL_FUNCTION_OBJECTSCD_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>

#undef CGAL_KD_TRACE
#undef CGAL_KD_TRACEN
#undef CGAL_KD_TRACEV
#define CGAL_KD_TRACE(t)  std::cerr << t
#define CGAL_KD_TRACEN(t) std::cerr << t << std::endl
#define CGAL_KD_TRACEV(t) std::cerr << #t << " = " << (t) << std::endl
 
namespace CGAL {

template <typename K>
class Compute_coordinateCd {
  typedef typename K::FT             FT;
  typedef typename K::Point_d        Point_d;
  public:
  typedef FT                         result_type;
  result_type 
    operator()(const Point_d& p, int i) const
  {
    return p.cartesian(i);
  }
};

template <typename K>
class Point_dimensionCd {
  typedef typename K::FT             FT;
  typedef typename K::Point_d        Point_d;
  public:
  typedef int                       result_type;
  result_type 
    operator()(const Point_d& p) const
  {
    return p.dimension();
  }
};

template <typename K>
class Less_coordinateCd {
  typedef typename K::FT             FT;
  typedef typename K::Point_d        Point_d;
  public:
  typedef bool                       result_type;
  result_type 
  operator()(const Point_d& p, const Point_d& q, int i) const
  {
    return p.cartesian(i)<q.cartesian(i);
  }
};

template <class R>
struct Lift_to_paraboloidCd {
typedef typename R::Point_d Point_d;
typedef typename R::FT FT;
typedef typename R::LA LA;

Point_d operator()(const Point_d& p) const
{ 
  int d = p.dimension();
  typename LA::Vector h(d+1);
  FT sum = 0;
  for (int i = 0; i<d; i++) {
    h[i] = p.cartesian(i);
    sum += h[i]*h[i];
  }
  h[d] = sum;
  return Point_d(d+1,h.begin(),h.end());
}
};

template <class R>
struct Project_along_d_axisCd {
typedef typename R::Point_d Point_d;
typedef typename R::FT FT;

Point_d operator()(const Point_d& p) const
{ return Point_d(p.dimension()-1, 
                 p.cartesian_begin(),p.cartesian_end()-1); }
};

template <class R>
struct MidpointCd {
typedef typename R::Point_d Point_d;
Point_d operator()(const Point_d& p, const Point_d& q) const
{ return Point_d(p + (q-p)/2); }
};

template <class R>
struct Center_of_sphereCd {
typedef typename R::Point_d Point_d;
typedef typename R::FT FT;
typedef typename R::LA LA;
template <class Forward_iterator>
Point_d operator()(Forward_iterator start, Forward_iterator end) const
{ CGAL_assertion(start!=end);
  int d = start->dimension();
  typename LA::Matrix M(d);
  typename LA::Vector b(d);
  Point_d pd = *start++;
  for (int i = 0; i < d; i++) { 
    // we set up the equation for p_i
    Point_d pi = *start++;
    b[i] = 0;
    for (int j = 0; j < d; j++) {
      M(i,j) = FT(2)*(pi.cartesian(j) - pd.cartesian(j));
      b[i] += (pi.cartesian(j) - pd.cartesian(j)) *
              (pi.cartesian(j) + pd.cartesian(j));
    }
  }
  FT D;
  typename LA::Vector x;
  LA::linear_solver(M,b,x,D);
  return Point_d(d,x.begin(),x.end());
}

}; // Center_of_sphereCd


template <class R>
struct Squared_distanceCd {
typedef typename R::Point_d Point_d;
typedef typename R::Vector_d Vector_d;
typedef typename R::FT FT;
FT operator()(const Point_d& p, const Point_d& q) const
{ Vector_d v = p-q; return v.squared_length(); }
};

template <class R>
struct Position_on_lineCd {
typedef typename R::Point_d Point_d;
typedef typename R::LA LA;
typedef typename R::FT FT;

bool operator()(const Point_d& p, const Point_d& s, const Point_d& t, 
     FT& l) const
{ int d = p.dimension(); 
  CGAL_assertion_msg((d==s.dimension())&&(d==t.dimension()&& d>0), 
  "position_along_line: argument dimensions disagree.");
  CGAL_assertion_msg((s!=t), 
  "Position_on_line_d: line defining points are equal.");
  FT lnum = (p.cartesian(0) - s.cartesian(0)); 
  FT lden = (t.cartesian(0) - s.cartesian(0)); 
  FT num(lnum), den(lden), lnum_i, lden_i;
  for (int i = 1; i < d; i++) {  
    lnum_i = (p.cartesian(i) - s.cartesian(i)); 
    lden_i = (t.cartesian(i) - s.cartesian(i)); 
    if (lnum*lden_i != lnum_i*lden) return false; 
    if (lden_i != FT(0)) { den = lden_i; num = lnum_i; }
  }
  l = num/den; return true; 
}
};

template <class R>
struct Barycentric_coordinatesCd {
typedef typename R::Point_d Point_d;
typedef typename R::LA LA;
typedef typename R::FT FT;

template <class ForwardIterator, class OutputIterator>
OutputIterator operator()(ForwardIterator first, ForwardIterator last, 
  const Point_d& p, OutputIterator result)
{ TUPLE_DIM_CHECK(first,last,Barycentric_coordinates_d);
  //int n = std::distance(first,last); //unused variable
  int d = p.dimension();
  typename R::Affine_rank_d affine_rank;
  CGAL_assertion(affine_rank(first,last)==d);
  std::vector< Point_d > V(first,last);
  typename LA::Matrix M(d+1,V.size());
  typename LA::Vector b(d+1), x;
  int i;
  for (i=0; i<d; ++i) {
    for (int j=0; j<V.size(); ++j) 
      M(i,j)=V[j].cartesian(i);
    b[i] = p.cartesian(i);
  }
  for (int j=0; j<V.size(); ++j) 
    M(d,j) = 1;
  b[d]=1;
  FT D;
  LA::linear_solver(M,b,x,D);
  for (i=0; i < x.dimension(); ++result, ++i) {
    *result= x[i];
  }
  return result;
}
};


template <class R>
struct OrientationCd { 
typedef typename R::Point_d Point_d;
typedef typename R::LA LA;

template <class ForwardIterator>
Orientation operator()(ForwardIterator first, ForwardIterator last)
{ TUPLE_DIM_CHECK(first,last,Orientation_d);
  int d = static_cast<int>(std::distance(first,last)); 
  // range contains d points of dimension d-1
  CGAL_assertion_msg(first->dimension() == d-1,
  "Orientation_d: needs first->dimension() + 1 many points.");
  typename LA::Matrix M(d); // quadratic
  for (int i = 0; i < d; ++first,++i) {
    for (int j = 0; j < d-1; ++j) 
      M(i,j) = first->cartesian(j);
    M(i,d-1) = 1;
  }
  int row_correction = ( (d % 2 == 0) ? -1 : +1 );
  // we invert the sign if the row number is even i.e. d is odd
  return Orientation(row_correction * LA::sign_of_determinant(M));
}
};

template <class R>
struct Side_of_oriented_sphereCd { 
typedef typename R::Point_d Point_d;
typedef typename R::LA LA;
typedef typename R::FT FT;

template <class ForwardIterator> 
Oriented_side operator()(ForwardIterator first, ForwardIterator last, 
                         const Point_d& x)
{ 
  TUPLE_DIM_CHECK(first,last,Side_of_oriented_sphere_d);
  int d = static_cast<int>(std::distance(first,last)); // |A| contains |d| points
  CGAL_assertion_msg((d-1 == first->dimension()), 
  "Side_of_oriented_sphere_d: needs first->dimension()+1 many input points.");
  typename LA::Matrix M(d + 1); 
  for (int i = 0; i < d; ++first, ++i) { 
    FT Sum = 0;
    M(i,0) = 1;
    for (int j = 0; j < d-1; j++) { 
      FT cj = first->cartesian(j);
      M(i,j + 1) = cj; Sum += cj*cj; 
    }
    M(i,d) = Sum; 
  }
  FT Sum = 0; 
  M(d,0) = 1; 
  for (int j = 0; j < d-1; j++) { 
    FT hj = x.cartesian(j);
    M(d,j + 1) = hj; Sum += hj*hj; 
  }
  M(d,d) = Sum;
  return - LA::sign_of_determinant(M);
}
};

template <class R>
struct Side_of_bounded_sphereCd { 
typedef typename R::Point_d Point_d;

template <class ForwardIterator> 
Bounded_side operator()(ForwardIterator first, ForwardIterator last, 
                        const Point_d& p)
{
  TUPLE_DIM_CHECK(first,last,region_of_sphere);
  typename R::Orientation_d _orientation;
  Orientation o = _orientation(first,last);
  CGAL_assertion_msg((o != 0), "Side_of_bounded_sphere_d: \
  A must be full dimensional.");
  typename R::Side_of_oriented_sphere_d _side_of_oriented_sphere;
  Oriented_side oside = _side_of_oriented_sphere(first,last,p);
  if (o == POSITIVE) {
    switch (oside) {
        case ON_POSITIVE_SIDE    :   return ON_BOUNDED_SIDE;
        case ON_ORIENTED_BOUNDARY:   return ON_BOUNDARY;
        case ON_NEGATIVE_SIDE    :   return ON_UNBOUNDED_SIDE;
    }       
  } else {
    switch (oside) {
        case ON_POSITIVE_SIDE    :   return ON_UNBOUNDED_SIDE;
        case ON_ORIENTED_BOUNDARY:   return ON_BOUNDARY;
        case ON_NEGATIVE_SIDE    :   return ON_BOUNDED_SIDE;
    }     
  }
  return ON_BOUNDARY; // never reached
}
};


template <class R>
struct Contained_in_simplexCd { 
typedef typename R::Point_d Point_d;
typedef typename R::LA LA;
typedef typename R::FT FT;

template <class ForwardIterator> 
bool operator()(ForwardIterator first, ForwardIterator last,
                const Point_d& p)
{
  TUPLE_DIM_CHECK(first,last,Contained_in_simplex_d);
  int k = static_cast<int>(std::distance(first,last)); // |A| contains |k| points
  int d = first->dimension(); 
  CGAL_assertion_code(
    typename R::Affinely_independent_d check_independence; )
  CGAL_assertion_msg(check_independence(first,last),
    "Contained_in_simplex_d: A not affinely independent.");
  CGAL_assertion(d==p.dimension());

  typename LA::Matrix M(d + 1,k); 
  typename LA::Vector b(d +1);
  for (int j = 0; j < k; ++first, ++j) {
    for (int i = 0; i < d; ++i) 
      M(i,j) = first->cartesian(i);
    M(d,j) = 1;
  }
  for (int i = 0; i < d; ++i) 
    b[i] = p.cartesian(i);
  b[d] = 1;

  FT D; 
  typename LA::Vector lambda; 
  if ( LA::linear_solver(M,b,lambda,D) ) { 
    for (int j = 0; j < k; j++) { 
      if (lambda[j] < FT(0)) return false;
    }
    return true;
  }
  return false; 
}
};

template <class R>
struct Contained_in_affine_hullCd { 
typedef typename R::Point_d Point_d;
typedef typename R::LA LA;

template <class ForwardIterator> 
bool operator()(ForwardIterator first, ForwardIterator last,
                const Point_d& p) 
{
  TUPLE_DIM_CHECK(first,last,Contained_in_affine_hull_d);
  int k = static_cast<int>(std::distance(first,last)); // |A| contains |k| points
  int d = first->dimension(); 
  typename LA::Matrix M(d + 1,k); 
  typename LA::Vector b(d + 1); 
  for (int j = 0; j < k; ++first, ++j) {
    for (int i = 0; i < d; ++i) 
      M(i,j) = first->cartesian(i);
    M(d,j) = 1;
  }
  for (int i = 0; i < d; ++i)
    b[i] = p.cartesian(i);
  b[d] = 1;
  return LA::is_solvable(M,b);
}
};


template <class R>
struct Affine_rankCd { 
typedef typename R::Point_d Point_d;
typedef typename R::Vector_d Vector_d;
typedef typename R::LA LA;

template <class ForwardIterator> 
int operator()(ForwardIterator first, ForwardIterator last) 
{
  TUPLE_DIM_CHECK(first,last,Affine_rank_d);
  int k = static_cast<int>(std::distance(first,last)); // |A| contains |k| points
  if (k == 0) return -1;
  if (k == 1) return 0; 
  int d = first->dimension();
  typename LA::Matrix M(d,--k);
  Point_d p0 = *first; ++first; // first points to second
  for (int j = 0; j < k; ++first, ++j) {
    Vector_d v = *first - p0;
    for (int i = 0; i < d; i++) 
      M(i,j) = v.cartesian(i); 
  }
  return LA::rank(M);
}
};

template <class R>
struct Affinely_independentCd { 
typedef typename R::Point_d Point_d;
typedef typename R::LA LA;

template <class ForwardIterator> 
bool operator()(ForwardIterator first, ForwardIterator last) 
{ typename R::Affine_rank_d rank;
  int n = static_cast<int>(std::distance(first,last));
  return rank(first,last) == n-1;
}
};


template <class R>
struct Compare_lexicographicallyCd {
typedef typename R::Point_d Point_d;
typedef typename R::Point_d PointD; //MSVC hack
Comparison_result operator()(const Point_d& p1, const Point_d& p2)
{ return PointD::cmp(p1,p2); }
};

template <class R>
struct Contained_in_linear_hullCd {
typedef typename R::LA LA;
typedef typename R::FT FT;
typedef typename R::Vector_d Vector_d;

template<class ForwardIterator>
bool operator()(
  ForwardIterator first, ForwardIterator last, const Vector_d& x) 
{ TUPLE_DIM_CHECK(first,last,Contained_in_linear_hull_d);
  int k = static_cast<int>(std::distance(first,last)); // |A| contains |k| vectors
  int d = first->dimension();
  typename LA::Matrix M(d,k);
  typename LA::Vector b(d); 
  for (int i = 0; i < d; i++) { 
     b[i] = x.cartesian(i); 
     for (int j = 0; j < k; j++) 
       M(i,j) = (first+j)->cartesian(i); 
  }
  return LA::is_solvable(M,b); 
}
};

template <class R>
struct Linear_rankCd {
typedef typename R::LA LA;
template <class ForwardIterator>
int operator()(ForwardIterator first, ForwardIterator last)
{ TUPLE_DIM_CHECK(first,last,linear_rank);
  int k = static_cast<int>(std::distance(first,last)); // k vectors
  int d = first->dimension(); 
  typename LA::Matrix M(d,k);
  for (int i = 0; i < d  ; i++)
     for (int j = 0; j < k; j++)  
       M(i,j) = (first + j)->cartesian(i);
  return LA::rank(M);
}
};

template <class R>
struct Linearly_independentCd {
typedef typename R::LA LA;
template <class ForwardIterator>
bool operator()(ForwardIterator first, ForwardIterator last)
{ typename R::Linear_rank_d rank;
  return rank(first,last) == static_cast<int>(std::distance(first,last));
}
};

template <class R>
struct Linear_baseCd {
typedef typename R::LA LA;
typedef typename R::FT FT;
typedef typename R::Vector_d Vector_d;
template <class ForwardIterator, class OutputIterator>
OutputIterator operator()(ForwardIterator first, ForwardIterator last,
  OutputIterator result)
{ TUPLE_DIM_CHECK(first,last,linear_base);
  int k = static_cast<int>(std::distance(first,last)); // k vectors
  int d = first->dimension();
  typename LA::Matrix M(d,k); 
  for (int j = 0; j < k; ++first, ++j)
    for (int i = 0; i < d; i++)
      M(i,j) = first->cartesian(i);

  std::vector<int> indcols;
  int r = LA::independent_columns(M,indcols);

  for (int l=0; l < r; l++) {
    typename LA::Vector v = M.column(indcols[l]);
    *result++ = Vector_d(d,v.begin(),v.end());
  }
  return result;
}
};



} //namespace CGAL
#endif //CGAL_FUNCTION_OBJECTSCD_H
