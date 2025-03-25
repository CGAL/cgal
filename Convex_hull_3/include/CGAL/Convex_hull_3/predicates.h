// Copyright (c) 2022 INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Samuel Hornus and Sébastien Loriot
//

#ifndef CGAL_CONVEX_HULL_3_PREDICATES_H
#define CGAL_CONVEX_HULL_3_PREDICATES_H

#include <CGAL/license/Convex_hull_3.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <vector>

namespace CGAL {

namespace Convex_hull_3 {

namespace predicates_impl
{

// spherical.h
//----------------------------

template <class Vector_3>
Vector_3 LInf_normalize(const Vector_3 &v){
  auto l=(max)(v.x(), (max)(v.y(),v.z()));
  return v/l;
}

template <class Vector_3>
struct SphericalPolygonElement {
  Vector_3 vertex_; // A vertex of the spherical polygon
  Vector_3 north_; // The north pole of the equatorial arc/edge leading OUT OF that vertex_ (arcs are oriented west-to-east, or CCW in a right-handed frame.
  // In the spherical polygon (v0, n0), (v1, n1), (v2, n2), ... we have
  // v1 = cross(n0, n1),  more generally: v_{i+1} = cross(n_i, n_{i+1})  and
  // n1 = cross(v1, v2),  more generally: n_i     = cross(v_i, v_{i+1}).
  SphericalPolygonElement(){}
  SphericalPolygonElement(const Vector_3 & n) : north_(LInf_normalize(n)) {}
  SphericalPolygonElement(const Vector_3 & v, const Vector_3 & n) : vertex_(v), north_(n) {}
};

template <class Vector_3>
struct SphericalPolygon : public std::vector<SphericalPolygonElement<Vector_3>> {

  typedef std::vector<SphericalPolygonElement<Vector_3>> Base;
  typedef typename Base::iterator iterator;
  typedef typename Base::const_iterator const_iterator;
  typedef typename Kernel_traits<Vector_3>::Kernel Kernel;
  typedef typename Kernel::FT NT;

  SphericalPolygon() {
    this->reserve(16);
  }

  Vector_3 averageDirection() const {
    // PRECONDITION : all northes are normalized.
    switch( this->size() ) {
      case 0 : return NULL_VECTOR; break;
      case 1 : return this->begin()->north_; break;
      case 2 :{   // The two vertex are opposite so we do not take their mean
                  // We want a direction with negative dot with bith north_
                  // We take the two perpendicular vector of resp north_0, and north_1 that lies on the same circle on them
                  // And we take a barycenter of them
                  Vector_3 perp1= LInf_normalize(cross_product((*this)[0].north_, (*this)[0].vertex_));
                  Vector_3 perp2= LInf_normalize(cross_product((*this)[1].north_, (*this)[1].vertex_));
                  return perp1 + perp2;
                  // return (*this)[0].north_ + (*this)[1].north_; //Old, need that both north_ was L2 normalized
               } break;
      default : {
                  Vector_3 avg;
                  for( const SphericalPolygonElement<Vector_3> & v : *this )
                    avg = avg + v.vertex_;
                  return avg;
                } break;
    }
  }

  void clip(const Vector_3& OrigVertex, SphericalPolygon<Vector_3> & result, bool doClean=true) const {
    // PRECONDITION : clipNorth, and all northes are normalized.
#define _ray_spherical_eps 1e-6f
    const int n = this->size();
    result.clear();
    switch( n ) {
      case 0 : break; // 0 means empty, so nothing to do
      case 1 : {
                 result = (*this);
                 Vector_3 clipNorth = LInf_normalize(OrigVertex);
                 NT dot = this->begin()->north_ * clipNorth;
                 if( dot < -0.99984769515 ) { // about one degree
                   // intersection of two almost opposite hemispheres ==> empty
                   result.clear();
                   break;
                 } else if( dot > 0.99984769515 ) {
                   break;
                 }
                 Vector_3 v(LInf_normalize(cross_product(clipNorth, this->begin()->north_)));
                 result.begin()->vertex_ = v;
                 result.emplace_back(-v, clipNorth);
                 break;
               }
      case 2 : {
                 result = (*this);
                 Vector_3 clipNorth = LInf_normalize(OrigVertex);
                 iterator next = result.begin();
                 iterator cur = next++;
                 NT vDot = this->begin()->vertex_ * clipNorth;
                 if( vDot >= _ray_spherical_eps ) {
                   // we'll get a triangle
                   next->vertex_ = LInf_normalize(cross_product(clipNorth, next->north_));
                   Vector_3 v( LInf_normalize(cross_product(cur->north_, clipNorth)));
                   result.emplace(next, v, clipNorth);
                 } else if( vDot <= - _ray_spherical_eps ) {
                   // we'll get a triangle
                   cur->vertex_ =  LInf_normalize(cross_product(clipNorth, cur->north_));
                   Vector_3 v( LInf_normalize(cross_product(next->north_, clipNorth)));
                   result.emplace_back(v, clipNorth);
                 } else {
                   // we keep a moon crescent
                   NT curTest(clipNorth *  cross_product(cur->north_, cur->vertex_));
                   Vector_3 nextTest( cross_product(next->north_, next->vertex_));
                   if( curTest > 0.0f ) {
                     if( (clipNorth * nextTest) <= 0.0f ) {
                       next->north_ = clipNorth;
                       cur->vertex_ =  LInf_normalize(cross_product(next->north_, cur->north_));
                       next->vertex_ = -cur->vertex_;
                     } else {
                       // the crescent is unchanged
                       //std::cerr << "kept a crescent\n";
                     }
                   } else {
                     if( (clipNorth * nextTest) > 0.0f ) {
                       cur->north_ = clipNorth;
                       next->vertex_ =  LInf_normalize(cross_product(cur->north_, next->north_));
                       cur->vertex_ = -next->vertex_;
                     } else {
                       //std::cerr << "killed a crescent\n";
                       result.clear();
                     }
                   }
                 }
                 break;
               }
      default : { // n >= 3
                  int nbKept(0);
                  const_iterator cur = this->begin();
                  Vector_3 clipNorth = LInf_normalize(OrigVertex);
                  NT nextDot, curDot = clipNorth * cur->vertex_;
                  while( cur != this->end() ) {
                    if( cur+1 == this->end() )
                      nextDot = clipNorth * this->begin()->vertex_;
                    else
                      nextDot = clipNorth * (cur+1)->vertex_;
                    if( curDot >= _ray_spherical_eps ) { // cur is "IN"
                      ++nbKept;
                      result.push_back(*cur);
                      if( nextDot <= -_ray_spherical_eps ) { // next is "OUT"
                        result.emplace_back( LInf_normalize(cross_product(cur->north_, clipNorth)), clipNorth);
                      }
                    } else if( curDot > -_ray_spherical_eps ) { // cur is "ON" the clipping plane
                      ++nbKept;
                      if ( nextDot <= -_ray_spherical_eps ) // next is "OUT"
                        result.emplace_back(cur->vertex_, clipNorth);
                      else
                        result.push_back(*cur);
                    } else { // cur is "OUT"
                      if ( nextDot >= _ray_spherical_eps ) { // next is "IN"
                        result.emplace_back( LInf_normalize(cross_product(clipNorth, cur->north_)), cur->north_);
                      }
                    }
                    curDot = nextDot;
                    ++cur;
                  }
                  if( (result.size() < 3/*too small*/) || ((nbKept == n)/*no change*/ && doClean) ) {
                    result.clear();
                  }
                  //if( nbKept == n ) {
                    //std::cerr << "**";
                  //}
                  break;
                }
    }
  }
};

//----------------------------

template <class Vector_3, class Convex>
bool differenceCoversZeroInDir(const Convex& A, const Convex& B, int & vA, int & vB, const Vector_3 & dir) {
  using Kernel= typename Kernel_traits<Vector_3>::Kernel;
  using NT= typename Kernel::FT;

  /* TODO
  Actual complexity is linear of the size of A and B
  With the graph of the convex_hull, we have O(sqrt(n)) by compute the product on adjacent vertices and
  progressed along the incrasing one.
  With some preprocess (adding temporary edge of the construction of the graph of the convex hull), this can be
  O(log n)
  */

  // difference above is: A - B
  NT maxOverA = Vector_3(ORIGIN, A[0]) * dir;
  NT minOverB = Vector_3(ORIGIN, B[0]) * dir;
  vA = vB = 0;
  const int na = A.size();
  for( int i = 1; i < na; ++i ) {
    NT tempA = Vector_3(ORIGIN, A[i]) * dir;
    if( tempA > maxOverA ) { maxOverA = tempA; vA = i; }
  }
  const int nb = B.size();
  for( int i = 1; i < nb; ++i ) {
    NT tempB = Vector_3(ORIGIN, B[i]) * dir;
    if( tempB < minOverB ) { minOverB = tempB; vB = i; }
  }
  return maxOverA >= minOverB;
}

template<typename Convex >
bool sphericalDisjoint(const Convex & a, const Convex & b, unsigned long INTER_MAX_ITER) {
  using Point_3 = std::remove_cv_t<typename std::iterator_traits<typename Convex::const_iterator>::value_type>;
  using Kernel= typename Kernel_traits<Point_3>::Kernel;
  using Vector_3= typename Kernel::Vector_3;

  SphericalPolygon<Vector_3> positiveBound, tempPoly;
  int vA, vB;
  Vector_3 dir(b[0] - a[0]);
  if( ! differenceCoversZeroInDir(a, b, vA, vB, dir) ) return true;
  positiveBound.clear();
  positiveBound.emplace_back(dir);
  positiveBound.clip(b[vB] - a[vA], tempPoly); positiveBound.swap(tempPoly);
  if( positiveBound.empty() ) return false;
  unsigned long planeStatPerPair = 0;
  do {
    if( ! differenceCoversZeroInDir(a, b, vA, vB, positiveBound.averageDirection()) ) return true;
    if(INTER_MAX_ITER!=0 && (++planeStatPerPair >= INTER_MAX_ITER))
    {
      return false;
    }
    positiveBound.clip(b[vB] - a[vA], tempPoly); positiveBound.swap(tempPoly);
    if( positiveBound.empty() ) return false;
  } while( true );
}

} // end of predicates_impl namespace

/**
* \ingroup PkgConvexHull3Predicates
*
* indicates if the convex hulls of two point sets intersect or not.
*
* @tparam PointRange1 a model of the concept `RandomAccessContainer`
* whose value type is a point type from a CGAL %Kernel
* @tparam PointRange2 a model of the concept `RandomAccessContainer`
* whose value type is the point type from a CGAL %Kernel
* @tparam NamedParameters_1 a sequence of \ref bgl_namedparameters "Named Parameters"
* @tparam NamedParameters_2 a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param r1 first range of points whose convex hull is considered in the do-intersect test
* @param r2 first range of points whose convex hull is considered in the do-intersect test
* @param np_1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
* @param np_2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of the range `r1` (`r2)}
*     \cgalParamType{a model of `ReadablePropertyMap` whose value types are the same for `np1` and `np2`}
*     \cgalParamDefault{`CGAL::Identity_property_map`}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{`np1` only}
*   \cgalParamNEnd
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{if not `0` (no limit), indicates the maximum number of iterations that the algorithm is allowed to do.
*                           If this value is not `0`, then an intersection might be reported even if the convex hulls does not intersect.
                            However, if the convex hulls are reported not to intersect, this is guaranteed.}
*     \cgalParamType{an positive integer convertible to `std::size_t`}
*     \cgalParamDefault{`0`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
*/
template <class PointRange1, class PointRange2,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters>
bool do_intersect(const PointRange1& r1, const PointRange2& r2,
                  const NamedParameters1 np1 = parameters::default_values(),
                  const NamedParameters2 np2 = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef Point_set_processing_3_np_helper<PointRange1, NamedParameters1> NP_helper1;
  typedef typename NP_helper1::Const_point_map PointMap1;
  typedef typename NP_helper1::Geom_traits Geom_traits;

  typedef Point_set_processing_3_np_helper<PointRange2, NamedParameters2> NP_helper2;
  typedef typename NP_helper2::Const_point_map PointMap2;

  typedef typename boost::property_traits<PointMap1>::value_type Point_3;
  CGAL_static_assertion((std::is_same<Point_3, typename boost::property_traits<PointMap2>::value_type>::value));

  PointMap1 point_map1 = NP_helper1::get_const_point_map(r1, np1);
  PointMap2 point_map2 = NP_helper2::get_const_point_map(r2, np2);

  // TODO: avoid doing a copy
  std::vector<typename Geom_traits::Point_3> a, b;
  a.reserve(r1.size());
  b.reserve(r2.size());
  for (const auto& p : r1)
    a.push_back(get(point_map1, p));
  for (const auto& p : r2)
    b.push_back(get(point_map2, p));

  unsigned int max_nb_iterations = choose_parameter(get_parameter(np1, internal_np::number_of_iterations), 0);

  return !predicates_impl::sphericalDisjoint(a, b, max_nb_iterations);
}

template <class PointRange1, class PointRange2,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
FT
#else
typename Point_set_processing_3_np_helper<PointRange1, NamedParameters1>::Geom_traits::FT
#endif
separation_distance(const PointRange1& r1, const PointRange2& r2,
                    const NamedParameters1 np1 = parameters::default_values(),
                    const NamedParameters2 np2 = parameters::default_values());



// TODO: add OBB dedicated code in OBB package (do-intersect for sure and check for separation_distance).

}} // CGAL::Convex_hull_3 namespace

#endif // CGAL_CONVEX_HULL_3_PREDICATES_H
