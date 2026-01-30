// Copyright (c) 2022 INRIA (France) and GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Samuel Hornus, Sébastien Loriot and Léo Valque
//

#ifndef CGAL_CONVEX_HULL_3_DO_INTERSECT_H
#define CGAL_CONVEX_HULL_3_DO_INTERSECT_H

#include <CGAL/license/Convex_hull_3.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Kernel_23/internal/Has_boolean_tags.h>

#include <CGAL/Convex_hull_hierarchy.h>
#include <CGAL/extreme_point_3.h>
#include <CGAL/intersections.h>

#include <CGAL/Convex_hull_3/internal/helpers.h>

#include <CGAL/Container_helper.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/property_map.h>

#include <boost/range/value_type.hpp>
#include <boost/range/reference.hpp>

#include <boost/graph/adjacency_list.hpp>

#include <vector>

namespace CGAL {

namespace Convex_hull_3 {

namespace internal {

template <class Vector_3>
Vector_3 LInf_normalize(const Vector_3 &v){
  if(v==NULL_VECTOR)
    return v;
  auto l=(max)(abs(v.x()), (max)(abs(v.y()),abs(v.z())));
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
  SphericalPolygonElement(const Vector_3 & n) : north_(n) {}
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
    // PRECONDITION : all vertices (that are unit vectors on the sphere) are normalized.
    switch( this->size() ) {
      case 0 : return Vector_3(1,0,0); break; // An arbitrary one
      case 1 : return this->begin()->north_; break;
      case 2 :{   // The two vertices are opposite so we do not take their mean
                  // We want a direction with negative dot with bith north_
                  // We take the two perpendicular vector of resp north_0, and north_1 that lies on the same circle on them
                  // And we take a barycenter of them
                  Vector_3 perp1= LInf_normalize(cross_product((*this)[0].north_, (*this)[0].vertex_));
                  Vector_3 perp2= LInf_normalize(cross_product((*this)[1].north_, (*this)[1].vertex_));
                  CGAL_assertion(is_positive((perp1+perp2) * (*this)[0].north_));
                  CGAL_assertion(is_positive((perp1+perp2) * (*this)[1].north_));
                  return (perp1 + perp2);
                  // return (*this)[0].north_ + (*this)[1].north_; //Old, need that both north_ was L2 normalized
               } break;
      default : {
                  Vector_3 avg(NULL_VECTOR);
                  for( const SphericalPolygonElement<Vector_3> & v : *this )
                    avg += v.vertex_;
                  return avg;
                } break;
    }
  }

  void clip(const Vector_3& clipNorth, SphericalPolygon<Vector_3> &result) const {
    const int n = this->size();
    result.clear();
    switch( n ) {
      case 0 : {
                  result.emplace_back(clipNorth);
                  break;
               }
      case 1 : {
                 result = (*this);
                //  NT dot = this->begin()->north_ * clipNorth;
                 Vector_3 v(LInf_normalize(cross_product(clipNorth, this->begin()->north_)));
                 if(v==NULL_VECTOR /* && is_negative(dot) */){
                  CGAL_assertion(is_negative(clipNorth * this->begin()->north_)); //Could not be two equal hemispheres
                   // intersection of two opposite hemispheres ==> empty
                   result.clear();
                   break;
                 }

                 result.begin()->vertex_ = v;
                 result.emplace_back(-v, clipNorth);
                 break;
               }
      case 2 : {
                 result = (*this);
                 iterator next = result.begin();
                 iterator cur = next++;
                 NT vDot = this->begin()->vertex_ * clipNorth;
                 if( is_positive(vDot) ) {
                   // we'll get a triangle
                   next->vertex_ = LInf_normalize(cross_product(clipNorth, next->north_));
                   Vector_3 v( LInf_normalize(cross_product(cur->north_, clipNorth)));
                   result.emplace(next, v, clipNorth);
                 } else if( is_negative(vDot) ) {
                   // we'll get a triangle
                   cur->vertex_ =  LInf_normalize(cross_product(clipNorth, cur->north_));
                   Vector_3 v( LInf_normalize(cross_product(next->north_, clipNorth)));
                   result.emplace_back(v, clipNorth);
                 } else {
                   // we keep a moon crescent
                   NT curTest(clipNorth *  cross_product(cur->north_, cur->vertex_));
                   Vector_3 nextTest( cross_product(next->north_, next->vertex_));
                   if( is_positive(curTest) ) {
                    //The clipNorth is not between the two previous ones
                    CGAL_assertion(!is_positive(clipNorth * nextTest));
                     next->north_ = clipNorth;
                     cur->vertex_ =  LInf_normalize(cross_product(next->north_, cur->north_));
                     next->vertex_ = -cur->vertex_;
                   } else {
                     if( is_positive(clipNorth * nextTest) ) {
                       cur->north_ = clipNorth;
                       next->vertex_ =  LInf_normalize(cross_product(cur->north_, next->north_));
                       cur->vertex_ = -next->vertex_;
                     } else {
                       //killed the crescent
                       result.clear();
                     }
                   }
                 }
                 break;
               }
      default : { // n >= 3
                  int nbKept(0);
                  const_iterator cur = this->begin();
                  NT nextDot, curDot = clipNorth * cur->vertex_;
                  while( cur != this->end() ) {
                    if( cur+1 == this->end() )
                      nextDot = clipNorth * this->begin()->vertex_;
                    else
                      nextDot = clipNorth * (cur+1)->vertex_;
                    if( is_positive(curDot) ) { // cur is "IN"
                      ++nbKept;
                      result.push_back(*cur);
                      if( is_negative(nextDot) ) { // next is "OUT"
                        result.emplace_back( LInf_normalize(cross_product(cur->north_, clipNorth)), clipNorth);
                      }
                    } else if( !is_negative(curDot) ) { // cur is "ON" the clipping plane
                      ++nbKept;
                      if ( is_negative(nextDot) ) // next is "OUT"
                        result.emplace_back(cur->vertex_, clipNorth);
                      else
                        result.push_back(*cur);
                    } else { // cur is "OUT"
                      if ( is_positive(nextDot) ) { // next is "IN"
                        result.emplace_back( LInf_normalize(cross_product(clipNorth, cur->north_)), cur->north_);
                      }
                    }
                    curDot = nextDot;
                    ++cur;
                  }
                  if( (result.size() < 3/*too small*/) || ((nbKept == n)/*no change*/) ) {
                    result.clear();
                  }
                  break;
                }
    }
  }
};

template<typename K>
struct Functor_do_intersect{

  template<typename Convex, typename NamedParameter1, typename NamedParameter2>
  bool operator()(const Convex &a, const Convex &b, const NamedParameter1 &np1, const NamedParameter2 &np2) const{
    using Vector_3= typename K::Vector_3;
    const int INTER_MAX_ITER=0;
    SphericalPolygon<Vector_3> positiveBound, tempPoly;
    positiveBound.clear();
    unsigned long planeStatPerPair = 0;
    do {
      Vector_3 dir = positiveBound.averageDirection();
      Vector_3 sp = extreme_point_3_wrapper(a, dir.direction(), np1) - extreme_point_3_wrapper(b, -dir.direction(), np2);
      if(sp==NULL_VECTOR) return true;
      if(is_negative(sp * dir)) return false;
      if(INTER_MAX_ITER!=0 && (++planeStatPerPair >= INTER_MAX_ITER)) return true;
      positiveBound.clip(-sp, tempPoly); positiveBound.swap(tempPoly);
    } while( !positiveBound.empty() );
    return true;
  }
};

} // end of internal namespace

//Do_intersect_traits, only used to filter the kernel
template<typename K,
         typename IK=K,
         typename Converter=Cartesian_converter<K, K>,
         bool Has_filtered_predicates_ = CGAL::internal::Has_filtered_predicates<K>::value>
struct Do_intersect_traits;

template<typename K, typename IK, typename Converter>
struct Do_intersect_traits<K, IK, Converter, false>{
  typedef internal::Functor_do_intersect<K> Do_intersect;
  Do_intersect do_intersect_object() const {
    return Do_intersect();
  }
};
template<typename K, typename Converter>
struct Do_intersect_traits<K, K, Converter, true> {
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Exact_kernel::Vector_3 EVector_3;
  typedef typename K::Approximate_kernel::Vector_3 FVector_3;

  typedef Cartesian_converter<K, K>  IdentityConverter;
  typedef typename K::C2E C2E;
  typedef typename K::C2F C2F;

  typedef Do_intersect_traits<typename K::Exact_kernel, K, C2E> Exact_traits;
  typedef Do_intersect_traits<typename K::Approximate_kernel, K, C2F> Filtering_traits;

  //The conversion are made lazely by Do_intersect, we thus use IdentityConverter in Filtered_predicate
  typedef Filtered_predicate<
              typename Exact_traits::Do_intersect,
              typename Filtering_traits::Do_intersect,
              IdentityConverter,
              IdentityConverter>  Do_intersect;

  Do_intersect do_intersect_object() const
  {
    typename Exact_traits::Do_intersect pe = Exact_traits().do_intersect_object();
    typename Filtering_traits::Do_intersect pf = Filtering_traits().do_intersect_object();

    return Do_intersect(pe, pf);
  }
};

/**
* \ingroup PkgConvexHull3Predicates
*
* \brief checks if the convex hulls intersect or not.
*
* each input can be provided as a range, a mesh, or as the specialized structure CGAL::Convex_hull_hierarchy.
* They are not required to use the same input type.
*
* @tparam Convex_1 is a model of `ConstRange` or a model of `VertexListGraph` and `AdjacencyGraph` or an instance of `CGAL::Convex_hull_hierarchy`
* @tparam Convex_1 same as Convex_2
* @tparam NamedParameters_1 a sequence of \ref bgl_namedparameters "Named Parameters"
* @tparam NamedParameters_2 a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param ch1 the first convex hull
* @param ch2 the second convex hull
* @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
* @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{when `ch1` (`ch2`) is a range, it is a property map associating points to its elements}
*     \cgalParamType{a model of `ReadablePropertyMap` whose value types are the same for `ch1` and `ch2`}
*     \cgalParamDefault{`CGAL::Identity_property_map`}
*     \cgalParamExtra{used only if `ch1` (`ch2`) is model of `ConstRange`}
*   \cgalParamNEnd
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{when `ch1` (`ch2`) is a mesh, it is a property map associating points to its vertices}
*     \cgalParamType{a model of `ReadablePropertyMap` whose value types are the same for `ch1` and `ch2`}
*     \cgalParamDefault{boost::get(CGAL::vertex_point, g)}
*     \cgalParamExtra{used only if `ch1` (`ch2`) is model of `VertexListGraph` and `AdjacencyGraph`.
*                     If this parameter is omitted, an internal property map for `CGAL::vertex_point_t` must be available in `Convex_1` (`Convex_2`).}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{An instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal kernel deduced from the point type of the input, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{`np1` only}
*   \cgalParamNEnd
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{if not `0` (no limit), indicates the maximum number of iterations performed by the algorithm.
*                           If this value is not `0`, then an intersection might be reported even if the convex hulls does not intersect.
*                           However, if the convex hulls are reported not to intersect, this is guaranteed.}
*     \cgalParamType{a positive integer convertible to `std::size_t`}
*     \cgalParamExtra{`np1` only}
*     \cgalParamDefault{`0`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \see `CGAL::Convex_hull_hierarchy`
*/
template <class Convex1, class Convex2,
          class NamedParameters_1 = parameters::Default_named_parameters,
          class NamedParameters_2 = parameters::Default_named_parameters>
bool do_intersect(const Convex1& ch1, const Convex2& ch2,
                  const NamedParameters_1& np1 = parameters::default_values(),
                  const NamedParameters_2& np2 = parameters::default_values()){
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;
  using GT= typename internal::GetGeomTraitsFromConvex<Convex1, NamedParameters_1>::type;
  // GT gt = choose_parameter<GT>(get_parameter(np1, internal_np::geom_traits));
  return Do_intersect_traits<GT>().do_intersect_object()(ch1, ch2, np1, np2);
}

}} // CGAL::Convex_hull_3 namespace

#endif // CGAL_CONVEX_HULL_3_DO_INTERSECT_H
