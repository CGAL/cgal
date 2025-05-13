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

#ifndef CGAL_CONVEX_HULL_3_PREDICATES_H
#define CGAL_CONVEX_HULL_3_PREDICATES_H

#include <CGAL/license/Convex_hull_3.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Kernel_23/internal/Has_boolean_tags.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Real_timer.h>

#include <vector>

int nb_visited_in_the_inner=0;
int nb_visited=0;

namespace CGAL {

namespace Convex_hull_3 {

namespace predicates_impl
{

// spherical.h
//----------------------------

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
    // PRECONDITION : all vertex are normalized.
    switch( this->size() ) {
      case 0 : return Vector_3(1,0,0); break; // An arbitrary one
      case 1 : return this->begin()->north_; break;
      case 2 :{   // The two vertex are opposite so we do not take their mean
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

template<class T, template<class> class U>
inline constexpr bool is_instance_of_v = std::false_type{};

template<template<class> class U, class V>
inline constexpr bool is_instance_of_v<U<V>,U> = std::true_type{};

template <class Value, class Converter, class Vector_3>
const Value extreme_point(const Surface_mesh<Value>& C, const Vector_3 &dir, const Converter& converter) {
  using Point_3= typename Kernel_traits<Vector_3>::Kernel::Point_3;
  using Convex= Surface_mesh<Value>;
  using FT= typename Kernel_traits<Vector_3>::Kernel::FT;

  if(C.vertices().size()<5){
    typename Convex::Vertex_index argmax=*C.vertices().begin();
    FT tmax=Vector_3(ORIGIN, converter(C.point(argmax)))*dir;
    for(auto vh=++(C.vertices().begin()); vh!=C.vertices().end(); ++vh){
      typename Convex::Vertex_index v=*vh;
      ++nb_visited;
      FT p=Vector_3(ORIGIN, converter(C.point(v)))*dir;
      if(compare(tmax, p)==SMALLER){
        tmax=p;
        argmax=v;
      }
    }
    return C.point(argmax);
  }

  typename Convex::Vertex_index argmax=*C.vertices().begin();
  nb_visited++;
  FT tmax= Vector_3(ORIGIN, converter(C.point(argmax)))*dir;
  bool is_local_max;
  do{
    is_local_max=true;
    for(auto v: vertices_around_target(argmax ,C)){
      FT p=Vector_3(ORIGIN, converter(C.point(v)))*dir;
      ++nb_visited;
      if(compare(tmax, p)==SMALLER){
        tmax=p;
        argmax=v;
        is_local_max=false; // repeat with the new vertex
        break;
      }
    }
  }while(!is_local_max);
  // Since convex, local maximum is a global maximum
  return C.point(argmax);
}

template <class Converter, class Value, class Vector_3>
Value extreme_point(const std::vector<Value>& C, const Vector_3 &dir, const Converter &converter) {
  using FT= typename Kernel_traits<Vector_3>::Kernel::FT;
  using Convex=std::vector<Value>;
  typename Convex::const_iterator argmax=C.begin();
  FT tmax= Vector_3(ORIGIN, converter(*argmax))*dir;
  for(typename Convex::const_iterator it=C.begin()+1; it!=C.end(); ++it){
    FT v=Vector_3(ORIGIN, converter(*it))*dir;
    if(compare(tmax, v)==SMALLER){
      tmax=v;
      argmax=it;
    }
  }
  return *argmax;
}

template<typename IK, typename OK, typename Converter>
struct Functor_do_intersect{
  typedef typename OK::Boolean  result_type;

  Converter c1;
  Converter c2;
  Functor_do_intersect():c1(Converter()), c2(Converter()){}
  Functor_do_intersect(Converter &&c1_, Converter &&c2_):c1(std::move(c1_)), c2(std::move(c2_)){}

  template<typename Convex>
  bool operator()(const Convex &a, const Convex &b, unsigned long INTER_MAX_ITER) const{
    using Point_3 = typename OK::Point_3;
    using Vector_3= typename OK::Vector_3;

    SphericalPolygon<Vector_3> positiveBound, tempPoly;
    positiveBound.clear();
    unsigned long planeStatPerPair = 0;
    do {
      Vector_3 dir = positiveBound.averageDirection();
      Vector_3 sp = c1(extreme_point(a, dir, c1)) - c2(extreme_point(b, -dir, c2));
      if(sp==NULL_VECTOR) return true;
      if(is_negative(sp * dir)) return false;
      if(INTER_MAX_ITER!=0 && (++planeStatPerPair >= INTER_MAX_ITER)) return true;
      positiveBound.clip(-sp, tempPoly); positiveBound.swap(tempPoly);
    } while( !positiveBound.empty() );
    return true;
  }
};

} // end of predicates_impl namespace

//Do_intersect_traits
template<typename K,
         typename IK=K,
         typename Converter=Cartesian_converter<K, K>,
         bool Has_filtered_predicates_ = CGAL::internal::Has_filtered_predicates<K>::value>
struct Do_intersect_traits;

template<typename K, typename IK, typename Converter>
struct Do_intersect_traits<K, IK, Converter, false>{
  typedef predicates_impl::Functor_do_intersect<IK, K, Converter> Do_intersect;
  Do_intersect do_intersect_object() const {
    return Do_intersect(Converter(), Converter());
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

template<typename K,
         typename PointMap,
         typename IK=K,
         typename Converter=Cartesian_converter<K, K>,
         bool Has_filtered_predicates_ = CGAL::internal::Has_filtered_predicates<K>::value>
struct Do_intersect_traits_with_point_maps{
  Do_intersect_traits_with_point_maps(const PointMap &map1_,const PointMap &map2_);
};

template<typename K, typename PointMap, typename IK, typename Converter>
struct Do_intersect_traits_with_point_maps<K, PointMap, IK, Converter, false>{
  const PointMap &map1;
  const PointMap &map2;

  Do_intersect_traits_with_point_maps(const PointMap &map1_,const PointMap &map2_):map1(map1_), map2(map2_){}

  struct PointMapConverter : Converter{
    const PointMap &map;

    PointMapConverter(const PointMap &map_):map(map_){}

    template<typename Vertex_index>
    typename K::Point_3 operator()(Vertex_index vi) const{
      return Converter()(map[vi]);
    }
  };

  typedef predicates_impl::Functor_do_intersect<IK, K, PointMapConverter> Do_intersect;
  Do_intersect do_intersect_object() const {
    return Do_intersect(PointMapConverter(map1), PointMapConverter(map2));
  }
};
template<typename K, typename PointMap, typename Converter>
struct Do_intersect_traits_with_point_maps<K, PointMap, K, Converter, true> {
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Exact_kernel::Vector_3 EVector_3;
  typedef typename K::Approximate_kernel::Vector_3 FVector_3;

  typedef Cartesian_converter<K, K>  IdentityConverter;
  typedef typename K::C2E C2E;
  typedef typename K::C2F C2F;

  typedef Do_intersect_traits_with_point_maps<typename K::Exact_kernel, PointMap, K, C2E> Exact_traits;
  typedef Do_intersect_traits_with_point_maps<typename K::Approximate_kernel, PointMap, K, C2F> Filtering_traits;

  typedef Filtered_predicate<
              typename Exact_traits::Do_intersect,
              typename Filtering_traits::Do_intersect,
              IdentityConverter,
              IdentityConverter>  Do_intersect;

  const PointMap &map1;
  const PointMap &map2;
  Do_intersect_traits_with_point_maps(const PointMap &map1_,const PointMap &map2_):map1(map1_), map2(map2_){}

  Do_intersect do_intersect_object() const
  {
    typename Exact_traits::Do_intersect pe = Exact_traits(map1, map2).do_intersect_object();
    typename Filtering_traits::Do_intersect pf = Filtering_traits(map1, map2).do_intersect_object();

    return Do_intersect(pe, pf);
  }
};

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
template <class Kernel, class Object1, class Object2, class Traits >
bool do_intersect(const Object1& obj1, const Object2& obj2, Traits traits)
{
  return traits.do_intersect_object()(obj1, obj2, 0);
}

template <class Kernel, class Object1, class Object2>
bool do_intersect(const Object1& obj1, const Object2& obj2)
{
  //Find kernel
  return Do_intersect_traits<Kernel>().do_intersect_object()(obj1, obj2, 0);
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
