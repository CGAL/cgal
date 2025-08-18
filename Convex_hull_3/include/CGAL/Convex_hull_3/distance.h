// Copyright (c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Léo Valque
//

#ifndef CGAL_CONVEX_HULL_3_DISTANCE_H
#define CGAL_CONVEX_HULL_3_DISTANCE_H

#include <CGAL/license/Convex_hull_3.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Distance_3/internal/squared_distance_utils_3.h>
#include <CGAL/Kernel_23/internal/Has_boolean_tags.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Convex_hull_hierarchy.h>
#include <CGAL/Convex_hull_3/intersections.h>
#include <CGAL/squared_distance_3.h>

#include <CGAL/IO/helpers.h>

#include <CGAL/Container_helper.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/property_map.h>

#include <boost/range/value_type.hpp>
#include <boost/range/reference.hpp>

#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Real_timer.h>

#include <vector>

namespace CGAL {
namespace Convex_hull_3 {

namespace predicates_impl {

// returns true iff pt is on the negative side of the plane defined by (ep0, ep1) and normal
template <class K>
inline bool
strictly_on_left_of_triangle_edge(const typename K::Point_3& pt,
                         const typename K::Vector_3& normal,
                         const typename K::Point_3& ep0,
                         const typename K::Point_3& ep1,
                         const K& k)
{
  typedef typename K::RT RT;
  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();

  const Vector_3 edge = vector(ep0, ep1);
  const Vector_3 diff = vector(ep0, pt);

  return (::CGAL::internal::wdot(::CGAL::internal::wcross(edge, normal, k), diff, k) < RT(0));
}

template <class K>
inline bool
strictly_on_left_of_triangle_edge(const typename K::Vector_3& pt,
                                  const typename K::Vector_3& normal,
                                  const typename K::Vector_3& edge,
                                  const K& k)
{
  auto cross=K().construct_cross_product_vector_3_object();
  auto dot=k.compute_scalar_product_3_object();
  auto orientation = K().orientation_3_object();

  // return (::CGAL::internal::wdot(::CGAL::internal::wcross(edge, normal, k), diff, k) < RT(0));
  return is_negative(dot(cross(edge, normal), pt));
}

template<typename K, typename Vector_3>
Vector_3 triangle_dir_to_origin(boost::container::small_vector<Vector_3, 4>& simplex){
  using Triangle_3=typename K::Triangle_3;
  using Plane_3=typename K::Plane_3;
  using FT=typename K::FT;

  const typename K::Point_3 O(0,0,0);

  auto cross=K().construct_cross_product_vector_3_object();
  auto dot=K().compute_scalar_product_3_object();
  auto orientation = K().orientation_3_object();

  Vector_3 &a = simplex[0];
  Vector_3 &b = simplex[1];
  Vector_3 &c = simplex[2];
  Vector_3 ab = b-a;
  Vector_3 ca = a-c;
  Vector_3 bc = c-b;

  Vector_3 n = cross(ab,-bc);
  // assert(n!=NULL_VECTOR);
  if(n==NULL_VECTOR){
    //We preserve the last found support
    std::swap(simplex[1],simplex[2]);
    simplex.pop_back();
    return -cross(cross(ab, a), ab);
  }

  // The origin cannot be on the side of ab, we abuse of this
  FT ratio_ab=dot(ab, -a)/ab.squared_length();
  FT ratio_ca=dot(ca, -c); ///ca.squared_length();
  //FT ratio_bc=dot(bc, -b)/bc.squared_length();
  FT ratio_cb=dot(-bc, -c);

  // bool on_left_of_ab = strictly_on_left_of_triangle_edge(-a, n, ab, K());
  bool on_left_of_ca = strictly_on_left_of_triangle_edge(-c, n, ca, K());
  bool on_left_of_bc = strictly_on_left_of_triangle_edge(-b, n, bc, K());

  if(/*!on_left_of_ab &&*/ !on_left_of_ca && !on_left_of_bc){
    if(is_negative(dot(n, -a))){
      // std::swap(simplex[1], simplex[2]);
      return -n;
    }
    std::swap(simplex[1], simplex[2]);
    return n;
  }

  // if(on_left_of_ab && ratio_ab>=0 && ratio_ab<=1){
  //   simplex.pop_back();
  //   return -cross(cross(ab, a), ab);
  // }
  if(on_left_of_ca && ratio_ca>=0 /*&& ratio_ca<=1*/){
    simplex[1] = simplex[2];
    simplex.pop_back();
    return -cross(cross(ca, c), ca);
  }
  if(on_left_of_bc && ratio_cb>=0/*ratio_bc>=0 && ratio_bc<=1*/){
    simplex[0] = simplex[2];
    simplex.pop_back();
    return -cross(cross(bc, b), bc);
  }

  // if(ratio_ab<0 && ratio_ca>1){
  //   simplex.pop_back();
  //   simplex.pop_back();
  //   return -a;
  // }
  // if(ratio_ab>1 && ratio_bc<0){
  //   simplex[0] = simplex[1];
  //   simplex.pop_back();
  //   simplex.pop_back();
  //   return -b;
  // }
  simplex[0] = simplex[2];
  simplex.pop_back();
  simplex.pop_back();
  return -c;
}

template<typename K, typename Vector_3>
Vector_3 tetrahedron_dir_to_origin(boost::container::small_vector<Vector_3, 4>& simplex) {
  using Triangle_3=typename K::Triangle_3;
  using Plane_3=typename K::Plane_3;
  using FT=typename K::FT;

  auto cp=K().construct_cross_product_vector_3_object();
  auto cross=K().construct_cross_product_vector_3_object();
  auto dot=K().compute_scalar_product_3_object();
  auto orientation = K().orientation_3_object();

  typename K::Point_3 O(0,0,0);

  // We use the fact that the origin do not be on the side of abc

  Vector_3 &a=simplex[0];
  Vector_3 &b=simplex[1];
  Vector_3 &c=simplex[2];
  Vector_3 &d=simplex[3];

  Vector_3 ab=b-a;
  Vector_3 ac=c-a;
  Vector_3 ad=d-a;

  // assert(orientation(ab, ac, ad) != NEGATIVE);
  if(orientation(ab, ac, ad) == NEGATIVE){
    std::swap(simplex[1], simplex[2]);
    // b=simplex[1];
    // c=simplex[2];
    std::swap(ab,ac);
  }

  bool ori_adb=orientation(a,d,b) == POSITIVE;
  bool ori_bdc=orientation(b,d,c) == POSITIVE;
  bool ori_cda=orientation(c,d,a) == POSITIVE;

  // std::cout << ori_adb << ori_bdc << ori_cda << std::endl;

  if(orientation(ab, ac, ad) == COPLANAR){
    ori_adb=!Triangle_3(O+a,O+b,O+d).is_degenerate();
    ori_bdc=!Triangle_3(O+c,O+b,O+d).is_degenerate();
    ori_cda=!Triangle_3(O+a,O+c,O+d).is_degenerate();
  }

  if(!ori_adb && !ori_bdc && !ori_cda)
    return NULL_VECTOR;

  if(ori_adb && Triangle_3(O+a,O+b,O+d).has_on(Plane_3(O+a,O+b,O+d).projection(O))){
    Vector_3 dir = cross(ab, ad);
    simplex[2]=simplex[3];
    simplex.pop_back();
    if(is_negative(dot(dir, -a))){
      std::swap(simplex[1], simplex[2]);
      return -dir;
    }
    return dir;
  }

  if(ori_cda && Triangle_3(O+a,O+c,O+d).has_on(Plane_3(O+a,O+c,O+d).projection(O))){
    Vector_3 dir = cross(ac, ad);
    simplex[1]=simplex[3];
    simplex.pop_back();
    if(is_negative(dot(dir, -a))){
      std::swap(simplex[1], simplex[2]);
      return -dir;
    }
    return dir;
  }

  if(ori_bdc && Triangle_3(O+b,O+d,O+c).has_on(Plane_3(O+b,O+d,O+c).projection(O))){
    Vector_3 dir = cross(d-b, d-c);
    simplex[0]=simplex[3];
    simplex.pop_back();
    if(is_negative(dot(dir, -b))){
      std::swap(simplex[1], simplex[2]);
      return -dir;
    }
    return dir;
  }


  if(ori_adb){
    simplex[2]=simplex[3];
    simplex.pop_back();
    return triangle_dir_to_origin<K>(simplex);
  }
  if(ori_bdc){
    simplex[0]=simplex[3];
    std::swap(simplex[0],simplex[2]);
    simplex.pop_back();
    return triangle_dir_to_origin<K>(simplex);
  }
  if(ori_cda){
    simplex[1]=simplex[3];
    std::swap(simplex[1],simplex[2]);
    simplex.pop_back();
    return triangle_dir_to_origin<K>(simplex);
  }

  assert(0);
}

template<typename K, typename Vector_3>
Vector_3 dir_to_origin(boost::container::small_vector<Vector_3, 4>& simplex) {
  using Triangle_3=typename K::Triangle_3;
  using Plane_3=typename K::Plane_3;
  using FT=typename K::FT;

  auto cross=K().construct_cross_product_vector_3_object();
  auto dot=K().compute_scalar_product_3_object();
  auto orientation = K().orientation_3_object();

  typename K::Point_3 O(0,0,0);

  if(simplex.size()==0)
    return Vector_3(1,0,0);

  if(simplex.size()==1)
    return -simplex[0];

  if(simplex.size()==2){

    Vector_3 &a = simplex[0];
    Vector_3 &b = simplex[1];
    Vector_3 ab=b-a;
    FT ratio=dot(ab, -a);
    if(is_negative(ratio)){
      simplex.pop_back();
      return -a;
    }
    if(ratio>ab.squared_length()){
      simplex[0]=simplex[1];
      simplex.pop_back();
      return -b;
    }
    Vector_3 dir= -cross(cross(ab, a), ab);
    return dir;
  }

  if(simplex.size()==3){
    return triangle_dir_to_origin<K>(simplex);
  }

  return tetrahedron_dir_to_origin<K>(simplex);

}

template<typename OK, typename K>
struct GJK_do_intersect{

  template<typename Convex, typename NamedParameters1, typename NamedParameters2>
  typename OK::FT operator()(const Convex &a, const Convex &b, const NamedParameters1 &np1, const NamedParameters2 &np2) const{
    using Point_3 = typename K::Point_3;
    using Vector_3= typename K::Vector_3;
    using Triangle_3= typename K::Triangle_3;
    using Plane_3= typename K::Plane_3;
    using FT = typename K::FT;
    const int INTER_MAX_ITER=0;
    boost::container::small_vector<Vector_3, 4> simplex;

    Cartesian_converter<K, OK> converter;

    auto cp=K().construct_cross_product_vector_3_object();
    auto dot=K().compute_scalar_product_3_object();

    FT min=100;
    Point_3 origin=Point_3(0,0,0);

    unsigned long planeStatPerPair = 0;
    do {
      // std::cout << "\nstep " << planeStatPerPair << std::endl;
      // for(auto v: simplex)
      //   std::cout << v << std::endl;

      Vector_3 dir=dir_to_origin<K>(simplex);
      if(dir==NULL_VECTOR) return 0;

      Vector_3 sp = extreme_point_3(a, dir.direction(), np1) - extreme_point_3(b, -dir.direction(), np2);
      if(sp==NULL_VECTOR) return 0;
      if(INTER_MAX_ITER!=0 && (++planeStatPerPair >= INTER_MAX_ITER)) return 0;

      bool closest_simplex=false;
      for(size_t i=0; i<simplex.size(); ++i)
        if(simplex[i]==sp)
          closest_simplex=true;

      if(closest_simplex)
      {
        FT dist;
        if(simplex.size()>=3)
          dist= squared_distance(Plane_3(origin+simplex[0],origin+simplex[1],origin+simplex[2]), origin);
        else if(simplex.size()==2)
          dist= squared_distance(typename K::Line_3(origin+simplex[0],origin+simplex[1]), origin);
        else
          dist= squared_distance(origin+simplex[0], origin);
        return converter(dist);
      }

      if(simplex.size()==4) simplex.pop_back();
      simplex.push_back(sp);

      assert(++planeStatPerPair<=30);
    } while( true );
    return true;
  }
};

} // end of predicates_impl namespace


#if DOXYGEN_RUNNING
/**
* \ingroup PkgConvexHull3Predicates
*
* provides a lower bound on the squared distance between the convex hulls of the two point sets.
*
* @tparam PointRange: is a model of `ConstRange`. The value type of its iterator is the key type of the named parameter `point_map`.
* @tparam NamedParameters_1 a sequence of \ref bgl_namedparameters "Named Parameters"
* @tparam NamedParameters_2 a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param r1 first point range
* @param r2 second point range
* @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
* @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of `r1` (`r2`)}
*     \cgalParamType{a model of `ReadablePropertyMap` whose value types are the same for `np1` and `np2`}
*     \cgalParamDefault{`CGAL::Identity_property_map`}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{An instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal kernel deduced from the point type, using `CGAL::Kernel_traits`}
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
*/
template <class PointRange,
          class NamedParameters_1 = parameters::Default_named_parameters,
          class NamedParameters_2 = parameters::Default_named_parameters>
FT separation_distance(const PointRange& r1, const PointRange& r2,
                       const NamedParameters_1& np1 = parameters::default_values(),
                       const NamedParameters_2& np2 = parameters::default_values());

/**
* \ingroup PkgConvexHull3Predicates
*
* provides a lower bound on the squared distance between the two convex graphs.
*
* @tparam AdjacencyGraph: is a model of `AdjacencyGraph`.
* @tparam NamedParameters_1 a sequence of \ref bgl_namedparameters "Named Parameters"
* @tparam NamedParameters_2 a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param g1 the first convex graph
* @param g2 the second convex graph
* @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
* @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* @warning The input graph must represent a convex object to guarantee a correct answer.
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `g1` (`g2`)}
*     \cgalParamType{a model of `ReadablePropertyMap` whose value types are the same for `np1` and `np2`}
*     \cgalParamDefault{boost::get(CGAL::vertex_point, g)}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t` must be available in ` AdjacencyGraph`.}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{An instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{`np1` only}
*   \cgalParamNEnd
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{if not `0` (no limit), indicates the maximum number of iterations performed by the algorithm.
                            if this value is not `0`, then the return value can be zero even if the convex hulls does not intersect.
*                           However, the value reported remains a lower bound of the distance between the convex.}
*     \cgalParamType{a positive integer convertible to `std::size_t`}
*     \cgalParamExtra{`np1` only}
*     \cgalParamDefault{`0`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
*/
template <class AdjacencyGraph,
          class NamedParameters_1 = parameters::Default_named_parameters,
          class NamedParameters_2 = parameters::Default_named_parameters>
FT separation_distance(const AdjacencyGraph& g1, const AdjacencyGraph& g2,
                       const NamedParameters_1& np1 = parameters::default_values(),
                       const NamedParameters_2& np2 = parameters::default_values());

/**
* \ingroup PkgConvexHull3Predicates
*
* provides a lower bound on the squared distance between the two convex hulls.
*
* @tparam PolygonMesh: is a model of `MutableFaceGraph`, more details in `CGAL::Convex_hull_hierarchy`
* @tparam NamedParameters_1 a sequence of \ref bgl_namedparameters "Named Parameters"
* @tparam NamedParameters_2 a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param ch1 the first convex hull
* @param ch2 the second convex hull
* @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
* @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `ch1` (`ch2`)}
*     \cgalParamType{a model of `ReadablePropertyMap` whose value types are the same for `np1` and `np2`}
*     \cgalParamDefault{boost::get(CGAL::vertex_point, g)}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t` must be available in `IncidenceGraph`.}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{An instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{`np1` only}
*   \cgalParamNEnd
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{if not `0` (no limit), indicates the maximum number of iterations performed by the algorithm.
                            If this value is not `0`, then the return value can be zero even if the convex hulls does not intersect.
*                           However, the value reported remains a lower bound of the distance between the convex.}
*     \cgalParamType{a positive integer convertible to `std::size_t`}
*     \cgalParamExtra{`np1` only}
*     \cgalParamDefault{`0`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
*/
template <class PolygonMesh,
          class NamedParameters_1 = parameters::Default_named_parameters,
          class NamedParameters_2 = parameters::Default_named_parameters>
FT separation_distance(const Convex_hull_hierarchy<PolygonMesh>& ch1, const Convex_hull_hierarchy<PolygonMesh>& ch2,
                       const NamedParameters_1& np1 = parameters::default_values(),
                       const NamedParameters_2& np2 = parameters::default_values());

#else
/**
* \ingroup PkgConvexHull3Predicates
*
* provides a lower bound on the distance between two convex sets. Consider the convex hull of point sets provide.
*
* @tparam Convex1: can be a model of the concept `Container`, `IncidenceGraph`, `Convex_hull_hierarchy` or any object 'M'
* such that exists a function 'extreme_point_3(M, Kernel::Vector_3, Converter)' returning a Kernel::Point_3
* @tparam Convex2: same as Convex1
* @tparam NamedParameters_1 a sequence of \ref bgl_namedparameters "Named Parameters"
* @tparam NamedParameters_2 a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param c1 first convex set considered in the do-intersect test
* @param c2 second convex set considered in the do-intersect test
* @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
* @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of `c1` (`c2`)}
*     \cgalParamType{a model of `ReadablePropertyMap` whose value types are the same for `np1` and `np2`}
*     \cgalParamDefault{`CGAL::Identity_property_map`}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{An instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{`np1` only}
*   \cgalParamNEnd
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{if not `0` (no limit), indicates the maximum number of iterations performed by the algorithm.
*                           If this value is not `0`, then the return value can be zero even if the convex hulls does not intersect.
*                           However, the value reported remains a lower bound of the distance between the convex.}
*     \cgalParamType{a positive integer convertible to `std::size_t`}
*     \cgalParamExtra{`np1` only}
*     \cgalParamDefault{`0`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \returns If the convex intersect, the return value is zero.
*
*/
template <class Convex1, class Convex2,
          class NamedParameters_1 = parameters::Default_named_parameters,
          class NamedParameters_2 = parameters::Default_named_parameters>
typename Point_set_processing_3_np_helper<Convex1, NamedParameters_1>::Geom_traits::FT
separation_distance(const Convex1& c1, const Convex2& c2,
                    const NamedParameters_1& np1 = parameters::default_values(),
                    const NamedParameters_2& np2 = parameters::default_values()){
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  using EPECK=Exact_predicates_exact_constructions_kernel;

  if constexpr(is_instance_of_v<Convex1, Convex_hull_hierarchy>){
    using GetGeomTraits = GetGeomTraits<typename Convex1::Mesh, NamedParameters_1>;
    using GT= typename GetGeomTraits::type;
    GT gt = choose_parameter<GT>(get_parameter(np1, internal_np::geom_traits));
    return predicates_impl::GJK_do_intersect<GT, EPECK>()(c1, c2, np1, np2);
  } else if constexpr(CGAL::IO::internal::is_Range_v<Convex1>){
    using NP_helper= Point_set_processing_3_np_helper<Convex1, NamedParameters_1>;
    using GT= typename NP_helper::Geom_traits;
    GT gt = NP_helper::get_geom_traits(c1, np1);
    return predicates_impl::GJK_do_intersect<GT, EPECK>()(c1, c2, np1, np2);
  } else {
    using GetGeomTraits = GetGeomTraits<Convex1, NamedParameters_1>;
    using GT= typename GetGeomTraits::type;
    GT gt = choose_parameter<GT>(get_parameter(np1, internal_np::geom_traits));
    return predicates_impl::GJK_do_intersect<GT, EPECK>()(c1, c2, np1, np2);
  }
  return true;
}

#endif

}} // CGAL::Convex_hull_3 namespace

#endif // CGAL_CONVEX_HULL_3_DISTANCE_H
