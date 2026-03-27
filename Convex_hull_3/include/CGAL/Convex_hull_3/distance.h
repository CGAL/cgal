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

#include <CGAL/Convex_hull_hierarchy_3.h>
#include <CGAL/Convex_hull_3/internal/helpers.h>
#include <CGAL/squared_distance_3.h>

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
namespace experimental {

/*
Note: The function is not stable with EPICK (possible infinite loops), It needs to be rewritten by moving all constructions inside predicates and by defining objects (Vector_3, Direction_3,...) are defined directly from the input points.
For instance, Vector_3 should be defined by two points and Direction_3 by 3 vectors (i.e 6 points) if the simplex size is 3 or 4.
Currently the function simply switches to EPECK and is considered experimental.
*/

namespace predicates_impl {

// returns true iff pt is on the negative side of the plane defined by (ep0, ep1) and normal
template <class K>
inline bool
strictly_on_left_of_triangle_edge(const typename K::Vector_3& pt,
                                  const typename K::Vector_3& normal,
                                  const typename K::Vector_3& edge,
                                  const K& k)
{
  const typename K::Construct_cross_product_vector_3 cross=k.construct_cross_product_vector_3_object();
  const typename K::Compute_scalar_product_3 dot=k.compute_scalar_product_3_object();

  return is_negative(dot(cross(edge, normal), pt));
}

/*
Compute the Voronoi cell containing the origin and return the direction of the simplex closest to it.
The origin necessarily lies below the edge formed by the first two vertices of the triangle,
we exploit this fact to simplify the computation.
*/
template<typename K, typename Vector_3>
Vector_3 triangle_dir_to_origin(boost::container::small_vector<Vector_3, 4>& simplex, const K &k){
  using FT=typename K::FT;

  const typename K::Construct_cross_product_vector_3 cross=k.construct_cross_product_vector_3_object();
  const typename K::Compute_scalar_product_3 dot=k.compute_scalar_product_3_object();

  Vector_3 &a = simplex[0];
  Vector_3 &b = simplex[1];
  Vector_3 &c = simplex[2];
  Vector_3 ab = b-a;
  Vector_3 ca = a-c;
  Vector_3 bc = c-b;

  Vector_3 n=cross(ca, bc);
  assert(n!=NULL_VECTOR);

  // FT ratio_ab=dot(ab, -a)/ab.squared_length();
  FT ratio_ca=dot(ca, -c);
  FT ratio_cb=dot(-bc, -c);

  bool on_left_of_ca = strictly_on_left_of_triangle_edge(-c, n, ca, k);
  bool on_left_of_bc = strictly_on_left_of_triangle_edge(-b, n, bc, k);

  // The origin cannot be on the side of ab, we abuse of this
  if(/*!on_left_of_ab &&*/ !on_left_of_ca && !on_left_of_bc){
    if(is_negative(dot(n, -a))){
      return -n;
    }
    std::swap(simplex[1], simplex[2]);
    return n;
  }

  if(on_left_of_ca && ratio_ca>=0 /*&& ratio_ca<=1*/){
    Vector_3 dir=-cross(cross(ca, c), ca);
    simplex[1] = simplex[2];
    simplex.pop_back();
    return dir;
  }
  if(on_left_of_bc && ratio_cb>=0 /*&& ratio_cb<=1*/){
    Vector_3 dir=-cross(cross(bc, b), bc);
    simplex[0] = simplex[2];
    simplex.pop_back();
    return dir;
  }

  simplex[0] = simplex[2];
  simplex.pop_back();
  simplex.pop_back();
  return -simplex[0];
}

/*
Compute the Voronoi cell containing the origin and return the direction of the simplex closest to it.
The origin necessarily lies below the triangle formed by the first three vertices of the tetrahedron,
we exploit this fact to simplify the computation.
*/
template<typename K, typename Vector_3>
Vector_3 tetrahedron_dir_to_origin(boost::container::small_vector<Vector_3, 4>& simplex, const K &k) {
  const typename K::Construct_cross_product_vector_3 cross=k.construct_cross_product_vector_3_object();
  const typename K::Compute_scalar_product_3 dot=k.compute_scalar_product_3_object();
  const typename K::Orientation_3 orientation = k.orientation_3_object();

  Vector_3 &a=simplex[0];
  Vector_3 &b=simplex[1];
  Vector_3 &c=simplex[2];
  Vector_3 &d=simplex[3];

  Vector_3 ad=d-a;
  Vector_3 bd=d-b;
  Vector_3 cd=d-c;

  assert(orientation(ad, bd, cd) == POSITIVE);

  assert(orientation(a,b,c) != POSITIVE); // The origin is below abc
  bool ori_adb=orientation(a,d,b) == POSITIVE;
  bool ori_bdc=orientation(b,d,c) == POSITIVE;
  bool ori_cda=orientation(c,d,a) == POSITIVE;

  //If the origin is inside the tetrahedron, no direction to return
  if(!ori_adb && !ori_bdc && !ori_cda)
    return NULL_VECTOR;

  //Compute the position of the origin with adb
  Vector_3 n_adb =cross(ad, bd);
  bool on_left_of_ad = strictly_on_left_of_triangle_edge(-d, n_adb, ad, k);
  bool on_right_of_bd = strictly_on_left_of_triangle_edge(-d, n_adb, -bd, k);

  if(ori_adb && !on_left_of_ad && !on_right_of_bd){
    // The origin is above the triangle adb
    simplex[2]=simplex[3];
    simplex.pop_back();
    return n_adb;
  }

  //Compute the position of the origin with cda
  Vector_3 n_cda =cross(cd, ad);
  bool on_left_of_cd = strictly_on_left_of_triangle_edge(-d, n_cda, cd, k);
  bool on_right_of_ad = strictly_on_left_of_triangle_edge(-d, n_cda, -ad, k);

  if(ori_cda && !on_left_of_cd && !on_right_of_ad){
    // The origin is above the triangle cda
    simplex[1]=simplex[3];
    simplex.pop_back();
    return n_cda;
  }

  //Compute the position of the origin with bdc
  Vector_3 n_bdc =cross(bd, cd);
  bool on_left_of_bd = strictly_on_left_of_triangle_edge(-d, n_bdc, bd, k);
  bool on_right_of_cd = strictly_on_left_of_triangle_edge(-d, n_bdc, -cd, k);

  if(ori_bdc && !on_left_of_bd && !on_right_of_cd){
    // The origin is above the triangle bdc
    simplex[0]=simplex[3];
    simplex.pop_back();
    return n_bdc;
  }

  // The origin is not above a face
  // Using the predicates already computed, we identify the edge below the origin
  if(on_left_of_ad && on_right_of_ad){
    if(is_negative(dot(ad, d))){
      //The origin is above the vertex d
      simplex[0]=simplex[3];
      simplex.pop_back();
      simplex.pop_back();
      simplex.pop_back();
      return -simplex[0];
    }
    Vector_3 dir=-cross(cross(ad, d), ad);
    simplex[1]=simplex[3];
    simplex.pop_back();
    simplex.pop_back();
    return dir;
  }

  if(on_left_of_bd && on_right_of_bd){
    if(is_negative(dot(bd, d))){
      //The origin is above the vertex d
      simplex[0]=simplex[3];
      simplex.pop_back();
      simplex.pop_back();
      simplex.pop_back();
      return -simplex[0];
    }
    Vector_3 dir=-cross(cross(bd, d), bd);
    simplex[0]=simplex[3];
    simplex.pop_back();
    simplex.pop_back();
    return dir;
  }

  assert(on_left_of_cd && on_right_of_cd);
  if(is_negative(dot(cd, d))){
    //The origin is above the vertex d
    simplex[0]=simplex[3];
    simplex.pop_back();
    simplex.pop_back();
    simplex.pop_back();
    return -simplex[0];
  }
  Vector_3 dir=-cross(cross(cd, d), cd);
  simplex[0]=simplex[2];
  simplex[1]=simplex[3];
  simplex.pop_back();
  simplex.pop_back();
  return dir;
}

template<typename K, typename Vector_3>
Vector_3 dir_to_origin(boost::container::small_vector<Vector_3, 4>& simplex, const K &k) {
  const typename K::Construct_cross_product_vector_3 cross=k.construct_cross_product_vector_3_object();
  const typename K::Compute_scalar_product_3 dot=k.compute_scalar_product_3_object();

  switch(simplex.size()){
  case 0:
    return Vector_3(1,0,0);
    break;
  case 1:
    return -simplex[0];
    break;
  case 2:
  {
    Vector_3 &a = simplex[0];
    Vector_3 &b = simplex[1];
    Vector_3 ab=b-a;
    typename K::FT ratio=dot(ab, -a);
    if(is_negative(ratio)){
      simplex.pop_back();
      return -a;
    }
    if(ratio>ab.squared_length()){
      Vector_3 dir=-b;
      simplex[0]=simplex[1];
      simplex.pop_back();
      return dir;
    }
    return -cross(cross(ab, a), ab);
    break;
  }
  case 3:
    return triangle_dir_to_origin(simplex, k);
    break;
  default: //4
    return tetrahedron_dir_to_origin(simplex, k);
    break;
  }
}

template<typename OK, typename K>
struct Separation_distance_functor{

  template<typename Convex, typename NamedParameters1, typename NamedParameters2>
  typename OK::FT operator()(const Convex &a, const Convex &b, const NamedParameters1 &np1, const NamedParameters2 &np2) const{
    using Vector_3= typename K::Vector_3;
    using FT = typename K::FT;

    const int INTER_MAX_ITER=0;
    boost::container::small_vector<Vector_3, 4> simplex;

    Cartesian_converter<K, OK> converter;

    //The Kernel used is EPECK and not the one of the named parameter
    const K &k=K();

    unsigned long planeStatPerPair = 0;
    do {
#ifdef CGAL_CONVEX_HULL_3_DISTANCE_VERBOSE
      std::cout << "\nIteration " << planeStatPerPair << std::endl;
      std::cout << "Simplex size: " << simplex.size() << ", points of the simplex:"
      for(const Vector_3 &v: simplex)
        std::cout << v << std::endl;
#endif

      Vector_3 dir=dir_to_origin<K>(simplex, k);
      if(dir==NULL_VECTOR) return 0;

      Vector_3 sp = Convex_hull_3::internal::extreme_point_3_wrapper(a, dir.direction(), np1) - Convex_hull_3::internal::extreme_point_3_wrapper(b, -dir.direction(), np2);
      if(sp==NULL_VECTOR) return 0;
      if(INTER_MAX_ITER!=0 && (++planeStatPerPair >= INTER_MAX_ITER)) return 0;

#ifdef CGAL_CONVEX_HULL_3_DISTANCE_VERBOSE
      std::cout << "Direction to origin: " << dir << std::endl;
      std::cout << "Support point: " << sp << std::endl;
#endif

      simplex.push_back(sp);
      //If the simplex is degenerated, the algorithm has terminate
      if( (simplex.size()==4 && orientation(ORIGIN+simplex[0],ORIGIN+simplex[1],ORIGIN+simplex[2],ORIGIN+simplex[3])==COPLANAR) ||
          (simplex.size()==3 && collinear(ORIGIN+simplex[0], ORIGIN+simplex[1], ORIGIN+simplex[2])) ||
          (simplex.size()==2 && (simplex[0]==simplex[1])) )
      {
        simplex.pop_back();
        break;
      }

      assert(++planeStatPerPair<=30);
    } while( true );
    FT dist;
    typename K::Point_3 ORIGIN(0,0,0);
    if(simplex.size()>=3)
      dist = squared_distance(typename K::Plane_3(ORIGIN+simplex[0], ORIGIN+simplex[1], ORIGIN+simplex[2]), ORIGIN);
    else if(simplex.size()==2)
      dist = squared_distance(typename K::Line_3(ORIGIN+simplex[0], ORIGIN+simplex[1]), ORIGIN);
    else
      dist = simplex[0].squared_length();
    return converter(dist);
  }
};

} // end of predicates_impl namespace


/**
* \ingroup PkgConvexHull3Intersections
*
* \brief compute the separation distance of two convex hulls returning zero if they intersect.
*
* Input hulls can be provided as a range of points or as a graph, and may be of different types.
* Furthermore, when many intersection queries use the same object, one should wrap that input in the class `CGAL::Convex_hull_hierarchy_3` as to construct an optimized view of the convex hull and accelerate intersection tests.
* This is especially true when the convex hull is made of a large number of vertices (see `CGAL::Convex_hull_hierarchy_3` for more details).
*
* @tparam Convex_1 is one of the following types:\n
*  - a model of `ConstRange`
*  - a model of `VertexListGraph` and `AdjacencyGraph`
*  - an instance of `CGAL::Convex_hull_hierarchy_3`
* @tparam Convex_2 same as `Convex_1`
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
*     \cgalParamExtra{used only if `ch1` (`ch2`) is model of `VertexListGraph` and `AdjacencyGraph`.}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t` must be available in `Convex_1` (`Convex_2`).}
*   \cgalParamNEnd
*
*   \cond SKIP_IN_MANUAL
*
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
*
* \endcond
*
* \cgalNamedParamsEnd
*
* \see `CGAL::Convex_hull_hierarchy_3`
*/
template <class Convex_1, class Convex_2,
          class NamedParameters_1 = parameters::Default_named_parameters,
          class NamedParameters_2 = parameters::Default_named_parameters>
typename internal::GetGeomTraitsFromConvex<Convex_1, NamedParameters_1>::type::FT
separation_distance(const Convex_1& c1, const Convex_2& c2,
                    const NamedParameters_1& np1 = parameters::default_values(),
                    const NamedParameters_2& np2 = parameters::default_values()){
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  //The function need exact computation to works correctly
  using EPECK=Exact_predicates_exact_constructions_kernel;
  using GT= typename internal::GetGeomTraitsFromConvex<Convex_1, NamedParameters_1>::type;
  // GT gt = choose_parameter<GT>(get_parameter(np1, internal_np::geom_traits));
  return predicates_impl::Separation_distance_functor<GT, EPECK>()(c1, c2, np1, np2);
}

} // end of experimental

}} // CGAL::Convex_hull_3 namespace

#endif // CGAL_CONVEX_HULL_3_DISTANCE_H
