// Copyright (c) 2001,2011  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//               : Amol Prakash <prakash@mpi-sb.mpg.de>
//               : Andreas Fabri

#ifndef CGAL_CONVEX_HULL_3_H
#define CGAL_CONVEX_HULL_3_H

#include <CGAL/license/Convex_hull_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Projection_traits_xz_3.h>
#include <CGAL/Projection_traits_yz_3.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/Convex_hull_2/ch_assertions.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/internal/Exact_type_selector.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_data_structure_2.h>
#include <CGAL/boost/graph/properties_Triangulation_data_structure_2.h>
#include <CGAL/Polyhedron_3_fwd.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/is_iterator.h>

#include <boost/bind.hpp>
#include <boost/next_prior.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <boost/graph/graph_traits.hpp>

#ifndef CGAL_CH_NO_POSTCONDITIONS
#include <CGAL/convexity_check_3.h>
#endif // CGAL_CH_NO_POSTCONDITIONS

#include <algorithm>
#include <iostream>
#include <list>
#include <memory>
#include <vector>
#include <type_traits>
#include <utility>

// first some internal stuff to avoid using a true Face_graph model for extreme_points_3
namespace CGAL {

// Forward declaration
template<class VertexPointMap,class Base_traits> class Extreme_points_traits_adapter_3;

namespace Convex_hull_3 {
namespace internal {

// wrapper used as a MutableFaceGraph to extract extreme points
template <class OutputIterator>
struct Output_iterator_wrapper
{
  OutputIterator out;
  Output_iterator_wrapper(OutputIterator out)
    : out(out)
  {}
};

template <class Point_3, class OutputIterator>
void add_isolated_points(const Point_3& point,
                         Convex_hull_3::internal::Output_iterator_wrapper<OutputIterator>& w)
{
    *w.out++ = point;
}

template <class Point_3, class PolygonMesh>
void add_isolated_points(const Point_3& point, PolygonMesh& P)
{
  put(get(CGAL::vertex_point, P), add_vertex(P), point);
}

template <class Point_3, class OutputIterator>
void copy_ch2_to_face_graph(const std::list<Point_3>& CH_2,
                            Convex_hull_3::internal::Output_iterator_wrapper<OutputIterator>& w)
{
  for(const Point_3& p : CH_2)
    *w.out++ = p;
}

} // namespace internal
} // namespace Convex_hull_3

template <class TDS, class OutputIterator>
void copy_face_graph(const TDS& tds, Convex_hull_3::internal::Output_iterator_wrapper<OutputIterator>& wrapper)
{
  typedef typename boost::graph_traits<TDS>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<TDS, boost::vertex_point_t >::const_type vpm = get(boost::vertex_point, tds);
  for(vertex_descriptor vh : vertices(tds))
  {
    *wrapper.out++ = get(vpm, vh);
  }
}

template <class OutputIterator>
void clear(Convex_hull_3::internal::Output_iterator_wrapper<OutputIterator>&)
{}

template <class Point, class OutputIterator>
void make_tetrahedron(const Point& p0, const Point& p1, const Point& p2, const Point& p3,
                      Convex_hull_3::internal::Output_iterator_wrapper<OutputIterator>& w)
{
  *w.out++ = p0;
  *w.out++ = p1;
  *w.out++ = p2;
  *w.out++ = p3;
}

} // CGAL

namespace boost {

// needed so that the regular make_tetrahedron of CGAL does not complain when tried to be instantiated
template <class OutputIterator>
struct graph_traits<CGAL::Convex_hull_3::internal::Output_iterator_wrapper<OutputIterator> >
{
  typedef void* halfedge_descriptor;
};

} // namespace boost

namespace CGAL {
namespace Convex_hull_3 {
namespace internal {

//struct to select the default traits class for computing convex hull
template< class Point_3,
          class PolygonMesh = Default,
          class Is_floating_point=typename boost::is_floating_point<typename Kernel_traits<Point_3>::Kernel::FT>::type,
          class Has_filtered_predicates_tag=typename Kernel_traits<Point_3>::Kernel::Has_filtered_predicates_tag >
struct Default_traits_for_Chull_3{
  typedef typename Kernel_traits<Point_3>::Kernel type;
};

//FT is a floating point type and Kernel is a filtered kernel
template <class Point_3, class PolygonMesh>
struct Default_traits_for_Chull_3<Point_3, PolygonMesh, boost::true_type,Tag_true>{
  typedef Convex_hull_traits_3< typename Kernel_traits<Point_3>::Kernel, PolygonMesh, Tag_true > type;
};

template <class Traits>
struct Default_polyhedron_for_Chull_3{
  typedef CGAL::Polyhedron_3<Traits> type;
};

template <class K, class P, class Tag>
struct Default_polyhedron_for_Chull_3<Convex_hull_traits_3<K, P, Tag> >{
  typedef typename  Convex_hull_traits_3<K, P, Tag>::Polygon_mesh type;
};

template <class T>
struct Is_cartesian_kernel
{
  typedef boost::false_type type;
};

template <class Kernel, class PolygonMesh>
struct Is_cartesian_kernel< Convex_hull_traits_3<Kernel, PolygonMesh, Tag_true> >
{
  // Rational here is that Tag_true can only be passed by us since it is not documented
  // so we can assume that Kernel is a CGAL Kernel
  typedef typename boost::is_same<typename Kernel::Kernel_tag, Cartesian_tag>::type type;
};

// Predicate internally used as a wrapper around has_on_positive_side
// We provide a partial specialization restricted to the case of CGAL Cartesian Kernels with inexact constructions below
//template <class Traits,class Tag_use_advanced_filtering=typename Use_advanced_filtering<Traits>::type >
template <class Traits, class Is_CK = typename Is_cartesian_kernel<Traits>::type >
class Is_on_positive_side_of_plane_3{
  typedef typename Traits::Point_3 Point_3;
  typename Traits::Plane_3 plane;
  typename Traits::Has_on_positive_side_3 has_on_positive_side;
public:
  typedef Protect_FPU_rounding<false> Protector;

  Is_on_positive_side_of_plane_3(const Traits& traits,const Point_3& p,const Point_3& q,const Point_3& r)
  :plane(traits.construct_plane_3_object()(p,q,r)),has_on_positive_side(traits.has_on_positive_side_3_object()) {}

  bool operator() (const Point_3& s) const
  {
    return has_on_positive_side(plane,s);
  }
};

template <class Base_traits, class VPM, class Is_CK>
class Is_on_positive_side_of_plane_3< Extreme_points_traits_adapter_3<VPM,Base_traits>, Is_CK>
  : public Is_on_positive_side_of_plane_3< Base_traits >
{
  typedef Extreme_points_traits_adapter_3<VPM, Base_traits> Traits;
  typedef Is_on_positive_side_of_plane_3< Base_traits > Base;
  typedef typename Traits::Point_3 Point_3;
  const Traits& m_traits;
public:
  typedef typename Base::Protector Protector;

  Is_on_positive_side_of_plane_3(const Traits& traits,
                                 const Point_3& p,const Point_3& q,const Point_3& r)
    : Base(static_cast<const Base_traits&>(traits), traits.get_point(p), traits.get_point(q), traits.get_point(r))
    , m_traits(traits)
  {}

  bool operator() (const Point_3& s) const
  {
    return static_cast<const Base*>(this)->operator()(m_traits.get_point(s));
  }
};

//This predicate uses copy of the code from the statically filtered version of
//Orientation_3. The rational is that the plane is a member of the functor
//so optimization are done to avoid doing several time operations on the plane.
//The main operator() first tries the static version of the predicate, then uses
//interval arithmetic (the protector must be created before using this predicate)
//and in case of failure, exact arithmetic is used.
template <class Kernel, class P>
class Is_on_positive_side_of_plane_3<Convex_hull_traits_3<Kernel, P, Tag_true>, boost::true_type >{
  typedef Simple_cartesian<CGAL::internal::Exact_field_selector<double>::Type>  Exact_K;
  typedef Simple_cartesian<Interval_nt_advanced >                               Approx_K;
  typedef Convex_hull_traits_3<Kernel, P, Tag_true>                             Traits;
  typedef typename Traits::Point_3                                              Point_3;

  Cartesian_converter<Kernel,Approx_K>                  to_AK;
  Cartesian_converter<Kernel,Exact_K>                   to_EK;

  template <typename K>
  struct Vector_plus_point {
    typename K::Vector_3 vector;
    typename K::Point_3  point;
  };

  const Point_3& p,q,r;
  mutable Vector_plus_point<Approx_K> ak_plane;
  mutable Vector_plus_point<Exact_K>* ek_plane_ptr;

  double m10,m20,m21,Maxx,Maxy,Maxz;

  static const int STATIC_FILTER_FAILURE = 555;

  //this function is a made from the statically filtered version of Orientation_3
  int static_filtered(double psx,double psy, double psz) const{

    // Then semi-static filter.
    double apsx = CGAL::abs(psx);
    double apsy = CGAL::abs(psy);
    double apsz = CGAL::abs(psz);

    double maxx = (Maxx < apsx)? apsx : Maxx;
    double maxy = (Maxy < apsy)? apsy : Maxy;
    double maxz = (Maxz < apsz)? apsz : Maxz;

    double det =  psx*m10 - m20*psy + m21*psz;

    // Sort maxx < maxy < maxz.
    if (maxx > maxz)
        std::swap(maxx, maxz);
    if (maxy > maxz)
        std::swap(maxy, maxz);
    else if (maxy < maxx)
        std::swap(maxx, maxy);

    // Protect against underflow in the computation of eps.
    if (maxx < 1e-97) /* cbrt(min_double/eps) */ {
      if (maxx == 0)
        return 0;
    }
    // Protect against overflow in the computation of det.
    else if (maxz < 1e102) /* cbrt(max_double [hadamard]/4) */ {
      double eps = 5.1107127829973299e-15 * maxx * maxy * maxz;
      if (det > eps)  return 1;
      if (det < -eps) return -1;
    }
    return STATIC_FILTER_FAILURE;
  }

public:
  typedef typename Interval_nt_advanced::Protector           Protector;

  Is_on_positive_side_of_plane_3(const Traits&,const Point_3& p_,const Point_3& q_,const Point_3& r_)
    : p(p_),q(q_),r(r_)
    , ak_plane()
    , ek_plane_ptr(nullptr)
  {
    ak_plane.vector =
      typename Approx_K::Vector_3(Interval_nt_advanced(0., std::numeric_limits<double>::infinity()),
                                  0., 0.);
    double pqx = q.x() - p.x();
    double pqy = q.y() - p.y();
    double pqz = q.z() - p.z();
    double prx = r.x() - p.x();
    double pry = r.y() - p.y();
    double prz = r.z() - p.z();


    m10 = pqy*prz - pry*pqz;
    m20 = pqx*prz - prx*pqz;
    m21 = pqx*pry - prx*pqy;

    double aprx = CGAL::abs(prx);
    double apry = CGAL::abs(pry);
    double aprz = CGAL::abs(prz);

    Maxx = CGAL::abs(pqx);
    if (Maxx < aprx) Maxx = aprx;
    Maxy = CGAL::abs(pqy);
    if (Maxy < apry) Maxy = apry;
    Maxz = CGAL::abs(pqz);
    if (Maxz < aprz) Maxz = aprz;
  }

  ~Is_on_positive_side_of_plane_3(){
    if (ek_plane_ptr!=nullptr) delete ek_plane_ptr;
  }

  bool operator() (const Point_3& s) const
  {
    double psx = s.x() - p.x();
    double psy = s.y() - p.y();
    double psz = s.z() - p.z();

    int static_res = static_filtered(psx,psy,psz);
    if (static_res != STATIC_FILTER_FAILURE)
      return static_res == 1;

    try{
      // infinity() is the sentinel for uninitialized `ak_plane`
      if (ak_plane.vector.x().sup() == std::numeric_limits<double>::infinity())
      {
        const typename Approx_K::Point_3 ap = to_AK(p);
        ak_plane.vector = cross_product(to_AK(q)-ap, to_AK(r)-ap);
        ak_plane.point = ap;
      }
      Uncertain<Sign> res =
        sign(scalar_product(to_AK(s) - ak_plane.point,
                            ak_plane.vector));
      if(is_certain(res)) {
        return (get_certain(res) == POSITIVE);
      }
    }
    catch (Uncertain_conversion_exception&){}
    if (ek_plane_ptr==nullptr) {
      const typename Exact_K::Point_3 ep = to_EK(p);
      ek_plane_ptr = new Vector_plus_point<Exact_K>;
      ek_plane_ptr->vector = cross_product(to_EK(q)-ep, to_EK(r)-ep);
      ek_plane_ptr->point = ep;
    }
    return sign(scalar_product(to_EK(s) - ek_plane_ptr->point,
                               ek_plane_ptr->vector)) == POSITIVE;
  }
};


BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Traits_has_typedef_Traits_xy_3,Traits_xy_3,false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Traits_has_typedef_Traits_yz_3,Traits_xy_3,false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Traits_has_typedef_Traits_xz_3,Traits_xy_3,false)

template <class T,bool has_projection_traits=
  Traits_has_typedef_Traits_xy_3<T>::value &&
  Traits_has_typedef_Traits_yz_3<T>::value &&
  Traits_has_typedef_Traits_xz_3<T>::value
>
struct Projection_traits{
    Projection_traits(const T&){}
  typedef typename Kernel_traits<typename T::Point_3>::Kernel K;
  typedef CGAL::Projection_traits_xy_3<K> Traits_xy_3;
  typedef CGAL::Projection_traits_yz_3<K> Traits_yz_3;
  typedef CGAL::Projection_traits_xz_3<K> Traits_xz_3;

    Traits_xy_3 construct_traits_xy_3_object()const
    {return Traits_xy_3();}
    Traits_yz_3 construct_traits_yz_3_object()const
    {return Traits_yz_3();}
    Traits_xz_3 construct_traits_xz_3_object()const
    {return Traits_xz_3();}
};

template <class T>
struct Projection_traits<T,true>{
  const T& traits;
  Projection_traits(const T& t):traits(t){}
  typedef typename T::Traits_xy_3 Traits_xy_3;
  typedef typename T::Traits_yz_3 Traits_yz_3;
  typedef typename T::Traits_xz_3 Traits_xz_3;
  Traits_xy_3 construct_traits_xy_3_object()const
  {return traits.construct_traits_xy_3_object();}
  Traits_yz_3 construct_traits_yz_3_object()const
  {return traits.construct_traits_yz_3_object();}
  Traits_xz_3 construct_traits_xz_3_object()const
  {return traits.construct_traits_xz_3_object();}
};

template <class Point_3, class PolygonMesh>
void copy_ch2_to_face_graph(const std::list<Point_3>& CH_2, PolygonMesh& P)
{
  typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type vpm
    = get(CGAL::vertex_point, P);
  typedef boost::graph_traits<PolygonMesh> Graph_traits;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  std::vector<vertex_descriptor> vertices;
  vertices.reserve(CH_2.size());
  for(const Point_3& p : CH_2){
    vertices.push_back(add_vertex(P));
    put(vpm, vertices.back(),p);
  }
  face_descriptor f = Euler::add_face(vertices, P);

  // Then triangulate that face
  const halfedge_descriptor he = halfedge(f, P);
  halfedge_descriptor other_he = next(next(he, P), P);
  for(std::size_t i = 3, end = vertices.size(); i < end; ++i) {
    const halfedge_descriptor next_he = next(other_he, P);
    Euler::split_face(other_he, he, P);
    other_he = next_he;
  }
}


template <class InputIterator, class Point_3, class PolygonMesh, class Traits>
void coplanar_3_hull(InputIterator first, InputIterator beyond,
                     const Point_3& p1, const Point_3& p2, const Point_3& p3,
                     PolygonMesh& P, const Traits&  traits )
{
  typedef typename Convex_hull_3::internal::Projection_traits<Traits> PTraits;
  typedef typename PTraits::Traits_xy_3 Traits_xy_3;
  typedef typename PTraits::Traits_yz_3 Traits_yz_3;
  typedef typename PTraits::Traits_xz_3 Traits_xz_3;

  PTraits ptraits(traits);

  std::list<Point_3> CH_2;

  Traits_xy_3 traits_xy = ptraits.construct_traits_xy_3_object();
  typename Traits_xy_3::Left_turn_2 left_turn_in_xy = traits_xy.left_turn_2_object();
  if ( left_turn_in_xy(p1,p2,p3) || left_turn_in_xy(p2,p1,p3) )
     convex_hull_points_2( first, beyond,
                           std::back_inserter(CH_2),
                           traits_xy );
  else{
    Traits_yz_3 traits_yz = ptraits.construct_traits_yz_3_object();
    typename Traits_yz_3::Left_turn_2 left_turn_in_yz = traits_yz.left_turn_2_object();
    if ( left_turn_in_yz(p1,p2,p3) || left_turn_in_yz(p2,p1,p3) )
       convex_hull_points_2( first, beyond,
                             std::back_inserter(CH_2),
                             traits_yz );
    else{
      Traits_xz_3 traits_xz = ptraits.construct_traits_xz_3_object();
      CGAL_assertion_code( typename Traits_xz_3::Left_turn_2 left_turn_in_xz = traits_xz.left_turn_2_object(); )
      CGAL_assertion( left_turn_in_xz(p1,p2,p3) || left_turn_in_xz(p2,p1,p3) );
         convex_hull_points_2( first, beyond,
                               std::back_inserter(CH_2),
                               traits_xz );
    }
  }
  copy_ch2_to_face_graph(CH_2, P);
}


//
// visible is the set of facets visible from point  and reachable from
// start_facet.
//
template <class TDS_2, class Traits>
void
find_visible_set(TDS_2& tds,
                 const typename Traits::Point_3& point,
                 typename TDS_2::Face_handle start,
                 std::list<typename TDS_2::Face_handle>& visible,
                 std::map<typename TDS_2::Vertex_handle, typename TDS_2::Edge>& outside,
                 const Traits& traits)
{
   typedef typename TDS_2::Face_handle Face_handle;
   typedef typename TDS_2::Vertex_handle Vertex_handle;

   std::vector<Vertex_handle> vertices;
   vertices.reserve(10);
   int VISITED=1, BORDER=2;
   visible.clear();
   typename std::list<Face_handle>::iterator  vis_it;
   visible.push_back(start);
   start->info() = VISITED;
   vertices.push_back(start->vertex(0));
   vertices.push_back(start->vertex(1));
   vertices.push_back(start->vertex(2));
   start->vertex(0)->info() = start->vertex(1)->info() = start->vertex(2)->info() = VISITED;

   for (vis_it = visible.begin(); vis_it != visible.end(); vis_it++)
   {
      // check all the neighbors of the current face to see if they have
      // already been visited or not and if not whether they are visible
      // or not.

      for(int i=0; i < 3; i++) {
        // the facet on the other side of the current halfedge
        Face_handle f = (*vis_it)->neighbor(i);
        // if haven't already seen this facet
        if (f->info() == 0) {
          f->info() = VISITED;
          Is_on_positive_side_of_plane_3<Traits> is_on_positive_side(
            traits,f->vertex(0)->point(),f->vertex(2)->point(),f->vertex(1)->point());
          int ind = f->index(*vis_it);
          if ( !is_on_positive_side(point) ){  // is visible
            visible.push_back(f);
            Vertex_handle vh = f->vertex(ind);
            if(vh->info() == 0){ vertices.push_back(vh); vh->info() = VISITED;}
          } else {
            f->info() = BORDER;
            f->vertex(TDS_2::cw(ind))->info() = BORDER;
            f->vertex(TDS_2::ccw(ind))->info() = BORDER;
            outside.insert(std::make_pair(f->vertex(TDS_2::cw(ind)),
                                          typename TDS_2::Edge(f,ind)));
          }
        } else if(f->info() == BORDER) {
          int ind = f->index(*vis_it);
          f->vertex(TDS_2::cw(ind))->info() = BORDER;
          f->vertex(TDS_2::ccw(ind))->info() = BORDER;
          outside.insert(std::make_pair(f->vertex(TDS_2::cw(ind)),
                                        typename TDS_2::Edge(f,ind)));
        }
      }
   }

   for(typename std::vector<Vertex_handle>::iterator vit =  vertices.begin();
       vit != vertices.end();
       ++vit){
     if((*vit)->info() != BORDER){
       tds.delete_vertex(*vit);
     } else {
       (*vit)->info() = 0;
     }
   }

}

// using a third template parameter for the point instead of getting it from
// the traits class as it should be is required by M$VC6
template <class Face_handle, class Traits, class Point>
typename std::list<Point>::iterator
farthest_outside_point(Face_handle f, std::list<Point>& outside_set,
                       const Traits& traits)
{

   typedef typename std::list<Point>::iterator Outside_set_iterator;
   CGAL_ch_assertion(!outside_set.empty());

   typename Traits::Plane_3 plane =
       traits.construct_plane_3_object()(f->vertex(0)->point(),
                                         f->vertex(1)->point(),
                                         f->vertex(2)->point());

   typename Traits::Less_signed_distance_to_plane_3 less_dist_to_plane =
            traits.less_signed_distance_to_plane_3_object();
   Outside_set_iterator farthest_it =
          std::max_element(outside_set.begin(),
                           outside_set.end(),
                           boost::bind(less_dist_to_plane, plane, _1, _2));
   return farthest_it;
}

template <class Face_handle, class Traits, class Point>
void
partition_outside_sets(const std::list<Face_handle>& new_facets,
                       std::list<Point>& vis_outside_set,
                       std::list<Face_handle>& pending_facets,
                       const Traits& traits)
{
  typename std::list<Face_handle>::const_iterator        f_list_it;
  typename std::list<Point>::iterator  point_it, to_splice;

  // walk through all the new facets and check each unassigned outside point
  // to see if it belongs to the outside set of this new facet.
  for (f_list_it = new_facets.begin(); (f_list_it != new_facets.end()) && (! vis_outside_set.empty());
        ++f_list_it)
  {
    Face_handle f = *f_list_it;
    Is_on_positive_side_of_plane_3<Traits> is_on_positive_side(
      traits,f->vertex(0)->point(),f->vertex(1)->point(),f->vertex(2)->point());
    std::list<Point>& point_list = f->points;

    for (point_it = vis_outside_set.begin();point_it != vis_outside_set.end();){
      if( is_on_positive_side(*point_it) ) {
        to_splice = point_it;
        ++point_it;
        point_list.splice(point_list.end(), vis_outside_set, to_splice);
      } else {
         ++point_it;
      }
    }
    if(! point_list.empty()){
      pending_facets.push_back(f);
      f->it = boost::prior(pending_facets.end());
    } else {
      f->it = pending_facets.end();
    }
  }


   for (; f_list_it != new_facets.end();++f_list_it)
    (*f_list_it)->it = pending_facets.end();
}



template <class TDS_2, class Traits>
void
ch_quickhull_3_scan(TDS_2& tds,
                    std::list<typename TDS_2::Face_handle>& pending_facets,
                    const Traits& traits)
{
  typedef typename TDS_2::Edge                            Edge;
  typedef typename TDS_2::Face_handle                     Face_handle;
  typedef typename TDS_2::Vertex_handle                   Vertex_handle;
  typedef typename Traits::Point_3                          Point_3;
  typedef std::list<Point_3>                              Outside_set;
  typedef typename std::list<Point_3>::iterator           Outside_set_iterator;
  typedef std::map<typename TDS_2::Vertex_handle, typename TDS_2::Edge> Border_edges;

  std::list<Face_handle>                     visible_set;
  typename std::list<Face_handle>::iterator  vis_set_it;
  Outside_set                                vis_outside_set;
  Border_edges                               border;

  while (!pending_facets.empty())
  {
     vis_outside_set.clear();

     Face_handle f_handle = pending_facets.front();

     Outside_set_iterator farthest_pt_it = farthest_outside_point(f_handle, f_handle->points, traits);
     Point_3 farthest_pt = *farthest_pt_it;
     f_handle->points.erase(farthest_pt_it);
     find_visible_set(tds, farthest_pt, f_handle, visible_set, border, traits);

     // for each visible facet
     for (vis_set_it = visible_set.begin(); vis_set_it != visible_set.end();
          vis_set_it++)
     {

        //   add its outside set to the global outside set list
       std::list<Point_3>& point_list = (*vis_set_it)->points;
       if(! point_list.empty()){
         vis_outside_set.splice(vis_outside_set.end(), point_list, point_list.begin(), point_list.end());
       }

       if((*vis_set_it)->it != pending_facets.end()){
         pending_facets.erase((*vis_set_it)->it);
       }
       (*vis_set_it)->info() = 0;
     }

     std::vector<Edge> edges;
     edges.reserve(border.size());
     typename Border_edges::iterator it = border.begin();
     Edge e = it->second;
     e.first->info() = 0;
     edges.push_back(e);
     border.erase(it);
     while(! border.empty()){
       it = border.find(e.first->vertex(TDS_2::ccw(e.second)));
       CGAL_ch_assertion(it != border.end());
       e = it->second;
       e.first->info() = 0;
       edges.push_back(e);
       border.erase(it);
     }

     // If we want to reuse the faces we must only pass |edges| many, and call delete_face for the others.
     // Also create facets if necessary
     std::ptrdiff_t diff = visible_set.size() - edges.size();
     if(diff < 0){
       for(int i = 0; i<-diff;i++){
         visible_set.push_back(tds.create_face());
       }
     } else {
       for(int i = 0; i<diff;i++){
         tds.delete_face(visible_set.back());
         visible_set.pop_back();
       }
     }
     Vertex_handle vh = tds.star_hole(edges.begin(), edges.end(), visible_set.begin(), visible_set.end());
     vh->point() = farthest_pt;
     vh->info() = 0;

     // now partition the set of outside set points among the new facets.

     partition_outside_sets(visible_set, vis_outside_set,
                            pending_facets, traits);

  }
}

template <class TDS_2, class Traits>
void non_coplanar_quickhull_3(std::list<typename Traits::Point_3>& points,
                              TDS_2& tds, const Traits& traits)
{
  typedef typename Traits::Point_3                        Point_3;

  typedef typename TDS_2::Face_handle                     Face_handle;
  typedef typename TDS_2::Face_iterator                     Face_iterator;
  typedef typename std::list<Point_3>::iterator           P3_iterator;

  std::list<Face_handle> pending_facets;

  typename Is_on_positive_side_of_plane_3<Traits>::Protector p;

  // for each facet, look at each unassigned point and decide if it belongs
  // to the outside set of this facet.
  for(Face_iterator fit = tds.faces_begin(); fit != tds.faces_end(); ++fit){
    Is_on_positive_side_of_plane_3<Traits> is_on_positive_side(
      traits,fit->vertex(0)->point(),fit->vertex(1)->point(),fit->vertex(2)->point() );
    for (P3_iterator point_it = points.begin() ; point_it != points.end(); )
    {
      if( is_on_positive_side(*point_it) ) {
        P3_iterator to_splice = point_it;
        ++point_it;
        fit->points.splice(fit->points.end(), points, to_splice);
      } else {
       ++point_it;
      }
    }
  }
  // add all the facets with non-empty outside sets to the set of facets for
  // further consideration
  for(Face_iterator fit = tds.faces_begin(); fit != tds.faces_end(); ++fit){
    if (! fit->points.empty()){
      pending_facets.push_back(fit);
      fit->it = boost::prior(pending_facets.end());
        } else {
      fit->it =  pending_facets.end();
    }
  }


  ch_quickhull_3_scan(tds, pending_facets, traits);

  //std::cout << "|V(tds)| = " << tds.number_of_vertices() << std::endl;
//  CGAL_ch_expensive_postcondition(all_points_inside(points.begin(),
//                                                    points.end(),P,traits));
//  CGAL_ch_postcondition(is_strongly_convex_3(P, traits));
}

template <class InputIterator, class PolygonMesh, class Traits>
void
ch_quickhull_face_graph(std::list<typename Traits::Point_3>& points,
                        InputIterator point1_it, InputIterator point2_it, InputIterator point3_it,
                        PolygonMesh& P,
                        const Traits& traits)
{
  typedef typename Traits::Point_3                            Point_3;
  typedef typename Traits::Plane_3                                Plane_3;
  typedef typename std::list<Point_3>::iterator           P3_iterator;

  typedef Triangulation_data_structure_2<
    Triangulation_vertex_base_with_info_2<int, GT3_for_CH3<Traits> >,
    Convex_hull_face_base_2<int, Traits> >                           Tds;

  typedef typename Tds::Vertex_handle                     Vertex_handle;
  typedef typename Tds::Face_handle                     Face_handle;

  // found three points that are not collinear, so construct the plane defined
  // by these points and then find a point that has maximum distance from this
  // plane.
  typename Traits::Construct_plane_3 construct_plane =
         traits.construct_plane_3_object();
  Plane_3 plane = construct_plane(*point3_it, *point2_it, *point1_it);
  typedef typename Traits::Less_signed_distance_to_plane_3      Dist_compare;
  Dist_compare compare_dist = traits.less_signed_distance_to_plane_3_object();

  typename Traits::Coplanar_3  coplanar = traits.coplanar_3_object();
  // find both min and max here since using signed distance.  If all points
  // are on the negative side of the plane, the max element will be on the
  // plane.
  std::pair<P3_iterator, P3_iterator> min_max;
  min_max = CGAL::min_max_element(points.begin(), points.end(),
                                  boost::bind(compare_dist, plane, _1, _2),
                                  boost::bind(compare_dist, plane, _1, _2));
  P3_iterator max_it;
  if (coplanar(*point1_it, *point2_it, *point3_it, *min_max.second))
  {
     max_it = min_max.first;
     // want the orientation of the points defining the plane to be positive
     // so have to reorder these points if all points were on negative side
     // of plane
     std::swap(*point1_it, *point3_it);
  }
  else
     max_it = min_max.second;

  // if the maximum distance point is on the plane then all are coplanar
  if (coplanar(*point1_it, *point2_it, *point3_it, *max_it)) {
     coplanar_3_hull(points.begin(), points.end(), *point1_it, *point2_it, *point3_it, P, traits);
  } else {
    Tds tds;
    Vertex_handle v0 = tds.create_vertex(); v0->set_point(*point1_it);
    Vertex_handle v1 = tds.create_vertex(); v1->set_point(*point2_it);
    Vertex_handle v2 = tds.create_vertex(); v2->set_point(*point3_it);
    Vertex_handle v3 = tds.create_vertex(); v3->set_point(*max_it);

    v0->info() = v1->info() = v2->info() = v3->info() = 0;
    Face_handle f0 = tds.create_face(v0,v1,v2);
    Face_handle f1 = tds.create_face(v3,v1,v0);
    Face_handle f2 = tds.create_face(v3,v2,v1);
    Face_handle f3 = tds.create_face(v3,v0,v2);
    tds.set_dimension(2);
    v0->set_face(f0);
    v1->set_face(f0);
    v2->set_face(f0);
    v3->set_face(f1);
    f0->set_neighbors(f2, f3, f1);
    f1->set_neighbors(f0, f3, f2);
    f2->set_neighbors(f0, f1, f3);
    f3->set_neighbors(f0, f2, f1);

    points.erase(point1_it);
    points.erase(point2_it);
    points.erase(point3_it);
    points.erase(max_it);
    if (!points.empty()){
      non_coplanar_quickhull_3(points, tds, traits);
      copy_face_graph(tds,P);
    }
    else{
      CGAL_assertion( traits.has_on_positive_side_3_object()(
            construct_plane(v2->point(),v1->point(),v0->point()),
            v3->point()) );
      make_tetrahedron(v0->point(),v1->point(),v3->point(),v2->point(),P);
    }
  }

}

} // namespace internal
} // namespace Convex_hull_3

template <class InputIterator, class Traits>
void
convex_hull_3(InputIterator first, InputIterator beyond,
              Object& ch_object,
              const Traits& traits)
{
  typedef typename Traits::Point_3                            Point_3;
  typedef std::list<Point_3>                              Point_3_list;
  typedef typename Point_3_list::iterator                 P3_iterator;
  typedef std::pair<P3_iterator,P3_iterator>              P3_iterator_pair;

  if (first == beyond)    // No point
    return;

  // If the first and last point are equal the collinearity test some lines below will always be true.
  Point_3_list points(first, beyond);
  std::size_t size = points.size();
  while((size > 1) && (points.front() == points.back())){
    points.pop_back();
    --size;
  }

  typename Traits::Collinear_3 collinear = traits.collinear_3_object();

  if ( size == 1 )                // 1 point
  {
      ch_object = make_object(*points.begin());
      return;
  }
  else if ( size == 2 )           // 2 points
  {
      typedef typename Traits::Segment_3                 Segment_3;
      typename Traits::Construct_segment_3 construct_segment =
             traits.construct_segment_3_object();
      Segment_3 seg = construct_segment(*points.begin(), *(++points.begin()));
      ch_object = make_object(seg);
      return;
  }
  else if ( ( size == 3 ) && (! collinear(*(points.begin()),
                                          *(++points.begin()),
                                          *(--points.end()) ) ) )           // 3 points
  {
      typedef typename Traits::Triangle_3                Triangle_3;
      typename Traits::Construct_triangle_3 construct_triangle =
             traits.construct_triangle_3_object();
      Triangle_3 tri = construct_triangle(*(points.begin()),
                                          *(++points.begin()),
                                          *(--points.end()));
      ch_object = make_object(tri);
      return;
  }

  // at least 4 points

  P3_iterator point1_it = points.begin();
  P3_iterator point2_it = points.begin();
  point2_it++;
  P3_iterator point3_it = points.end();
  point3_it--;

  // find three that are not collinear
  while (point2_it != points.end() &&
         collinear(*point1_it,*point2_it,*point3_it))
    point2_it++;


  // all are collinear, so the answer is a segment
  if (point2_it == points.end())
  {
     typedef typename Traits::Less_distance_to_point_3      Less_dist;

     Less_dist less_dist = traits.less_distance_to_point_3_object();
     P3_iterator_pair endpoints =
      min_max_element(points.begin(), points.end(),
                      boost::bind(less_dist, *points.begin(), _1, _2),
                      boost::bind(less_dist, *points.begin(), _1, _2));

     typename Traits::Construct_segment_3 construct_segment =
            traits.construct_segment_3_object();
     typedef typename Traits::Segment_3                 Segment_3;

     Segment_3 seg = construct_segment(*endpoints.first, *endpoints.second);
     ch_object = make_object(seg);
     return;
  }

  // result will be a polyhedron
  typedef typename Convex_hull_3::internal::Default_polyhedron_for_Chull_3<Traits>::type Polyhedron;
  Polyhedron P;

  P3_iterator minx, maxx, miny, it;
  minx = maxx = miny = it = points.begin();
  ++it;
  for(; it != points.end(); ++it){
    if(it->x() < minx->x()) minx = it;
    if(it->x() > maxx->x()) maxx = it;
    if(it->y() < miny->y()) miny = it;
  }
  if(! collinear(*minx, *maxx, *miny) ){
    Convex_hull_3::internal::ch_quickhull_face_graph(points, minx, maxx, miny, P, traits);
  } else {
    Convex_hull_3::internal::ch_quickhull_face_graph(points, point1_it, point2_it, point3_it, P, traits);
  }
  CGAL_assertion(num_vertices(P)>=3);
  typename boost::graph_traits<Polyhedron>::vertex_iterator b,e;
  boost::tie(b,e) = vertices(P);
  if (num_vertices(P) == 3){
    typename boost::property_map<Polyhedron, vertex_point_t>::type vpmap  = get(CGAL::vertex_point, P);
    typedef typename Traits::Triangle_3                Triangle_3;
    typename Traits::Construct_triangle_3 construct_triangle =
           traits.construct_triangle_3_object();
    typedef typename Traits::Point_3 Point_3;
    Point_3 p = get(vpmap, *b); ++b;
    Point_3 q = get(vpmap, *b); ++b;
    Point_3 r = get(vpmap, *b);
    Triangle_3 tri = construct_triangle(p,q,r);
    ch_object = make_object(tri);
  }
  else
    ch_object = make_object(P);
}


template <class InputIterator>
void convex_hull_3(InputIterator first, InputIterator beyond,
                   Object& ch_object)
{
   typedef typename std::iterator_traits<InputIterator>::value_type Point_3;
   typedef typename Convex_hull_3::internal::Default_traits_for_Chull_3<Point_3>::type Traits;
   convex_hull_3(first, beyond, ch_object, Traits());
}

template <class InputIterator, class PolygonMesh, class Traits>
void convex_hull_3(InputIterator first, InputIterator beyond,
                   PolygonMesh& polyhedron,
                   const Traits& traits)
{
  typedef typename Traits::Point_3                Point_3;
  typedef std::list<Point_3>                      Point_3_list;
  typedef typename Point_3_list::iterator         P3_iterator;
  typedef std::pair<P3_iterator,P3_iterator>      P3_iterator_pair;

  if(first == beyond)
    return;

  Point_3_list points(first, beyond);
  clear(polyhedron);

  typename Traits::Collinear_3 collinear = traits.collinear_3_object();
  typename Traits::Equal_3 equal = traits.equal_3_object();
  P3_iterator point1_it = points.begin();
  P3_iterator point2_it = points.begin();
  point2_it++;

  while (point2_it != points.end() && equal(*point1_it,*point2_it))
    ++point2_it;

  // if there is only one point or all points are equal
  if(point2_it == points.end())
  {
    Convex_hull_3::internal::add_isolated_points(*point1_it, polyhedron);
    return;
  }

  P3_iterator point3_it = point2_it;
  ++point3_it;

  while (point3_it != points.end() && collinear(*point1_it,*point2_it,*point3_it))
    ++point3_it;

  if (point3_it == points.end())
  {
    // only 2 points or all points are collinear
    typedef typename Traits::Less_distance_to_point_3 Less_dist;
    Less_dist less_dist = traits.less_distance_to_point_3_object();
    P3_iterator_pair endpoints =
    min_max_element(points.begin(), points.end(),
                    boost::bind(less_dist, *points.begin(), _1, _2),
                    boost::bind(less_dist, *points.begin(), _1, _2));
    Convex_hull_3::internal::add_isolated_points(*endpoints.first, polyhedron);
    Convex_hull_3::internal::add_isolated_points(*endpoints.second, polyhedron);
    return;
  }

  Convex_hull_3::internal::ch_quickhull_face_graph(points, point1_it, point2_it, point3_it,
    polyhedron, traits);
}

template <class InputIterator, class PolygonMesh>
void convex_hull_3(InputIterator first, InputIterator beyond,
                   PolygonMesh& polyhedron,
                   // workaround to avoid ambiguity with next overload.
                   typename std::enable_if<CGAL::is_iterator<InputIterator>::value>::type* = 0)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Point_3;
  typedef typename Convex_hull_3::internal::Default_traits_for_Chull_3<Point_3, PolygonMesh>::type Traits;
  convex_hull_3(first, beyond, polyhedron, Traits());
}

template <class VertexListGraph, class PolygonMesh, class NamedParameters>
void convex_hull_3(const VertexListGraph& g,
                   PolygonMesh& pm,
                   const NamedParameters& np)
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  typedef typename GetVertexPointMap<VertexListGraph, NamedParameters>::const_type Vpmap;
  typedef CGAL::Property_map_to_unary_function<Vpmap> Vpmap_fct;
  Vpmap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(boost::vertex_point, g));

  Vpmap_fct v2p(vpm);
  convex_hull_3(boost::make_transform_iterator(vertices(g).begin(), v2p),
                boost::make_transform_iterator(vertices(g).end(), v2p), pm);
}

template <class VertexListGraph, class PolygonMesh>
void convex_hull_3(const VertexListGraph& g,
                   PolygonMesh& pm)
{
  convex_hull_3(g,pm,CGAL::parameters::all_default());
}

template <class InputRange, class OutputIterator, class Traits>
OutputIterator
extreme_points_3(const InputRange& range,
                 OutputIterator out,
                 const Traits& traits)
{
  Convex_hull_3::internal::Output_iterator_wrapper<OutputIterator> wrapper(out);
  convex_hull_3(range.begin(), range.end(), wrapper, traits);
  return out;
}

template <class InputRange, class OutputIterator>
OutputIterator
extreme_points_3(const InputRange& range, OutputIterator out)
{
  typedef typename InputRange::const_iterator Iterator_type;
  typedef typename std::iterator_traits<Iterator_type>::value_type Point_3;
  typedef typename Convex_hull_3::internal::Default_traits_for_Chull_3<Point_3>::type Traits;

  return extreme_points_3(range, out, Traits());
}

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_CONVEX_HULL_3_H
