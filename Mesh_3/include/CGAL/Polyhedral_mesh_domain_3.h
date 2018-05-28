// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// Copyright (c) 2011 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_POLYHEDRAL_MESH_DOMAIN_3_H
#define CGAL_POLYHEDRAL_MESH_DOMAIN_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/boost/parameter.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>

#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <sstream>

#include <CGAL/Default.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Mesh_3/Profile_counter.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>

#include <boost/optional.hpp>
#include <boost/none.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/mpl/or.hpp>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/tuple.h>
#include <boost/format.hpp>
#include <boost/variant.hpp>
#include <boost/math/special_functions/round.hpp>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/enumerable_thread_specific.h>
#endif

// To handle I/O for Surface_patch_index if that is a pair of `int` (the
// default)
#include <CGAL/internal/Mesh_3/Handle_IO_for_pair_of_int.h>

namespace CGAL {

namespace Mesh_3 {

namespace details {

inline
double
max_length(const Bbox_3& b)
{
  return (std::max)(b.xmax()-b.xmin(),
                    (std::max)(b.ymax()-b.ymin(),b.zmax()-b.zmin()) );
}

// -----------------------------------
// Index_generator
// Don't use boost::variant if types are the same type
// -----------------------------------
template < typename Subdomain_index, typename Surface_patch_index >
struct Index_generator
{
  typedef boost::variant<Subdomain_index,Surface_patch_index> Index;
  typedef Index                                         type;
};

template < typename T >
struct Index_generator<T, T>
{
  typedef T       Index;
  typedef Index   type;
};

// -----------------------------------
// Geometric traits generator
// -----------------------------------
template < typename Gt,
           typename Use_exact_intersection_construction_tag >
struct IGT_generator {};

template < typename Gt >
struct IGT_generator<Gt,CGAL::Tag_true>
{
#ifdef CGAL_MESH_3_NEW_ROBUST_INTERSECTION_TRAITS
  typedef CGAL::Mesh_3::Robust_intersection_traits_3_new<Gt> type;
#else // NOT CGAL_MESH_3_NEW_ROBUST_INTERSECTION_TRAITS
  typedef CGAL::Mesh_3::Robust_intersection_traits_3<Gt> type;
#endif // NOT CGAL_MESH_3_NEW_ROBUST_INTERSECTION_TRAITS
  typedef type Type;
};

template < typename Gt >
struct IGT_generator<Gt,CGAL::Tag_false>
{
  typedef Gt type;
  typedef type Type;
};

}  // end namespace details

}  // end namespace Mesh_3


namespace internal { namespace Mesh_3 {

template <typename Polyhedron_type,
          bool = CGAL::graph_has_property<Polyhedron_type,
                                           CGAL::face_index_t>::value>
class Get_face_index_pmap {
public:
  typedef typename boost::property_map<Polyhedron_type,
                                       CGAL::face_index_t>::const_type Pmap;
  Get_face_index_pmap(const Polyhedron_type&) {}
  Pmap operator()(const Polyhedron_type& polyhedron) {
    return get(CGAL::face_index, polyhedron);
  }
};

template <typename Polyhedron_type>
class Get_face_index_pmap<Polyhedron_type, false> {
  typedef typename boost::graph_traits<Polyhedron_type>::face_descriptor
                                                              face_descriptor;
  typedef std::map<face_descriptor, int> Map;
public:
  Get_face_index_pmap(const Polyhedron_type& polyhedron) {
    int id = 0;
    BOOST_FOREACH(face_descriptor f, faces(polyhedron))
    {
      face_ids[f] = id++;
    }
  }
  typedef boost::associative_property_map<Map> Pmap;

  Pmap operator()(const Polyhedron_type&) {
    return Pmap(face_ids);
  }
private:
  Map face_ids;
};

} // end namespace Mesh_3
} // end namespace internal

/**
 * @class Polyhedral_mesh_domain_3
 *
 *
 */
template<class Polyhedron,/*FaceGraph*/
         class IGT_,
         class = CGAL::Default,
         class Patch_id_ = void,
         class Use_exact_intersection_construction_tag = CGAL::Tag_true>
class Polyhedral_mesh_domain_3
{
public:
  typedef typename Mesh_3::details::IGT_generator<
    IGT_,Use_exact_intersection_construction_tag>::type IGT;


  typedef Patch_id_ Patch_id;

  /// Geometric object types
  typedef typename IGT::Point_3    Point_3;
  typedef typename IGT::Segment_3  Segment_3;
  typedef typename IGT::Ray_3      Ray_3;
  typedef typename IGT::Line_3     Line_3;
  typedef typename IGT::Vector_3   Vector_3;
  typedef typename IGT::Sphere_3   Sphere_3;

  //-------------------------------------------------------
  // Index Types
  //-------------------------------------------------------
  /// Type of indexes for cells of the input complex
  typedef int Subdomain_index;
  typedef boost::optional<Subdomain_index> Subdomain;

  /// Type of indexes for surface patch of the input complex
  typedef typename boost::property_map<Polyhedron,
                                       face_patch_id_t<Patch_id>
                                       >::type            Face_patch_id_pmap;
  typedef typename boost::property_traits<
    Face_patch_id_pmap>::value_type                       Surface_patch_index;
  typedef boost::optional<Surface_patch_index>            Surface_patch;
  /// Type of indexes to characterize the lowest dimensional face of the input
  /// complex on which a vertex lie
  typedef typename Mesh_3::details::Index_generator<
    Subdomain_index, Surface_patch_index>::type           Index;

  typedef CGAL::cpp11::tuple<Point_3,Index,int> Intersection;


  typedef typename IGT::FT         FT;

  // Kernel_traits compatibility
  typedef IGT R;

  BOOST_MPL_HAS_XXX_TRAIT_DEF(HalfedgeDS)

  template <typename P>
  struct Primitive {
    typedef typename boost::graph_traits<P>::face_descriptor face_descriptor_t;

    face_descriptor_t face_descriptor;
    const P* graph;

    typedef Triangle_from_face_descriptor_map<P>   Triangle_pmap;
    typedef One_point_from_face_descriptor_map<P>  Point_pmap;

    typedef typename Triangle_pmap::value_type              Datum;
    typedef typename Point_pmap::value_type                 Point;
    typedef Primitive                                       Id;

    Primitive() : face_descriptor(), graph() {}

    template <typename IndexIterator>
    Primitive(IndexIterator it, const P& polyhedron)
      : face_descriptor(*it), graph(&polyhedron)
    {}

    bool operator==(const Primitive& other) const {
      return other.face_descriptor == this->face_descriptor
        && other.graph == this->graph;
    }

    Datum datum() const {
      return get(Triangle_pmap(graph), face_descriptor);
    }

    Point reference_point() const {
      return get(Point_pmap(graph), face_descriptor);
    }

    Id id() const {
      return *this;
    }
  }; // Primitive template, for face graphs

  template <typename P, bool is_a_CGAL_Polyhedron_3 = has_HalfedgeDS<P>::value>
  struct Primitive_type {
    typedef Primitive<P> type;

    static
    typename IGT_::Triangle_3 datum(const typename type::Id primitive_id) {
      CGAL::Triangle_from_face_descriptor_map<P> pmap(primitive_id.graph);
      return get(pmap, primitive_id.face_descriptor);
    }

    static Surface_patch_index get_index(const typename type::Id primitive_id) {
      return get(get(face_patch_id_t<Patch_id>(),
                     *primitive_id.graph),
                 primitive_id.face_descriptor);
    }
  }; // Primitive_type (for non-Polyhedron_3)

  template <typename P> struct Primitive_type<P, true> {
    typedef AABB_face_graph_triangle_primitive<P > type;

    static
    typename IGT_::Triangle_3 datum(const typename type::Id face_handle) {
      typedef typename IGT_::Point_3 Point;
      const Point& a = face_handle->halfedge()->vertex()->point();
      const Point& b = face_handle->halfedge()->next()->vertex()->point();
      const Point& c = face_handle->halfedge()->next()->next()->vertex()->point();
      return typename IGT_::Triangle_3(a,b,c);
    }

    static Surface_patch_index get_index(const typename type::Id face_handle,
                                         Tag_false)
    {
      typename boost::property_map<P, face_patch_id_t<Patch_id> >::type pmap;
      return get(pmap, face_handle);
    }

    static Surface_patch_index get_index(const typename type::Id,
                                         Tag_true)
    {
      return Surface_patch_index(0,1);
    }

    static Surface_patch_index get_index(const typename type::Id face_handle)
    {
      namespace m = boost::mpl;
      return get_index(face_handle,
                       Boolean_tag<m::or_<boost::is_same<Patch_id, void>,
                                          boost::is_same<Patch_id, Tag_false>
                       >::value>());
    }
  }; // Primitive_type specialized for CGAL::Polyehdron_3

public:
  typedef typename Primitive_type<Polyhedron>::type       Ins_fctor_primitive;
  typedef CGAL::AABB_traits<IGT, Ins_fctor_primitive>     Ins_fctor_traits;
  typedef CGAL::AABB_tree<Ins_fctor_traits>               Ins_fctor_AABB_tree;

  typedef Side_of_triangle_mesh<Polyhedron,
                                IGT,
                                Default,
                                Ins_fctor_AABB_tree>      Inside_functor;
  typedef typename Inside_functor::AABB_tree              AABB_tree_;
  BOOST_STATIC_ASSERT((boost::is_same<AABB_tree_, Ins_fctor_AABB_tree>::value));
  typedef typename AABB_tree_::AABB_traits                AABB_traits;
  typedef typename AABB_tree_::Primitive                  AABB_primitive;
  typedef typename AABB_tree_::Primitive_id               AABB_primitive_id;
  typedef typename AABB_traits::Bounding_box              Bounding_box;

public:

  /// Default constructor
  Polyhedral_mesh_domain_3(CGAL::Random* p_rng = NULL)
    : tree_()
    , bounding_tree_(&tree_)
    , p_rng_(p_rng)
  {
  }

  /**
   * @brief Constructor. Contruction from a polyhedral surface
   * @param polyhedron the polyhedron describing the polyhedral surface
   */
  Polyhedral_mesh_domain_3(const Polyhedron& p,
                           CGAL::Random* p_rng = NULL)
    : tree_()
    , bounding_tree_(&tree_) // the bounding tree is tree_
    , p_rng_(p_rng)
  {
    this->add_primitives(p);
    if(! is_triangle_mesh(p)) {
      std::cerr << "Your input polyhedron must be triangulated!\n";
      CGAL_error_msg("Your input polyhedron must be triangulated!");
    }
    this->build();
  }

  Polyhedral_mesh_domain_3(const Polyhedron& p,
                           const Polyhedron& bounding_polyhedron,
                           CGAL::Random* p_rng = NULL)
    : tree_()
    , bounding_tree_(new AABB_tree_)
    , p_rng_(p_rng)
  {
    this->add_primitives(p);
    this->add_primitives(bounding_polyhedron);
    if(!bounding_polyhedron.empty()) {
      this->add_primitives_to_bounding_tree(bounding_polyhedron);
    } else {
      this->set_surface_only();
    }
    this->build();
  }

  /**
   * Constructor.
   *
   * Constructor from a sequence of polyhedral surfaces, and a bounding
   * polyhedral surface.
   *
   * @param InputPolyhedraPtrIterator must an iterator of a sequence of
   * pointers to polyhedra
   *
   * @param bounding_polyhedron reference to the bounding surface
   */
  template <typename InputPolyhedraPtrIterator>
  Polyhedral_mesh_domain_3(InputPolyhedraPtrIterator begin,
                           InputPolyhedraPtrIterator end,
                           const Polyhedron& bounding_polyhedron,
                           CGAL::Random* p_rng = NULL)
    : p_rng_(p_rng)
    , delete_rng_(false)
  {
    if(begin != end) {
      for(; begin != end; ++begin) {
        this->add_primitives(**begin);
      }
      this->add_primitives(bounding_polyhedron);
    }
    if(!bounding_polyhedron.empty()) {
      this->add_primitives_to_bounding_tree(bounding_polyhedron);
    } else {
      this->set_surface_only();
    }
    this->build();
  }

  /**
   * Constructor.
   *
   * Constructor from a sequence of polyhedral surfaces, without bounding
   * surface. The domain will always answer false to "is_in_domain"
   * queries.
   *
   * @param InputPolyhedraPtrIterator must an iterator of a sequence of
   * pointers to polyhedra
   */
  template <typename InputPolyhedraPtrIterator>
  Polyhedral_mesh_domain_3(InputPolyhedraPtrIterator begin,
                           InputPolyhedraPtrIterator end,
                           CGAL::Random* p_rng = NULL)
    : p_rng_(p_rng)
  {
    if(begin != end) {
      for(; begin != end; ++begin) {
        this->add_primitives(**begin);
      }
      tree_.build();
    }
    bounding_tree_ = 0;
  }

  /// Destructor
  ~Polyhedral_mesh_domain_3() {
    if(bounding_tree_ != 0 && bounding_tree_ != &tree_) {
      delete bounding_tree_;
    }
  }

  void set_surface_only() {
    bounding_tree_ = 0;
  }

  /**
   * Constructs  a set of \ccc{n} points on the surface, and output them to
   *  the output iterator \ccc{pts} whose value type is required to be
   *  \ccc{std::pair<Points_3, Index>}.
   */
  struct Construct_initial_points
  {
    Construct_initial_points(const Polyhedral_mesh_domain_3& domain)
      : r_domain_(domain) {}

    template<class OutputIterator>
    OutputIterator operator()(OutputIterator pts, const int n = 8) const;

  private:
    const Polyhedral_mesh_domain_3& r_domain_;
  };

  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }


  /**
   * Returns a bounding box of the domain
   */
  Bbox_3 bbox() const {
    return tree_.bbox();
  }


  /**
   * Returns true if point~\ccc{p} is in the domain. If \ccc{p} is in the
   *  domain, the parameter index is set to the index of the subdomain
   *  including $p$. It is set to the default value otherwise.
   */
  struct Is_in_domain
  {
    Is_in_domain(const Polyhedral_mesh_domain_3& domain)
      : r_domain_(domain) {}

    Subdomain operator()(const Point_3& p) const;
  private:
    const Polyhedral_mesh_domain_3& r_domain_;
  };

  Is_in_domain is_in_domain_object() const { return Is_in_domain(*this); }

  Point_3 project_on_surface(const Point_3& p) const
  {
    return tree_.closest_point(p);
  }

  /// Allowed query types
  typedef boost::mpl::vector<Segment_3, Ray_3, Line_3> Allowed_query_types;

  /**
   * Returns true is the element \ccc{type} intersect properly any of the
   * surface patches describing the either the domain boundary or some
   * subdomain boundary.
   * \ccc{Type} is either \ccc{Segment_3}, \ccc{Ray_3} or \ccc{Line_3}.
   * Parameter index is set to the index of the intersected surface patch
   * if \ccc{true} is returned and to the default \ccc{Surface_patch_index}
   * value otherwise.
   */
  struct Do_intersect_surface
  {
    Do_intersect_surface(const Polyhedral_mesh_domain_3& domain)
      : r_domain_(domain) {}

    template <typename Query>
    typename boost::enable_if<typename boost::mpl::contains<Allowed_query_types,
                                                            Query>::type,
                              Surface_patch>::type
    operator()(const Query& q) const
    {
      CGAL_MESH_3_PROFILER(std::string("Mesh_3 profiler: ") + std::string(CGAL_PRETTY_FUNCTION));

      boost::optional<AABB_primitive_id> primitive_id = r_domain_.tree_.any_intersected_primitive(q);
      if ( primitive_id )
      {
        r_domain_.cache_primitive(q, *primitive_id);
        return Surface_patch(r_domain_.make_surface_index(*primitive_id));
      } else {
        return Surface_patch();
      }
    }

  private:
    const Polyhedral_mesh_domain_3& r_domain_;
  };

  Do_intersect_surface do_intersect_surface_object() const
  {
    return Do_intersect_surface(*this);
  }

  /**
   * Returns a point in the intersection of the primitive \ccc{type}
   * with some boundary surface.
   * \ccc{Type1} is either \ccc{Segment_3}, \ccc{Ray_3} or \ccc{Line_3}.
   * The integer \ccc{dimension} is set to the dimension of the lowest
   * dimensional face in the input complex containing the returned point, and
   * \ccc{index} is set to the index to be stored at a mesh vertex lying
   * on this face.
   */
  struct Construct_intersection
  {
    Construct_intersection(const Polyhedral_mesh_domain_3& domain)
      : r_domain_(domain) {}

    template <typename Query>
    typename boost::enable_if<typename boost::mpl::contains<Allowed_query_types,
                                                            Query>::type,
                              Intersection>::type
    operator()(const Query& q) const
    {
      CGAL_MESH_3_PROFILER(std::string("Mesh_3 profiler: ") + std::string(CGAL_PRETTY_FUNCTION));
      typedef typename AABB_tree_::template Intersection_and_primitive_id<Query>::Type
        Intersection_and_primitive_id;
      typedef boost::optional<Intersection_and_primitive_id> AABB_intersection;
      typedef Point_3 Bare_point;

      AABB_intersection intersection;
#ifndef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
      if(r_domain_.query_is_cached(q))
      {
        const AABB_primitive_id primitive_id = r_domain_.cached_primitive_id();
        typename cpp11::result_of<
          typename IGT::Intersect_3(typename Primitive::Datum, Query)>::type o
            = IGT().intersect_3_object()(Primitive(primitive_id).datum(),q);
        intersection = o ?
          Intersection_and_primitive_id(*o, primitive_id) :
          AABB_intersection();
      } else
#endif // not CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
      {
#ifndef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
        CGAL_precondition(r_domain_.do_intersect_surface_object()(q)
                          != boost::none);
#endif // NOT CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3

        intersection = r_domain_.tree_.any_intersection(q);
      }
      if ( intersection )
      {
        // Get primitive
        AABB_primitive_id primitive_id = intersection->second;

        // intersection may be either a point or a segment
#if CGAL_INTERSECTION_VERSION > 1
        if ( const Bare_point* p_intersect_pt =
             boost::get<Bare_point>( &(intersection->first) ) )
#else
        if ( const Bare_point* p_intersect_pt =
             object_cast<Bare_point>( &(intersection->first) ) )
#endif
        {
          return Intersection(*p_intersect_pt,
                              r_domain_.index_from_surface_patch_index(
                                r_domain_.make_surface_index(primitive_id)),
                              2);
        }
#if CGAL_INTERSECTION_VERSION > 1
        else if ( const Segment_3* p_intersect_seg =
                  boost::get<Segment_3>(&(intersection->first)))
#else
        else if ( const Segment_3* p_intersect_seg =
                  object_cast<Segment_3>(&(intersection->first)))
#endif
        {
          CGAL_MESH_3_PROFILER("Mesh_3 profiler: Intersection is a segment");

          return Intersection(p_intersect_seg->source(),
                              r_domain_.index_from_surface_patch_index(
                                r_domain_.make_surface_index(primitive_id)),
                              2);
        }
        else {
#ifndef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
          std::stringstream stream;
          stream.precision(17);
          set_pretty_mode(stream);
          stream <<
            "Mesh_3 error : AABB_tree any_intersection result is "
            "not a point nor a segment\n";
          if(intersection->first.empty()) {
            stream <<  "The intersection is empty!";
          } else {
            stream <<  "The intersection typeinfo name is ";
            stream <<  intersection->first.type().name();
          }
          stream << "\nThe query was: ";
          stream << q << std::endl;
          stream << "The intersecting primitive in the AABB tree was: "
                 << AABB_primitive(intersection->second).datum() << std::endl;
          CGAL_error_msg(stream.str().c_str());
#endif // not CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
        }
      }

      // Should not happen
      // unless CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3 is defined
      return Intersection();
    }

  private:
    const Polyhedral_mesh_domain_3& r_domain_;
  };

  Construct_intersection construct_intersection_object() const
  {
    return Construct_intersection(*this);
  }


  /**
   * Returns the index to be stored in a vertex lying on the surface identified
   * by \c index.
   */
  Index index_from_surface_patch_index(const Surface_patch_index& index) const
  { return Index(index); }

  /**
   * Returns the index to be stored in a vertex lying in the subdomain
   * identified by \c index.
   */
  Index index_from_subdomain_index(const Subdomain_index& index) const
  { return Index(index); }

  /**
   * Returns the \c Surface_patch_index of the surface patch
   * where lies a vertex with dimension 2 and index \c index.
   */
  Surface_patch_index surface_patch_index(const Index& index) const
  { return boost::get<Surface_patch_index>(index); }

  /**
   * Returns the index of the subdomain containing a vertex
   *  with dimension 3 and index \c index.
   */
  Subdomain_index subdomain_index(const Index& index) const
  { return boost::get<Subdomain_index>(index); }

  // -----------------------------------
  // Backward Compatibility
  // -----------------------------------
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index   Surface_index;

  Index index_from_surface_index(const Surface_index& index) const
  { return index_from_surface_patch_index(index); }

  Surface_index surface_index(const Index& index) const
  { return surface_patch_index(index); }
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  // -----------------------------------
  // End backward Compatibility
  // -----------------------------------

public:
  Surface_patch_index make_surface_index(
    const AABB_primitive_id& primitive_id = AABB_primitive_id() ) const
  {
    return Primitive_type<Polyhedron>::get_index(primitive_id);
  }

  // Undocumented function, used to implement a sizing field that
  // computes lfs using this AABB tree. That avoids to rebuild the same
  // tree.
  typedef AABB_tree_ AABB_tree;
  const AABB_tree& aabb_tree() const {
    return tree_;
  }

  const AABB_tree* bounding_aabb_tree_ptr() const {
    return bounding_tree_;
  }

protected:
  void add_primitives(const Polyhedron& p)
  {
    tree_.insert(faces(p).begin(), faces(p).end(), p);
  }

  void add_primitives_to_bounding_tree(const Polyhedron& p)
  {
    if(bounding_tree_ == &tree_ || bounding_tree_ == 0) {
      bounding_tree_ = new AABB_tree_;
    }
    bounding_tree_->insert(faces(p).begin(), faces(p).end(), p);
  }

  void build() {
    tree_.build();
    CGAL_assertion(!tree_.empty());
    if(bounding_tree_ != &tree_ && bounding_tree_ != 0) {
      bounding_tree_->build();
      CGAL_assertion(!bounding_tree_->empty());
    }
  }

private:
  /// The AABB tree: intersection detection and more
  AABB_tree_ tree_;

  AABB_tree_* bounding_tree_;

  // cache queries and intersected primitive
  typedef typename boost::make_variant_over<Allowed_query_types>::type Cached_query;
  struct Query_cache
  {
    Query_cache() : has_cache(false) {}
    bool has_cache;
    Cached_query cached_query;
    AABB_primitive_id cached_primitive_id;
  };
#ifdef CGAL_LINKED_WITH_TBB
  mutable tbb::enumerable_thread_specific<Query_cache> query_cache;
#else
  mutable Query_cache query_cache;
#endif

  //random number generator for Construct_initial_points
  CGAL::Random* p_rng_;
  bool delete_rng_;

public:

  template <typename Query>
  void cache_primitive(const Query& q,
                       const AABB_primitive_id id) const
  {
#ifdef CGAL_LINKED_WITH_TBB
    Query_cache &qc = query_cache.local();
    qc.cached_query = Cached_query(q);
    qc.has_cache = true;
    qc.cached_primitive_id = id;
#else
    query_cache.cached_query = Cached_query(q);
    query_cache.has_cache = true;
    query_cache.cached_primitive_id = id;
#endif
  }

  template <typename Query>
  bool query_is_cached(const Query& q) const {
#ifdef CGAL_LINKED_WITH_TBB
    Query_cache &qc = query_cache.local();
    return qc.has_cache && (qc.cached_query == Cached_query(q));
#else
    return query_cache.has_cache 
      && (query_cache.cached_query == Cached_query(q));
#endif
  }

  AABB_primitive_id cached_primitive_id() const {
#ifdef CGAL_LINKED_WITH_TBB
    return query_cache.local().cached_primitive_id;
#else
    return query_cache.cached_primitive_id;
#endif
  }

  void set_random_generator(CGAL::Random* p_rng)
  {
    p_rng_ = p_rng;
  }

private:
  // Disabled copy constructor & assignment operator
  typedef Polyhedral_mesh_domain_3 Self;
  Polyhedral_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Polyhedral_mesh_domain_3





template<typename P_, typename IGT_, typename TA,
         typename Tag, typename E_tag_>
template<class OutputIterator>
OutputIterator
Polyhedral_mesh_domain_3<P_,IGT_,TA,Tag,E_tag_>::
Construct_initial_points::operator()(OutputIterator pts,
                                     const int n) const
{
  typename IGT::Construct_ray_3 ray = IGT().construct_ray_3_object();
  typename IGT::Construct_vector_3 vector = IGT().construct_vector_3_object();

  const Bounding_box bbox = r_domain_.tree_.bbox();
  const Point_3 center( FT( (bbox.xmin() + bbox.xmax()) / 2),
                        FT( (bbox.ymin() + bbox.ymax()) / 2),
                        FT( (bbox.zmin() + bbox.zmax()) / 2) );

  CGAL::Random& rng = *(r_domain_.p_rng_ != 0 ?
                        r_domain_.p_rng_ :
                        new Random(0));
  Random_points_on_sphere_3<Point_3> random_point(1., rng);

  int i = n;
# ifdef CGAL_MESH_3_VERBOSE
  std::cerr << "construct initial points:" << std::endl;
# endif
  // Point construction by ray shooting from the center of the enclosing bbox
  while ( i > 0 )
  {
    const Ray_3 ray_shot = ray(center, vector(CGAL::ORIGIN,*random_point));

#ifdef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
    Intersection intersection = r_domain_.construct_intersection_object()(ray_shot);
    if(CGAL::cpp0x::get<2>(intersection) != 0) {
#else
    if(r_domain_.do_intersect_surface_object()(ray_shot)) {
      Intersection intersection = r_domain_.construct_intersection_object()(ray_shot);
#endif
      *pts++ = std::make_pair(CGAL::cpp0x::get<0>(intersection),
                              CGAL::cpp0x::get<1>(intersection));

      --i;

#ifdef CGAL_MESH_3_VERBOSE
      std::cerr << boost::format("\r             \r"
                                 "%1%/%2% initial point(s) found...")
        % (n - i)
        % n;
# endif
    }
    ++random_point;
  }

#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << std::endl;
#endif
  if(r_domain_.p_rng_ == 0) delete &rng;
  return pts;
}


template<typename P_, typename IGT_, typename TA,
         typename Tag, typename E_tag_>
typename Polyhedral_mesh_domain_3<P_,IGT_,TA,Tag,E_tag_>::Subdomain
Polyhedral_mesh_domain_3<P_,IGT_,TA,Tag,E_tag_>::
Is_in_domain::operator()(const Point_3& p) const
{
  if(r_domain_.bounding_tree_ == 0) return Subdomain();

  Inside_functor inside_functor(*(r_domain_.bounding_tree_));
  Bounded_side side = inside_functor(p);

  if(side == CGAL::ON_UNBOUNDED_SIDE) { return Subdomain(); }
  else { return Subdomain(Subdomain_index(1)); } // case ON_BOUNDARY && ON_BOUNDED_SIDE
}




}  // end namespace CGAL

#include <CGAL/enable_warnings.h>
  
#endif // POLYHEDRAL_MESH_TRAITS_3_H_
