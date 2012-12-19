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

#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Mesh_3/Triangle_accessor_primitive.h>
#include <CGAL/Triangle_accessor_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/Mesh_3/Creator_weighted_point_3.h>

#include <boost/optional.hpp>
#include <boost/none.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/contains.hpp>
#include <CGAL/tuple.h>
#include <boost/format.hpp>
#include <boost/variant.hpp>

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
// Surface_patch_index_generator
// To use patch_id enclosed in AABB_primitives or not
// -----------------------------------
template < typename Subdomain_index, typename Tag >
struct Surface_patch_index_generator
{
  typedef std::pair<Subdomain_index,Subdomain_index>  Surface_patch_index;
  typedef Surface_patch_index                         type;
  
  template < typename Primitive_id >
  Surface_patch_index operator()(const Primitive_id&)
  { return Surface_patch_index(0,1); }
};
  
template < typename Subdomain_index >
struct Surface_patch_index_generator<Subdomain_index, CGAL::Tag_true>
{
  typedef Subdomain_index       Surface_patch_index;
  typedef Surface_patch_index   type;

  template < typename Primitive_id >
  Surface_patch_index operator()(const Primitive_id& primitive_id)
  { return primitive_id->patch_id(); }
};


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
  typedef CGAL::Mesh_3::Robust_intersection_traits_3<Gt> type;
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


/**
 * @class Polyhedral_mesh_domain_3
 *
 *
 */
template<class Polyhedron,
         class IGT_,
         class TriangleAccessor=Triangle_accessor_3<Polyhedron,IGT_>,
         class Use_patch_id_tag=Tag_false,
         class Use_exact_intersection_construction_tag = CGAL::Tag_true>
class Polyhedral_mesh_domain_3
{
  typedef typename Mesh_3::details::IGT_generator<
    IGT_,Use_exact_intersection_construction_tag>::type IGT;
  
public:
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
  typedef typename Mesh_3::details::Surface_patch_index_generator<
    Subdomain_index,Use_patch_id_tag>::type               Surface_patch_index;
  typedef boost::optional<Surface_patch_index>            Surface_patch;
  /// Type of indexes to characterize the lowest dimensional face of the input
  /// complex on which a vertex lie
  typedef typename Mesh_3::details::Index_generator<
    Subdomain_index, Surface_patch_index>::type           Index;

  typedef CGAL::cpp11::tuple<Point_3,Index,int> Intersection;


  typedef typename IGT::FT         FT;

  // Kernel_traits compatibility
  typedef IGT R;

private:
  typedef Mesh_3::Triangle_accessor_primitive<
    TriangleAccessor, IGT>                              AABB_primitive;
  typedef class AABB_traits<IGT,AABB_primitive>         AABB_traits;
  typedef class AABB_tree<AABB_traits>                  AABB_tree_;
private:
  typedef typename AABB_tree_::Primitive_id              AABB_primitive_id;
  typedef typename AABB_traits::Bounding_box            Bounding_box;
  
public:

  /// Default constructor
  Polyhedral_mesh_domain_3()
    : tree_()
    , bounding_tree_(&tree_) {}
  
  /**
   * @brief Constructor. Contruction from a polyhedral surface
   * @param polyhedron the polyhedron describing the polyhedral surface
   */
  Polyhedral_mesh_domain_3(const Polyhedron& p)
    : tree_(TriangleAccessor().triangles_begin(p),
            TriangleAccessor().triangles_end(p)),
      bounding_tree_(&tree_) // the bounding tree is tree_
  { 
    if(!p.is_pure_triangle()) {
      std::cerr << "Your input polyhedron must be triangulated!\n";
      CGAL_error_msg("Your input polyhedron must be triangulated!");
    }
  }

  Polyhedral_mesh_domain_3(const Polyhedron& p,
                           const Polyhedron& bounding_polyhedron)
    : tree_(TriangleAccessor().triangles_begin(p),
            TriangleAccessor().triangles_end(p))
    , bounding_tree_(new AABB_tree_(TriangleAccessor().triangles_begin(bounding_polyhedron),
                                    TriangleAccessor().triangles_end(bounding_polyhedron)))
  { 
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
                           const Polyhedron& bounding_polyhedron)
  {
    if(begin != end) { 
      for(; begin != end; ++begin) {
        tree_.insert(TriangleAccessor().triangles_begin(**begin),
                     TriangleAccessor().triangles_end(**begin));
      }
      tree_.build();
      bounding_tree_ = 
        new AABB_tree_(TriangleAccessor().triangles_begin(bounding_polyhedron),
                       TriangleAccessor().triangles_end(bounding_polyhedron));
    }
    else {
      tree_.rebuild(TriangleAccessor().triangles_begin(bounding_polyhedron),
                    TriangleAccessor().triangles_end(bounding_polyhedron));
      bounding_tree_ = &tree_;
    }
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
                           InputPolyhedraPtrIterator end)
  {
    if(begin != end) {
      for(; begin != end; ++begin) {
        tree_.insert(TriangleAccessor().triangles_begin(**begin),
                     TriangleAccessor().triangles_end(**begin));
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
      // Check first the bounding_tree
      boost::optional<AABB_primitive_id> primitive_id =
        r_domain_.bounding_tree_ == 0 ? 
        boost::none :
        r_domain_.bounding_tree_->any_intersected_primitive(q);
      
      if ( primitive_id )
      { 
        return Surface_patch(r_domain_.make_surface_index(*primitive_id));
      }
      // then the other trees
      else if ( r_domain_.bounding_tree_ != &r_domain_.tree_ )
      {
        primitive_id = r_domain_.tree_.any_intersected_primitive(q);
        if ( primitive_id )
        { 
          return Surface_patch(r_domain_.make_surface_index(*primitive_id));
        }
      }
      
      return Surface_patch();
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
      typedef boost::optional<
        typename AABB_tree_::Object_and_primitive_id> AABB_intersection;
      typedef Point_3 Bare_point;

      CGAL_precondition(r_domain_.do_intersect_surface_object()(q));

      AABB_intersection intersection = 
        r_domain_.bounding_tree_ == 0 ?
        AABB_intersection() :
        r_domain_.bounding_tree_->any_intersection(q);

      if(! intersection && 
         r_domain_.bounding_tree_ != &r_domain_.tree_) {
        intersection = r_domain_.tree_.any_intersection(q);
      }
      if ( intersection )
      {
        // Get primitive
        AABB_primitive_id primitive_id = (*intersection).second;
        
        // intersection may be either a point or a segment
        if ( const Bare_point* p_intersect_pt =
                              object_cast<Bare_point>(&(*intersection).first) )
        {
          return Intersection(*p_intersect_pt,
                              r_domain_.index_from_surface_patch_index(
                                r_domain_.make_surface_index(primitive_id)),
                              2);
        }
        else if ( const Segment_3* p_intersect_seg =
                              object_cast<Segment_3>(&(*intersection).first) )
        {
          return Intersection(p_intersect_seg->source(),
                              r_domain_.index_from_surface_patch_index(
                                r_domain_.make_surface_index(primitive_id)),
                              2);
        }
        else
          CGAL_error_msg("Mesh_3 error : AABB_tree any_intersection result is "
                         "not a point nor a segment");
      }

      // Should not happen
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
    Mesh_3::details::Surface_patch_index_generator<Subdomain_index,Use_patch_id_tag>
      generator;
    
    return generator(primitive_id);
  }

  // Undocumented function, used to implement a sizing field that
  // computes lfs using this AABB tree. That avoids to rebuild the same
  // tree.
  typedef AABB_tree_ AABB_tree;
  const AABB_tree& aabb_tree() const {
    return tree_;
  }

protected:
  void add_primitives(const Polyhedron& p)
  {
    tree_.insert(TriangleAccessor().triangles_begin(p),
                 TriangleAccessor().triangles_end(p));
    
    tree_.build();
  }

private:
  /// The AABB tree: intersection detection and more
  AABB_tree_ tree_;

  AABB_tree_* bounding_tree_;

private:
  // Disabled copy constructor & assignment operator
  typedef Polyhedral_mesh_domain_3 Self;
  Polyhedral_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Polyhedral_mesh_domain_3





template<typename P_, typename IGT_, typename TA, typename Tag, typename E_tag_>
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
  
  Random_points_on_sphere_3<Point_3> random_point(1.);

  int i = n;
#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << "construct initial points:" << std::endl;
#endif
  // Point construction by ray shooting from the center of the enclosing bbox
  while ( i > 0 )
  {
    const Ray_3 ray_shot = ray(center, vector(CGAL::ORIGIN,*random_point));
    if(r_domain_.do_intersect_surface_object()(ray_shot)) {
      Intersection intersection = r_domain_.construct_intersection_object()(ray_shot);
      *pts++ = std::make_pair(CGAL::cpp11::get<0>(intersection),
                              CGAL::cpp11::get<1>(intersection));
        
      --i;
        
#ifdef CGAL_MESH_3_VERBOSE
      std::cerr << boost::format("\r             \r"
                                 "%1%/%2% initial point(s) found...")
        % (n - i)
        % n;
#endif
    }

    ++random_point;
  }
  
#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << std::endl;
#endif
  return pts;
}


template<typename P_, typename IGT_, typename TA, typename Tag, typename E_tag_>
typename Polyhedral_mesh_domain_3<P_,IGT_,TA,Tag,E_tag_>::Subdomain
Polyhedral_mesh_domain_3<P_,IGT_,TA,Tag,E_tag_>::
Is_in_domain::operator()(const Point_3& p) const
{
  if(r_domain_.bounding_tree_ == 0) return Subdomain();
  const Bounding_box bbox = r_domain_.bounding_tree_->bbox();

  if(   p.x() < bbox.xmin() || p.x() > bbox.xmax()
     || p.y() < bbox.ymin() || p.y() > bbox.ymax()
     || p.z() < bbox.zmin() || p.z() > bbox.zmax() )
  {
    return Subdomain();
  }
  
  // Shoot ray
  typename IGT::Construct_ray_3 ray = IGT().construct_ray_3_object();
  typename IGT::Construct_vector_3 vector = IGT().construct_vector_3_object();
  
  Random_points_on_sphere_3<Point_3> random_point(1.);

  const Ray_3 ray_shot = ray(p, vector(CGAL::ORIGIN,*random_point));

  if ( (r_domain_.bounding_tree_->number_of_intersected_primitives(ray_shot)&1) == 1 )
    return Subdomain(Subdomain_index(1));
  else
    return Subdomain();
}




}  // end namespace CGAL


#endif // POLYHEDRAL_MESH_TRAITS_3_H_
