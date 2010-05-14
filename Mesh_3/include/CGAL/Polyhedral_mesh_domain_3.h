// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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

#include <CGAL/Mesh_3/Triangle_accessor_primitive.h>
#include <CGAL/Triangle_accessor_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/Mesh_3/Creator_weighted_point_3.h>

#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/contains.hpp>
#include <CGAL/tuple.h>
#include <boost/format.hpp>

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

}  // end namespace details

}  // end namespace Mesh_3


/**
 * @class Polyhedral_mesh_domain_3
 *
 *
 */
template<class Polyhedron,
         class IGT,
         class TriangleAccessor=Triangle_accessor_3<Polyhedron, IGT> >
class Polyhedral_mesh_domain_3
{
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
  typedef std::pair<Subdomain_index, Subdomain_index> Surface_index;
  typedef boost::optional<Surface_index> Surface_patch;
  /// Type of indexes to characterize the lowest dimensional face of the input
  /// complex on which a vertex lie
  typedef boost::variant<Subdomain_index, Surface_index> Index;

  typedef CGAL::cpp0x::tuple<Point_3,Index,int> Intersection;


  typedef typename IGT::FT         FT;

  // Kernel_traits compatibility
  typedef IGT R;


  /**
   * @brief Constructor. Contruction from a polyhedral surface
   * @param polyhedron the polyhedron describing the polyhedral surface
   */
  Polyhedral_mesh_domain_3(const Polyhedron& p)
    : tree_(TriangleAccessor().triangles_begin(p),
            TriangleAccessor().triangles_end(p))
  {}

  /// Destructor
  ~Polyhedral_mesh_domain_3() {}


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
   * if \ccc{true} is returned and to the default \ccc{Surface_index}
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
      if ( r_domain_.tree_.do_intersect(q) )
        return Surface_patch(r_domain_.make_surface_index());
      else
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
        typename AABB_tree::Object_and_primitive_id> AABB_intersection;
      typedef Point_3 Bare_point;

      CGAL_precondition(r_domain_.do_intersect_surface_object()(q));

      AABB_intersection intersection = r_domain_.tree_.any_intersection(q);
      if ( intersection )
      {
        // intersection may be either a point or a segment
        if ( const Bare_point* p_intersect_pt =
                              object_cast<Bare_point>(&(*intersection).first) )
        {
          return Intersection(*p_intersect_pt,
                              r_domain_.index_from_surface_index(
                                                r_domain_.make_surface_index()),
                              2);
        }
        else if ( const Segment_3* p_intersect_seg =
                              object_cast<Segment_3>(&(*intersection).first) )
        {
          return Intersection(p_intersect_seg->source(),
                              r_domain_.index_from_surface_index(
                                                r_domain_.make_surface_index()),
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
  Index index_from_surface_index(const Surface_index& index) const
  { return Index(index); }

  /**
   * Returns the index to be stored in a vertex lying in the subdomain
   * identified by \c index.
   */
  Index index_from_subdomain_index(const Subdomain_index& index) const
  { return Index(index); }

  /**
   * Returns the \c Surface_index of the surface patch
   * where lies a vertex with dimension 2 and index \c index.
   */
  Surface_index surface_index(const Index& index) const
  { return boost::get<Surface_index>(index); }

  /**
   * Returns the index of the subdomain containing a vertex
   *  with dimension 3 and index \c index.
   */
  Subdomain_index subdomain_index(const Index& index) const
  { return boost::get<Subdomain_index>(index); }

public:
  Surface_index make_surface_index() const
  {
    return Surface_index(0,1);
  }
  
private:
  typedef Mesh_3::Triangle_accessor_primitive<TriangleAccessor, IGT>
                                                                AABB_primitive;
  typedef class AABB_traits<IGT,AABB_primitive> AABB_traits;
  typedef class AABB_tree<AABB_traits> AABB_tree;
  typedef typename AABB_traits::Bounding_box Bounding_box;



private:
  /// The AABB tree: intersection detection and more
  AABB_tree tree_;

private:
  // Disabled copy constructor & assignment operator
  typedef Polyhedral_mesh_domain_3<Polyhedron, IGT, TriangleAccessor> Self;
  Polyhedral_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Polyhedral_mesh_domain_3





template<typename P_, typename IGT, typename TA>
template<class OutputIterator>
OutputIterator
Polyhedral_mesh_domain_3<P_,IGT,TA>::Construct_initial_points::operator()(
                                                    OutputIterator pts,
                                                    const int n) const
{
  typedef boost::optional<typename AABB_tree::Object_and_primitive_id>
                                                            AABB_intersection;
  
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
    
    AABB_intersection intersection = r_domain_.tree_.any_intersection(ray_shot);
    if ( intersection )
    {
      const Point_3* p_pt = object_cast<Point_3>(&(*intersection).first);
      if ( NULL != p_pt )
      {
        *pts++ = std::make_pair(*p_pt,
          r_domain_.index_from_surface_index(r_domain_.make_surface_index()));
        
        --i;
        
#ifdef CGAL_MESH_3_VERBOSE
        std::cerr << boost::format("\r             \r"
                                   "%1%/%2% initial point(s) found...")
        % (n - i)
        % n;
#endif
      }
    }

    ++random_point;
  }
  
#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << std::endl;
#endif
  return pts;
}


template<typename P_, typename IGT, typename TA>
typename Polyhedral_mesh_domain_3<P_,IGT,TA>::Subdomain
Polyhedral_mesh_domain_3<P_,IGT,TA>::Is_in_domain::operator()(
                                                    const Point_3& p) const
{
  const Bounding_box bbox = r_domain_.tree_.bbox();

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

  if ( (r_domain_.tree_.number_of_intersected_primitives(ray_shot)&1) == 1 )
    return Subdomain(Subdomain_index(1));
  else
    return Subdomain();
}




}  // end namespace CGAL


#endif // POLYHEDRAL_MESH_TRAITS_3_H_
