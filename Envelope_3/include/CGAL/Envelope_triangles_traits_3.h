// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $Source: /CVSROOT/CGAL/Packages/Envelope_3/include/CGAL/Envelope_triangles_traits_3.h,v $
// $Revision$ $Date$
// $Name:  $
//
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>

/*! \file CGAL/Envelope_triangles_traits_3.h
 * \brief Model for CGAL's EnvelopeTraits_3 concept.
 * \endlink
 */

#ifndef ENVELOPE_TRIANGLES_TRAITS_3_H
#define ENVELOPE_TRIANGLES_TRAITS_3_H

#include <CGAL/Handle_for.h>
#include <CGAL/Object.h>
#include <CGAL/enum.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Timer.h>
#include <CGAL/Envelope_base.h>
#include <CGAL/Arr_segment_traits_2.h>

#include <vector>

//#define CGAL_DEBUG_ENVELOPE_3_TRIANGLES_TRAITS

//#define PRINT_TRAINGLES_INTERSECTION_STATS

// this cache doesn't help much
#define CGAL_ENV_TRIANGLES_TRAITS_CACHE_POINT_ON

CGAL_BEGIN_NAMESPACE

template <class Kernel_> class Envelope_triangle_3;
 
template <class Kernel_>
class Envelope_triangles_traits_3 : public Arr_segment_traits_2<Kernel_>
{
public:
  typedef Arr_segment_traits_2<Kernel_>             Traits_2;
  typedef typename Traits_2::Point_2                Point_2;
  typedef typename Traits_2::Curve_2                Curve_2;
  typedef typename Traits_2::X_monotone_curve_2     X_monotone_curve_2;

  typedef Kernel_                                   Kernel;
  typedef Envelope_triangles_traits_3<Kernel>       Self;

  typedef typename Kernel::Point_3                  Point_3;

  /*!
   * \class Representation of a 3d triangle with cached data.
   */
  class _Triangle_cached_3 //: public Ref_counted
  {
  public:

    typedef typename Kernel::Plane_3               Plane_3;
    typedef typename Kernel::Triangle_3            Triangle_3;
    typedef typename Kernel::Point_3               Point_3;

  protected:

    Plane_3 pl;          // The plane that supports the triangle.
    Point_3 vertices[3]; // The vertices of the triangle.
    bool    is_vert;     // Is this a vertical triangle.
    bool    is_degen;    // Is the triangle degenerate (a single point/segment).

  public:

    /*!
     * Default constructor.
     */
    _Triangle_cached_3() :
      is_vert(false),
      is_degen(true)
    {}

    /*!
     * Constructor from a triangle.
     * \param tri The triangle.
     */
    _Triangle_cached_3(const Triangle_3 & tri)
    {
      Kernel   kernel;

      typename Kernel::Construct_vertex_3
        construct_vertex = kernel.construct_vertex_3_object();

      vertices[0] = construct_vertex(tri, 0);
      vertices[1] = construct_vertex(tri, 1);
      vertices[2] = construct_vertex(tri, 2);

      is_degen = kernel.is_degenerate_3_object()(tri);

      if (! is_degen)
      {
        pl = kernel.construct_plane_3_object()(vertices[0],
                                               vertices[1], vertices[2]);
        Self self;
        is_vert = kernel.collinear_2_object()(self.project(vertices[0]),
                                              self.project(vertices[1]),
                                              self.project(vertices[2]));
      }
    }

    /*!
     * Construct a triangle from three end-points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \param p3 The third point.
     */
    _Triangle_cached_3(const Point_3 &p1, const Point_3 &p2,
                       const Point_3 &p3)
    {
      vertices[0] = p1;
      vertices[1] = p2;
      vertices[2] = p3;

      Kernel   kernel;

      is_degen = kernel.collinear_3_object()(vertices[0], vertices[1], vertices[2]);

      if (! is_degen)
      {
        pl = kernel.construct_plane_3_object()(vertices[0], vertices[1], vertices[2]);
        Self self;
        is_vert = kernel.collinear_2_object()(self.project(vertices[0]),
                                              self.project(vertices[1]),
                                              self.project(vertices[2]));
      }
    }

    /*!
     * Construct a triangle from 3 end-points on a supporting plane.
     * \param supp_plane The supporting plane.
     * \param p1 The first point.
     * \param p2 The second point.
     * \param p3 The third point.
     * \pre The 3 endpoints are not the collinear and all lie on the given plane.

     */
    _Triangle_cached_3(const Plane_3& supp_plane,
                       const Point_3 &p1,
                       const Point_3 &p2,
                       const Point_3 &p3) :
      pl(supp_plane)
    {
      Kernel   kernel;

      CGAL_precondition(kernel.has_on_3_object() (pl, p1) &&
                        kernel.has_on_3_object() (pl, p2) &&
                        kernel.has_on_3_object() (pl, p3));      
      CGAL_precondition(!kernel.collinear_3_object()(p1, p2, p3));

      vertices[0] = p1;
      vertices[1] = p2;
      vertices[2] = p3;

      is_degen = false;

      Self self;
      is_vert = kernel.collinear_2_object()(self.project(vertices[0]),
                                            self.project(vertices[1]),
                                            self.project(vertices[2]));
    }

    /*!
     * Assignment operator.
     * \param seg the source segment to copy from
     */
    const _Triangle_cached_3& operator=(const Triangle_3 &tri)
    {

      Kernel   kernel;

      typename Kernel_::Construct_vertex_3
        construct_vertex = kernel.construct_vertex_3_object();

      vertices[0] = construct_vertex(tri, 0);
      vertices[1] = construct_vertex(tri, 1);
      vertices[2] = construct_vertex(tri, 2);

      is_degen = kernel.is_degenerate_3_object()(tri);

      if (! is_degen)
      {
        pl = kernel.construct_plane_3_object()(vertices[0],
                                               vertices[1], vertices[2]);
        Self self;
        is_vert = kernel.collinear_2_object()(self.project(vertices[0]),
                                              self.project(vertices[1]),
                                              self.project(vertices[2]));
      }

      return (*this);
    }

    /*!
     * Get the ith endpoint.
     */
    const Point_3& vertex(unsigned int i) const
    {
      if (i > 2)
        i = i % 3;
      CGAL_assertion(i>=0 && i<=2);
      return vertices[i];
    }

    /*!
     * Get the supporting plane.
     */
    const Plane_3& plane() const
    {
      CGAL_precondition (!is_degen);
      return (pl);
    }

    /*!
     * Check if the triangel is vertical.
     */
    bool is_vertical() const
    {
      CGAL_precondition (!is_degen);
      return (is_vert);
    }

    /*!
     * Check if the triangle is degenerate.
     */
    bool is_degenerate() const
    {
      return (is_degen);
    }

  };

public:
  // types for EnvelopeTraits_3 concept
  //! type of xy-monotone surfaces
  typedef Envelope_triangle_3<Kernel>               Xy_monotone_surface_3;
  //! type of surfaces
  typedef Xy_monotone_surface_3                     Surface_3;

  // we have a collision between the Kernel's Intersect_2 and the one
  // from the segment traits
  typedef typename Traits_2::Intersect_2            Intersect_2;

protected:
  typedef typename Kernel::FT                       FT;
  typedef typename Kernel::Triangle_2               Triangle_2;
  typedef typename Kernel::Segment_2                Segment_2;

  typedef typename Kernel::Point_3                  Point_3;
  typedef typename Kernel::Segment_3                Segment_3;
  typedef typename Kernel::Triangle_3               Triangle_3;
  typedef typename Kernel::Plane_3                  Plane_3;

  typedef typename Kernel::Assign_2                 Assign_2;
  typedef typename Kernel::Construct_vertex_2       Construct_vertex_2;

  typedef typename Kernel::Assign_3                 Assign_3;
  typedef typename Kernel::Intersect_3              Intersect_3;
  typedef typename Kernel::Construct_vertex_3       Construct_vertex_3;


  typedef typename Kernel::Line_2                   Line_2;
  typedef typename Kernel::Direction_2              Direction_2;

  typedef typename Kernel::Line_3                   Line_3;
  typedef typename Kernel::Direction_3              Direction_3;

  #ifdef CGAL_ENV_TRIANGLES_TRAITS_CACHE_POINT_ON
    // caching the computation of a surface 3d point from a 2d point
    typedef std::pair<Xy_monotone_surface_3, Point_2> Surface_point_pair;
    struct Less_surface_point_pair
    {
      bool operator() (const Surface_point_pair& sp1,
                       const Surface_point_pair& sp2) const
      {
        // Compare the pairs of IDs lexicographically.
        return (sp1.first < sp2.first ||
                (sp1.first == sp2.first && Kernel().less_xy_2_object()(sp1.second,sp2.second)));
      }
    };
    typedef std::map<Surface_point_pair, Point_3,
                     Less_surface_point_pair>          Surface_point_cache;
  #endif

  typedef std::pair<Curve_2, Intersection_type>        Intersection_curve;
public:

  /***************************************************************************/
  // EnvelopeTraits_3 functors
  /***************************************************************************/

  /*!\brief
   * Subdivide the given surface into envelope relevant xy-monotone 
   * parts, and insert them into the output iterator.
   * 
   * The iterator value-type is Xy_monotone_surface_3
   */
  class Construct_envelope_xy_monotone_parts_3
  {
  protected:
    const Self *parent;
  public:

    Construct_envelope_xy_monotone_parts_3(const Self* p)
      : parent(p)
    {}
    // create xy-monotone surfaces from a general surface
    // return a past-the-end iterator
    template <class OutputIterator>
    OutputIterator operator()(const Surface_3& surf, OutputIterator o) const
    {
      parent->total_timer.start();
      // a triangle is already xy-monotone
      *o++ = surf;
      parent->total_timer.stop();
      return o;
    }
  };

  /*! Get a Construct_envelope_xy_monotone_parts_3 functor object. */
  Construct_envelope_xy_monotone_parts_3
  construct_envelope_xy_monotone_parts_3_object() const
  {
    return Construct_envelope_xy_monotone_parts_3(this);
  }

  /*!\brief
   * Insert all 2D curves, which form the boundary of the vertical
   * projection of the surface onto the xy-plane, into the output iterator.
   * The iterator value-type is Curve_2.
   */
  class Construct_projected_boundary_curves_2
  {
  protected:
    const Self *parent;
  public:

    Construct_projected_boundary_curves_2(const Self* p)
      : parent(p)
    {}

    // insert into the OutputIterator all the (2d) curves of the boundary of
    // the vertical projection of the surface on the xy-plane
    // the OutputIterator value type is Curve_2
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& surf,
                              OutputIterator o) const
    {
      parent->total_timer.start();
      parent->pboundary_timer.start();
      
      const Point_3 &a1 = surf.vertex(0),
                     a2 = surf.vertex(1),
                     a3 = surf.vertex(2);
      Point_2 b1 = parent->project(a1),
              b2 = parent->project(a2),
              b3 = parent->project(a3);

      if (!surf.is_vertical())
      {
        *o++ = Curve_2(b1, b2);
        *o++ = Curve_2(b2, b3);
        *o++ = Curve_2(b3, b1);
      }
      else
      {
        // only 2 curves in the output
        Kernel k;
        if (k.collinear_are_ordered_along_line_2_object()(b1, b2, b3))
        {
          *o++ = Curve_2(b1, b2);
          *o++ = Curve_2(b2, b3);
        }
        else if (k.collinear_are_ordered_along_line_2_object()(b1, b3, b2))
        {
          *o++ = Curve_2(b1, b3);
          *o++ = Curve_2(b3, b2);
        }
        else
        {
          *o++ = Curve_2(b2, b1);
          *o++ = Curve_2(b1, b3);
        }
      }
      parent->total_timer.stop();
      parent->pboundary_timer.stop();
      return o;
    }  
  };  
  




  /*! Get a Construct_projected_boundary_curves_2 functor object. */
  Construct_projected_boundary_curves_2
  construct_projected_boundary_curves_2_object() const
  {
    return Construct_projected_boundary_curves_2(this);
  }

  /*!\brief
   * Insert all the 2D projections (onto the xy-plane) of the 
   * intersection objects between s1 and s2 into the output iterator.
   *
   * The iterator value-type is Object. An Object may be:
   * 1. A pair<Curve_2,Intersection_type>, where the intersection 
   * type is an enumeration that can take the values
   * {Transversal, Tangency, Unknown}.
   * 2. A Point_2 instance (in degenerate cases).
   */
  class Construct_projected_intersections_2
  {
  protected:
    const Self *parent;
  public:

    Construct_projected_intersections_2(const Self* p)
      : parent(p)
    {}
    
    // insert into OutputIterator all the (2d) projections on the xy plane of
    // the intersection objects between the 2 surfaces
    // the data type of OutputIterator is Object
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s1,
                              const Xy_monotone_surface_3& s2,
                              OutputIterator o) const
    {
      parent->total_timer.start();
      parent->intersection_timer.start();
      Kernel k;
      if (!k.do_intersect_3_object()(s1, s2))
      {
        parent->total_timer.stop();
        parent->intersection_timer.stop();
        return o;
      }
        
      Object inter_obj = parent->intersection(s1,s2);
      if (inter_obj.is_empty())
      {
        parent->total_timer.stop();
        parent->intersection_timer.stop();
        return o;
      }

      Point_3 point;
      Segment_3 curve;
      if (k.assign_3_object()(point, inter_obj))
        *o++ = make_object(parent->project(point));
      else if (k.assign_3_object()(curve, inter_obj))
      {
        Curve_2 projected_cv = parent->project(curve);
        if (projected_cv.is_degenerate())
          *o++ = make_object(projected_cv.left());
        else
        {
          Intersection_curve inter_cv(projected_cv, TRANSVERSAL);
          *o++ = make_object(inter_cv);
        }
      }
      else
      {
        // if we get here the triangles may be overlapping, and the
        // intersection is a polygon
        // it is important to the current algorithm only if the surfaces
        // are vertical
        if (s1.is_vertical())
        {
          CGAL_assertion(s2.is_vertical());
          std::vector<Point_3> poly;
          CGAL_assertion(k.assign_3_object()(poly, inter_obj));
          k.assign_3_object()(poly, inter_obj);
          unsigned int n = poly.size();
          // project the points
          std::vector<Point_2> proj_poly(n);
          for(unsigned int i=0; i<n; ++i)
            proj_poly[i] = parent->project(poly[i]);
          // return the projected intersections
          for(unsigned int i=0; i<n; ++i)
          {
            // deal with segment between points i and i+1 
	    // (mod proj_poly.size())
            if (proj_poly[i] == proj_poly[(i+1) % n])
              *o++ = make_object(proj_poly[i]);
            else
            {
              Curve_2 cv(proj_poly[i], proj_poly[(i+1) % n]);
              Intersection_curve inter_cv(cv, UNKNOWN);
              *o++ = make_object(inter_cv);
            }
          }
        }
//        // it is irrelevant to the current d&c algorithm (redundant data)
//        // but we return it anyway, so we can only invoke compare_left over
//        // intersection curves
      }
      
      parent->total_timer.stop();
      parent->intersection_timer.stop();
      return o;
    }  
  };  

  /*! Get a Construct_projected_intersections_2 functor object. */
  Construct_projected_intersections_2
  construct_projected_intersections_2_object() const
  {
    return Construct_projected_intersections_2(this);
  }



  /*!\brief
   * Check if the surface s1 is closer/equally distanced/farther 
   * from the envelope with respect to s2 at the xy-coordinates of p/c.
   */
  class Compare_distance_to_envelope_3
  {
  protected:
    const Self *parent;

  public:

    Compare_distance_to_envelope_3(const Self* p)
      : parent(p)
    {}

    // check which of the surfaces is closer to the envelope at the xy 
    // coordinates of point
    // (i.e. lower if computing the lower envelope, or upper if computing 
    // the upper envelope)
    // precondition: the surfaces are defined in point
    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& surf1,
                                 const Xy_monotone_surface_3& surf2) const
    {
      bool use_timer = true;
      if (parent->total_timer.is_running())
        use_timer = false;

      if (use_timer)
      {
        parent->total_timer.start();
        parent->compare_timer.start();

      }
      
      #ifdef CGAL_DEBUG_ENVELOPE_3_TRIANGLES_TRAITS
        std::cout << "traits: point= " << p << " surf1= " << surf1 
		  << " surf2= " << surf2 << std::endl;
      #endif

      // we compute the points on the planes, and then compare their z 
      // coordinates
      const Plane_3& plane1 = surf1.plane();
      const Plane_3& plane2 = surf2.plane();

      // if the 2 triangles have the same supporting plane, and they are not 
      // vertical, then they have the same z coordinate over this point
      if ((plane1 == plane2 || plane1 == plane2.opposite()) &&
          !surf1.is_vertical())
      {
        if (use_timer)
        {
          parent->total_timer.stop();
          parent->compare_timer.stop();
        }
        return EQUAL;
      }

      // these should contain the points on the surfaces that we need to 
      // compare for the envelope
      Point_3   ip1, ip2;
      bool ip1_found = false, ip2_found = false;

      #ifdef CGAL_ENV_TRIANGLES_TRAITS_CACHE_POINT_ON
        // first try the cache:

        typename Surface_point_cache::iterator  cache_iter;
//        Surface_point_pair  r1(const_cast<Xy_monotone_surface_3>(surf1), p);
//        Surface_point_pair  pair2(const_cast<Xy_monotone_surface_3>(surf2), p);
        Surface_point_pair                      spair1(surf1, p);
        Surface_point_pair                      spair2(surf2, p);

        cache_iter = parent->point_on_cache.find(spair1);
        if (cache_iter != (parent->point_on_cache).end())
        {
          ip1 = (*cache_iter).second;
          ip1_found = true;
          (parent->cache_hits)++;
        }

        cache_iter = parent->point_on_cache.find(spair2);
        if (cache_iter != (parent->point_on_cache).end())
        {
          ip2 = (*cache_iter).second;
          ip2_found = true;
          (parent->cache_hits)++;
        }
      #endif
      
      Kernel k;
      if (!ip1_found || !ip2_found)
      {
        parent->compute_z_at_xy_timer.start();
        // should calculate at least one point

        // Compute the intersetion between the vertical line and the given 
	// surfaces
        if (!ip1_found)
        {
          ip1 = parent->envelope_point_of_surface(p, surf1);
          #ifdef CGAL_ENV_TRIANGLES_TRAITS_CACHE_POINT_ON
            // update the cache
            (parent->point_on_cache)[spair1] = ip1;
          #endif
        }

        if (!ip2_found)
        {   
          ip2 = parent->envelope_point_of_surface(p, surf2);
          #ifdef CGAL_ENV_TRIANGLES_TRAITS_CACHE_POINT_ON
            // update the cache
            (parent->point_on_cache)[spair2] = ip2;
          #endif
        }
        parent->compute_z_at_xy_timer.stop();
      }

      
      if (use_timer)
      {
        parent->total_timer.stop();
        parent->compare_timer.stop();
      }
      #ifdef CGAL_DEBUG_ENVELOPE_3_TRIANGLES_TRAITS
        std::cout << "in compare on point, compare the points " << ip1 << " and "
                  << ip2 << std::endl;
      #endif
      // the answer changes when we compute lower/upper envelope
      if (parent->get_envelope_type() == LOWER)
        return k.compare_z_3_object()(ip1, ip2);
      else
        return k.compare_z_3_object()(ip2, ip1);
    }

    // check which of the surfaces is closer to the envelope at the xy 
    // coordinates of cv
    // (i.e. lower if computing the lower envelope, or upper if computing the
    // upper envelope)
    // precondition: the surfaces are defined in all points of cv, 
    //               and the answer is the same for each of these points
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& surf1,
                                 const Xy_monotone_surface_3& surf2) const
    {
      bool use_timer = true;
      if (parent->total_timer.is_running())
        use_timer = false;

      if (use_timer)
      {
        parent->total_timer.start();

        parent->compare_on_cv_timer.start();
      }

      // first try the endpoints, if cannot be sure, use the mid point
      Comparison_result res;
      res = parent->compare_distance_to_envelope_3_object()(cv.left(), 
                                                            surf1, surf2);

      if (res == EQUAL)
      {
        res = parent->compare_distance_to_envelope_3_object()(cv.right(), 
                                                              surf1, surf2);
        if (res == EQUAL)
        {
          Point_2 mid = parent->construct_middle_point(cv);
          res = parent->compare_distance_to_envelope_3_object()(mid, 
                                                                surf1, surf2);
        }
      }
      if (use_timer)
      {
        parent->compare_on_cv_timer.stop();
        parent->total_timer.stop();
      }
      return res;
    }
  
  };
   
  /*! Get a Compare_distance_to_envelope_3 functor object. */
  Compare_distance_to_envelope_3 
  compare_distance_to_envelope_3_object() const
  {
    return Compare_distance_to_envelope_3(this);
  }

  /*!\brief 
   * Check if the surface s1 is closer/equally distanced/farther 
   * from the envelope with
   * respect to s2 immediately above the curve c. 
   */
  class Compare_distance_to_envelope_above_3
  {
  protected:
    const Self *parent;

  public:

    Compare_distance_to_envelope_above_3(const Self* p)
      : parent(p)
    {}
    
    // check which of the surfaces is closer to the envelope on the points 
    // above the curve cv
    // (i.e. lower if computing the lower envelope, or upper if computing the
    // upper envelope)
    // precondition: the surfaces are defined above cv (to the left of cv, 
    //               if cv is directed from min point to max point)
    //               the choise between surf1 and surf2 for the envelope is 
    //               the same for every point in the infinitesimal region 
    //               above cv 
    //               the surfaces are EQUAL over the curve cv
    Comparison_result
    operator()(const X_monotone_curve_2& cv,
               const Xy_monotone_surface_3& surf1,
               const Xy_monotone_surface_3& surf2) const
    {
      parent->total_timer.start();
      parent->compare_side_timer.start();

      #ifdef CGAL_DEBUG_ENVELOPE_3_TRIANGLES_TRAITS
        std::cout << "compare left from : " << cv << std::endl;
      #endif

      // a vertical surface cannot be defined in the infinitesimal region above
      // a curve
      CGAL_precondition(!surf1.is_vertical());
      CGAL_precondition(!surf2.is_vertical());

      CGAL_precondition(parent->compare_distance_to_envelope_3_object()
                              (cv, surf1, surf2) == EQUAL);
      CGAL_precondition(parent->compare_distance_to_envelope_3_object()
                              (cv.source(), surf1, surf2) == EQUAL);
      CGAL_precondition(parent->compare_distance_to_envelope_3_object()
                              (cv.target(), surf1, surf2) == EQUAL);

      
      if (parent->do_overlap(surf1, surf2))
      {
        parent->total_timer.stop();
        parent->compare_side_timer.stop();
        return EQUAL;
      }

      #ifdef CGAL_DEBUG_ENVELOPE_3_TRIANGLES_TRAITS
        std::cout << "have the case where intersect over the curve" 
		  << std::endl;
      #endif

      // now we must have 2 different non-vertical planes:
 	    // plane1: a1*x + b1*y + c1*z + d1 = 0  , c1 != 0
 	    // plane2: a2*x + b2*y + c2*z + d2 = 0  , c2 != 0

      const Plane_3& plane1 = surf1.plane();
      const Plane_3& plane2 = surf2.plane();

      FT a1 = plane1.a(), b1 = plane1.b(), c1 = plane1.c();
      FT a2 = plane2.a(), b2 = plane2.b(), c2 = plane2.c();

 	    // our line is a3*x + b3*y + c3 = 0
 	    // it is assumed that the planes intersect over this line
      const Line_2& line = cv.line(); 
      FT a3 = line.a(), b3 = line.b(), c3 = line.c();

      // if the line was parallel to the y-axis (i.e x = const),
      // then it was enough to compare dz/dx of both planes
 	    // for general line, we change coordinates to (v, w), preserving
 	    // orientation, so the line is the w-axis in the new coordinates
 	    // (i.e v = const).
 	    //
 	    // ( v )  =  A ( x )    where A = (  a3  b3 )
 	    //   w           y                  -b3  a3
 	    //
 	    // so v =  a3*x + b3*y
 	    //    w = -b3*x + a3*y
 	    // preserving orientation since detA = a3^2 +b3^2 > 0
 	    //
 	    // We compute the planes equations in the new coordinates
 	    // and compare dz/dv
 	    //
 	    // ( x )  =  A^(-1) ( v )    where A^(-1) = ( a3  -b3 ) * detA^(-1)
 	    //   y                w                       b3   a3
 	    // so x = (a3*v - b3*w)*(1/detA)
 	    //    y = (b3*v + a3*w)*(1/detA)
 	    // plane1 ==> (a1a3 + b1b3)v + (b1a3 - a1b3)w + (c1z + d1)*detA = 0
 	    // plane2 ==> (a2a3 + b2b3)v + (b2a3 - a2b3)w + (c2z + d2)*detA = 0
 	    //
 	    // dz/dv(1) = (-a1a3 - b1b3) / c1*detA
 	    // dz/dv(2) = (-a2a3 - b2b3) / c2*detA
 	    // since detA>0 we can omit it.
 	    //
      Sign s1 = CGAL_NTS sign((a2*a3+b2*b3)/c2-(a1*a3+b1*b3)/c1);

 	    // We only need to make sure that w is in the correct direction
 	    // (going from down to up)
 	    // the original segment endpoints p1=(x1,y1) and p2=(x2,y2)
 	    // are transformed to (v1,w1) and (v2,w2), so we need that w2 > w1
 	    // (otherwise the result should be multiplied by -1)

      const Point_2& p1 = cv.left();
      const Point_2& p2 = cv.right();
      FT x1 = p1.x(), y1 = p1.y(), x2 = p2.x(), y2 = p2.y();

      Sign s2 = CGAL_NTS sign(-b3*x1+a3*y1-(-b3*x2+a3*y2));
      parent->total_timer.stop();
      parent->compare_side_timer.stop();
      // the answer is reversed when computing upper envelope vs. lower envelope
      if (parent->get_envelope_type() == LOWER)
        return Comparison_result(s1 * s2);
      else
        return Comparison_result(s1 * s2 * -1);      
    }  
  };


  /*! Get a Compare_distance_to_envelope_above_3 functor object. */
  Compare_distance_to_envelope_above_3
  compare_distance_to_envelope_above_3_object() const
  {
    return Compare_distance_to_envelope_above_3(this);
  }

  /*!\brief 
   * Check if the surface s1 is closer/equally distanced/farther 
   * from the envelope with
   * respect to s2 immediately below the curve c. 
   */
  class Compare_distance_to_envelope_below_3
  {
  protected:
    const Self *parent;

  public:

    Compare_distance_to_envelope_below_3(const Self* p)
      : parent(p)
    {}
    
    Comparison_result
    operator()(const X_monotone_curve_2& cv,
               const Xy_monotone_surface_3& surf1,
               const Xy_monotone_surface_3& surf2) const
    {
      Comparison_result left_res = parent->compare_distance_to_envelope_above_3_object()(cv, surf1, surf2);
      if (left_res == LARGER)
        return SMALLER;
      else if (left_res == SMALLER)
        return LARGER;
      else
        return EQUAL;
    }  
  };

  /*! Get a Compare_distance_to_envelope_below_3 functor object. */
  Compare_distance_to_envelope_below_3
  compare_distance_to_envelope_below_3_object() const
  {
    return Compare_distance_to_envelope_below_3(this);
  }

  /***************************************************************************/

  // checks if xy-monotone surface is vertical
  class Is_vertical_3
  {
  public:

    bool operator()(const Xy_monotone_surface_3& s) const
    {
      return s.is_vertical();
    }
  };

  /*! Get a Is_vertical_3 functor object. */
  Is_vertical_3 is_vertical_3_object() const
  {
    return Is_vertical_3();
  }
  
  /***************************************************************************/

  // public method needed for testing

  // checks if point is in the xy-range of surf
  class Is_defined_over
  {
  public:
    // checks if point is in the xy-range of surf
    bool operator()(const Point_2& point, 
		    const Xy_monotone_surface_3& surf) const

    {
      Kernel k;
      Self parent;

      // project the surface on the plane
      Triangle_2 boundary = parent.project(surf);

      // if surface is not vertical (i.e. boundary is not degenerate)
      // check if the projected point is inside the projected boundary
      if (!k.is_degenerate_2_object()(boundary))
        return (!k.has_on_unbounded_side_2_object()(boundary, point));

      // if surface is vertical, we check if the point is collinear
      // with the projected vertices, and on one of the projected segments
      // of the boundary
      Point_2 v1 = k.construct_vertex_2_object()(boundary, 0);
      Point_2 v2 = k.construct_vertex_2_object()(boundary, 1);
      Point_2 v3 = k.construct_vertex_2_object()(boundary, 2);

      if (!k.collinear_2_object()(v1, v2, point))
        return false;

      // enough to check 2 edges, because the 3rd is part of their union
      return (k.collinear_are_ordered_along_line_2_object()(v1, point, v2) ||
              k.collinear_are_ordered_along_line_2_object()(v2, point, v3));

    }
  };

  /*! Get a Is_defined_over functor object. */
  Is_defined_over is_defined_over_object() const
  {
    return Is_defined_over();
  }


  Curve_2 project(const Segment_3& segment_3) const
  {
    Kernel k;
    Construct_vertex_3 vertex_on = k.construct_vertex_3_object();

    const Point_3&  end1 = vertex_on(segment_3, 0),
                    end2 = vertex_on(segment_3, 1);
    Point_2 projected_end1(end1.x(), end1.y()),
            projected_end2(end2.x(), end2.y());
    return Curve_2(projected_end1, projected_end2);
  }

  Point_2 project(const Point_3& obj) const
  {
    return Point_2(obj.x(), obj.y());
  }
  
  Triangle_2 project(const Xy_monotone_surface_3& triangle_3) const
  {
    const Point_3&  end1 = triangle_3.vertex(0),
                    end2 = triangle_3.vertex(1),
                    end3 = triangle_3.vertex(2);
    Point_2 projected_end1(end1.x(), end1.y()),
            projected_end2(end2.x(), end2.y()),

            projected_end3(end3.x(), end3.y());
    return Triangle_2(projected_end1, projected_end2, projected_end3);
  }

  // triangles overlap if they lie on the same plane and intersect on it.
  bool do_overlap(const Xy_monotone_surface_3& s1, 
		  const Xy_monotone_surface_3& s2) const
  {
    Kernel k;
    if (!k.do_intersect_3_object()(s1, s2))
      return false;

    // check if they are coplanar
    Point_3 a1 = s1.vertex(0),
            b1 = s1.vertex(1),
            c1 = s1.vertex(2);
    Point_3 a2 = s2.vertex(0),
            b2 = s2.vertex(1),
            c2 = s2.vertex(2);
    bool b = k.coplanar_3_object()(a1, b1, c1, a2);
    if (!b) return false;

    b = k.coplanar_3_object()(a1, b1, c1, b2);
    if (!b) return false;

    b = k.coplanar_3_object()(a1, b1, c1, c2);
    return b;    
  }

  // intersect two 3D-triangles
  // the result can be a triangle, a segment, a polygon or a point
  Object intersection(const Xy_monotone_surface_3& s1, 
                      const Xy_monotone_surface_3& s2) const
  {
    CGAL_precondition( !s1.is_degenerate() );
    CGAL_precondition( !s2.is_degenerate() );
    Kernel k;

    // first, try to intersect the bounding boxes of the triangles,
    // eficiently return empty object when the triangles are faraway
    if (!CGAL::do_overlap(s1.bbox(), s2.bbox()))
      return Object();

    // if both triangle lie on the same plane, should calculate the intersection
    // as in the case of dimention 2
    Plane_3 p1 = s1.plane();
    Plane_3 p2 = s2.plane();
    if  (p1 == p2 || p1 == p2.opposite())
    {
      if (s1.is_vertical())
      {
        inter_overlap_timer.start();
        Object res = intersection_on_plane_3(p1, s1, s2);
        inter_overlap_timer.stop();
        return res;
      }
      else
        // we don't care about overlap, because they are not passed to the
        // algorithm anyway, so we save the costly computation
        return Object();
    }

    // calculate intersection between a triangle and the other triangle's 
    // supporting plane
    // if there is no intersection - then the triangles have no intersection 
    // between them.
    inter_plane_tri_timer.start();
    Object inter_obj = intersection(p1, s2);
    inter_plane_tri_timer.stop();
    if (inter_obj.is_empty())
      return Object();

    // otherwise, if the intersection in a point, we should check if it lies
    // inside the first triangle
    Assign_3 assign_obj = k.assign_3_object();
    Point_3 inter_point;
    if (assign_obj(inter_point, inter_obj))
    {
      inter_on_plane_tri_pt_timer.start();
      Object res = intersection_on_plane_3(p1, s1, inter_point);
      inter_on_plane_tri_pt_timer.stop();
      return res;
    }
    else
    {
      // if the intersection is a segment, we check the intersection of the
      // other plane-triangle pair
      Segment_3 inter_seg;
      CGAL_assertion(assign_obj(inter_seg, inter_obj));
      assign_obj(inter_seg, inter_obj);

      inter_plane_tri_timer.start();
      inter_obj = intersection(p2, s1);
      inter_plane_tri_timer.stop();

      // if there is no intersection - then the triangles have no intersection 
      // between them.
      if (inter_obj.is_empty())
	return Object();
      
      if (assign_obj(inter_point, inter_obj))
      {
	// if the intersection is a point, which lies on the segment, 
	// than it is the result,
	// otherwise, empty result
	 if (k.has_on_3_object()(inter_seg, inter_point))
	   return make_object(inter_point);
	 else
	   return Object();
      }
      else
      {
	// both plane-triangle intersections are segments, which are collinear,
	// and lie on the line which is the intersection of the two supporting
	// planes
        Segment_3 inter_seg2;
	CGAL_assertion(assign_obj(inter_seg2, inter_obj));
	assign_obj(inter_seg2, inter_obj);
	
	inter_on_plane_seg_seg_timer.start();

	Point_3 min1 = k.construct_min_vertex_3_object()(inter_seg),
  	        max1 = k.construct_max_vertex_3_object()(inter_seg); 
	Point_3 min2 = k.construct_min_vertex_3_object()(inter_seg2),
	        max2 = k.construct_max_vertex_3_object()(inter_seg2); 

	CGAL_assertion((k.collinear_3_object()(min1, min2, max1) &&
                  k.collinear_3_object()(min1, max2, max1)));
 
	// we need to find the overlapping part, if exists
	Point_3 min, max;
	if (k.less_xyz_3_object()(min1, min2))
	  min = min2;
	else
	  min = min1;
	if (k.less_xyz_3_object()(max1, max2))
	  max = max1;
	else
	  max = max2;
	
	Object res;
	Comparison_result comp_res = k.compare_xyz_3_object()(min, max);
	if (comp_res == EQUAL)
	  res = make_object(min);
	else if (comp_res == SMALLER)
	  res = make_object(Segment_3(min, max));
	// else - empty result
	inter_on_plane_seg_seg_timer.stop();
	return res;
      }
    }
  }

  // intersect two 3D-triangles - old and slow method!
  // the result can be a triangle, a segment, a polygon or a point
  Object intersection_old(const Xy_monotone_surface_3& s1, 
                          const Xy_monotone_surface_3& s2) const
  {
    CGAL_precondition( !s1.is_degenerate() );
    CGAL_precondition( !s2.is_degenerate() );
    Kernel k;

    // first, try to intersect the bounding boxes of the triangles,
    // eficiently return empty object when the triangles are faraway
    if (!CGAL::do_overlap(s1.bbox(), s2.bbox()))
      return Object();

    // if both triangle lie on the same plane, should calculate the intersection
    // as in the case of dimention 2
    Plane_3 p1 = s1.plane();
    Plane_3 p2 = s2.plane();
    if  (p1 == p2 || p1 == p2.opposite())
    {
      inter_overlap_timer.start();
      Object res = intersection_on_plane_3(p1, s1, s2);
      inter_overlap_timer.stop();
      return res;
    }

    // calculate intersection between a triangle and the other triangle's supporting plane
    // if there is no intersection - then the triangles have no intersection between them.
    inter_plane_tri_timer.start();
    Object inter_obj = intersection(p1, s2);
    inter_plane_tri_timer.stop();
    if (inter_obj.is_empty())
      return Object();

    // otherwise, we should intersect the result (which is a point or a segment only)
    // with the triangle (on the same plane in 3d)
    Assign_3 assign_obj = k.assign_3_object();
    Point_3 inter_point;
    if (assign_obj(inter_point, inter_obj))
    {
      inter_on_plane_tri_pt_timer.start();
      Object res = intersection_on_plane_3(p1, s1, inter_point);
      inter_on_plane_tri_pt_timer.stop();
      return res;
    }
    else
    {
      Segment_3 inter_seg;
      CGAL_assertion(assign_obj(inter_seg, inter_obj));
      assign_obj(inter_seg, inter_obj);
      inter_on_plane_tri_seg_timer.start();
      Object res = intersection_on_plane_3(p1, s1, inter_seg);
      inter_on_plane_tri_seg_timer.stop();
      return res;
    }
  }

  // calculate intersection between triangle & point on the same plane plane
  Object intersection_on_plane_3(const Plane_3& plane,
                                 const Xy_monotone_surface_3& triangle,
                                 const Point_3& point) const
  {
    Kernel k;
    CGAL_precondition( !triangle.is_degenerate() );
    CGAL_precondition( !k.is_degenerate_3_object()(plane) );
    CGAL_precondition( triangle.plane() == plane ||
                       triangle.plane() == plane.opposite());
    CGAL_precondition( k.has_on_3_object()(plane, point) );

    // if the point is inside the triangle, then the point is the intersection
    // otherwise there is no intersection
    if (k.has_on_3_object()(triangle, point))
      return make_object(point);
    else
      return Object();
  }

  // calculate intersection between triangle & segment on the same plane plane
  Object intersection_on_plane_3(const Plane_3& plane,
                                 const Xy_monotone_surface_3& triangle,
                                 const Segment_3& segment) const
  {
    Kernel k;
    CGAL_precondition( !triangle.is_degenerate() );
    CGAL_precondition( !k.is_degenerate_3_object()(plane) );
    CGAL_precondition( triangle.plane() == plane ||
                       triangle.plane() == plane.opposite());
    CGAL_precondition( k.has_on_3_object()
                          (plane, k.construct_vertex_3_object()(segment, 0)) );
    CGAL_precondition( k.has_on_3_object()
                          (plane, k.construct_vertex_3_object()(segment, 1)) );

    Construct_vertex_3 vertex_on_3 = k.construct_vertex_3_object();
    
    // for simplicity, we transform the triangle & segment to the xy-plane,
    // compute the intersection there, and transform it back to the 3d plane.
    Point_2 src_t = plane.to_2d(vertex_on_3(segment, 0)),
            tar_t = plane.to_2d(vertex_on_3(segment, 1));
    Segment_2 segment_t(src_t, tar_t);

    Point_2 v1 = plane.to_2d(triangle.vertex(0)),
            v2 = plane.to_2d(triangle.vertex(1)),
            v3 = plane.to_2d(triangle.vertex(2));
    Triangle_2 triangle_t(v1, v2, v3);

	  Object inter_obj = k.intersect_2_object()(triangle_t, segment_t);
	  Assign_2 assign_2 = k.assign_2_object();
	  if (inter_obj.is_empty())
		  return inter_obj;

  	Point_2 inter_point;
  	if (assign_2(inter_point, inter_obj))
  	   return make_object(plane.to_3d(inter_point));
  	else
  	{
	    // intersection is a segment
      Segment_2 inter_segment;
      CGAL_assertion(assign_2(inter_segment, inter_obj));
	    assign_2(inter_segment, inter_obj);
	    return make_object(Segment_3(plane.to_3d(k.construct_vertex_2_object()(inter_segment, 0)),
                                 plane.to_3d(k.construct_vertex_2_object()(inter_segment, 1))));


	  }

	  return Object();
  }

  // calculate intersection between 2 triangle on the same plane plane
  Object intersection_on_plane_3(const Plane_3& plane,
                                 const Xy_monotone_surface_3& tri1,
                                 const Xy_monotone_surface_3& tri2) const
  {
    Kernel k;
    CGAL_precondition( !tri1.is_degenerate() );
    CGAL_precondition( !tri2.is_degenerate() );
    CGAL_precondition( !k.is_degenerate_3_object()(plane) );
    CGAL_precondition( tri1.plane() == plane ||
                       tri1.plane() == plane.opposite());
    CGAL_precondition( tri2.plane() == plane ||
                       tri2.plane() == plane.opposite());

  	// for simplicity, we transform the triangles to the xy-plane,
  	// compute the intersection there, and transform it back to the 3d plane.

    Point_2 v1 = plane.to_2d(tri1.vertex(0)),
            v2 = plane.to_2d(tri1.vertex(1)),
            v3 = plane.to_2d(tri1.vertex(2));
    Triangle_2 triangle1_t(v1, v2, v3);

  	Point_2 u1 = plane.to_2d(tri2.vertex(0)),
            u2 = plane.to_2d(tri2.vertex(1)),
            u3 = plane.to_2d(tri2.vertex(2));
  	Triangle_2 triangle2_t(u1, u2, u3);

  	Object inter_obj = k.intersect_2_object()(triangle1_t, triangle2_t);
  	Assign_2 assign_2 = k.assign_2_object();
  	if (inter_obj.is_empty())
  		return inter_obj;

  	Point_2 inter_point;
    Segment_2 inter_segment;
    Triangle_2 inter_tri;

   if (assign_2(inter_point, inter_obj))
  	   return make_object(plane.to_3d(inter_point));
  	else if (assign_2(inter_segment, inter_obj))
  	   return make_object(Segment_3(
                  plane.to_3d(k.construct_vertex_2_object()(inter_segment, 0)),
									plane.to_3d(k.construct_vertex_2_object()(inter_segment, 1))));

  	else if (assign_2(inter_tri, inter_obj))
  	   return make_object(Triangle_3(
                  plane.to_3d(k.construct_vertex_2_object()(inter_tri, 0)),
  								plane.to_3d(k.construct_vertex_2_object()(inter_tri, 1)),
  								plane.to_3d(k.construct_vertex_2_object()(inter_tri, 2))));
  	else
  	{
  	   // intersection is a polygon, given as a vector of points
  		std::vector<Point_2> inter_poly;
  		assign_2(inter_poly, inter_obj);
  		std::vector<Point_3> poly(inter_poly.size());

      for(unsigned int i=0; i<inter_poly.size(); ++i)
  			poly[i] = plane.to_3d(inter_poly[i]);

  		return make_object(poly);
  	}
  	return Object();

  }

  // calculate the intersection between a (non degenerate) triangle
  // and a (non degenerate) plane in 3d
  // the result object can be empty, a point, a segment or the original
  // triangle
  Object intersection(const Plane_3& pl, 
		                  const Xy_monotone_surface_3& tri) const
  {
    Kernel k;
    CGAL_precondition( !tri.is_degenerate() );
    CGAL_precondition( !k.is_degenerate_3_object()(pl) );

    // first, check for all 3 vertices of tri on which side of pl they lie on
    int points_on_plane[3]; // contains the indices of vertices that lie on pl
    int points_on_positive[3]; // contains the indices of vertices that lie on the positive side of pl
    int points_on_negative[3]; // contains the indices of vertices that lie on the negative side of pl

    int n_points_on_plane = 0;
    int n_points_on_positive = 0;
    int n_points_on_negative = 0;

    Oriented_side side;
    for (int i=0; i<3; ++i)
    {
      side = pl.oriented_side(tri.vertex(i));
      if (side == ON_NEGATIVE_SIDE)
        points_on_negative[n_points_on_negative++] = i;
      else if (side == ON_POSITIVE_SIDE)
        points_on_positive[n_points_on_positive++] = i;
      else
        points_on_plane[n_points_on_plane++] = i;
    }

    assert(n_points_on_plane + n_points_on_positive + n_points_on_negative == 3);

    // if all vertices of tri lie on the same size (positive/negative) of pl,
    // there is no intersection
    if (n_points_on_positive == 3 || n_points_on_negative == 3)
      return Object();

    // if all vertices of tri lie on pl then we return tri
    if (n_points_on_plane == 3)
       return make_object(tri);

    // if 2 vertices lie on pl, then return the segment between them
    if (n_points_on_plane == 2)
    {
      int point_idx1 = points_on_plane[0], point_idx2 = points_on_plane[1];
      return make_object(Segment_3(tri.vertex(point_idx1),tri.vertex(point_idx2)));
    }

    // if only 1 lie on pl, should check the segment opposite to it on tri
    if (n_points_on_plane == 1)
    {
      int point_on_plane_idx = points_on_plane[0];

      // if the other 2 vertices are on the same side of pl,
      // then the answer is just this vertex
      if (n_points_on_negative == 2 || n_points_on_positive == 2)
        return make_object(tri.vertex(point_on_plane_idx));

      // now it is known that one vertex is on pl, and the segment of tri
      // opposite to it should intersect pl

      // the segment of tri opposite of tri[point_on_plane_idx]
      Segment_3 tri_segment(tri.vertex(point_on_plane_idx+1),
                            tri.vertex(point_on_plane_idx+2));

      Object inter_result = k.intersect_3_object()(pl, tri_segment);
      Point_3 inter_point;
      CGAL_assertion( k.assign_3_object()(inter_point, inter_result) );
      k.assign_3_object()(inter_point, inter_result);

      // create the resulting segment
      // (between tri[point_on_plane_idx] and inter_point)
      return make_object(Segment_3(tri.vertex(point_on_plane_idx), inter_point));

    }

    assert( n_points_on_plane == 0 );
    assert( n_points_on_positive + n_points_on_negative == 3 );
    assert( n_points_on_positive != 0 );
    assert( n_points_on_negative != 0 );

    // now it known that there is an intersection between 2 segments of tri
    // and pl, it is also known which segments are those.
    Point_3 inter_points[2];
    int pos_it, neg_it, n_inter_points = 0;
    for(pos_it = 0; pos_it < n_points_on_positive; ++pos_it)
      for(neg_it = 0; neg_it < n_points_on_negative; ++neg_it)
      {
        Segment_3 seg(tri.vertex(points_on_positive[pos_it]),
                      tri.vertex(points_on_negative[neg_it]));
        Object inter_result = k.intersect_3_object()(pl, seg);
        Point_3 inter_point;
        // the result of the intersection must be a point
        CGAL_assertion( k.assign_3_object()(inter_point, inter_result) );
        k.assign_3_object()(inter_point, inter_result);
        inter_points[n_inter_points++] = inter_point;
      }

    assert( n_inter_points == 2 );
    return make_object(Segment_3(inter_points[0], inter_points[1]));
  }

  // compare the value of s1 in p1 to the value of s2 in p2
  Comparison_result
  compare_z(const Point_2& p1,
            const Xy_monotone_surface_3& s1,
            const Point_2& p2,
            const Xy_monotone_surface_3& s2)
  {
    CGAL_precondition(is_defined_over_object()(p1, s1));
    CGAL_precondition(is_defined_over_object()(p2, s2));

    Point_3 v1 = envelope_point_of_surface(p1, s1);
    Point_3 v2 = envelope_point_of_surface(p2, s2);
    Kernel k;
    return k.compare_z_3_object()(v1, v2);
  }
  
  // find the envelope point of the surface over the given point
  // precondition: the surface is defined in point
  Point_3
  envelope_point_of_surface(const Point_2& p,
                            const Xy_monotone_surface_3& s) const
  {
    CGAL_precondition(is_defined_over_object()(p, s));

    // Construct a vertical line passing through point
    Kernel k;
    Point_3 point(p.x(), p.y(), 0);
    Direction_3 dir (0, 0, 1);
    Line_3      vl = k.construct_line_3_object() (point, dir);

    // Compute the intersetion between the vertical line and the given surfaces
    if (s.is_vertical())
      return envelope_point_of_vertical_surface(point, s);
    else
    {
      const Plane_3& plane = s.plane();
      Object    res = k.intersect_3_object()(plane, vl);
      CGAL_assertion(!res.is_empty());
      Point_3 ip;
      CGAL_assertion(k.assign_3_object()(ip, res));
      k.assign_3_object()(ip, res);

      return ip;
    }
  }

  // find the envelope point of the surface over the given point
  // precondition: the surface is defined in point and is vertical
  Point_3
  envelope_point_of_vertical_surface(const Point_3& point,
                                     const Xy_monotone_surface_3& surf) const
  {
    Kernel k;
    CGAL_precondition(surf.is_vertical());
    CGAL_precondition(is_defined_over_object()(project(point), surf));

    const Plane_3& plane = surf.plane();


    // Construct a vertical line passing through point
    Direction_3 dir (0, 0, 1);
    Line_3      vl = k.construct_line_3_object() (point, dir);
    // we need 2 points on this line, to be transformed to 2d,
    // and preserve the direction of the envelope
    Point_3 vl_point1 = k.construct_point_on_3_object()(vl, 0),
            vl_point2 = k.construct_point_on_3_object()(vl, 1);

    // the surface and the line are on the same plane(plane),
    // so we transform them to the xy-plane, compute the intersecting segment
    // and transform it back to plane.
    Point_3 v1 = k.construct_vertex_3_object()(surf, 0);
    Point_3 v2 = k.construct_vertex_3_object()(surf, 1);
    Point_3 v3 = k.construct_vertex_3_object()(surf, 2);

    Point_2 t1 = plane.to_2d(v1);
    Point_2 t2 = plane.to_2d(v2);
    Point_2 t3 = plane.to_2d(v3);

    Point_2 tvl_point1 = plane.to_2d(vl_point1);
    Point_2 tvl_point2 = plane.to_2d(vl_point2);
    Line_2 l(tvl_point1, tvl_point2);

    Triangle_2 tri(t1, t2, t3);
    Object inter_obj = k.intersect_2_object()(tri, l);
    Segment_2 inter_seg;
    Point_2 inter_point;
    bool is_inter_point = k.assign_2_object()(inter_point, inter_obj);
    if (is_inter_point)
       return plane.to_3d(inter_point);

    // intersection must be a segment because the transformation keeps
    // combinatorics
    CGAL_assertion(k.assign_2_object()(inter_seg, inter_obj));
    k.assign_2_object()(inter_seg, inter_obj);

    // now find the vertex of inter_seg which is the transformed envelope
    // point


    #ifdef CGAL_DEBUG_ENVELOPE_3_TRIANGLES_TRAITS
      std::cout << "in envelope_point_of_vertical_surface - segment is " 
                << plane.to_3d(k.construct_min_vertex_2_object()(inter_seg))
                << " --> "
                << plane.to_3d(k.construct_max_vertex_2_object()(inter_seg))
                << endl;

    #endif
    
    Point_2 result;
    Segment_3 seg_3(plane.to_3d(k.construct_vertex_2_object()(inter_seg, 0)),
                    plane.to_3d(k.construct_vertex_2_object()(inter_seg, 1)));
    if (get_envelope_type() == LOWER)
      return k.construct_min_vertex_3_object()(seg_3);
    else
      return k.construct_max_vertex_3_object()(seg_3);
    
  }

//  // return true if the surface is ortogonal to the xy-plane
//  bool surface_is_vertical(const Xy_monotone_surface_3& s) const
//  {
//     CGAL_precondition(!s.is_degenerate());
//     return project(s).is_degenerate();
//  }

  Point_2 construct_middle_point(const Point_2& p1, const Point_2& p2) const
  {
    Kernel k;
    return k.construct_midpoint_2_object()(p1, p2);
  }

  Point_2 construct_middle_point(const X_monotone_curve_2& cv) const
  {
    Kernel k;
//    Construct_vertex_2 construct_vertex = k.construct_vertex_2_object();
//    return k.construct_midpoint_2_object()(construct_vertex(cv, 0), construct_vertex(cv, 1));
    return k.construct_midpoint_2_object()(cv.source(), cv.target());
  }

  /***************************************************************************/
  // for vertical decomposition
  /***************************************************************************/
  
  class Construct_vertical_2
  {
  public:
    X_monotone_curve_2 operator()(const Point_2& p1, const Point_2& p2) const
    {
      return X_monotone_curve_2(p1, p2);

    }
  };

  /*! Get a Construct_vertical_2 functor object. */
  Construct_vertical_2 construct_vertical_2_object() const
  {
    return Construct_vertical_2();
  }
 

  Point_2 vertical_ray_shoot_2 (const Point_2& pt,
                                const X_monotone_curve_2& cv) const
  {
    CGAL_precondition(!cv.is_vertical());
    Kernel k;
    // If the curve contains pt, return it.
    if (k.has_on_2_object() (cv, pt))
      return (pt);

    // Construct a vertical line passing through pt.
    typename Kernel::Direction_2  dir (0, 1);
    typename Kernel::Line_2        vl = k.construct_line_2_object() (pt, dir);

    // Compute the intersetion between the vertical line and the given curve.
    Object    res = k.intersect_2_object()(cv, vl);
    Point_2   ip;
    bool      ray_shoot_successful = k.assign_2_object()(ip, res);

    if (! ray_shoot_successful)
      CGAL_assertion (ray_shoot_successful);

    return (ip);
  }

  Envelope_triangles_traits_3() : type(LOWER), cache_hits(0)
  {}

  Envelope_triangles_traits_3(Envelope_type t) : type(t), cache_hits(0)
  {}

  Envelope_type get_envelope_type() const
  {
    return type;
  }


  
  virtual ~Envelope_triangles_traits_3() {}

  void print_times()
  {
    std::cout << total_timer.intervals()
              << " number of calls to traits interface functions" << std::endl;
    std::cout << "total time: " << total_timer.time() << " seconds"
              << std::endl << std::endl;



    std::cout << "projected boundary        #: "
              << pboundary_timer.intervals()
              << " total time: " << pboundary_timer.time()
              << " average time: " 
              << pboundary_timer.time()/pboundary_timer.intervals()
              << std::endl;
    std::cout << "projected intersections   #: "
              << intersection_timer.intervals()
              << " total time: " << intersection_timer.time()
              << " average time: " 
              << intersection_timer.time()/intersection_timer.intervals()
              << std::endl;              
    std::cout << "compare_distance on point #: " 
              << compare_timer.intervals()
              << " total time: " << compare_timer.time()
              << " average time: " 
              << compare_timer.time()/compare_timer.intervals()
              << std::endl;
//    std::cout << "total time spent on computing z at xy: " 
//              << compute_z_at_xy_timer.time()
//              << " seconds, hits# = " << cache_hits << std::endl;
    std::cout << "compare_distance on curve #: "
              << compare_on_cv_timer.intervals()
              << " total time: " << compare_on_cv_timer.time()
              << " average time: " 
              << compare_on_cv_timer.time()/compare_on_cv_timer.intervals()
              << std::endl;
    std::cout << "compare_distance on side  #: "
              << compare_side_timer.intervals()
              << " total time: " << compare_side_timer.time()
              << " average time: " 
              << compare_side_timer.time()/compare_side_timer.intervals()
              << std::endl;
    std::cout << std::endl;              

#ifdef PRINT_TRAINGLES_INTERSECTION_STATS
    std::cout << "Projected intersection calculation statistics: " << std::endl
              << "inter plane & triangle: " << inter_plane_tri_timer.time()
              << std::endl
              << "overlap: " << inter_overlap_timer.time() << std::endl
              << "inter on plane triangle & point: " 
              << inter_on_plane_tri_pt_timer.time() << std::endl
              << "inter on plane triangle & segment: " 
              << inter_on_plane_tri_seg_timer.time() << std::endl
              << "inter two segments: " 
              << inter_on_plane_seg_seg_timer.time() << std::endl;
    std::cout << std::endl;        
#endif      
              
  }

  /*! reset statistics */
  void reset()
  {
    total_timer.reset();
    intersection_timer.reset();
    pboundary_timer.reset();
    compare_timer.reset();
    compute_z_at_xy_timer.reset();
    compare_on_cv_timer.reset();
    compare_side_timer.reset();

    inter_plane_tri_timer.reset();
    inter_overlap_timer.reset();
    inter_on_plane_tri_pt_timer.reset();
    inter_on_plane_tri_seg_timer.reset();
    inter_on_plane_seg_seg_timer.reset();

    cache_hits = 0;
  }
  
protected:
  Envelope_type type;

public:
  // measure times:
  // measure the total time for all interface methods
  mutable Timer total_timer;
  // measure the total time for the intersection method
  mutable Timer intersection_timer;
  // measure total time for projected boundary method
  mutable Timer pboundary_timer;
  // measure the total time for the compare_distance_to_envelope method
  mutable Timer compare_timer;
  mutable Timer compute_z_at_xy_timer;
  mutable Timer compare_on_cv_timer;

  // measure the total time for the compare_distance_to_envelope_side method
  mutable Timer compare_side_timer;

  // measure intersection calculation parts, to see where is the weak part
  mutable Timer inter_plane_tri_timer;
  mutable Timer inter_overlap_timer;
  mutable Timer inter_on_plane_tri_pt_timer;
  mutable Timer inter_on_plane_tri_seg_timer;
  mutable Timer inter_on_plane_seg_seg_timer;

  #ifdef CGAL_ENV_TRIANGLES_TRAITS_CACHE_POINT_ON
    mutable Surface_point_cache point_on_cache;
  #endif
  mutable int     cache_hits;

};


/*!
 * \class A representation of a triangle, as used by the 
 * Envelope_triangles_traits_3 traits-class.
 */
template <class Kernel_>
class Envelope_triangle_3 :
    public Handle_for<typename Envelope_triangles_traits_3<Kernel_>::_Triangle_cached_3>
{
  typedef Kernel_                                                  Kernel;
  typedef typename Kernel::Triangle_3                              Triangle_3;
  typedef typename Kernel::Point_3                                 Point_3;
  typedef typename Kernel::Plane_3                                 Plane_3;

  typedef typename Envelope_triangles_traits_3<Kernel>::_Triangle_cached_3
                                                                   Rep_type;

public:

  /*!
   * Default constructor.
   */
  Envelope_triangle_3() :
    Handle_for<Rep_type>(Rep_type())
  {}

  /*!
   * Constructor from a "kernel" triangle.
   * \param seg The segment.
   */
  Envelope_triangle_3(const Triangle_3& tri) :
    Handle_for<Rep_type>(Rep_type(tri))
  {}

  /*!
   * Construct a triangle from 3 end-points.
   * \param p1 The first point.
   * \param p2 The second point.
   * \param p3 The third point.
   */
    Envelope_triangle_3(const Point_3 &p1, const Point_3 &p2,
                        const Point_3 &p3) :
      Handle_for<Rep_type>(Rep_type(p1, p2, p3))
  {}

  /*!
   * Construct a triangle from a plane and 3 end-points.
   * \param pl The supporting plane.
   * \param p1 The first point.
   * \param p2 The second point.
   * \param p3 The third point.
   * \pre All points must be on the supporting plane.
   */
  Envelope_triangle_3(const Plane_3& pl,
                      const Point_3 &p1, const Point_3 &p2,
                      const Point_3 &p3) :
    Handle_for<Rep_type>(Rep_type(pl, p1, p2, p3))

  {}

  /*!
   * Cast to a triangle.
   */
  operator Triangle_3() const
  {
    return (Triangle_3(ptr()->vertex(0), ptr()->vertex(1), ptr()->vertex(2)));
  }

  /*!
   * Create a bounding box for the triangle.
   */
  Bbox_3 bbox() const
  {
    Triangle_3 tri(ptr()->vertex(0), ptr()->vertex(1), ptr()->vertex(2));
    return (tri.bbox());
  }

  /*!
   * Get the ith endpoint.
   */
  const Point_3& vertex(unsigned int i) const
  {
    return ptr()->vertex(i);
  }

  /*!
   * Get the supporting plane.
   */
  const Plane_3& plane() const
  {
    return ptr()->plane();
  }

  /*!
   * Check if the triangel is vertical.
   */
  bool is_vertical() const
  {
    return ptr()->is_vertical();
  }

  /*!
   * Check if the triangle is degenerate.
   */
  bool is_degenerate() const
  {
    return ptr()->is_degenerate();
  }

};

template <class Kernel>
bool
operator<(const Envelope_triangle_3<Kernel> &a,
          const Envelope_triangle_3<Kernel> &b)
{
  return (a.id() < b.id());
}
template <class Kernel>
bool
operator==(const Envelope_triangle_3<Kernel> &a,
           const Envelope_triangle_3<Kernel> &b)
{
  return (a.id() == b.id());
}

/*!
 * Exporter for the triangle class used by the traits-class.
 */
template <class Kernel, class OutputStream>
OutputStream& operator<< (OutputStream& os, const Envelope_triangle_3<Kernel>& tri)
{
  os << static_cast<typename Kernel::Triangle_3>(tri);
  return (os);
}

/*!
 * Importer for the triangle class used by the traits-class.
 */
template <class Kernel, class InputStream>
InputStream& operator>> (InputStream& is, Envelope_triangle_3<Kernel>& tri)
{
  typename Kernel::Triangle_3   kernel_tri;
  is >> kernel_tri;
  tri = kernel_tri;
  return (is);
}

CGAL_END_NAMESPACE

#endif 
