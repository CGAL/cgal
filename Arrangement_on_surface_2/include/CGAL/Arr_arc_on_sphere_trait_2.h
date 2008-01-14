// Copyright (c) 2008 INRIA Sophia-Antipolis
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
// $URL: 
// $Id: 
// 
//
// Author(s)     : Sebastien Loriot         <sebastien.loriot@sophia.inria.fr>

#ifndef CGAL_ARR_ARC_ON_SPHERE_TRAITS_2_H
#define CGAL_ARR_ARC_ON_SPHERE_TRAITS_2_H


/*! \file
 * A class that handles circular arcs embedded on spheres suitable
 * as a geometry traits class for the arrangement on surface package.
 */

#include <CGAL/tags.h>
#include <CGAL/intersections.h>
#if defined(CGAL_ARR_PLANE)
#include <CGAL/Arr_geometry_traits/Arr_plane_3.h>
#endif
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/Spherical_kernel_3.h>


#include <fstream>

CGAL_BEGIN_NAMESPACE

/*! A traits class-template for constructing and maintaining arcs of great
 * circles embedded on spheres. It is parameterized from a (linear) geometry
 * kernel, which it also derives from
 */
template <typename T_Kernel>
class Arr_arc_on_sphere_traits_2 : public T_Kernel {
  friend class Arr_x_monotone_great_circular_arc_on_sphere_3<T_Kernel>;
  friend class Arr_geodesic_arc_on_sphere_3<T_Kernel>;
  friend class Arr_extended_direction_3<T_Kernel>;

public:
  typedef T_Kernel                              Kernel;
  
  // Category tags:
  typedef Tag_true                              Has_left_category;
  typedef Tag_true                              Has_merge_category;
  typedef Arr_bounded_boundary_tag              Boundary_category;

  /*! Default constructor */
  Arr_arc_on_sphere_traits_2(){}

protected:
  typedef typename Kernel::Vector_3             Vector_3;
  typedef typename Kernel::Line_3               Line_3;
      
  typedef typename Kernel::Plane_3              Plane_3;
  typedef typename Kernel::Point_3              Point_3;

    
public:
  
  typedef CGAL::Algebraic_kernel_for_spheres_2_3<typename Kernel::FT>       CGAL_AK;
  typedef CGAL::Spherical_kernel_3<Kernel,CGAL_AK>                          CGAL_SK;
  
  // Traits objects
  typedef Circular_arc_point_on_reference_sphere_3<Kernel>              Point_2;
  typedef Circular_arc_on_reference_sphere_3<Kernel>                 X_monotone_curve_2;
  typedef Circular_arc_on_reference_sphere_3<Kernel>                 Curve_2;
  typedef unsigned int                                  Multiplicity;

  CGAL_SK Spherical_kernel;

public:
  /// \name Basic functor definitions
  //@{

  /*! A functor that compares the x-coordinates of two directional points */
  class Compare_x_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_x_2(const Traits * traits) : m_traits(traits) {}

    friend class Arr_arc_on_sphere_traits_2<Kernel>;
    
  public:    
    /*! Compare the x-coordinates of two directional points.
     * \param p1 the first directional point.
     * \param p2 the second directional point.
     * \return SMALLER - x(p1) < x(p2);
     *         EQUAL   - x(p1) = x(p2);
     *         LARGER  - x(p1) > x(p2).
     * \pre p1 does not lie on the closed discontinuity arc.
     * \pre p2 does not lie on the closed discontinuity arc.
     */
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      CGAL_precondition(p1.is_no_boundary());
      CGAL_precondition(p2.is_no_boundary());

      return m_traits->Spherical_kernel.compare_theta_3_object()(p1, p2);
    }
  };

  /*! Obtain a Compare_x_2 function object */
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(this); }

  /*! A functor that compares two directional points lexigoraphically:
   * by x, then by y.
   */
  class Compare_xy_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_xy_2(const Traits * traits) : m_traits(traits) {}

    friend class Arr_arc_on_sphere_traits_2<Kernel>;
    
  public:
    /*! Compare two directional points lexigoraphically: by x, then by y.
     * \param p1 the first enpoint directional point.
     * \param p2 the second endpoint directional point.
     * \return SMALLER - x(p1) < x(p2);
     *         SMALLER - x(p1) = x(p2) and y(p1) < y(p2);
     *         EQUAL   - x(p1) = x(p2) and y(p1) = y(p2);
     *         LARGER  - x(p1) = x(p2) and y(p1) > y(p2);
     *         LARGER  - x(p1) > x(p2).
     * \pre p1 is not on the discontinuity arc.
     * \pre p2 is not on the discontinuity arc.
     */
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      CGAL_precondition(p1.is_no_boundary());
      CGAL_precondition(p2.is_no_boundary());
      
      return m_traits->Spherical_kernel.compare_theta_z_3_object()(p1, p2);
    }
  };

  /*! Obtain a Compare_xy_2 function object */
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(this); }

  /*! A functor that obtain the left endpoint of an x-monotone arc */
  class Construct_min_vertex_2 {
  public:
    /*! Obtain the left endpoint of and arc.
     * \param xc the arc.
     * \return the left endpoint.
     */
    const Point_2 & operator()(const X_monotone_curve_2 & xc) const
    {
      CGAL::Comparison_result res=SK.compare_theta_z_3_object()(xc.source(),xc.target());
      if (res == CGAL::LARGER)
        return xc.target()
      return xc.source(); 
    }
  };

  /*! Obtain a Construct_min_vertex_2 function object */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(); }

  /*! A functor that obtain the right endpoint of an x-monotone arc */
  class Construct_max_vertex_2 {
  public:
    /*! Obtain the right endpoint of an arc.
     * \param xc the arc.
     * \return the right endpoint.
     */
    const Point_2 & operator()(const X_monotone_curve_2 & xc) const
    {
      CGAL::Comparison_result res=SK.compare_theta_z_3_object()(xc.source(),xc.target());
      if (res == CGAL::LARGER)
        return xc.source()
      return xc.target(); 
    }
  };

  /*! Obtain a Construct_max_vertex_2 function object */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(); }

  /*! A functor that checks whether an x-monotone arc is a vertical */
  class Is_vertical_2 {
  public:
    /*! Check whether a given arc is a vertical.
     * \param xc the curve.
     * \return true if the curve is a vertical spherical_arc; false otherwise.
     * \pre the arc is not degenerate (consists of a single point)
     */
    bool operator()(const X_monotone_curve_2 & xc) const
    {
      CGAL_precondition(!xc.is_degenerate());
      return ( CGAL::certainly(xc.supporting_circle().supporting_plane().c() == Kernel::FT(0)) );
    }
  };

  /*! Obtain an Is_vertical_2 function object */
  Is_vertical_2 is_vertical_2_object() const { return Is_vertical_2(); }

  /*! A functor that compares the y-coordinates of a directional point and
   * an x-monotone arc at the point x-coordinate
   */
  class Compare_y_at_x_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_at_x_2(const Traits * traits) : m_traits(traits) {}

    friend class Arr_arc_on_sphere_traits_2<Kernel>;
    
  public:
    /*! Return the location of the given point with respect to the input arc.
     * \param xc the curve.
     * \param p the point.
     * \return SMALLER - y(p) < xc(x(p)), i.e. the point is below the curve;
     *         EQUAL   - p lies on the curve.
     *         LARGER  - y(p) > xc(x(p)), i.e. the point is above the curve;
     * \pre p is not a singularity point.
     * \pre p is in the x-range of xc.
     */
    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & xc) const
    { 
      if (this->is_vertical_2_object()(xc)){
        Point_2* pt1=& (xc.source());
        Point_2* pt2=& (xc.target());
        CGAL::Compare_z_3<Spherical_kernel> cmp=m_traits->SK.compare_z_3_object()
        if (cmp(*pt1,*pt2)==CGAL::LARGER)
          std::swap(pt1,pt2);
        CGAL::Comparison_result res=cmp(p,*pt1);
        if(res!=CGAL::LARGER)
          return res;
        res=cmp(p,*pt2);
        if (res!=CGAL::SMALLER)
          return res;
        return CGAL::EQUAL;
      }
      return CGAL::opposite(m_traits->SK.compare_z_at_theta_3_object()(xc,p));
    }
  };

  /*! Obtain a Compare_y_at_x_2 function object */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(this); }

  /*! A functor that compares the y-coordinates of two x-monotone arcs
   * immediately to the left of their intersection directional point.
   */
  class Compare_y_at_x_left_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_at_x_left_2(const Traits * traits) : m_traits(traits) {}

    friend class Arr_arc_on_sphere_traits_2<Kernel>;
    
  public:
    /*! Compare the y value of two x-monotone curves immediately to the left
     * of their intersection point.
     * \param xc1 the first curve.
     * \param xc2 the second curve.
     * \param p the intersection point.
     * \return the relative position of xc1 with respect to xc2 immdiately to
     *         the left of p: SMALLER, EQUAL, or LARGER.
     * \pre the point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \pre the arcs are not degenerate
     */
    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 const X_monotone_curve_2 & xc2,
                                 const Point_2 & p) const
    {
      
      // If Both arcs are vertical, they overlap:
      Is_vertical_2 pred=m_traits->is_vertical_2_object();
      if ( pred(xc1) && pred(xc2) ) return EQUAL;
      if ( pred(xc1) ) return SMALLER;
      if ( pred(xc2) ) return LARGER;
      
      //TAG_TODO
      CGAL::Hcircle_type t1=????????????????;
      CGAL::Hcircle_type t2=????????????????;
      
      //if they overlap
      if ( t1==t2 && xc1.supporting_circle()==xc2.supporting_circle() ) return EQUAL;
      
      CGAL::Compare_z_3 cmp=m_traits->SK.compare_z_3_object();
      
      if (cmp(xc1.supporting_circle().extremal_point_z(),p)==CGAL::EQUAL){
        if (cmp(xc2.supporting_circle().extremal_point_z(),p)==CGAL::EQUAL){
          //p is the extremal point of both circles
          if (t1==t2)
            if (t1=CGAL::UPPER)
              return ( xc1.supporting_circle().squared_radius() < xc2.supporting_circle().squared_radius() )?CGAL::SMALLER:CGAL::LARGER;
            else
              return ( xc1.supporting_circle().squared_radius() > xc2.supporting_circle().squared_radius() )?CGAL::SMALLER:CGAL::LARGER;
          return (t1==CGAL::UPPER)?CGAL::LARGER:CGAL::SMALLER;
        }
        //~ return m_traits->SK.compare_z_at_theta_3_object()(p,xc2);
        return (t1==CGAL::UPPER)?CGAL::LARGER:CGAL::SMALLER;
      }
      if (cmp(xc2.supporting_circle().extremal_point_z(),p)==CGAL::EQUAL)
        return (t2==CGAL::LOWER)?CGAL::LARGER:CGAL::SMALLER;
        //~ return CGAL::opposite(m_traits->SK.compare_z_at_theta_3_object()(p,xc1));
      
      return m_traits->SK.compare_z_to_left_3_object()(xc1,xc2,p);
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 function object */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(this); }
  
  /*! A functor that compares the y-coordinates of two x-monotone arcs
   * immediately to the right of their intersection directional point.
   */
  class Compare_y_at_x_right_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_at_x_right_2(const Traits * traits) : m_traits(traits) {}

    friend class Arr_arc_on_sphere_traits_2<Kernel>;
    
  public:
    /*! Compare the y value of two x-monotone curves immediately to the right
     * of their intersection point.
     * \param xc1 the first curve.
     * \param xc2 the second curve.
     * \param p the intersection point.
     * \return the relative position of xc1 with respect to xc2 immdiately to
     *         the right of p: SMALLER, EQUAL, or LARGER.
     * \pre the point p lies on both curves, and both of them must also be
     *      defined to its right (lexicographically).
     * \pre the arcs are not degenerate
     */
    Comparison_result operator()(const X_monotone_curve_2 & xc1,
                                 const X_monotone_curve_2 & xc2,
                                 const Point_2 & p) const
    {
      
      // If Both arcs are vertical, they overlap:
      Is_vertical_2 pred=m_traits->is_vertical_2_object();
      if ( pred(xc1) && pred(xc2) ) return EQUAL;
      if ( pred(xc1) ) return LARGER;
      if ( pred(xc2) ) return SMALLER;
      
      //TAG_TODO
      CGAL::Hcircle_type t1=????????????????;
      CGAL::Hcircle_type t2=????????????????;
      
      //if they overlap
      if ( t1==t2 && xc1.supporting_circle()==xc2.supporting_circle() ) return EQUAL;      
      
      
      CGAL::Compare_z_3 cmp=m_traits->SK.compare_z_3_object();
      
      if (cmp(xc1.supporting_circle().extremal_point_z(),p)==CGAL::EQUAL){
        if (cmp(p,xc2.supporting_circle().extremal_point_z())==CGAL::EQUAL){
          //p is the extremal point of both circles
          if (t1==t2)
            if (t1==CGAL::UPPER)
              return ( xc1.supporting_circle().squared_radius() < xc2.supporting_circle().squared_radius() )?CGAL::SMALLER:CGAL::LARGER;
            else
              return ( xc1.supporting_circle().squared_radius() > xc2.supporting_circle().squared_radius() )?CGAL::SMALLER:CGAL::LARGER;
          return (t1==CGAL::UPPER)?CGAL::LARGER:CGAL::SMALLER;
        }
        return (t1==CGAL::UPPER)?CGAL::LARGER:CGAL::SMALLER;
        //~ return m_traits->SK.compare_z_at_theta_3_object()(p,xc2);
      }
      if (cmp(xc2.supporting_circle().extremal_point_z(),p)==CGAL::EQUAL)
        return (t2==CGAL::LOWER)?CGAL::LARGER:CGAL::SMALLER;
        //~ return CGAL::opposite(m_traits->SK.compare_z_at_theta_3_object()(p,xc1));
      
      typename Spherical_kernel::V_compare_theta_tgt cmpt(p);
      
      if (cmpt.sign_of_delta(xc1,xc2)==0)
        return m_traits->SK.compare_z_to_left_3_object()(xc1,xc2,p);
      return CGAL::opposite(m_traits->SK.compare_z_to_left_3_object()(xc1,xc2,p));
    }
  };

  /*! Obtain a Compare_y_at_x_right_2 function object */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(this); }

  /*! A functor that checks whether two directional points and two x-monotone
   * arcs are identical.
   */
  class Equal_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Equal_2(const Traits * traits) : m_traits(traits) {}

    friend class Arr_arc_on_sphere_traits_2<Kernel>;
    
  public:
    /*! Determines whether the two x-monotone curves are the same (have the
     * same graph).
     * \param xc1 the first curve.
     * \param xc2 the second curve.
     * \return true if the two curves are the same; false otherwise.
     */
    bool operator()(const X_monotone_curve_2 & xc1,
                    const X_monotone_curve_2 & xc2) const
    {
      return m_traits->SK.equal_3_object(xc1,xc2);
    }

    /*! Determines whether the two points are the same.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return true if the two point are the same; false otherwise.
     */
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      const Kernel * kernel = m_traits;
      return kernel->equal_3_object()(p1, p2);
    }
  };

  /*! Obtain an Equal_2 function object */
  Equal_2 equal_2_object() const { return Equal_2(this); }
  //@}

  /// \name Functor definitions to handle boundaries
  //@{

  /*! A function object that obtains the parameter space of a geometric
   * entity along the x-axis
   */
  class Parameter_space_in_x_2 {
  public:
    /*! Obtains the parameter space at the end of an arc along the x-axis .
     * Note that if the arc end coincides with a pole, then unless the arc
     * coincides with the identification arc, the arc end is considered to
     * be approaching the boundary, but not on the boundary.
     * If the arc coincides with the identification arc, it is assumed to
     * be smaller than any other object.
     * \param xcv the arc
     * \param ce the arc end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the arc xcv.
     *   ARR_LEFT_BOUNDARY  - the arc approaches the identification arc from
     *                        the right at the arc left end.
     *   ARR_INTERIOR       - the arc does not approache the identification arc.
     *   ARR_RIGHT_BOUNDARY - the arc approaches the identification arc from
     *                        the left at the arc right end.
     * \pre xcv does not coincide with the vertical identification arc.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xcv,
                                   Arr_curve_end ce) const
    {
      if (xcv.is_vertical()) {
        CGAL_precondition(!xcv.is_on_boundary());
        return ARR_INTERIOR;
      }

      Point_2 p = (ce == ARR_MIN_END) ?
        (*this(m_traits.construct_min_vertex_2_object(xcv))
        :*this(m_traits.construct_max_vertex_2_object(xcv)));
      return (*this)(p);
    }

    /*! Obtains the parameter space at a point along the x-axis.
     * \param p the point.
     * \return the parameter space at p.
     * \pre p does not lie on the vertical identification curve.
     */
    Arr_parameter_space operator()(const Point_2& p) const
    {
      if (p.get_hq()!=0.5)
        return ARR_LEFT_BOUNDARY;
      if (p.get_hq()!=8.5 )
        return ARR_RIGHT_BOUNDARY;
      return ARR_INTERIOR;
    }
  };

  /*! Obtain a Parameter_space_in_x_2 function object */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { return Parameter_space_in_x_2(); }
  
  /*! A function object that obtains the parameter space of a geometric
   * entity along the y-axis
   */
  class Parameter_space_in_y_2 {
  public:
    /*! Obtains the parameter space at the end of an arc along the y-axis .
     * Note that if the arc end coincides with a pole, then unless the arc
     * coincides with the identification arc, the arc end is considered to
     * be approaching the boundary, but not on the boundary.
     * If the arc coincides with the identification arc, it is assumed to
     * be smaller than any other object.
     * \param xcv the arc
     * \param ce the arc end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the arc xcv.
     *   ARR_BOTTOM_BOUNDARY  - the arc approaches the south pole at the arc
     *                          left end.
     *   ARR_INTERIOR         - the arc does not approache a contraction point.
     *   ARR_TOP_BOUNDARY     - the arc approaches the north pole at the arc
     *                          right end.
     * There are no horizontal identification arcs!
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xcv,
                                   Arr_curve_end ce) const
    {
      Point_2 p = (ce == ARR_MIN_END) ?
        (*this(m_traits.construct_min_vertex_2_object(xcv))
        :*this(m_traits.construct_max_vertex_2_object(xcv)));
      return (*this)(p);
    }

    /*! Obtains the parameter space at a point along the y-axis.
     * \param p the point.
     * \return the parameter space at p.
     * \pre p does not lie on the horizontal identification curve.
     * There are no horizontal identification arcs!
     */
    Arr_parameter_space operator()(const Point_2& p) const
    {
      if (CGAL::certainly(CGAL::square(p.z()-p.reference_sphere().center().z())-p.reference_sphere().squared_radius()==0)){
        //TAG_TODO
        //#warning cannot use this since z coord of polar and bipolar are designed for the sorting
        return (CGAL::certainly(CGAL::Sign(p.z()-p.reference_sphere().center().z())> 0))?ARR_TOP_BOUNDARY:ARR_BOTTOM_BOUNDARY;
      }
      return ARR_INTERIOR;
    }
  };

  /*! Obtain a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(); }

  /*! A functor that compares the x-coordinates of arc ends near the
   * boundary of the parameter space.
   */
  class Compare_x_near_boundary_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_x_near_boundary_2(const Traits * traits) : m_traits(traits) {}

    friend class Arr_arc_on_sphere_traits_2<Kernel>;
    
  public:
    /*! Compare the x-coordinate of a direction with the x-coordinate of an
     * arc end near the boundary.
     * \param p the point direction.
     * \param xcv the arc, the endpoint of which is compared.
     * \param ce the arc-end indicator -
     *            ARR_MIN_END - the minimal end of xc or
     *            ARR_MAX_END - the maximal end of xc.
     * \return the comparison result:
     *         SMALLER - x(p) < x(xc, ce);
     *         EQUAL   - x(p) = x(xc, ce);
     *         LARGER  - x(p) > x(xc, ce).     
     * \pre p lies in the interior of the parameter space.
     * \pre the ce end of the arc xcv lies on a boundary.
     * \pre xcv does not coincide with the vertical identification curve.
     */
    Comparison_result operator()(const Point_2 & point,
                                 const X_monotone_curve_2 & xcv,
                                 Arr_curve_end ce) const
    {
      //TAG_TODO
      //#here add the case of polar and bipolar circles
      const Point_2 & p2 = (ce == ARR_MIN_END) ? m_traits->SK.construct_min_vertex_2_object(xcv) : m_traits->SK.construct_max_vertex_2_object(xcv);
      return m_traits->SK.compare_x_2_object(point,p2);
    }

    /*! Compare the x-coordinates of 2 arc ends near the boundary of the
     * parameter space.
     * \param xcv1 the first arc.
     * \param ce1 the first arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv1 or
     *            ARR_MAX_END - the maximal end of xcv1.
     * \param xcv2 the second arc.
     * \param ce2 the second arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv2 or
     *            ARR_MAX_END - the maximal end of xcv2.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce1) < x(xcv2, ce2);
     *         EQUAL   - x(xcv1, ce1) = x(xcv2, ce2);
     *         LARGER  - x(xcv1, ce1) > x(xcv2, ce2).
     * \pre the ce1 end of the arc xcv1 lies on a boundary.
     * \pre the ce2 end of the arc xcv2 lies on a boundary.
     * \pre xcv1 does not coincide with the vertical identification curve.
     * \pre xcv2 does not coincide with the vertical identification curve.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce2) const
    {
      //TAG_TODO
      //~ #here add the case of polar and bipolar circles
      const Point_2 & p1 = (ce1 == ARR_MIN_END) ? m_traits->SK.construct_min_vertex_2_object(xcv1) : m_traits->SK.construct_max_vertex_2_object(xcv1);
      const Point_2 & p2 = (ce2 == ARR_MIN_END) ? m_traits->SK.construct_min_vertex_2_object(xcv2) : m_traits->SK.construct_max_vertex_2_object(xcv2);
      return m_traits->SK.compare_x_2_object(p1,p2);
    }
  };

  /*! Obtain a Compare_x_near_boundary_2 function object */
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const
  { return Compare_x_near_boundary_2(this); }
    

  /*! A functor that compares the y-coordinates of arc ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_near_boundary_2(const Traits * traits) : m_traits(traits) {}

    friend class Arr_arc_on_sphere_traits_2<Kernel>;
    
  public:
    /*! Compare the y-coordinates of 2 curves at their ends near the boundary
     * of the parameter space.
     * \param xcv1 the first arc.
     * \param xcv2 the second arc.
     * \param ce the arc end indicator.
     * \return the second comparison result.
     * \pre the ce ends of the arcs xcv1 and xcv2 lie either on the left
     * boundary or on the right boundary of the parameter space.
     * There is no horizontal identification curve!
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce) const
    {
      const Point_2& p1=(ce==ARR_MIN_END)?m_traits->SK.construct_min_vertex_2_object(xcv1):m_traits->SK.construct_max_vertex_2_object(xcv1);
      const Point_2& p2=(ce==ARR_MIN_END)?m_traits->SK.construct_min_vertex_2_object(xcv2):m_traits->SK.construct_max_vertex_2_object(xcv2);
      return m_traits->SK.compare_z_3_object()(p1,p2);
    }
  };

  /*! Obtain a Compare_y_near_boundary_2 function object */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(this); }
  
  /*! A functor that indicates whether a geometric object lies on the
   * horizontal identification arc.
   */
  class Is_on_x_identification_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;
    
  public:
    /*! Determine whether a point lies on the horizontal identification arc.
     * \param p the point.
     * \return a Boolean indicating whether p lies on the vertical
     * identification arc.
     */
    bool operator()(const Point_2 & p) const
    {
      return false;
    }

    /*! Determine whether an arc coincides with the horizontal identification
     * arc.
     * \param xcv the arc.
     * \return a Boolean indicating whether xcv coincides with the vertical
     * identification arc.
     */
    bool operator()(const X_monotone_curve_2 & xcv) const
    {
      return false;
    }
  };
  
  /*! Obtain a Is_on_x_identification_2 function object */
  Is_on_x_identification_2 is_on_x_identification_2_object() const
  { return Is_on_x_identification_2(); }
  
  /*! A functor that indicates whether a geometric object lies on the
   * vertical identification arc.
   */
  class Is_on_y_identification_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;
    
  public:
    /*! Determine whether a point lies on the vertical identification arc.
     * \param p the point.
     * \return a Boolean indicating whether p lies on the vertical
     * identification arc.
     */
    bool operator()(const Point_2 & p) const
    {
      return ( p.get_hq()==0.5 || p.get_hq()==8.5 );
    }

    /*! Determine whether an arc coincides with the vertical identification
     * arc.
     * \param xcv the arc.
     * \return a Boolean indicating whether xcv coincides with the vertical
     * identification arc.
     */
    bool operator()(const X_monotone_curve_2 & xcv) const
    {
      //TAG_TODO      #problem of bipolar circle!!!!!!!!!
      // If the curve is not vertical, it cannot coincide with the ident. arc:
      if (!m_traits->is_vertical_2_object(xcv)) return false;

      if ( *this(xcv.source()) && *this(xcv.target()) )
        return true;
      return false;
    }
  };
  
  /*! Obtain a Is_on_y_identification_2 function object */
  Is_on_y_identification_2 is_on_y_identification_2_object() const
  { return Is_on_y_identification_2(); }
  
  /*! A functor that compares the x-coordinate of two given points
   * that lie on the horizontal identification arc.
   * The parameter space does not contain a horizontal identification arc.
   */
  class Compare_x_on_identification_2 {
  public:
    /*! Compare the x-coordinate of two given points that lie on the
     * horizontal identification arc.
     * \param p1 the first point.
     * \param p2 the second point.
     * There is no horizontal identification arc!
     */
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      CGAL_error_msg("There is no horizontal identification arc!");
      return SMALLER;
    }
  };

  /*! Obtain a Compare_x_on_identification_2 function object */
  Compare_x_on_identification_2 compare_x_on_identification_2_object() const
  { return Compare_x_on_identification_2(); }

  /*! A functor that compares the y-coordinate of two given points
   * that lie on the vertical identification arc.
   */
  class Compare_y_on_identification_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_on_identification_2(const Traits * traits) : m_traits(traits) {}

    friend class Arr_arc_on_sphere_traits_2<Kernel>;
    
  public:
    /*! Compare the y-coordinate of two given points that lie on the vertical
     * identification curve.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return SMALLER - p1 is lexicographically smaller than p2;
     *         EQUAL   - p1 and p2 coincides;
     *         LARGER  - p1 is lexicographically larger than p2;
     * \pre p1 lies on the vertical identification arc.
     * \pre p2 lies on the vertical identification arc.
     */
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      return m_traits->compare_y(p1, p2);
    }
  };

  /*! Obtain a Compare_y_on_identification_2 function object */
  Compare_y_on_identification_2 compare_y_on_identification_2_object() const
  { return Compare_y_on_identification_2(this); }
  //@}

  /// \name Functor definitions for supporting intersections.
  //@{

  /*! A functor that divides an arc into x-monotone arcs. That are, arcs that
   * do not cross the identification arc.
   */
  class Make_x_monotone_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Make_x_monotone_2(const Traits * traits) : m_traits(traits) {}

    friend class Arr_arc_on_sphere_traits_2<Kernel>;
    
  public:    
    /*! Cut the given curve into x-monotone subcurves and insert them into the
     * given output iterator. As spherical_arcs are always x_monotone, only one
     * object will be contained in the iterator.
     * \param xc the curve.
     * \param oi the output iterator, whose value-type is Object. The output
     *           object is a wrapper of either an X_monotone_curve_2, or - in
     *           case the input spherical_arc is degenerate - a Point_2 object.
     * \return the past-the-end iterator.
     */
    template<typename OutputIterator>
    OutputIterator operator()(const Curve_2 & c, OutputIterator oi) const
    {
      if (c.is_degenerate()) {
        // The spherical_arc is a degenerate point - wrap it with an object:
        *oi++ = make_object(c.source());
        return oi;
      }
      std::list<CGAL::Object> L;
      m_traits->SK.make_theta_monotonic_3_object()(c,std::back_inserter(L));
      if (L.size==0)
        return oi;
      typename Spherical_kernel::Circular_arc_on_reference_sphere_3 arc;
      if (CGAL::assign(arc,*(L.first())) )
        if (!m_traits->is_on_y_identification_2_object()(arc)){
          //TAG_TODO #arc does not cross else cut it
          *oi++=*(L.first());
        }
      if (L.size()>1 && CGAL::assign(arc,*(++L.first())))
        if (!m_traits->is_on_y_identification_2_object()(arc))
          *oi++=*(++L.first());
      if (L.size()>2 && CGAL::assign(arc,*(++(++L.first()))))
        if (!m_traits->is_on_y_identification_2_object()(arc))
          *oi++=*(++(++L.first()));        
      return oi;
    }
  };

  /*! Obtain a Make_x_monotone_2 function object */
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(this); }

  /*! A functor that splits an x-monotone arc at a directional point. */
  class Split_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Split_2(const Traits * traits) : m_traits(traits) {}

    friend class Arr_arc_on_sphere_traits_2<Kernel>;
    
  public:
    /*! Split a given x-monotone curve at a given point into two sub-curves.
     * \param xc the curve to split
     * \param p the split point.
     * \param xc1 (output) the left resulting subcurve. p is its right
     * endpoint.
     * \param xc2 (output) the right resulting subcurve. p is its left
     * endpoint.
     * \pre p lies on xc but is not one of its endpoints.
     * \pre xc is not degenerate
     */
    void operator()(const X_monotone_curve_2 & xc, const Point_2 & p,
                    X_monotone_curve_2 & xc1, X_monotone_curve_2 & xc2) const
    {
      const Point_2 & source = xc.source();
      const Point_2 & target = xc.target();
      CGAL_precondition_code(const Kernel * kernel = m_traits);
      CGAL_precondition_code
        (typename Kernel::Equal_3 equal_3 = kernel->equal_3_object());
      CGAL_precondition(!equal_3(p, source));
      CGAL_precondition(!equal_3(p, target));
      xc1=X_monotone_curve_2(xc.supporting_circle(),xc.source(),p);
      xc2=X_monotone_curve_2(xc.supporting_circle(),p,xc.target());
      return;
    }
  };

  /*! Obtain a Split_2 function object */
  Split_2 split_2_object() const { return Split_2(this); }

  /*! A functor that computes intersections between x-monotone arcs. */
  class Intersect_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Intersect_2(const Traits * traits) : m_traits(traits) {}

    friend class Arr_arc_on_sphere_traits_2<Kernel>;
    
  public:
    /*! Find the intersections of the two given curves and insert them into the
     * given output iterator. As two spherical_arcs may itersect only once,
     * only a single intersection will be contained in the iterator.
     * \param xc1 the first curve.
     * \param xc2 the second curve.
     * \param oi the output iterator.
     * \return the past-the-end iterator.
     * \pre xc1 and xc2 are not degenerate
     */
    template<typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2 & xc1,
                              const X_monotone_curve_2 & xc2,
                              OutputIterator oi) const
    {
      CGAL_precondition(!xc1.is_degenerate());
      CGAL_precondition(!xc2.is_degenerate());

      std::list<CGAL::Object> L;
      m_traits->SK.intersect_3_object()(static_cast<typename Spherical_kernel::Circular_arc_3>(xc1),
                                        static_cast<typename Spherical_kernel::Circular_arc_3>(xc2),
                                        std::back_inserter(L));
      
      if (L.size()==0)
        return oi;
      
      std::pair<typename Spherical_kernel::Circular_arc_point_3,unsigned> P;
      if (CGAL::assign(P,*(L.first()) )
          && !m_traits->is_on_y_identification_2_object()(P.first) 
          && CGAL::certainly(CGAL::square(P.first.z()-xc1.reference_sphere().center().z())-xc1.reference_sphere().squared_radius()==0) ){
        *oi++=CGAL::make_object(std::make_pair(typename Spherical_kernel::Circular_arc_point_on_reference_sphere_3(CGAL::hquadrant(P),P.first),P.second()));
        if (L.size()==1) return oi;
      }
      if (L.size()==2)
        if (CGAL::assign(P,*(++L.first()) )
            && !m_traits->is_on_y_identification_2_object()(P.first) 
            && CGAL::certainly(CGAL::square(P.first.z()-xc1.reference_sphere().center().z())-xc1.reference_sphere().squared_radius()==0) ){
          *oi++=CGAL::make_object(std::make_pair(typename Spherical_kernel::Circular_arc_point_on_reference_sphere_3(CGAL::hquadrant(P),P.first),P.second()));
          return oi;
        }
      typename Spherical_kernel::Circular_arc_3 arc;
      if (CGAL::assign(arc,L.first())){
        //TAG_TODO
        //~ #*oi++=CGAL::make_object( X_monotone_curve_2(,,) );
      }
      Point_2 pt;
      if (CGAL::assign(pt,*(L.first())) 
        && !m_traits->is_on_y_identification_2_object()(pt.first) 
        && CGAL::certainly(CGAL::square(pt.first.z()-xc1.reference_sphere().center().z())-xc1.reference_sphere().squared_radius()==0))
          *oi++=CGAL::make_object(pt);
     return oi; 
    }
  };

  /*! Obtain an Intersect_2 function object */
  Intersect_2 intersect_2_object() const { return Intersect_2(this); }

  /*! A functor that tests whether two x-monotone arcs can be merged. */
  class Are_mergeable_2 {
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Are_mergeable_2(const Traits * traits) : m_traits(traits) {}

    friend class Arr_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Check whether it is possible to merge two given x-monotone curves.
     * \param xc1 the first curve.
     * \param xc2 the second curve.
     * \return true if the two arcs are mergeable; false otherwise.
     * Two arcs are mergeable if:
     * 1. they are supported by the same plane, and
     * 2. share a common endpoint that is not on the identification arc
     */
    bool operator()(const X_monotone_curve_2 & xc1,
                    const X_monotone_curve_2 & xc2) const
    {
      if (xc1.supporting_circle()!=xc2.supporting_circle())
        return false;
      
      if (xc1.is_empty() || xc2.is_empty()) return true;
      if (xc1.is_full() && xc2.is_full()) return false;
        
      if (CGAL::Sign(
        xc1.supporting_circle().supporting_plane().orthogonal_vector()*
        xc2.supporting_circle().supporting_plane().orthogonal_vector()) > 0
      ){
        if ( m_traits->SK.equal_3_object()(xc1.source(),xc2.target()) 
          || m_traits->SK.equal_3_object()(xc2.source(),xc1.target()) 
        )
          return true;
        return false;
      }
      
      if ( m_traits->SK.equal_3_object()(xc1.target(),xc2.target()) 
        || m_traits->SK.equal_3_object()(xc2.source(),xc1.source()))
          return true;
        return false;
    }
  };

  /*! Obtain an Are_mergeable_2 function object */
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(this); }

  /*! A functor that merges two x-monotone arcs into one */
  class Merge_2 {
  protected:
    typedef Arr_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits * m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Traits * traits) : m_traits(traits) {}

    friend class Arr_arc_on_sphere_traits_2<Kernel>;
    
  public:
    /*! Merge two given x-monotone curves into a single curve (spherical_arc).
     * \param xc1 the first curve.
     * \param xc2 the second curve.
     * \param xc Output: the merged curve.
     * \pre the two curves are mergeable. That is, they are supported by the
     *      same plane or oposite planes and share a common endpoint that is
     *      not on the discontinuity arc.
     */
    void operator()(const X_monotone_curve_2 & xc1,
                    const X_monotone_curve_2 & xc2,
                    X_monotone_curve_2 & xc) const
    {
      if (xc1.is_degenerate() || xc1.is_empty()) {
        xc = xc2;
        return;
      }

      if (xc2.is_degenerate() || xc2.is_empty()) {
        xc = xc1;
        return;
      }
        
      if (CGAL::Sign(
        xc1.supporting_circle().supporting_plane().orthogonal_vector()*
        xc2.supporting_circle().supporting_plane().orthogonal_vector()) > 0
      ){
        if ( m_traits->SK.equal_3_object()(xc1.source(),xc2.target()) ){
          xc=X_monotone_curve_2(xc1.supporting_circle(),xc2.source,xc1.target());
          return;
        }
        if(m_traits->SK.equal_3_object()(xc2.source(),xc1.target()) ){
          xc=X_monotone_curve_2(xc1.supporting_circle(),xc1.source,xc2.target());
          return;          
        }
      }
      
      if ( m_traits->SK.equal_3_object()(xc1.target(),xc2.target()) ){
        xc=X_monotone_curve_2(xc1.supporting_circle(),xc1.source,xc2.source());
        return;        
      }
      
      if( m_traits->SK.equal_3_object()(xc2.source(),xc1.source())){
        xc=X_monotone_curve_2(xc1.supporting_circle(),xc2.target(),xc1.target());
        return;
      }
    }
  };

  /*! Obtain a Merge_2 function object */
  Merge_2 merge_2_object() const { return Merge_2(this); }
  //@}

  /// \name Functor definitions for the landmarks point-location strategy.
  //@{
  typedef double                          Approximate_number_type;

  class Approximate_2 {
  public:

    /*! Return an approximation of a point coordinate.
     * \param p the exact point.
     * \param i the coordinate index (either 0 or 1).
     * \pre i is either 0 or 1.
     * \return an approximation of p's x-coordinate (if i == 0), or an 
     *         approximation of p's y-coordinate (if i == 1).
     */
    Approximate_number_type operator()(const Point_2 & p, int i) const
    {
      CGAL_precondition(i == 0 || i == 1);
      return (i == 0) ? p.get_theta_approx() : CGAL::to_double(p.z());
    }
  };

  /*! Obtain an Approximate_2 function object */
  Approximate_2 approximate_2_object() const { return Approximate_2(); }

  /// \name Functor definitions for the Boolean set-operation traits.
  //@{
  class Compare_endpoints_xy_2 {
  public:

    /*!
     * Compare the endpoints of an $x$-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param xc the curve.
     * \return SMALLER if the curve is directed right;
     *         LARGER if the curve is directed left.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xc)
    { return m_traits->SK.compare_theta_z_3_object(xc.source(),xc.target); }
  };

  /*! Obtain a Compare_endpoints_xy_2 function object */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(); }

  class Construct_opposite_2 {
  public:
    /*! Construct an opposite x-monotone (with swapped source and target).
     * \param xc the curve.
     * \return the opposite curve.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2 & xc)
    { return X_monotone_curve_2(c.supporting_circle(),xc.target(),xc.source()); }
  };

  /*! Obtain a Construct_opposite_2 function object */
  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(); }

};

CGAL_END_NAMESPACE

#endif
