// Copyright (c) 2006,2007,2009,2010,2011,2013 Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein        <wein@post.tau.ac.il>
//                 Efi Fogel       <efif@post.tau.ac.il>

#ifndef CGAL_ARR_BATCHED_POINT_LOCATION_TRAITS_2_H
#define CGAL_ARR_BATCHED_POINT_LOCATION_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*!
 * Definition of the Arr_batched_point_location_traits_2<Arrangement> class.
 */

#include <CGAL/Arr_tags.h>

namespace CGAL {

/*! \class
 * A traits-class decorator for the use of the batched point-location process.
 */
template <typename Arrangement_>
class Arr_batched_point_location_traits_2
{
public:

  typedef Arrangement_                                  Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2     Base_traits_2;

  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;

  typedef typename Base_traits_2::X_monotone_curve_2   Base_x_monotone_curve_2;
  typedef typename Base_traits_2::Point_2              Base_point_2;

  typedef typename Base_traits_2::Construct_min_vertex_2
                                                    Base_construct_min_vertex_2;
  typedef typename Base_traits_2::Construct_max_vertex_2
                                                    Base_construct_max_vertex_2;
  typedef typename Base_traits_2::Compare_x_2       Base_compare_x_2;
  typedef typename Base_traits_2::Compare_xy_2      Base_compare_xy_2;
  typedef typename Base_traits_2::Compare_y_at_x_2  Base_compare_y_at_x_2;
  typedef typename Base_traits_2::Compare_y_at_x_right_2
                                                    Base_compare_y_at_x_right_2;
  typedef typename Base_traits_2::Equal_2           Base_equal_2;
  typedef typename Base_traits_2::Is_vertical_2     Base_is_vertical_2;

  typedef typename Base_traits_2::Has_do_intersect_category
                                                    Has_do_intersect_category;

  typedef typename internal::Arr_complete_left_side_category< Base_traits_2 >::Category
                                                    Left_side_category;
  typedef typename internal::Arr_complete_bottom_side_category< Base_traits_2 >::Category
                                                    Bottom_side_category;
  typedef typename internal::Arr_complete_top_side_category< Base_traits_2 >::Category
                                                   Top_side_category;
  typedef typename internal::Arr_complete_right_side_category< Base_traits_2 >::Category
                                                    Right_side_category;

  /* Overlay is implemented as sweep-line visitor. The sweep-line algorithm
   * never uses Compare_y_at_x_left_2, and it never performs merging of curves.
   * Thus, AreMergeable_2 and Merge_2 are not needed either.
   */
  typedef Tag_false                                 Has_left_category;
  typedef Tag_false                                 Has_merge_category;

protected:

  const Base_traits_2*    m_base_traits;

public:

  /*! Constructor. */
  Arr_batched_point_location_traits_2 (const Base_traits_2& tr) :
    m_base_traits (&tr)
  {}

  /*! \class
   * Nested extension of the x-monotone curve type.
   */
  class Ex_x_monotone_curve_2
  {
  public:

    typedef Base_x_monotone_curve_2    Base;

  protected:

    Base_x_monotone_curve_2 m_base_xcv;  // The base x-monotone curve.
    Halfedge_const_handle   m_he;        // The corresponding arrangement edge.

  public:

    Ex_x_monotone_curve_2 ():
      m_base_xcv(),
      m_he()
    {}

    Ex_x_monotone_curve_2 (const Base& xcv):
      m_base_xcv(xcv),
      m_he()
    {}

    Ex_x_monotone_curve_2 (const Base& xcv, Halfedge_const_handle he) :
      m_base_xcv(xcv),
      m_he(he)
    {
      CGAL_precondition (he->direction() == ARR_RIGHT_TO_LEFT);
    }

    Halfedge_const_handle halfedge_handle() const
    {
      return (m_he);
    }

    const Base& base () const
    {
      return (m_base_xcv);
    }

    Base& base ()
    {
      return (m_base_xcv);
    }

    operator const Base&() const
    {
      return (m_base_xcv);
    }

    operator Base&()
    {
      return (m_base_xcv);
    }

  };

  /*! \class
   * Nested extension of the point type.
   */
  class Ex_point_2
  {
  public:

    typedef  Base_point_2            Base;

  protected:

    Base                   m_base_pt; // The base point.
    Vertex_const_handle    m_v;       // The corresponding arrangement vertex.

  public:

    Ex_point_2 ():
      m_base_pt(),
      m_v()
    {}

    Ex_point_2 (const Base& pt):
      m_base_pt (pt),
      m_v()
    {}

    Ex_point_2 (const Base& pt, Vertex_const_handle v):
      m_base_pt (pt),
      m_v (v)
    {}

    const Base& base() const
    {
      return (m_base_pt);
    }

    Base& base()
    {
      return (m_base_pt);
    }

    operator const Base&() const
    {
      return (m_base_pt);
    }

    operator Base&()
    {
      return (m_base_pt);
    }

    Vertex_const_handle vertex_handle() const
    {
      return m_v;
    }
  };

  typedef Ex_x_monotone_curve_2                     X_monotone_curve_2;
  typedef Ex_point_2                                Point_2;

  // For debugging purposes:
  friend std::ostream& operator<< (std::ostream& os,
                                   const X_monotone_curve_2& xcv)
  {
    os << xcv.base();
    return (os);
  }

  // For debugging purposes:
  friend std::ostream& operator<< (std::ostream& os,
                                   const Point_2& pt)
  {
    os << pt.base();
    return (os);
  }

  /*! A functor that obtains the left endpoint of an x-monotone curve. */
  class Construct_min_vertex_2 {
  protected:
    //! The base operator.
    Base_construct_min_vertex_2 m_base_min_v;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Construct_min_vertex_2 (const Base_construct_min_vertex_2& base):
        m_base_min_v(base)
    {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param xcv The curve.
     * \return The left endpoint.
     */
    Point_2 operator() (const X_monotone_curve_2 & xcv)
    {
      // Note that the halfedge associated with the curve is always directed
      // from right to left, so its target is the leftmost vertex.
      Vertex_const_handle vh = xcv.halfedge_handle()->target();
      return (Point_2 (m_base_min_v (xcv.base()), vh));
    }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return
      Construct_min_vertex_2 (m_base_traits->construct_min_vertex_2_object());
  }

  /*! A functor that obtains the right endpoint of an x-monotone curve. */
  class Construct_max_vertex_2 {
  protected:
    //! The base operator.
    Base_construct_max_vertex_2 m_base_max_v;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Construct_max_vertex_2 (const Base_construct_max_vertex_2& base):
      m_base_max_v(base)
    {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    /*!
     * Get the right endpoint of the x-monotone curve .
     * \param xcv The curve.
     * \return The right endpoint.
     */
    Point_2 operator() (const X_monotone_curve_2 & xcv)
    {
      // Note that the halfedge associated with the curve is always directed
      // from right to left, so its source is the rightmost vertex.
      Vertex_const_handle vh = xcv.halfedge_handle()->source();
      return (Point_2 (m_base_max_v (xcv.base()), vh));
    }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return
      Construct_max_vertex_2 (m_base_traits->construct_max_vertex_2_object());
  }

  /*! A functor that compares two points lexigoraphically: by x, then by y. */
  class Compare_xy_2 {
  protected:
    //! The base operator.
    Base_compare_xy_2    m_base_cmp_xy;

    Vertex_const_handle  invalid_v;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_xy_2(const Base_compare_xy_2& base) :
      m_base_cmp_xy(base),
      invalid_v()
    {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param xcv The curve.
     * \return The left endpoint.
     */
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      if (p1.vertex_handle() == p2.vertex_handle() &&
          p1.vertex_handle() != invalid_v)
        return (EQUAL);

      return (m_base_cmp_xy (p1.base(), p2.base()));
    }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const
  {
    return Compare_xy_2(m_base_traits->compare_xy_2_object());
  }

  /*! A functor that compares the y-coordinates of a point and an
   * x-monotone curve at the point x-coordinate.
   */
  class Compare_y_at_x_2 {
  protected:
    //! The base operator.
    Base_compare_y_at_x_2 m_base_cmp_y_at_x;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_at_x_2(const Base_compare_y_at_x_2& base) :
      m_base_cmp_y_at_x(base)
    {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    Comparison_result operator() (const Point_2& p,
                                  const X_monotone_curve_2& xcv) const
    {
      return (m_base_cmp_y_at_x (p.base(), xcv.base()));
    }
  };

  /*! Obtain a Compare_y_at_x_2 function object. */
  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return (Compare_y_at_x_2 (m_base_traits->compare_y_at_x_2_object()));
  }

  /*! A functor that compares compares the y-coordinates of two x-monotone
   * curves immediately to the right of their intersection point.
   */
  class Compare_y_at_x_right_2 {
  protected:
    //! The base operator.
    Base_compare_y_at_x_right_2 m_base_cmp_y_at_x_right;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_at_x_right_2(const Base_compare_y_at_x_right_2& base) :
      m_base_cmp_y_at_x_right(base)
    {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    Comparison_result operator() (const X_monotone_curve_2& xcv1,
                                  const X_monotone_curve_2& xcv2,
                                  const Point_2& p) const
    {
      return (m_base_cmp_y_at_x_right(xcv1.base(), xcv2.base(), p.base()));
    }
  };

  /*! Obtain a Compare_y_at_x_right_2 function object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object () const
  {
    return (Compare_y_at_x_right_2
	    (m_base_traits->compare_y_at_x_right_2_object()));
  }

  /*! A functor that checks whether two points and two x-monotone curves are
   * identical.
   */
  class Equal_2 {
  protected:
    //! The base operator.
    Base_equal_2           m_base_eq;

    Vertex_const_handle    invalid_v;
    Halfedge_const_handle  invalid_he;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Equal_2(const Base_equal_2& base) :
      m_base_eq(base),
      invalid_v(),
      invalid_he()
    {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    /*! Check if two curves are the same. */
    bool operator() (const X_monotone_curve_2& xcv1,
		     const X_monotone_curve_2& xcv2) const
    {
      if (xcv1.halfedge_handle() == xcv2.halfedge_handle() &&
          xcv1.halfedge_handle() != invalid_he)
        return (true);

      return (m_base_eq(xcv1.base(), xcv2.base()));
    }

    /*! Check if the two points are the same. */
    bool operator() (const Point_2& p1, const Point_2& p2) const
    {
      if (p1.vertex_handle() == p2.vertex_handle() &&
          p1.vertex_handle() != invalid_v)
        return (true);

      return (m_base_eq(p1.base(), p2.base()));
    }
  };

  /*! Obtain a Equal_2 function object. */
  Equal_2 equal_2_object () const
  {
    return (Equal_2 (m_base_traits->equal_2_object()));
  }

  /*! A functor that compares the x-coordinates of two points */
  class Compare_x_2 {
  protected:
    //! The base operator.
    Base_compare_x_2 m_base_cmp_x;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_2(const Base_compare_x_2& base) : m_base_cmp_x(base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      return (m_base_cmp_x(p1.base(), p2.base()));
    }
  };

  /*! Obtain a Compare_x_2 function object. */
  Compare_x_2 compare_x_2_object () const
  {
    return (Compare_x_2 (m_base_traits->compare_x_2_object()));
  }

  /*! A functor that checks whether a given x-monotone curve is vertical. */
  class Is_vertical_2 {
  protected:
    //! The base operator.
    Base_is_vertical_2 m_base_is_vert;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Is_vertical_2(const Base_is_vertical_2& base) : m_base_is_vert(base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    bool operator() (const X_monotone_curve_2& xcv) const
    {
      return (m_base_is_vert(xcv.base()));
    }
  };

  /*! Obtain a Is_vertical_2 function object. */
  Is_vertical_2 is_vertical_2_object() const
  {
    return (Is_vertical_2(m_base_traits->is_vertical_2_object()));
  }


  // left-right

  /*! A functor that determines whether an endpoint of an x-monotone curve lies
   * on a boundary of the parameter space along the x axis.
   */
  class Parameter_space_in_x_2 {
  protected:
    //! The base traits.
    const Base_traits_2      *m_base;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Parameter_space_in_x_2 (const Base_traits_2 *tr) : m_base (tr) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    Arr_parameter_space operator() (const X_monotone_curve_2& xcv,
                                    Arr_curve_end ce) const
    {
      return m_base->parameter_space_in_x_2_object() (xcv.base(), ce);
    }

    Arr_parameter_space operator() (const Point_2 & p) const
    {
      return m_base->parameter_space_in_x_2_object() (p.base());
    }

    Arr_parameter_space operator() (const X_monotone_curve_2 & xcv) const
    {
      return m_base->parameter_space_in_x_2_object() (xcv.base());
    }
  };

  /*! Obtain a Parameter_space_in_x_2 function object */
  Parameter_space_in_x_2 parameter_space_in_x_2_object () const
  {
    return Parameter_space_in_x_2 (m_base_traits);
  }


  /*! A function object that determines whether an x-monotone curve or a
   * point coincide with the vertical identification curve.
   */
  class Is_on_x_identification_2 {
  protected:
    //! The base traits.
    const Base_traits_2      *m_base;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Is_on_x_identification_2 (const Base_traits_2* tr) : m_base (tr) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    bool operator() (const Point_2 & p) const
    {
      return m_base->is_on_x_identification_2_object() (p.base());
    }

    bool operator() (const X_monotone_curve_2 & xcv) const
    {
      return m_base->is_on_x_identification_2_object() (xcv.base());
    }
  };

  /*! Obtain a Is_on_x_identification_2 function object */
  Is_on_x_identification_2 is_on_x_identification_2_object () const
  {
    return Is_on_x_identification_2 (m_base_traits);
  }

  /*! A functor that compares the y-coordinate of two given points
   * that lie on the vertical identification curve.
   */
  class Compare_y_on_boundary_2 {
  protected:
    //! The base traits.
    const Base_traits_2 * m_base;

    /*! Constructor.
     * \param tr The base traits class. It must be passed, to handle
     *           non stateless traits (e.g., it stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_on_boundary_2(const Base_traits_2 * tr) : m_base(tr) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      return m_base->compare_y_on_boundary_2_object()(p1.base(), p2.base());
    }
    Comparison_result operator() (const Point_2 & pt,
                                  const X_monotone_curve_2& xcv, Arr_curve_end ce) const
    {
      return m_base->compare_y_on_boundary_2_object()(pt.base(), xcv.base(), ce);
    }
    Comparison_result operator() (const X_monotone_curve_2& xcv1, Arr_curve_end ce1,
                                  const X_monotone_curve_2& xcv2, Arr_curve_end ce2) const
    {
      return m_base->compare_y_on_boundary_2_object()(xcv1.base(), ce1, xcv2.base(), ce2);
    }
  };

  /*! Obtain a Compare_y_on_boundary_2 functor object. */
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object () const
  {
    return Compare_y_on_boundary_2(m_base_traits);
  }

  /*! A function object that compares the y-coordinates of curve ends near the
   * boundary of the parameter space
   */
  class Compare_y_near_boundary_2 {
  protected:
    //! The base traits.
    const Base_traits_2 * m_base;

    /*! Constructor.
     * \param tr The base traits class. It must be passed, to handle
     *           non stateless traits (e.g., it stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_near_boundary_2(const Base_traits_2 * tr) : m_base(tr) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce) const
    {
      // If the traits class does not support open curves, we just
      // return EQUAL, as this comparison will not be invoked anyway.
      return m_base->compare_y_near_boundary_2_object()(xcv1.base(),
                                                        xcv2.base(), ce);
    }
  };

  /*! Obtain a Compare_y_near_boundary_2 functor object. */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object () const
  {
    return Compare_y_near_boundary_2(m_base_traits);
  }

  // bottom-top

  /*! A functor that determines whether an endpoint of an x-monotone arc lies
   * on a boundary of the parameter space along the y axis.
   */
  class Parameter_space_in_y_2 {
  protected:
    //! The base traits.
    const Base_traits_2      *m_base;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Parameter_space_in_y_2(const Base_traits_2 *tr) : m_base (tr) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    Arr_parameter_space operator() (const X_monotone_curve_2& xcv,
                                    Arr_curve_end ce) const
    {
      return m_base->parameter_space_in_y_2_object() (xcv.base(), ce);
    }

    Arr_parameter_space operator() (const Point_2 & p) const
    {
      return m_base->parameter_space_in_y_2_object()(p.base());
    }

    Arr_parameter_space operator() (const X_monotone_curve_2 & xcv) const
    {
      return m_base->parameter_space_in_y_2_object()(xcv.base());
    }

  };

  /*! Obtain a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object () const
  {
    return Parameter_space_in_y_2 (m_base_traits);
  }

  /*! A function object that determines whether an x-monotone curve or a
   * point coincide with the horizontal identification curve.
   */
  class Is_on_y_identification_2 {
  protected:
    //! The base traits.
    const Base_traits_2      *m_base;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Is_on_y_identification_2 (const Base_traits_2* tr) : m_base (tr) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    bool operator() (const Point_2 & p) const
    {
      return m_base->is_on_y_identification_2_object() (p.base());
    }

    bool operator() (const X_monotone_curve_2 & xcv) const
    {
      return m_base->is_on_y_identification_2_object() (xcv.base());
    }
  };

  /*! Obtain a Is_on_y_identification_2 function object */
  Is_on_y_identification_2 is_on_y_identification_2_object () const
  {
    return Is_on_y_identification_2 (m_base_traits);
  }

  /*! A functor that compares the x-limits of curve ends on the
   * boundary of the parameter space.
   */
  class Compare_x_at_limit_2 {
  protected:
    //! The base traits.
    const Base_traits_2 * m_base;

    /*! Constructor.
     * \param tr The base traits class. It must be passed, to handle
     *           non stateless traits (e.g., it stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_at_limit_2(const Base_traits_2 * tr) : m_base(tr) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    Comparison_result operator() (const Point_2& p,
                                  const X_monotone_curve_2 & xcv,
                                  Arr_curve_end ce) const
    {
      return m_base->compare_x_at_limit_2_object()(p.base(),
                                                   xcv.base(), ce);
    }

    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce2) const
    {
      return m_base->compare_x_at_limit_2_object()(xcv1.base(), ce1,
                                                   xcv2.base(), ce2);
    }
  };

  /*! Obtain a Compare_x_at_limit_2 function object. */
  Compare_x_at_limit_2 compare_x_at_limit_2_object () const
  {
    return Compare_x_at_limit_2(m_base_traits);
  }

  /*! A functor that compares the x-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_x_near_limit_2 {
  protected:
    //! The base traits.
    const Base_traits_2 * m_base;

    /*! Constructor.
     * \param tr The base traits class. It must be passed, to handle
     *           non stateless traits (e.g., it stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_near_limit_2(const Base_traits_2 * tr) : m_base(tr) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce) const
    {
      return m_base->compare_x_near_limit_2_object()(xcv1.base(),
                                                     xcv2.base(),
                                                     ce);
    }
  };

  /*! Obtain a Compare_x_near_limit_2 function object. */
  Compare_x_near_limit_2 compare_x_near_limit_2_object () const
  {
    return Compare_x_near_limit_2(m_base_traits);
  }

  /*! A functor that compares the x-coordinate of two given points
   * that lie on the horizontal identification curve.
   */
  class Compare_x_on_boundary_2 {
  protected:
    //! The base traits.
    const Base_traits_2 * m_base;

    /*! Constructor.
     * \param tr The base traits class. It must be passed, to handle
     *           non stateless traits (e.g., it stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_on_boundary_2(const Base_traits_2 * tr) : m_base(tr) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      return m_base->compare_x_on_boundary_2_object()(p1.base(), p2.base());
    }

    Comparison_result operator()(const Point_2 & p1, const X_monotone_curve_2& xcv, Arr_curve_end ce) const
    {
      return m_base->compare_x_on_boundary_2_object()(p1.base(), xcv.base(), ce);
    }

    Comparison_result operator()( const X_monotone_curve_2& xcv1, Arr_curve_end ce1,
                                  const X_monotone_curve_2& xcv2, Arr_curve_end ce2) const
    {
      return m_base->compare_x_on_boundary_2_object()(xcv1.base(), ce1, xcv2.base(), ce2);
    }
  };

  /*! Obtain a Compare_x_on_boundary_2 functor object. */
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object () const
  {
    return Compare_x_on_boundary_2(m_base_traits);
  }

  /*! A functor that compares the x-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_x_near_boundary_2 {
  protected:
    //! The base traits.
    const Base_traits_2 * m_base;

    /*! Constructor.
     * \param tr The base traits class. It must be passed, to handle
     *           non stateless traits (e.g., it stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_near_boundary_2(const Base_traits_2 * tr) : m_base(tr) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_batched_point_location_traits_2<Arrangement_2>;

  public:
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce) const
    {
      return m_base->compare_x_near_boundary_2_object()(xcv1.base(),
                                                        xcv2.base(),
							ce);
    }
  };

  /*! Obtain a Compare_x_near_boundary_2 function object. */
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object () const
  {
    return Compare_x_near_boundary_2(m_base_traits);
  }

};

} //namespace CGAL

#endif
