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
// $URL$
// $Id$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_ARR_BATCHED_POINT_LOCATION_TRAITS_2_H
#define CGAL_ARR_BATCHED_POINT_LOCATION_TRAITS_2_H

/*!
 * Definition of the Arr_batched_point_location_traits_2<Arrangement> class.
 */

CGAL_BEGIN_NAMESPACE

/*! \class
 * A traits-class decorator for the use of the batched point-location process.
 */
template <class Arrangement_>
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
  typedef typename Base_traits_2::Compare_x_2      Base_compare_x_2;
  typedef typename Base_traits_2::Compare_xy_2     Base_compare_xy_2;
  typedef typename Base_traits_2::Compare_y_at_x_2 Base_compare_y_at_x_2;
  typedef typename Base_traits_2::Compare_y_at_x_right_2 
                                                   Base_compare_y_at_x_right_2;
  typedef typename Base_traits_2::Equal_2          Base_equal_2;
  typedef typename Base_traits_2::Is_vertical_2    Base_is_vertical_2;

  typedef typename Base_traits_2::Has_boundary_category
                                                   Base_has_boundary_category;
  typedef Tag_true                                 Has_boundary_category;
  typedef Tag_false                                Has_left_category;

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
    
    Base_x_monotone_curve_2 m_base_cv;  // The base x-monotone curve.
    Halfedge_const_handle   m_he;       // The corresponding arrangement edge.

  public:

    Ex_x_monotone_curve_2 ():
      m_base_cv(),
      m_he()
    {}

    Ex_x_monotone_curve_2 (const Base& cv):
      m_base_cv(cv),
      m_he()
    {}

    Ex_x_monotone_curve_2 (const Base& cv, Halfedge_const_handle he) :
      m_base_cv(cv),
      m_he(he)
    {
      CGAL_precondition (he->direction() == RIGHT_TO_LEFT);
    }

    Halfedge_const_handle halfedge_handle() const
    {
      return (m_he);
    }

    const Base& base () const
    {
      return (m_base_cv);
    }

    Base& base ()
    {
      return (m_base_cv);
    }

    operator const Base&() const
    {
      return (m_base_cv);
    }

    operator Base&()
    {
      return (m_base_cv);
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


  /*! \class
   * The Boundary_in_x_2 functor.
   */
  class Boundary_in_x_2
  {
  private:

    const Base_traits_2      *m_base;

  public:

    Boundary_in_x_2 (const Base_traits_2 *tr) :
      m_base (tr)
    {}

    Boundary_type operator() (const X_monotone_curve_2& cv,
                              Curve_end ind) const
    {
      return _boundary_in_x_imp (cv, ind, Base_has_boundary_category());
    }

  private:

    Boundary_type _boundary_in_x_imp (const X_monotone_curve_2& cv,
                                      Curve_end ind,
                                      Tag_true) const
    {
      return (m_base->boundary_in_x_2_object() (cv.base(), ind));
    }

    Boundary_type _boundary_in_x_imp (const X_monotone_curve_2& , Curve_end ,
                                      Tag_false) const
    {
      return (NO_BOUNDARY);
    }
  };

  Boundary_in_x_2 boundary_in_x_2_object () const
  {
    return (Boundary_in_x_2 (m_base_traits));
  }

  /*! \class
   * The Boundary_in_y_2 functor.
   */
  class Boundary_in_y_2
  {
  private:

    const Base_traits_2      *m_base;

  public:

    Boundary_in_y_2 (const Base_traits_2 *tr) :
      m_base (tr)
    {}

    Boundary_type operator() (const X_monotone_curve_2& cv,
                              Curve_end ind) const
    {
      return _boundary_in_y_imp (cv, ind, Base_has_boundary_category());
    }

  private:

    Boundary_type _boundary_in_y_imp (const X_monotone_curve_2& cv,
                                      Curve_end ind,
                                      Tag_true) const
    {
      return (m_base->boundary_in_y_2_object() (cv.base(), ind));
    }

    Boundary_type _boundary_in_y_imp (const X_monotone_curve_2& , Curve_end ,
                                      Tag_false) const
    {
      return (NO_BOUNDARY);
    }
  };

  Boundary_in_y_2 boundary_in_y_2_object () const
  {
    return (Boundary_in_y_2 (m_base_traits));
  }

  /*! \class
   * The Construct_min_vertex_2 functor.
   */
  class Construct_min_vertex_2
  {
  private:
    Base_construct_min_vertex_2 m_base_min_v;

  public:

    Construct_min_vertex_2(const Base_construct_min_vertex_2& base):
        m_base_min_v(base)
    {}

    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    Point_2 operator() (const X_monotone_curve_2 & cv) 
    {
      // Note that the halfedge associated with the curve is always directed
      // from right to left, so its target is the leftmost vertex.
      Vertex_const_handle vh = cv.halfedge_handle()->target();
      return (Point_2 (m_base_min_v (cv.base()), vh));
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return 
      (Construct_min_vertex_2(m_base_traits->construct_min_vertex_2_object()));
  }

  /*! \class
   * The Construct_max_vertex_2 functor.
   */
  class Construct_max_vertex_2
  {
  private:
    Base_construct_max_vertex_2 m_base_max_v;

  public:

    Construct_max_vertex_2 (const Base_construct_max_vertex_2& base):
        m_base_max_v(base)
    {}

    /*!
     * Get the right endpoint of the x-monotone curve .
     * \param cv The curve.
     * \return The right endpoint.
     */
    Point_2 operator() (const X_monotone_curve_2 & cv) 
    {
      // Note that the halfedge associated with the curve is always directed
      // from right to left, so its source is the rightmost vertex.
      Vertex_const_handle vh = cv.halfedge_handle()->source();
      return (Point_2 (m_base_max_v (cv.base()), vh));
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return
      (Construct_max_vertex_2(m_base_traits->construct_max_vertex_2_object()));
  }

  /*! \class
   * The Comapre_xy_2 functor.
   */
  class Compare_xy_2
  {
  private:
    Base_compare_xy_2    m_base_cmp_xy;
    Vertex_const_handle  invalid_v;

  public:

    Compare_xy_2(const Base_compare_xy_2& base):
      m_base_cmp_xy(base),
      invalid_v()
    {}

    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
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

  /*! Get a Construct_min_vertex_2 functor object. */
  Compare_xy_2 compare_xy_2_object () 
  {
    return Compare_xy_2(m_base_traits->compare_xy_2_object());
  }

  /*! \class
   * The Comapre_y_at_x_2 functor.
   */
  class Compare_y_at_x_2
  {
  private:
    Base_compare_y_at_x_2 m_base_cmp_y_at_x;

  public:
    Compare_y_at_x_2(const Base_compare_y_at_x_2& base):
        m_base_cmp_y_at_x(base)
    {}

    Comparison_result operator() (const Point_2& p,
                                  const X_monotone_curve_2& cv) const
    {
      return (m_base_cmp_y_at_x (p.base(), cv.base()));
    }

    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2, 
                                  Curve_end ind) const
    {
      // The function is implemented based on the Has_boundary category.
      // If the traits class does not support unbounded curves, we just
      // return EQUAL, as this comparison will not be invoked anyway.
      return _comp_y_at_infinity_imp (cv1, cv2, ind, 
                                      Base_has_boundary_category());
    }

  private:

    Comparison_result _comp_y_at_infinity_imp (const X_monotone_curve_2& cv1,
                                               const X_monotone_curve_2& cv2, 
                                               Curve_end ind,
                                               Tag_true) const
    {
      return (m_base_cmp_y_at_x (cv1.base(), cv2.base(), ind));
    }

    Comparison_result _comp_y_at_infinity_imp (const X_monotone_curve_2& ,
                                               const X_monotone_curve_2& , 
                                               Curve_end ,
                                               Tag_false) const
    {
      return (EQUAL);
    }
  };

  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return (Compare_y_at_x_2 (m_base_traits->compare_y_at_x_2_object()));
  }

  /*! \class
   * The Comapre_y_at_x_right_2 functor.
   */
  class Compare_y_at_x_right_2
  {
  private:
    Base_compare_y_at_x_right_2 m_base_cmp_y_at_x_right;

  public:
    Compare_y_at_x_right_2(const Base_compare_y_at_x_right_2& base):
        m_base_cmp_y_at_x_right(base)
    {}

    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& p) const
    {
      return (m_base_cmp_y_at_x_right(cv1.base(),
                                      cv2.base(),
                                      p.base()));
    }
  };

  Compare_y_at_x_right_2 compare_y_at_x_right_2_object () const
  {
    return (Compare_y_at_x_right_2
	    (m_base_traits->compare_y_at_x_right_2_object()));
  }

  /*! \class
   * The Equal_2 functor.
   */
  class Equal_2
  {
  private:
    Base_equal_2           m_base_eq;
    Vertex_const_handle    invalid_v;
    Halfedge_const_handle  invalid_he;

  public:
    
    Equal_2(const Base_equal_2& base):
      m_base_eq(base),
      invalid_v(),
      invalid_he()
    {}

    /*! Check if two curves are the same. */
    bool operator() (const X_monotone_curve_2& cv1,
		     const X_monotone_curve_2& cv2) const
    {
      if (cv1.halfedge_handle() == cv2.halfedge_handle() &&
          cv1.halfedge_handle() != invalid_he)
        return (true);

      return (m_base_eq(cv1.base(), cv2.base()));
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

  Equal_2 equal_2_object () const
  {
    return (Equal_2 (m_base_traits->equal_2_object()));
  }

  /*! \class
   * The Comapre_x_2 functor.
   */
  class Compare_x_2
  {
  private:
    Base_compare_x_2 m_base_cmp_x;

  public:
    Compare_x_2(const Base_compare_x_2& base):
        m_base_cmp_x(base)
    {}

    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      return (m_base_cmp_x(p1.base(), p2.base()));
    }

    Comparison_result operator() (const Point_2& p,
                                  const X_monotone_curve_2& cv,
                                  Curve_end ind) const
    {
      return (_compare_point_curve_imp (p, cv, ind,
                                        Base_has_boundary_category()));
    }

    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  Curve_end ind1,
                                  const X_monotone_curve_2& cv2,
                                  Curve_end ind2) const
    {
      return (_compare_curves_imp (cv1, ind1, cv2, ind2,
                                   Base_has_boundary_category()));
    }

  private:

    Comparison_result _compare_point_curve_imp (const Point_2& p,
                                                const X_monotone_curve_2& cv,
                                                Curve_end ind,
                                                Tag_true) const
    {
      return (m_base_cmp_x (p.base(), cv.base(), ind));
    }

    Comparison_result _compare_point_curve_imp (const Point_2& ,
                                                const X_monotone_curve_2& ,
                                                Curve_end ,
                                                Tag_false) const
    {
      return (EQUAL);
    }

    Comparison_result _compare_curves_imp (const X_monotone_curve_2& cv1, 
                                           Curve_end ind1,
                                           const X_monotone_curve_2& cv2,
                                           Curve_end ind2,
                                           Tag_true) const
    {
      return (m_base_cmp_x (cv1.base(), ind1, cv2.base(), ind2));
    }

    Comparison_result _compare_curves_imp (const X_monotone_curve_2& ,
                                           Curve_end ,
                                           const X_monotone_curve_2& , 
                                           Curve_end ,
                                           Tag_false) const
    {
      return (EQUAL);
    }

  };

  Compare_x_2 compare_x_2_object () const
  {
    return (Compare_x_2 (m_base_traits->compare_x_2_object()));
  }

  /*! \class
   * The Is_vertical_2 functor.
   */
  class Is_vertical_2
  {
  private:
    Base_is_vertical_2 m_base_is_vert;

  public:
    Is_vertical_2(const Base_is_vertical_2& base):
        m_base_is_vert(base)
    {}

     bool operator() (const X_monotone_curve_2& cv) const
    {
      return (m_base_is_vert(cv.base()));
    }
  };

  Is_vertical_2 is_vertical_2_object() const
  {
    return (Is_vertical_2(m_base_traits->is_vertical_2_object()));
  }
};


CGAL_END_NAMESPACE

#endif
