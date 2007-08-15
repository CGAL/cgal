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

#ifndef CGAL_ARR_OVERLAY_TRAITS_2_H
#define CGAL_ARR_OVERLAY_TRAITS_2_H

/*!
 * Defintion of the Arr_overlay_traits_2 class-template.
 */

#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A meta-traits class that stores a red or a blue halfedge handle with every
 * x-monotone curve, and a red or blue vertex handle with each point. This
 * information is used to speed up the overlay of a red arrangement and a blue
 * arrangement one on top of the other.
 */
template <class Traits_,
          class ArrangementRed_,
          class ArrangementBlue_>
class Arr_overlay_traits_2 : public Traits_
{
public:

  typedef Traits_                                   Traits_2;
  typedef ArrangementRed_                           Arrangement_red_2;
  typedef ArrangementBlue_                          Arrangement_blue_2;

  typedef typename Arrangement_red_2::Face_const_handle 
                                                    Face_handle_red;
  typedef typename Arrangement_blue_2::Face_const_handle
                                                    Face_handle_blue;

  typedef typename Arrangement_red_2::Halfedge_const_handle 
                                                    Halfedge_handle_red;
  typedef typename Arrangement_blue_2::Halfedge_const_handle
                                                    Halfedge_handle_blue;

  typedef typename Arrangement_red_2::Vertex_const_handle 
                                                    Vertex_handle_red;
  typedef typename Arrangement_blue_2::Vertex_const_handle
                                                    Vertex_handle_blue;

  typedef typename Traits_2::X_monotone_curve_2     Base_x_monotone_curve_2;
  typedef typename Traits_2::Point_2                Base_point_2;
  typedef typename Traits_2::Intersect_2            Base_intersect_2;
  typedef typename Traits_2::Split_2                Base_split_2;
  typedef typename Traits_2::Construct_min_vertex_2 
                                                   Base_construct_min_vertex_2;
  typedef typename Traits_2::Construct_max_vertex_2
                                                   Base_construct_max_vertex_2;
  typedef typename Traits_2::Compare_xy_2          Base_compare_xy_2;
  typedef typename Traits_2::Compare_y_at_x_2      Base_compare_y_at_x_2;
  typedef typename Traits_2::Compare_y_at_x_right_2 
                                                   Base_compare_y_at_x_right_2;
  typedef typename Traits_2::Compare_x_2           Base_compare_x_2;
  typedef typename Traits_2::Equal_2               Base_equal_2;
  typedef typename Traits_2::Has_boundary_category Base_has_boundary_category;
  
  typedef Tag_true                                    Has_boundary_category;
  typedef Tag_false                                   Has_left_category;
  typedef Tag_false                                   Has_merge_category;

  // The color of a feature.
  enum Color
  {
    RED,          // From the "red" arrangement.
    BLUE,         // From the "blue" arrangement.
    RB_OVERLAP    // Red-blue overlap.
  };

private:

  Traits_2    *m_base_traits;        // The base traits object.

public:

  /*! Default constructor. */
  Arr_overlay_traits_2()
  {}

  /*! Constructor from a base traits class. */
  Arr_overlay_traits_2 (Traits_2& base_tr) :
    m_base_traits (&base_tr)
  {}

  /*! \class
   * Nested extension of the x-monotone curve type.
   */
  class Ex_x_monotone_curve_2   
  {
  public:

    typedef  Base_x_monotone_curve_2     Base;

  protected:
    
    Base                  m_base_cv;                // The base curve.
    Halfedge_handle_red   m_red_halfedge_handle;    // The red halfedge.
    Halfedge_handle_blue  m_blue_halfedge_handle;   // The blue halfedge.

  public:

    /*! Default constructor. */
    Ex_x_monotone_curve_2() : 
      m_base_cv(),
      m_red_halfedge_handle(),
      m_blue_halfedge_handle()
    {}

    /*! Constructor from a curve. */
    Ex_x_monotone_curve_2 (const Base& cv):
      m_base_cv(cv),
      m_red_halfedge_handle(),
      m_blue_halfedge_handle()
    {}

    /*! Constructor from a curve and halfedge handles. */
    Ex_x_monotone_curve_2 (const Base& cv,
                           Halfedge_handle_red  he_r,
                           Halfedge_handle_blue he_b): 
      m_base_cv(cv),
      m_red_halfedge_handle (he_r),
      m_blue_halfedge_handle (he_b)
    {
      CGAL_precondition (he_r == Halfedge_handle_red() ||
                         he_r->direction() == RIGHT_TO_LEFT);
      CGAL_precondition (he_b == Halfedge_handle_blue() ||
                         he_b->direction() == RIGHT_TO_LEFT);
    }

    /*! Get the base curve (const version). */
    const Base& base () const
    {
      return (m_base_cv);
    }

    /*! Get the base curve (non-const version). */
    Base& base ()
    {
      return (m_base_cv);
    }

    /*! Casting to a base curve (const version). */
    operator const Base& () const
    {
      return (m_base_cv);
    }

    /*! Casting to a base curve (const version). */
    operator Base&()
    {
      return (m_base_cv);
    }

    /*! Get the red halfedge handle. */
    Halfedge_handle_red red_halfedge_handle() const
    {
      return (m_red_halfedge_handle);
    }
    
    /*! Get the blue halfedge handle. */
    Halfedge_handle_blue blue_halfedge_handle() const
    {
      return (m_blue_halfedge_handle);
    }

    /*! Set the red halfedge handle. */
    void set_red_halfedge_handle (Halfedge_handle_red he_r)
    {
      CGAL_precondition (he_r == Halfedge_handle_red() ||
                         he_r->direction() == RIGHT_TO_LEFT);

      m_red_halfedge_handle = he_r;
    }

    /*! Set the blue halfedge handle. */
    void set_blue_halfedge_handle (Halfedge_handle_blue he_b)
    {
      CGAL_precondition (he_b == Halfedge_handle_blue() ||
                         he_b->direction() == RIGHT_TO_LEFT);

      m_blue_halfedge_handle = he_b;
    }

    /*! Get the color of the subcurve. */
    Color color () const
    {
      Halfedge_handle_red     null_red_he;
      Halfedge_handle_blue    null_blue_he;

      if (m_red_halfedge_handle != null_red_he &&
          m_blue_halfedge_handle == null_blue_he)
        return (RED);

      if (m_blue_halfedge_handle != null_blue_he &&
          m_red_halfedge_handle == null_red_he)
        return (BLUE);
      
      // Overlap, return the RB_OVERLAP color:
      CGAL_assertion (m_red_halfedge_handle != null_red_he && 
                      m_blue_halfedge_handle != null_blue_he);
      return (RB_OVERLAP);
    }
    
  }; // nested class Ex_x_monotone_curve_2 - END
  
  typedef Ex_x_monotone_curve_2                     X_monotone_curve_2;

#ifdef VERBOSE
  // For debugging purposes:
  friend std::ostream& operator<< (std::ostream& os,
                                   const X_monotone_curve_2& cv)
  {
    os << cv.base();
    return (os);
  }
#endif
 
  /*! \class
   * Nested extension of the point type.
   */
  class Ex_point_2 
  {
  public:

    typedef Base_point_2    Base;

  protected:

    Base          m_base_pt;        // The base point.
    Object        m_red_obj;        // The "red" object.
    Object        m_blue_obj;       // The "blue" object.

  public:

    /*! Default constructor. */
    Ex_point_2() :
      m_base_pt(),
      m_red_obj(),
      m_blue_obj()
    {}

    /*! Constructor from a base point. */
    Ex_point_2 (const Base& pt) :
      m_base_pt (pt),
      m_red_obj(),
      m_blue_obj()
    {}

    /*! Constructor from a base point with red and blue objects. */
    Ex_point_2 (const Base& pt, const Object& obj_r, const Object& obj_b) :
      m_base_pt (pt),
      m_red_obj (obj_r),
      m_blue_obj (obj_b)
    {}

    /*! Get the base point (const version). */
    const Base& base () const
    {
      return (m_base_pt);
    }

    /*! Get the base point (non-const version). */
    Base& base ()
    {
      return (m_base_pt);
    }

    /*! Casting to a base point (const version). */
    operator const Base& () const
    {
      return (m_base_pt);
    }

    /*! Casting to a base point (non-const version). */
    operator Base& ()
    {
      return (m_base_pt);
    }

    /*! Get the red object. */
    const Object& red_object() const  
    { 
      return (m_red_obj);  
    }

    /*! Get the blue object. */
    const Object& blue_object() const
    {
      return (m_blue_obj); 
    }

    /*! Check if the red object is empty. */
    bool is_red_object_empty () const
    {
      return (m_red_obj.is_empty());
    }

    /*! Check if the blue object is empty. */
    bool is_blue_object_empty () const
    {
      return (m_blue_obj.is_empty());
    }

    /*! Set the red object. */
    void set_red_object (const Object& obj_r)
    {
      m_red_obj = obj_r;
    }

    /*! Set the blue object. */
    void set_blue_object (const Object& obj_b)
    {
      m_blue_obj = obj_b;
    }

  };

  typedef Ex_point_2                                Point_2;

#ifdef VERBOSE
  // For debugging purposes:
  friend std::ostream& operator<< (std::ostream& os,
                                   const Point_2& pt)
  {
    os << pt.base();
    return (os);
  }
#endif

  /*! \class
   * The Intersect_2 functor.
   */ 
  class Intersect_2
  {
  private:

    Traits_2         *m_base_tr;

  public:
   
    /*! Constructor. */
    Intersect_2 (Traits_2 *base) :
      m_base_tr (base)
    {}

    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi)
    {
      // In case the curves originate from the same arrangement, they are
      // obviously interior-disjoint.
      if (cv1.color() == cv2.color())
        return (oi);
      
      if (cv1.color() == RB_OVERLAP || cv2.color() == RB_OVERLAP)
        return (oi);

      // Compute the intersection points between the curves. Note that if
      // cv1 and cv2 are subcruves of x-monotone curves that had intersected
      // before the current point on the status line, we may get a filter
      // failure if we send the subcurve whose left endpoint is to the left
      // of the other curve - this is because their previously computed
      // intersection point p may be equal to the this left endpoint. As many
      // traits classes start by computing the intersection between the
      // supporting curves and then check whether the result is in the x-range
      // of both subcurves, this will result in a filter failure. However, if
      // we send cv1 first, then p is obviusly not in its x-range and there is
      // no filter failure.
      //
      //              / cv1
      //             /
      //            /
      //       ----+--
      //          /
      //         /
      //      p +------------- cv2
      //              ^
      //              |
      //              status line
      //
      // Note that we do not bother with curves whose left ends are unbounded,
      // since such curved did not intersect before.
      const std::pair<Base_point_2, unsigned int>   *base_ipt;
      const Base_x_monotone_curve_2                 *overlap_cv;
      bool                                           send_cv1_first = true;
      OutputIterator                                 oi_end;

      Boundary_in_x_2      inf_in_x (m_base_tr);
      Boundary_in_y_2      inf_in_y (m_base_tr);
      const Boundary_type  bx1 = inf_in_x (cv1, MIN_END);
      const Boundary_type  by1 = inf_in_y (cv1, MIN_END);
      const Boundary_type  bx2 = inf_in_x (cv2, MIN_END);
      const Boundary_type  by2 = inf_in_y (cv2, MIN_END);

      if (bx1 == NO_BOUNDARY && by1 == NO_BOUNDARY &&
          bx2 == NO_BOUNDARY && by2 == NO_BOUNDARY)
      {
        send_cv1_first =
          (m_base_tr->compare_xy_2_object()
           (m_base_tr->construct_min_vertex_2_object()(cv1.base()),
            m_base_tr->construct_min_vertex_2_object()(cv2.base())) == LARGER);
      }

      if (send_cv1_first)
        oi_end = m_base_tr->intersect_2_object()(cv1.base(),
                                                 cv2.base(), oi);
      else
        oi_end = m_base_tr->intersect_2_object()(cv2.base(),
                                                 cv1.base(), oi);

      // Convert objects that are associated with Base_x_monotone_curve_2 to
      // the exteneded X_monotone_curve_2. 
      for (; oi != oi_end; ++oi)
      {
        base_ipt = object_cast<std::pair<Base_point_2, unsigned int> >(&(*oi));

        if (base_ipt != NULL)
        {
          // We have an red-blue intersection point, so we attach the
          // intersecting red and blue halfedges to it.
          Object   red_obj;
          Object   blue_obj;

          if (cv1.color() == RED)
          {
            CGAL_assertion (cv2.color() == BLUE);
            red_obj = CGAL::make_object (cv1.red_halfedge_handle());
            blue_obj = CGAL::make_object (cv2.blue_halfedge_handle());
          }
          else
          {
            CGAL_assertion(cv2.color() == RED &&
                           cv1.color() == BLUE);
            red_obj = CGAL::make_object(cv2.red_halfedge_handle());
            blue_obj = CGAL::make_object(cv1.blue_halfedge_handle());
          }

          // Create the extended point and add the multiplicity.
          Point_2   ex_point (base_ipt->first,
                              red_obj, blue_obj);
          *oi = CGAL::make_object(std::make_pair (ex_point, 
                                                  base_ipt->second));
        }
        else
        {
          overlap_cv = object_cast<Base_x_monotone_curve_2> (&(*oi));
          CGAL_assertion (overlap_cv != NULL);

          // We have a red-blue overlap, so we mark the curve accordingly.
          Halfedge_handle_red        red_he;
          Halfedge_handle_blue       blue_he;
          
          if (cv1.color() == RED)
          {
            red_he = cv1.red_halfedge_handle();
            
            // Overlap can occur only between curves from a different color.
            CGAL_assertion(cv2.color() == BLUE);
            blue_he = cv2.blue_halfedge_handle();
          }
          else
          {
            CGAL_assertion(cv1.color() == BLUE &&
                           cv2.color() == RED);

            red_he = cv2.red_halfedge_handle();
            blue_he = cv1.blue_halfedge_handle();
          }

          *oi = CGAL::make_object (X_monotone_curve_2 (*overlap_cv,
                                                       red_he, blue_he));
        }
      }

      // Return the past-the-end iterator.
      return (oi_end);
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const
  {
    return Intersect_2(m_base_traits); 
  }


  /*! \class
   * The Split_2 functor.
   */
  class Split_2
  {
  private:
    Base_split_2    m_base_split;

  public:

    /*! Constructor. */
    Split_2 (const Base_split_2& base) :
        m_base_split (base)
    {}

    void operator() (const X_monotone_curve_2& cv, const Point_2 & p,
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2)
    {
      // Split the base curve.
      m_base_split(cv.base(),
                   p.base(),
                   c1.base(),
                   c2.base());

      // Duplicate the halfedge data to the resulting curves.
      c1.set_red_halfedge_handle (cv.red_halfedge_handle());
      c1.set_blue_halfedge_handle (cv.blue_halfedge_handle());

      c2.set_red_halfedge_handle (cv.red_halfedge_handle());
      c2.set_blue_halfedge_handle (cv.blue_halfedge_handle());
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object () const
  {
    return Split_2(m_base_traits->split_2_object());
  }

  /*! \class
   * The Construct_min_vertex_2 functor.
   */
  class Construct_min_vertex_2
  {
  private:
    Base_construct_min_vertex_2  m_base_min_v;
    Base_equal_2                 m_base_equal;

  public:

    Construct_min_vertex_2 (const Base_construct_min_vertex_2& base_min_v,
                            const Base_equal_2& base_equal) :
      m_base_min_v (base_min_v),
      m_base_equal (base_equal)
    {}

    Point_2 operator() (const X_monotone_curve_2& cv) 
    {
      // Create the objects that wrap the arrangement vertex.
      // Note that the halfedges associated with the curves are always
      // directed from right to left, so their target is the smaller end.
      const Base_point_2&   base_p = m_base_min_v (cv.base());
      Object                obj_red, obj_blue;

      if (cv.color() == RED)
      {
        obj_red = CGAL::make_object (cv.red_halfedge_handle()->target());
      }
      else if (cv.color() == BLUE)
      {
        obj_blue = CGAL::make_object (cv.blue_halfedge_handle()->target());
      }
      else
      {
        CGAL_assertion (cv.color() == RB_OVERLAP);

        if (! cv.red_halfedge_handle()->target()->is_at_infinity() &&
            m_base_equal (base_p,
                          cv.red_halfedge_handle()->target()->point()))
        {
          obj_red = CGAL::make_object (cv.red_halfedge_handle()->target());
        }

        if (! cv.blue_halfedge_handle()->target()->is_at_infinity() &&
            m_base_equal (base_p,
                          cv.blue_halfedge_handle()->target()->point()))
        {
          obj_blue = CGAL::make_object (cv.blue_halfedge_handle()->target());
        }
      }

      return (Point_2 (base_p, obj_red, obj_blue));
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return 
      (Construct_min_vertex_2 (m_base_traits->construct_min_vertex_2_object(),
                               m_base_traits->equal_2_object()));
  }

  /*! \class
   * The Construct_max_vertex_2 functor.
   */
  class Construct_max_vertex_2
  {
  private:
    Base_construct_max_vertex_2  m_base_max_v;
    Base_equal_2                 m_base_equal;

  public:

    Construct_max_vertex_2 (const Base_construct_max_vertex_2& base_max_v,
                            const Base_equal_2& base_equal) :
      m_base_max_v (base_max_v),
      m_base_equal (base_equal)
    {}

    Point_2 operator() (const X_monotone_curve_2 & cv) const
    {
      // Create the objects that wrap the arrangement vertex.
      // Note that the halfedges associated with the curves are always
      // directed from right to left, so their target is the smaller end.
      const Base_point_2&   base_p = m_base_max_v (cv.base());
      Object                obj_red, obj_blue;

      if(cv.color() == RED)
      {
        obj_red = CGAL::make_object (cv.red_halfedge_handle()->source());
      }
      else if(cv.color() == BLUE)
      {
        obj_blue = CGAL::make_object (cv.blue_halfedge_handle()->source());
      }
      else
      {
        CGAL_assertion(cv.color() == RB_OVERLAP);

        if (! cv.red_halfedge_handle()->source()->is_at_infinity() &&
            m_base_equal (base_p,
                          cv.red_halfedge_handle()->source()->point()))
        {
          obj_red = CGAL::make_object (cv.red_halfedge_handle()->source());
        }

        if (! cv.blue_halfedge_handle()->source()->is_at_infinity() &&
            m_base_equal (base_p,
                          cv.blue_halfedge_handle()->source()->point()))
        {
          obj_blue = CGAL::make_object (cv.blue_halfedge_handle()->source());
        }
      }

      return (Point_2 (base_p, obj_red, obj_blue));
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return
      (Construct_max_vertex_2 (m_base_traits->construct_max_vertex_2_object(),
                               m_base_traits->equal_2_object()));
  }

  /*! \class
   * The Comapre_xy_2 functor.
   */
  class Compare_xy_2
  {
  private:
    Base_compare_xy_2 m_base_cmp_xy;

  public:

    Compare_xy_2(const Base_compare_xy_2& base):
        m_base_cmp_xy(base)
    {}

    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      // Check if there wither points represent red or blue vertices.
      Vertex_handle_red  vr1;
      Vertex_handle_red  vr2;

      Vertex_handle_blue vb1;
      Vertex_handle_blue vb2;

      const bool         assign_v1_red  = CGAL::assign (vr1, p1.red_object());
      const bool         assign_v2_red  = CGAL::assign (vr2, p2.red_object());
      const bool         assign_v1_blue = CGAL::assign (vb1, p1.blue_object());
      const bool         assign_v2_blue = CGAL::assign (vb2, p2.blue_object());

      if ((assign_v1_red && assign_v1_blue) ||
          (assign_v2_red && assign_v2_blue))
      {
        // In case of an overlapping vertex, just perform the comparison.
        return (m_base_cmp_xy (p1.base(), p2.base()));
      }

      // If both points are vertices of the same color, avoid the geometric
      // comparison if they refer to the same vertex.
      if ((assign_v1_red && assign_v2_red && vr1 == vr2) ||
          (assign_v1_blue && assign_v2_blue && vb1 == vb2))
      {
        return (EQUAL);
      }

      return (m_base_cmp_xy (p1.base(), p2.base()));
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const
  {
    return (Compare_xy_2(m_base_traits->compare_xy_2_object()));
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
      return (m_base_cmp_x (p1.base(), p2.base()));
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
                                           Curve_end,
                                           const X_monotone_curve_2& , 
                                           Curve_end ind2,
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

    Comparison_result operator() (const Point_2 & p,
                                  const X_monotone_curve_2 & cv) const
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
      return _comp_y_at_boundary_imp (cv1, cv2, ind, 
                                      Base_has_boundary_category());
    }

  private:

    Comparison_result _comp_y_at_boundary_imp (const X_monotone_curve_2& cv1,
                                               const X_monotone_curve_2& cv2, 
                                               Curve_end ind,
                                               Tag_true) const
    {
      return (m_base_cmp_y_at_x (cv1.base(), cv2.base(), ind));
    }

    Comparison_result _comp_y_at_boundary_imp (const X_monotone_curve_2& ,
                                               const X_monotone_curve_2& , 
                                               Curve_end ,
                                               Tag_false) const
    {
      return (EQUAL);
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return (Compare_y_at_x_2(m_base_traits->compare_y_at_x_2_object()));
  }

  /*! \class
   * The Comapre_y_at_x_right_2 functor.
   */
  class Compare_y_at_x_right_2
  {
  private:
    Base_compare_y_at_x_right_2    m_base_cmp_y_at_x_right;

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

  /*! Get a Construct_min_vertex_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object () const
  {
    return 
      (Compare_y_at_x_right_2(m_base_traits->compare_y_at_x_right_2_object()));
  }

  /*!
   * The Boundary_in_x_2 functor.
   */
  class Boundary_in_x_2
  {
  private:
    const Traits_2      *m_base;

  public:
       
    Boundary_in_x_2 (const Traits_2 *tr) :
      m_base (tr)
    {}
 
    Boundary_type operator() (const X_monotone_curve_2& cv,
                              Curve_end ind) const
    {
      return _boundary_in_x_imp (cv, ind, Base_has_boundary_category());
    }

  private:
    
    Boundary_type _boundary_in_x_imp (const X_monotone_curve_2& cv,
                                      Curve_end ind, Tag_true) const
    {
      return (m_base->boundary_in_x_2_object() (cv.base(), ind));
    }

    Boundary_type _boundary_in_x_imp (const X_monotone_curve_2& cv,
                                      Curve_end ind, Tag_false) const
    {
      return (NO_BOUNDARY);
    }
  };

  /*! Get an Boundary_in_x_2 functor object. */
  Boundary_in_x_2 boundary_in_x_2_object () const
  {
    return (Boundary_in_x_2 (m_base_traits));
  } 

  /*!
   * The Boundary_in_y_2 functor.
   */
  class Boundary_in_y_2
  {
  private:
    const Traits_2      *m_base;

  public:
       
    Boundary_in_y_2 (const Traits_2 *tr) :
      m_base (tr)
    {}
   
    Boundary_type operator() (const X_monotone_curve_2& cv,
                              Curve_end ind) const
    {
      return _boundary_in_y_imp(cv, ind, Base_has_boundary_category());
    }

  private:

    Boundary_type _boundary_in_y_imp(const X_monotone_curve_2& cv,
                                     Curve_end ind, Tag_true) const
    {
      return (m_base->boundary_in_y_2_object() (cv.base(), ind));
    }

    Boundary_type _boundary_in_y_imp(const X_monotone_curve_2& cv,
                                     Curve_end ind, Tag_false) const
    {
      return (NO_BOUNDARY);
    }
  };

  /*! Get an Boundary_in_x_2 functor object. */
  Boundary_in_y_2 boundary_in_y_2_object () const
  {
    return (Boundary_in_y_2 (m_base_traits));
  } 
};

CGAL_END_NAMESPACE

#endif
