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
// $Date$
// 
//
// Author(s)     : Ron Wein             <wein@post.tau.ac.il>
//                 Efi Fogel            <efif@post.tau.ac.il>
//                 (based on old version by Iddo Hanniel
//                                          Eyal Flato
//                                          Oren Nechushtan
//                                          Efi Fogel
//                                          Ron Wein
//                                          Idit Haran)
#ifndef CGAL_ARR_TRAITS_ADAPTOR_2_H
#define CGAL_ARR_TRAITS_ADAPTOR_2_H

/*! \file
 * Definitions of the adaptor classes for the arrangement traits class.
 */

#include <CGAL/config.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_enums.h>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A traits-class adaptor that extends the basic traits-class interface.
 */
template <class ArrangementBasicTraits_>
class Arr_traits_basic_adaptor_2 : public ArrangementBasicTraits_
{
public:

  // Traits-class geometric types.
  typedef ArrangementBasicTraits_               Base;
  typedef Arr_traits_basic_adaptor_2<Base>      Self;
  typedef typename Base::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename Base::Point_2                Point_2;

  // Tags.
  typedef typename Base::Has_left_category      Has_left_category;
  typedef typename Base::Has_boundary_category  Has_boundary_category;

  /// \name Construction.
  //@{
  /*! Default constructor. */
  Arr_traits_basic_adaptor_2 () :
    Base()
  {}

  /*! Constructor from a base-traits class. */
  Arr_traits_basic_adaptor_2 (const Base& traits) :
    Base (traits)
  {}
  //@}

  // Inherited functors:
  typedef typename Base::Compare_xy_2           Compare_xy_2;
  typedef typename Base::Construct_min_vertex_2 Construct_min_vertex_2;
  typedef typename Base::Construct_max_vertex_2 Construct_max_vertex_2;
  typedef typename Base::Is_vertical_2          Is_vertical_2;
  typedef typename Base::Compare_y_at_x_right_2 Compare_y_at_x_right_2;
  typedef typename Base::Equal_2                Equal_2;

  /// \name Overriden functors.
  //@{
  class Compare_x_2
  {
  public:

    /*!
     * Constructor
     * \param base The base traits class. It must be passed, to handle the
     *             case it is not stateless (e.g., it stores data).
     */
    Compare_x_2(const Base * base) :
      m_base(base)
    {}
    
    /*!
     * Redefining the mandatory comparison operator.
     */
    Comparison_result operator() (const Point_2& p1,
                                  const Point_2& p2) const
    {
      return (m_base->compare_x_2_object() (p1, p2));
    }

    /*!
     * Compare the relative x-positions of a vertical curve and another given
     * curve end.
     * \param p A reference point; we refer to a vertical line incident to p.
     * \param cv The compared curve.
     * \param ind MIN_END if we refer to cv's minimal end; 
     *            MAX_END if we refer to its maximal end.
     * \pre cv's relevant end has a boundary condition in y.
     * \return SMALLER if p lies to the left of cv;
     *         LARGER if p lies to the right cv;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator() (const Point_2& p,
                                  const X_monotone_curve_2& cv,
                                  Curve_end ind) const
    {
      return (_compare_point_curve_imp (p, cv, ind,
                                        Has_boundary_category()));
    }

    /*!
     * Compare the relative x-positions of two curve ends.
     * \param cv1 The first curve.
     * \param ind1 MIN_END if we refer to cv1's minimal end; 
     *             MAX_END if we refer to its maximal end.
     * \param cv2 The second curve.
     * \param ind2 MIN_END if we refer to cv2's minimal end; 
     *             MAX_END if we refer to its maximal end.
     * \pre Both curve ends have a boundary condition in y.
     * \return SMALLER if cv1 lies to the left of cv2;
     *         LARGER if cv1 lies to the right cv2;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  Curve_end ind1,
                                  const X_monotone_curve_2& cv2,
                                  Curve_end ind2) const
    {
      return (_compare_curves_imp (cv1, ind1, cv2, ind2,
                                   Has_boundary_category()));
    }

  private:

    const Base    *m_base;           // The base traits.
    
    /*!
     * Implementation of the operator() in case the Has_boundary tag is true.
     */
    Comparison_result _compare_point_curve_imp (const Point_2& p,
                                                const X_monotone_curve_2& cv,
                                                Curve_end ind,
                                                Tag_true) const
    {
      return (m_base->compare_x_2_object() (p, cv, ind));
    }

    /*!
     * Implementation of the operator() in case the Has_boundary tag is false.
     */
    Comparison_result _compare_point_curve_imp (const Point_2& ,
                                                const X_monotone_curve_2& ,
                                                Curve_end ,
                                                Tag_false) const
    {
      return (EQUAL);
    }

    /*!
     * Implementation of the operator() in case the Has_boundary tag is true.
     */
    Comparison_result _compare_curves_imp (const X_monotone_curve_2& cv1, 
                                           Curve_end ind1,
                                           const X_monotone_curve_2& cv2,
                                           Curve_end ind2,
                                           Tag_true ) const
    {
      return (m_base->compare_x_2_object() (cv1, ind1, cv2, ind2));
    }

    /*!
     * Implementation of the operator() in case the Has_boundary tag is false.
     */
    Comparison_result _compare_curves_imp (const X_monotone_curve_2& ,
                                           Curve_end ,
                                           const X_monotone_curve_2& , 
                                           Curve_end ,
                                           Tag_false ) const
    {
      return (EQUAL);
    }

  };

  /*! Get a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object () const
  {
    return Compare_x_2(this);
  }


  class Compare_y_at_x_2
  {
  public:

    /*!
     * Constructor
     * \param base The base traits class. It must be passed, to handle the
     *             case it is not stateless (e.g., it stores data)
     */
    Compare_y_at_x_2 (const Base *base) :
      m_base (base)
    {}

    /*!
     * Redefining the mandatory comparison operator.
     */
    Comparison_result operator() (const Point_2& p,
                                  const X_monotone_curve_2& cv) const
    {
      return (m_base->compare_y_at_x_2_object() (p, cv));
    }

    /*!
     * Compare the relative y-positions of two curve ends.
     * \param cv1 The first curve.
     * \param ind1 The relevant end of cv1.
     * \param cv2 The second curve.
     * \param ind2 The relevant end of cv2.
     * \pre Both curve ends have a boundary condition in x.
     * \return SMALLER if cv1 lies below cv2;
     *         LARGER if cv1 lies above cv2;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  Curve_end ind1,
                                  const X_monotone_curve_2& cv2, 
                                  Curve_end ind2) const
    {
      // The function is implemented based on the Has_boundary category.
      // If the traits class does not support unbounded curves, we just
      // return EQUAL, as this comparison will not be invoked anyway.
      return _comp_y_at_infinity_imp (cv1, ind1, cv2, ind2, 
                                      Has_boundary_category());
    }

    /*!
     * Compare the relative y-positions of two curve ends.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param ind The relevant end of cv1 and cv2.
     * \pre Both curve ends have a boundary condition in x.
     * \return SMALLER if cv1 lies below cv2;
     *         LARGER if cv1 lies above cv2;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2, 
                                  Curve_end ind) const
    {
      // The function is implemented based on the Has_boundary category.
      // If the traits class does not support unbounded curves, we just
      // return EQUAL, as this comparison will not be invoked anyway.
      return _comp_y_at_infinity_imp (cv1, cv2, ind, 
                                      Has_boundary_category());
    }

  private:


    const Base    *m_base;            // The base traits.

    /*!
     * Implementation of the operator() in case the Has_boundary tag is true.
     */
    Comparison_result _comp_y_at_infinity_imp (const X_monotone_curve_2& cv1,
                                               const X_monotone_curve_2& cv2, 
                                               Curve_end ind,
                                               Tag_true) const
    {
      return (m_base->compare_y_at_x_2_object() (cv1, cv2, ind));
    }

    /*!
     * Implementation of the operator() in case the Has_boundary tag is true.
     */
    Comparison_result _comp_y_at_infinity_imp (const X_monotone_curve_2& cv1,
                                               Curve_end ind1,
                                               const X_monotone_curve_2& cv2, 
                                               Curve_end ind2,
                                               Tag_true) const
    {
      return (m_base->compare_y_at_x_2_object() (cv1, ind1, cv2, ind2));
    }
    
    /*!
     * Implementation of the operator() in case the Has_boundary tag is false.
     */
    Comparison_result _comp_y_at_infinity_imp (const X_monotone_curve_2& ,
                                               Curve_end ,
                                               const X_monotone_curve_2& , 
                                               Curve_end ,
                                               Tag_false) const
    {
      return (EQUAL);
    }

    /*!
     * Implementation of the operator() in case the Has_boundary tag is false.
     */
    Comparison_result _comp_y_at_infinity_imp (const X_monotone_curve_2& ,
                                               const X_monotone_curve_2& , 
                                               Curve_end ,
                                               Tag_false) const
    {
      return (EQUAL);
    }
  };

  /*! Get a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return Compare_y_at_x_2(this);
  }


  class Compare_y_at_x_left_2
  {
  public:
    
    /*!
     * Constructor
     * \param self The traits class. It must be passed, to handle the case
     *             it is not stateless (e.g., it stores data).
     */
    Compare_y_at_x_left_2 (const Self *self) :
      m_self (self)
    {}

    /*!
     * Compare two curves immediately to the left of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The query point.
     * \pre The two curves intersect at p, and they are defined to its left.
     * \return SMALLER if cv1 lies below cv2 to the left of q;
     *         LARGER if cv1 lies above cv2 to the left of q;
     *         EQUAL in case of an overlap to the left of q.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& p) const 
    {
      // The function is implemented based on the Has_left category. If the 
      // category indicates that the "left" version is available, it calls the
      // function with same name defined in the base class. Otherwise, it 
      // uses other predicates to provide this comparison.
      return _compare_y_at_x_left_imp (cv1, cv2, p, Has_left_category());
    }

  private:

    const Self    *m_self;       // The traits.

    /*!
     * Implementation of the operator() in case the HasLeft tag is true.
     */
    Comparison_result _compare_y_at_x_left_imp (const X_monotone_curve_2& cv1,
                                                const X_monotone_curve_2& cv2,
                                                const Point_2& p,
                                                Tag_true) const
    {
      const Base   *base = m_self;
      return (base->compare_y_at_x_left_2_object() (cv1, cv2, p));
    }

    /*!
     * Implementation of the operator() in case the HasLeft tag is false.
     */
    Comparison_result _compare_y_at_x_left_imp (const X_monotone_curve_2& cv1,
                                                const X_monotone_curve_2& cv2,
                                                const Point_2& p,
                                                Tag_false) const
    {
      Boundary_in_x_2         boundary_x = m_self->boundary_in_x_2_object();
      Boundary_in_y_2         boundary_y = m_self->boundary_in_y_2_object();
      Construct_min_vertex_2  min_vertex =
                                    m_self->construct_min_vertex_2_object();
      Equal_2                 equal = m_self->equal_2_object();

      // Check if the left ends of the curves are bounded endpoints.
      const Boundary_type  bx1 = boundary_x (cv1, MIN_END);
      const Boundary_type  by1 = (bx1 != NO_BOUNDARY ? 
                                  NO_BOUNDARY : boundary_y (cv1, MIN_END));
      const bool           has_left1 = (bx1 == NO_BOUNDARY &&
                                        by1 == NO_BOUNDARY);

      const Boundary_type  bx2 = boundary_x (cv2, MIN_END);
      const Boundary_type  by2 = (bx2 != NO_BOUNDARY ? 
                                  NO_BOUNDARY : boundary_y (cv2, MIN_END));
      const bool           has_left2 = (bx2 == NO_BOUNDARY &&
                                        by2 == NO_BOUNDARY);

      CGAL_assertion (CGAL::sign (bx1) != POSITIVE &&
                      CGAL::sign (bx2) != POSITIVE);

      // Make sure that p lies on both curves, and that both are defined to its
      // right (so their right endpoint is lexicographically larger than p).
      CGAL_precondition_code (
        Compare_xy_2       compare_xy = m_self->compare_xy_2_object();
        Compare_y_at_x_2   compare_y_at_x = m_self->compare_y_at_x_2_object();
      );

      CGAL_precondition (compare_y_at_x (p, cv1) == EQUAL &&
                         compare_y_at_x (p, cv2) == EQUAL);

      CGAL_precondition ((! has_left1 ||
                          compare_xy(min_vertex (cv1), p) == SMALLER) &&
                         (! has_left2 ||
                          compare_xy(min_vertex (cv2), p) == SMALLER));

      // If one of the curves is vertical, it is below the other one.
      Is_vertical_2       is_vertical = m_self->is_vertical_2_object();

      if (is_vertical(cv1))
      {
        if (is_vertical (cv2))
          // Both are vertical:
          return (EQUAL);
        else
          return (SMALLER);
      }
      else if (is_vertical (cv2))
      {
        return (LARGER);
      }

      // Perform the comparison based on the existance of bounded left
      // endpoints.
      if (has_left1 && has_left2)
      {
        // Get the left endpoints of cv1 and cv2.
        Point_2        left1 = min_vertex(cv1);
        Point_2        left2 = min_vertex(cv2);
          
        if (equal (left1, left2))
        {
          // The two curves have a common left endpoint:
          // Compare them to the right of this point.
          return (m_self->compare_y_at_x_right_2_object() (cv1, cv2, left1));
        }    
      }

      // We now that the curves do not share a common endpoint, and we can
      // compare their relative y-position (which does not change to the left
      // of the given point p).
      return (m_self->compare_y_position_2_object() (cv1, cv2));      
    }
  };

  /*! Get a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object () const
  {
    return Compare_y_at_x_left_2(this);
  }


  class Boundary_in_x_2
  {
  public:

    /*!
     * Constructor.
     * \param self The traits class. It must be passed, to handle the case
     *             it is not stateless (e.g., it stores data).
     */
    Boundary_in_x_2 (const Base * base) :
      m_base (base)
    {}

    /*!
     * Get the boundary condition of the given curve end in x.
     * \param cv The curve.
     * \param ind MIN_END if we refer to cv's minimal end,
     *            MAX_END if we refer to its maximal end.
     * \return The boundary condition of the curve end.
     */
    Boundary_type operator() (const X_monotone_curve_2& cv,
                              Curve_end ind) const
    {
      // The function is implemented based on the Has_boundary category.
      // If the traits class does not support boundary conditions, we just
      // return NO_BOUNDARY.
      return _boundary_in_x_imp (cv, ind, Has_boundary_category());
    }

  private:

    const Base    *m_base;            // The base traits.

    /*!
     * Implementation of the operator() in case the Has_boundary tag is true.
     */
    Boundary_type _boundary_in_x_imp (const X_monotone_curve_2& cv,
                                      Curve_end ind,
                                      Tag_true) const
    {
      return (m_base->boundary_in_x_2_object() (cv, ind));
    }

    /*!
     * Implementation of the operator() in case the Has_boundary tag is false.
     */
    Boundary_type _boundary_in_x_imp (const X_monotone_curve_2& , Curve_end ,
                                      Tag_false) const
    {
      return (NO_BOUNDARY);
    }
  };

  /*! Get an Boundary_in_x_2 functor object. */
  Boundary_in_x_2 boundary_in_x_2_object () const
  {
    return Boundary_in_x_2(this);
  }


  class Boundary_in_y_2
  {
  public:
    
    /*!
     * Constructor.
     * \param self The traits class. It must be passed, to handle the case
     *             it is not stateless (e.g., it stores data).
     */
    Boundary_in_y_2 (const Base * base) :
      m_base (base)
    {}

    /*!
     * Get the boundary condition of the given curve end in y.
     * \param cv The curve.
     * \param ind MIN_END if we refer to cv's minimal end,
     *            MAX_END if we refer to its maximal end.
     * \return The boundary condition of the curve end.
     */
    Boundary_type operator() (const X_monotone_curve_2& cv,
                              Curve_end ind) const
    {
      // The function is implemented based on the Has_boundary category.
      // If the traits class does not support boundary conditions, we just
      // return NO_BOUNDARY.
      return _boundary_in_y_imp (cv, ind, Has_boundary_category());
    }

  private:

    const Base    *m_base;          // The base traits.

    /*!
     * Implementation of the operator() in case the Has_boundary tag is true.
     */
    Boundary_type _boundary_in_y_imp (const X_monotone_curve_2& cv,
                                      Curve_end ind,
                                      Tag_true) const
    {
      return (m_base->boundary_in_y_2_object() (cv, ind));
    }

    /*!
     * Implementation of the operator() in case the Has_boundary tag is false.
     */
    Boundary_type _boundary_in_y_imp (const X_monotone_curve_2& , Curve_end ,
                                      Tag_false) const
    {
      return (NO_BOUNDARY);
    }
  };

  /*! Get an Boundary_in_y_2 functor object. */
  Boundary_in_y_2 boundary_in_y_2_object () const
  {
    return Boundary_in_y_2(this);
  }
  //@}

  /// \name Additional auxiliary functors.
  //@{
  class Is_in_x_range_2
  {
  public:

    /*!
     * Constructor.
     * \param self The traits class. It must be passed, to handle the case
     *             it is not stateless (e.g., it stores data)
     */
    Is_in_x_range_2(const Self * self) :
      m_self (self)
    {}

    /*!
     * Check whether the given point is in the x-range of the given x-monotone
     * curve.
     * \param cv The x-monotone curve.
     * \param p The point.
     * \return true if x(cv_left) <= x(p) <= x(cv_right), false otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv, const Point_2& p) const
    {
      Boundary_in_x_2     boundary_x = m_self->boundary_in_x_2_object();
      Boundary_in_y_2     boundary_y = m_self->boundary_in_y_2_object();
      Compare_x_2         compare_x =  m_self->compare_x_2_object();

      // Compare p to the position of the left end of the curve.
      // Note that if the left end of cv lies at x boundary, p is obviously to
      // its right.
      Boundary_type       bx, by;
      Comparison_result   res;

      bx = boundary_x (cv, MIN_END);

      if (bx == NO_BOUNDARY)
      {
        by = boundary_y (cv, MIN_END);

        if (by == NO_BOUNDARY)
        {
          // The left endpoint of cv is a normal point.
          res = compare_x (p, m_self->construct_min_vertex_2_object() (cv));
        }
        else
        {
          // The left end of cv lies at y boundary.
          res = compare_x (p, cv, MIN_END);
        }

        if (res == SMALLER)
          return (false);         // p is to the left of the x-range.
        else if (res == EQUAL)
          return (true);
      }

      // If necessary, compare p to the right end of the curve.
      // Note that if this end lies at x boundary, p is obviously to its left.
      bx = boundary_x (cv, MAX_END);

      if (bx != NO_BOUNDARY)
        return (true);
      
      by = boundary_y (cv, MAX_END);

      if (by == NO_BOUNDARY)
      {
        // The right endpoint of cv is a normal point.
        res = compare_x (p, m_self->construct_max_vertex_2_object() (cv));
      }
      else
      {
        // The right end of cv lies at y boundary:
        res = compare_x (p, cv, MAX_END);
      }

      return (res != LARGER);
    }

    /*!
     * Check whether the x-ranges of the given x-monotone curves overlap.
     * \param cv1 The first x-monotone curve.
     * \param cv2 The second x-monotone curve.
     * \return (true) if there is an overlap in the x-ranges of the given
     *         curves; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2) const
    {
      Boundary_in_x_2         boundary_x = m_self->boundary_in_x_2_object();
      Boundary_in_y_2         boundary_y = m_self->boundary_in_y_2_object();
      Compare_x_2             compare_x = m_self->compare_x_2_object();
      Construct_min_vertex_2  min_vertex =
                                    m_self->construct_min_vertex_2_object();
      Construct_max_vertex_2  max_vertex =
                                    m_self->construct_max_vertex_2_object();

      // Locate the rightmost of the two left endpoints of the two curves.
      // Note that we guard for curve ends with boundary conditions.
      Boundary_type             bx1, by1;
      Boundary_type             bx2, by2;
      const X_monotone_curve_2 *cv_left;
      Boundary_type             by_left;
      Comparison_result         res;

      bx1 = boundary_x (cv1, MIN_END);
      bx2 = boundary_x (cv2, MIN_END);

      if (bx1 != NO_BOUNDARY)
      {
        // If both curves are defined at x boundary, they obviously overlap in
        // their x-ranges.
        if (bx2 != NO_BOUNDARY)
          return (true);

        // As cv2 is not defined at x boundary, take its left end as the
        // rightmost of the two left curve ends.
        cv_left = &cv2;
        by_left = boundary_y (cv2, MIN_END);
      }
      else if (bx2 != NO_BOUNDARY)
      {
        // As cv1 is not defined at x boundary, take its left end as the
        // rightmost of the two left curve ends.
        cv_left = &cv1;
        by_left = boundary_y (cv1, MIN_END);
      }
      else
      {
        // Compare the (finite) x-coordinates of the two left ends.
        // We take special care of the case of boundaries in y.
        by1 = boundary_y (cv1, MIN_END);
        by2 = boundary_y (cv2, MIN_END);

        if (by1 == NO_BOUNDARY)
        {
          if (by2 == NO_BOUNDARY)
          {
            res = compare_x (min_vertex (cv1), min_vertex (cv2));
          }
          else
          {
            res = compare_x (min_vertex (cv1), cv2, MIN_END);
          }
        }
        else
        {
          if (by2 == NO_BOUNDARY)
          {
            res = compare_x (min_vertex (cv2), cv1, MIN_END);
            if (res != EQUAL)
              res = (res == SMALLER) ? LARGER : SMALLER;
          }
          else
          {
            res = compare_x (cv1, MIN_END, cv2, MIN_END);
          }
        }

        if (res == LARGER)
        {
          cv_left = &cv1;
          by_left = by1;
        }
        else
        {
          cv_left = &cv2;
          by_left = by2;
        }
      }
      
      // Locate the leftmost of the two right endpoints of the two curves.
      // Note that we guard for curve ends with boundary conditions.
      const X_monotone_curve_2 *cv_right;
      Boundary_type             by_right;

      bx1 = boundary_x (cv1, MAX_END);
      bx2 = boundary_x (cv2, MAX_END);

      if (bx1 != NO_BOUNDARY)
      {
        // If both curves are defined at x boundary, they obviously overlap in
        // their x-ranges.
        if (bx2 != NO_BOUNDARY)
          return (true);

        // As cv2 is not defined at x boundary, take its right end as the
        // leftmost of the two right curve ends.
        cv_right = &cv2;
        by_right = boundary_y (cv2, MAX_END);
      }
      else if (bx2 != NO_BOUNDARY)
      {
        // As cv1 is not defined at x boundary, take its right end as the
        // leftmost of the two right curve ends.
        cv_right = &cv1;
        by_right = boundary_y (cv1, MAX_END);
      }
      else
      {
        // Compare the (finite) x-coordinates of the two right ends.
        // We take special care of the case of boundaries in y.
        by1 = boundary_y (cv1, MAX_END);
        by2 = boundary_y (cv2, MAX_END);

        if (by1 == NO_BOUNDARY)
        {
          if (by2 == NO_BOUNDARY)
          {
            res = compare_x (max_vertex (cv1), max_vertex (cv2));
          }
          else
          {
            res = compare_x (max_vertex (cv1), cv2, MAX_END);
           }
        }
        else
        {
          if (by2 == NO_BOUNDARY)
          {
            res = compare_x (max_vertex (cv2), cv1, MAX_END);
            if (res != EQUAL)
              res = (res == SMALLER) ? LARGER : SMALLER;
          }
          else
          {
            res = compare_x (cv1, MAX_END, cv2, MAX_END);
          }
        }

        if (res == SMALLER)
        {
          cv_right = &cv1;
          by_right = by1;
        }
        else
        {
          cv_right = &cv2;
          by_right = by2;
        }
      }
          
      // Now compare the (finite) x-coordiates of the left end of cv_left and
      // the right end of cv_right.
      if (by_left == NO_BOUNDARY)
      {
        if (by_right == NO_BOUNDARY)
        {
          res = compare_x (min_vertex (*cv_left), max_vertex (*cv_right));
        }
        else
        {
          res = compare_x (min_vertex (*cv_left), *cv_right, MAX_END);
        }
      }
      else
      {
        if (by_right == NO_BOUNDARY)
        {
          res = compare_x (max_vertex (*cv_right), *cv_left, MIN_END);
          if (res != EQUAL)
            res = (res == SMALLER) ? LARGER : SMALLER;
        }
        else
        {
          res = compare_x (*cv_left, MIN_END, *cv_right, MAX_END);
        }
      }

      // The two curves overlap in their x-range if and only if the left end
      // of cv_left is not to the right if the right end of cv_right.
      return (res != LARGER);
    }

  private:

    const Self    *m_self;         // The traits.
  };

  /*! Get an Is_in_x_range_2 functor object. */
  Is_in_x_range_2 is_in_x_range_2_object () const
  {
    return Is_in_x_range_2(this);
  }


  class Compare_y_position_2
  {
  public:

    /*!
     * Constructor
     * \param self The traits class. It must be passed, to handle the case
     *             it is not stateless (e.g., it stores data).
     */
    Compare_y_position_2 (const Self *self) :
      m_self (self)
    {}

    /*!
     * Get the relative of two x-monotone curves with overlapping x-ranges
     * that are disjoint in their interiors.
     * \param cv1 The first x-monotone curve.
     * \param cv2 The second x-monotone curve.
     * \pre The x-ranges of the two curves overlap.
     * \return SMALLER if cv1 lies below cv2;
     *         LARGER if cv1 lies above cv2;
     *         EQUAL in case the common x-range is a single point.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2) const
    {
      CGAL_precondition_code (
        Is_in_x_range_2  is_in_x_range = m_self->is_in_x_range_2_object();
      );
      CGAL_precondition (is_in_x_range (cv1, cv2));

      Boundary_in_x_2         boundary_x = m_self->boundary_in_x_2_object();
      Boundary_in_y_2         boundary_y = m_self->boundary_in_y_2_object();
      Compare_y_at_x_2        compare_y_at_x =
                                          m_self->compare_y_at_x_2_object();
      Construct_min_vertex_2  min_vertex =
                                    m_self->construct_min_vertex_2_object();
 
      // First check whether any of the curves is defined at x boundary.
      const Boundary_type     bx1 = boundary_x (cv1, MIN_END);
      const Boundary_type     bx2 = boundary_x (cv2, MIN_END);
      Comparison_result       res;

      CGAL_assertion (CGAL::sign (bx1) != POSITIVE &&
                      CGAL::sign (bx2) != POSITIVE);

      if (bx1 != NO_BOUNDARY)
      {
        if (bx2 != NO_BOUNDARY)
        {
          // Compare the relative position of the curves at x boundary.
          return (compare_y_at_x (cv1, cv2, MIN_END));
        }
        
        // Check if the left end of cv2 lies at y boundary.
        const Boundary_type    by2 = boundary_y (cv2, MIN_END);
        const CGAL::Sign       sign_by2 = CGAL::sign (by2);

        if (sign_by2 == NEGATIVE)
          return (LARGER);          // cv2 is obviously below cv1.
        else if (sign_by2 == POSITIVE)
          return (SMALLER);         // cv2 is obviously above cv1.
          
        // Compare the position of the left end of cv2 (which is a normal
        // point) to cv1.
        res = compare_y_at_x (min_vertex (cv2), cv1);

        // Swap the result.
        if (res == EQUAL)
          return (EQUAL);

        return ((res == SMALLER) ? LARGER : SMALLER);
      }
      else if (bx2 != NO_BOUNDARY)
      {
        // Check if the left end of cv1 lies at y boundary.
        const Boundary_type    by1 = boundary_y (cv1, MIN_END);
        const CGAL::Sign       sign_by1 = CGAL::sign (by1);

        if (sign_by1 == NEGATIVE)
          return (SMALLER);         // cv1 is obviously below cv2.
        else if (sign_by1 == POSITIVE)
          return (LARGER);          // cv1 is obviously above cv2.

        // Compare the position of the left end of cv1 (which is a normal
        // point) to cv2.
        res = compare_y_at_x (min_vertex (cv1), cv2);

        return (res);
      }
      
      // Check if the left curve end lies at y = +/- oo.
      const Boundary_type     by1 = boundary_y (cv1, MIN_END);
      const CGAL::Sign        sign_by1 = CGAL::sign (by1);
      const Boundary_type     by2 = boundary_y (cv2, MIN_END);
      const CGAL::Sign        sign_by2 = CGAL::sign (by2);
      Comparison_result       l_res;

      if (by1 != NO_BOUNDARY)
      {
        if (by2 != NO_BOUNDARY)
        {
          // The curve ends have boundary conditions with oposite signs in y,
          // we readily know their relative position (recall that they do not
          // instersect).
          if (sign_by1 == NEGATIVE && sign_by2 == POSITIVE)
            return (SMALLER);
          else if (sign_by1 == POSITIVE && sign_by2 == NEGATIVE)
            return (LARGER);

          // Both curves have vertical asymptotes with the same sign in y.
          // Check which asymptote is the rightmost. Note that in this case
          // the vertical asymptotes cannot be equal.
          l_res = m_self->compare_x_2_object() (cv1, MIN_END,
                                                cv2, MIN_END);
          CGAL_assertion (l_res != EQUAL);

          if (sign_by1 == POSITIVE)
            return (l_res);
          else
            return ((l_res == SMALLER) ? LARGER : SMALLER);
        }

        // cv1 has a vertical asymptote and cv2 has a normal left endpoint.
        // Compare the x-positions of this endpoint and the asymptote.
        const Point_2&  left2 = min_vertex(cv2);
        
        l_res = m_self->compare_x_2_object() (left2, cv1, MIN_END);

        if (l_res == LARGER)
        {
          // left2 lies in the x-range of cv1, so it is safe to compare:
          res = compare_y_at_x (left2, cv1);

          // Swap the result.
          if (res == EQUAL)
            return (EQUAL);

          return ((res == SMALLER) ? LARGER : SMALLER);
        }
        else
        {
          if (sign_by1 == NEGATIVE)
            return (SMALLER);          // cv1 is obviously below cv2.
          else
            return (LARGER);           // cv2 is obviously above cv1.
        }
      }
      else if (by2 != NO_BOUNDARY)
      {
        // cv2 has a vertical asymptote and cv1 has a normal left endpoint.
        // Compare the x-positions of this endpoint and the asymptote.
        const Point_2&  left1 = min_vertex(cv1);
        
        l_res = m_self->compare_x_2_object() (left1, cv2, MIN_END);

        if (l_res == LARGER)
        {
          // left1 lies in the x-range of cv2, so it is safe to compare:
          return (compare_y_at_x (left1, cv2));
        }
        else
        {
          if (sign_by2 == NEGATIVE)
            return (LARGER);           // cv2 is obviously below cv1.
          else
            return (SMALLER);          // cv1 is obviously above cv2.
        }
      }

      // In this case we compare two normal points.
      Compare_xy_2            compare_xy = m_self->compare_xy_2_object();
      Compare_y_at_x_right_2  compare_y_at_x_right =
                                 m_self->compare_y_at_x_right_2_object();

      // Get the left endpoints of cv1 and cv2.
      const Point_2&  left1 = min_vertex(cv1);
      const Point_2&  left2 = min_vertex(cv2);

      // Locate the rightmost point of left1 and left2 and compare its position
      // to the other curve.
      l_res = compare_xy (left1, left2);

      if (l_res != SMALLER)
      {
        // left1 is in the x-range of cv2:
        res = compare_y_at_x (left1, cv2);

        if (res == EQUAL)
        {
          // The two curves intersect at left1. If both curves are defined to
          // the right of the reference point, we can compare them to its 
          // right. Otherwise, their share a common endpoint (which is the only
          // overlap in their x-ranges) and are really equal.
          if (l_res == EQUAL)
            res = compare_y_at_x_right (cv1, cv2, left1);
        }

        return (res);
      }
      else
      {
        // left2 is in the x-range of cv1:
        res = compare_y_at_x (left2, cv1);

        if (res == EQUAL)
        {
          // The two curves share a common endpoint (which is the only overlap
          // in their x-ranges) and are really equal.
          return (EQUAL);
        }

        // Swap the result:
        return ((res == SMALLER) ? LARGER : SMALLER);
      }
    }

  private:

    const Self    *m_self;            // The traits.
  };

  /*! Get a Compare_y_position_2 functor object. */
  Compare_y_position_2 compare_y_position_2_object () const
  {
    return Compare_y_position_2(this);
  }


  class Is_between_cw_2
  {
  public:

    /*!
     * Constructor.
     * \param self The traits class. It must be passed, to handle the case
     *             it is not stateless (e.g., it stores data).
     */
    Is_between_cw_2 (const Self *self) :
      m_self (self)
    {}

    /*!
     * Check whether the given query curve is encountered when rotating the
     * first curve in a clockwise direction around a given point until reaching
     * the second curve.
     * \param cv The query curve.
     * \param cv_to_right Is cv directed from left to right (that is, the
     *                    common vertex is cv's left endpoint).
     * \param cv1 The first curve.
     * \param cv1_to_right Is cv1 directed from left to right.
     * \param cv2 The second curve.
     * \param cv2_to_right Is cv2 directed from left to right.
     * \param p The point around which we rotate cv1.
     * \param cv_equal_cv1 Output: does cv equal cv1.
     * \param cv_equal_cv2 Output: does cv equal cv2.
     * \pre p is an end-point of all three curves.
     * \return (true) if cv is between cv1 and cv2; (false) otherwise.
     *         If cv overlaps cv1 or cv2 the result is always (false).
     *         If cv1 and cv2 overlap, the result is (true), unless cv 
     *         also overlaps them.
     */
    bool operator() (const X_monotone_curve_2& cv, bool cv_to_right,
                     const X_monotone_curve_2& cv1, bool cv1_to_right,
                     const X_monotone_curve_2& cv2, bool cv2_to_right,
                     const Point_2& p,
                     bool& cv_equal_cv1, 
                     bool& cv_equal_cv2) const
    {
      Compare_y_at_x_left_2   compare_y_at_x_left =
                                     m_self->compare_y_at_x_left_2_object();
      Compare_y_at_x_right_2  compare_y_at_x_right =
                                    m_self->compare_y_at_x_right_2_object();

      // Initialize output flags.
      cv_equal_cv1 = false;
      cv_equal_cv2 = false;
    
      // Take care of the general 4 cases:
      Comparison_result  l_res, r_res;
      Comparison_result  res1, res2;
    
      if (!cv1_to_right && !cv2_to_right)
      {
        // Case 1: Both cv1 and cv2 are defined to the left of p.
        l_res = compare_y_at_x_left (cv1, cv2, p);
      
        if (l_res == LARGER)
        {
          // Case 1(a) : cv1 is above cv2.
          if (!cv_to_right)
          {
            res1 = compare_y_at_x_left (cv1, cv, p);
            res2 = compare_y_at_x_left (cv2, cv, p);
          
            if (res1 == EQUAL)
              cv_equal_cv1 = true;
            if (res2 == EQUAL)
              cv_equal_cv2 = true;
          
            return (res1 == SMALLER || res2 == LARGER);
          }
          return (true);
        }
        else if (l_res == SMALLER)
        {
          // Case 1(b): cv1 is below cv2.
          if (!cv_to_right)
          {
            res1 = compare_y_at_x_left (cv1, cv, p);
            res2 = compare_y_at_x_left (cv2, cv, p);
          
            if (res1 == EQUAL)
              cv_equal_cv1 = true;
            if (res2 == EQUAL)
              cv_equal_cv2 = true;
          
            return (res1 == SMALLER && res2  == LARGER);
          }
          return (false);
        }
        else
        {
          // Overlapping segments.
          if (!cv_to_right)
          {
            res1 = compare_y_at_x_left (cv1, cv, p);
            if (res1 == EQUAL)
            {
              cv_equal_cv1 = true;
              cv_equal_cv2 = true;
              return (false);
            }
            return (true);
          }
          return (true);
        }
      }
      
      if (cv1_to_right && cv2_to_right)
      {
        // Case 2: Both cv1 and cv2 are defined to the right of p.
        r_res = compare_y_at_x_right (cv1, cv2, p);

        if (r_res == LARGER)
        {
          // Case 2(a) : cv1 is above cv2.
          if (cv_to_right)
          {
            res1 = compare_y_at_x_right (cv1, cv, p);
            res2 = compare_y_at_x_right (cv2, cv, p);

            if (res1 == EQUAL)
              cv_equal_cv1 = true;
            if (res2 == EQUAL)
              cv_equal_cv2 = true;

            return (res1 == LARGER && res2 == SMALLER);
          }
          return (false);
        }
        else if (r_res == SMALLER)
        {
          // Case 2(b): cv1 is below cv2.
          if (cv_to_right)
          {
            res1 = compare_y_at_x_right (cv1, cv, p);
            res2 = compare_y_at_x_right (cv2, cv, p);

            if (res1 == EQUAL)
              cv_equal_cv1 = true;
            if (res2 == EQUAL)
              cv_equal_cv2 = true;

            return (res1 == LARGER || res2 == SMALLER);
          }
          return (true);
        }
        else
        {
          // Overlapping segments.
          if (cv_to_right)
          {
            res1 = compare_y_at_x_right (cv1, cv, p);
          
            if (res1 == EQUAL)
            {
              cv_equal_cv1 = true;
              cv_equal_cv2 = true;             
              return (false);
            }
            return (true);
          }
          return (true);
        }
      }

      if (!cv1_to_right && cv2_to_right)
      {
        // Case 3: cv1 is defined to the left of p, and cv2 to its right.
        if (!cv_to_right)
        {
          res1 = compare_y_at_x_left (cv1, cv, p);

          if (res1 == EQUAL)
            cv_equal_cv1 = true;
        
          return (res1 == SMALLER);
        }
        else
        {
          res2 = compare_y_at_x_right (cv2, cv, p);

          if (res2 == EQUAL)
            cv_equal_cv2 = true;

          return (res2 == SMALLER);
        }
      }

      CGAL_assertion (cv1_to_right && !cv2_to_right);

      // Case 4: cv1 is defined to the right of p, and cv2 to its left.
      if (cv_to_right)
      {
	res1 = compare_y_at_x_right (cv1, cv, p);
	
	if (res1 == EQUAL)
	  cv_equal_cv1 = true;
        
	return (res1  == LARGER);
      }
      else
      {
        res2 = compare_y_at_x_left (cv2, cv, p);
        
        if (res2 == EQUAL)
          cv_equal_cv2 = true;
	
        return (res2 == LARGER);
      }
    }
    
  private:

    const Self    *m_self;        // The traits.
  };

  /*! Get an Is_between_cw_2 functor object. */
  Is_between_cw_2 is_between_cw_2_object () const
  {
    return Is_between_cw_2(this);
  }


  class Compare_cw_around_point_2
  {
  public:

    /*!
     * Constructor.
     * \param self The traits class. It must be passed, to handle the case
     *             it is not stateless (e.g., it stores data)
     */
    Compare_cw_around_point_2 (const Self *self) :
      m_self (self)
    {}
    
    /*!
     * Compare the two interior disjoint x-monotone curves in a clockwise
     * order around their common endpoint.
     * \param cv1 The first curve.
     * \param cv1_to_right Is cv1 directed from left to right.
     * \param cv2 The second curve.
     * \param cv2_to_right Is cv2 directed from left to right.
     * \param p The common endpoint.
     * \param from_top (true) if we start from 12 o'clock, 
     *                 (false) if we start from 6 o'clock.
     * \pre The point p is an endpoint of both curves.
     * \return SMALLER if we encounter cv1 before cv2;
     *         LARGER if we encounter cv2 before cv1;
     *         EQUAL otherwise.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  bool cv1_to_right,
                                  const X_monotone_curve_2& cv2,
                                  bool cv2_to_right,
                                  const Point_2& p,
                                  bool from_top = true) const
    {
      // Act according to where cv1 and cv2 lie.
      if (!cv1_to_right && !cv2_to_right)
      {
        // Both are defined to the left of p, and we encounter cv1 before
        // cv2 if it is below cv2:
        return (m_self->compare_y_at_x_left_2_object() (cv1, cv2, p));
      }
      
      if (cv1_to_right && cv2_to_right)
      {
        // Both are defined to the right of p, and we encounter cv1 before
        // cv2 if it is above cv2. We therefore reverse the order of the
        // curves when we invoke compare_y_at_x_right:
        return (m_self->compare_y_at_x_right_2_object() (cv2, cv1, p));
      }
      
      if (!cv1_to_right && cv2_to_right)
      {
        // If we start from the top, we encounter the right curve (which
        // is cv2) first. If we start from the bottom, we encounter cv1 first.
        return (from_top ? LARGER : SMALLER);
      }

      CGAL_assertion (cv1_to_right && !cv2_to_right);

      // If we start from the top, we encounter the right curve (which
      // is cv1) first. If we start from the bottom, we encounter cv2 first.
      return (from_top ? SMALLER : LARGER);
    }

  private:
    
    const Self    *m_self;       // The traits.
  };

  /*! Get a Compare_cw_around_point_2 functor object. */
  Compare_cw_around_point_2 compare_cw_around_point_2_object () const
  {
    return Compare_cw_around_point_2(this);
  }
  //@}
};

/*! \class
 * A traits-class adaptor that extends the basic traits-class interface.
 */
template <class ArrangementTraits_>
class Arr_traits_adaptor_2 :
  public Arr_traits_basic_adaptor_2<ArrangementTraits_>
{
public:

  // Traits-class geometric types.
  typedef ArrangementTraits_                             Base_traits_2;
  typedef Arr_traits_basic_adaptor_2<ArrangementTraits_> Base;

  typedef typename Base_traits_2::Curve_2                Curve_2;
  typedef typename Base::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Base::Point_2                         Point_2;

  // Tags.
  typedef typename Base::Has_left_category               Has_left_category;
  typedef typename Base::Has_boundary_category           Has_boundary_category;
  typedef typename Base::Has_merge_category              Has_merge_category;

  /// \name Construction.
  //@{
  /*! Default constructor. */
  Arr_traits_adaptor_2 () :
    Base()
  {}

  /*! Constructor from a base-traits class. */
  Arr_traits_adaptor_2 (const Base_traits_2& traits) :
    Base (traits)
  {}
  //@}

  // Inherited functors:
  typedef typename Base::Compare_x_2            Compare_x_2;
  typedef typename Base::Compare_xy_2           Compare_xy_2;
  typedef typename Base::Construct_min_vertex_2 Construct_min_vertex_2;
  typedef typename Base::Construct_max_vertex_2 Construct_max_vertex_2;
  typedef typename Base::Is_vertical_2          Is_vertical_2;
  typedef typename Base::Compare_y_at_x_2       Compare_y_at_x_2;
  typedef typename Base::Compare_y_at_x_right_2 Compare_y_at_x_right_2;
  typedef typename Base::Compare_y_at_x_left_2  Compare_y_at_x_left_2;
  typedef typename Base::Equal_2                Equal_2;

  // Note that the basic adaptor does not have to support these functors:
  typedef typename Base_traits_2::Make_x_monotone_2  Make_x_monotone_2;
  typedef typename Base_traits_2::Split_2            Split_2; 
  typedef typename Base_traits_2::Intersect_2        Intersect_2;

  /// \name Overriden functors.
  //@{

  class Are_mergeable_2
  {
  public:

    /*!
     * Constructor.
     * \param base The base traits class. It must be passed, to handle the
     *             case it is not stateless (e.g., it stores data).
     */
    Are_mergeable_2 (const Base *base) :
      m_base (base)
    {}

    /*!
     * Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are mergeable - if they are supported
     *         by the same line and share a common endpoint; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2) const
    {
      // The function is implemented based on the Has_merge category.
      return (_are_mergeable_imp (cv1, cv2, Has_merge_category()));
    }

  private:
    
    const Base    *m_base;         // The base traits.

    /*!
     * Implementation of the operator() in case the Has_merge tag is true.
     */
    bool _are_mergeable_imp (const X_monotone_curve_2& cv1,
           const X_monotone_curve_2& cv2,
           Tag_true) const
    {
      return (m_base->are_mergeable_2_object() (cv1, cv2));      
    }

    /*!
     * Implementation of the operator() in case the Has_merge tag is false.
     */
    bool _are_mergeable_imp (const X_monotone_curve_2& ,
           const X_monotone_curve_2& ,
           Tag_false) const
    {
      // Curve merging is not supported:
      return (false);
    }
  };

  /*! Get an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object () const
  {
    return Are_mergeable_2(this);
  }


  class Merge_2
  {
  public:

    /*!
     * Constructor.
     * \param base The base traits class. It must be passed, to handle the
     *             case it is not stateless (e.g., it stores data).
     */
    Merge_2 (const Base *base) :
      m_base (base)
    {}

    /*!
     * Merge two given x-monotone curves into a single curve.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable, that is they are supported by the
     *      curve line and share a common endpoint.
     */
    void operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2,
                     X_monotone_curve_2& c) const
    {
      // The function is implemented based on the Has_merge category.
      _merge_imp (cv1, cv2, c, Has_merge_category());
    }

  private:


    const Base    *m_base;            // The base traits.

    /*!
     * Implementation of the operator() in case the HasMerge tag is true.
     */
    void _merge_imp (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2,
                     X_monotone_curve_2& c,
                     Tag_true) const
    {
      return (m_base->merge_2_object() (cv1, cv2, c));      
    }

    /*!
     * Implementation of the operator() in case the HasMerge tag is false.
     */
    void _merge_imp (const X_monotone_curve_2& ,
                     const X_monotone_curve_2& ,
                     X_monotone_curve_2& ,
                     Tag_false) const
    {
      // This function should never be called!
      CGAL_assertion_msg (false,
                          "Merging curves is not supported.");
      return;
    }
  };

  /*! Get a Merge_2 functor object. */
  Merge_2 merge_2_object () const
  {
    return Merge_2(this);
  }
  //@}

};

CGAL_END_NAMESPACE

#endif
