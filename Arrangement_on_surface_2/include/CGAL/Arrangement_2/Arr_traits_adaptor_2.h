// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// $Date$
//
// Author(s): Ron Wein             <wein@post.tau.ac.il>s
//            Efi Fogel            <efif@post.tau.ac.il>
//            Eric Berberich       <eric@mpi-inf.mpg.de>
//            (based on old version by Iddo Hanniel
//                                     Eyal Flato
//                                     Oren Nechushtan
//                                     Efi Fogel
//                                     Ron Wein
//                                     Idit Haran)

#ifndef CGAL_ARR_TRAITS_ADAPTOR_2_H
#define CGAL_ARR_TRAITS_ADAPTOR_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 * Definitions of the adaptor classes for the arrangement traits class.
 */

#include <list>

#include <CGAL/config.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2_dispatching.h>

namespace CGAL {

/*! \class
 * A traits-class adaptor that extends the basic traits-class interface.
 */
template <typename ArrangementBasicTraits_>
class Arr_traits_basic_adaptor_2 : public ArrangementBasicTraits_ {
public:
  // Traits-class geometric types.
  typedef ArrangementBasicTraits_                   Base;
  typedef Arr_traits_basic_adaptor_2<Base>          Self;
  typedef typename Base::X_monotone_curve_2         X_monotone_curve_2;
  typedef typename Base::Point_2                    Point_2;
  typedef typename Base::Multiplicity               Multiplicity;

  // Categories
  typedef typename Base::Has_left_category          Has_left_category;
  typedef typename Base::Has_do_intersect_category  Has_do_intersect_category;

  typedef typename internal::Arr_complete_left_side_category< Base >::Category
                                                    Left_side_category;
  typedef typename internal::Arr_complete_bottom_side_category< Base >::Category
                                                    Bottom_side_category;
  typedef typename internal::Arr_complete_top_side_category< Base >::Category
                                                    Top_side_category;
  typedef typename internal::Arr_complete_right_side_category< Base >::Category
                                                    Right_side_category;

protected:
  // All sides
  typedef typename Arr_are_all_sides_oblivious_tag<Left_side_category,
                                                   Bottom_side_category,
                                                   Top_side_category,
                                                   Right_side_category>::result
    Are_all_sides_oblivious_category;

  // left-right dispatch
  typedef CGAL::internal::Arr_left_right_implementation_dispatch<
    Left_side_category, Right_side_category > LR;

  typedef typename LR::Parameter_space_in_x_2_curve_end_tag
    Psx_2_curve_end_tag;
  typedef typename LR::Parameter_space_in_x_2_curve_tag        Psx_2_curve_tag;
  typedef typename LR::Parameter_space_in_x_2_point_tag        Psx_2_point_tag;
  typedef typename LR::Is_on_y_identification_2_curve_tag      Ioyi_2_curve_tag;
  typedef typename LR::Is_on_y_identification_2_point_tag      Ioyi_2_point_tag;
  typedef typename LR::Compare_y_on_boundary_2_points_tag
    Cmp_y_ob_2_points_tag;
  typedef typename LR::Compare_y_near_boundary_2_curve_ends_tag
    Cmp_y_nb_2_curve_ends_tag;

  // bottom-top dispatch
  typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
    Bottom_side_category, Top_side_category > BT;

  typedef typename BT::Parameter_space_in_y_2_curve_end_tag
    Psy_2_curve_end_tag;
  typedef typename BT::Parameter_space_in_y_2_curve_tag        Psy_2_curve_tag;
  typedef typename BT::Parameter_space_in_y_2_point_tag        Psy_2_point_tag;
  typedef typename BT::Is_on_x_identification_2_curve_tag      Ioxi_2_curve_tag;
  typedef typename BT::Is_on_x_identification_2_point_tag      Ioxi_2_point_tag;

  typedef typename BT::Compare_x_at_limit_2_point_curve_end_tag
    Cmp_x_al_2_point_curve_end_tag;
  typedef typename BT::Compare_x_at_limit_2_curve_ends_tag
    Cmp_x_al_2_curve_ends_tag;
  typedef typename BT::Compare_x_near_limit_2_curve_ends_tag
    Cmp_x_nl_2_curve_ends_tag;

  typedef typename BT::Compare_x_on_boundary_2_points_tag
    Cmp_x_ob_2_points_tag;
  typedef typename BT::Compare_x_on_boundary_2_point_curve_end_tag
    Cmp_x_ob_2_point_curve_end_tag;
  typedef typename BT::Compare_x_on_boundary_2_curve_ends_tag
    Cmp_x_ob_2_curve_ends_tag;
  typedef typename BT::Compare_x_near_boundary_2_curve_ends_tag
    Cmp_x_nb_2_curve_ends_tag;

public:

  /// \name Construction.
  //@{
  /*! Default constructor. */
  Arr_traits_basic_adaptor_2() : Base() {}

  /*! Constructor from a base-traits class. */
  Arr_traits_basic_adaptor_2(const Base& traits) : Base(traits) {}
  //@}

  // Inherited functors:
  typedef typename Base::Compare_x_2            Compare_x_2;
  typedef typename Base::Compare_xy_2           Compare_xy_2;
  typedef typename Base::Compare_y_at_x_2       Compare_y_at_x_2;
  typedef typename Base::Construct_min_vertex_2 Construct_min_vertex_2;
  typedef typename Base::Construct_max_vertex_2 Construct_max_vertex_2;
  typedef typename Base::Is_vertical_2          Is_vertical_2;
  typedef typename Base::Compare_y_at_x_right_2 Compare_y_at_x_right_2;
  typedef typename Base::Equal_2                Equal_2;

  /// \name Overriden functors for bounded boundaries.
  //@{

  /*! A functor that compares the y-coordinates of two x-monotone curves
   * immediately to the left of their intersection point.
   */
  class Compare_y_at_x_left_2 {
  public:
    /*!
     * Compare two curves immediately to the left of their intersection point.
     * \param xcv1 The first curve.
     * \param xcv2 The second curve.
     * \param p The query point.
     * \pre The two curves intersect at p, and they are defined to its left.
     * \return SMALLER if xcv1 lies below xcv2 to the left of q;
     *         LARGER if xcv1 lies above xcv2 to the left of q;
     *         EQUAL in case of an overlap to the left of q.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 const Point_2& p) const
    {
      // The function is implemented based on the Has_left category. If the
      // category indicates that the "left" version is available, it calls the
      // function with same name defined in the base class. Otherwise, it
      // uses other predicates to provide this comparison.
      return _compare_y_at_x_left_imp(xcv1, xcv2, p, Has_left_category());
    }

  protected:
    //! The base traits.
    const Self* m_self;

    /*! Constructor.
     * \param tr The base traits class. It must be passed, to handle non
     *           stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_at_x_left_2(const Self* self) : m_self(self) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

    /*!
     * Implementation of the operator() in case the HasLeft tag is true.
     */
    Comparison_result _compare_y_at_x_left_imp(const X_monotone_curve_2& xcv1,
                                               const X_monotone_curve_2& xcv2,
                                               const Point_2& p,
                                               Tag_true) const
    {
      const Base* base = m_self;
      return (base->compare_y_at_x_left_2_object()(xcv1, xcv2, p));
    }

    /*!
     * Implementation of the operator() in case the HasLeft tag is false.
     */
    Comparison_result
    _compare_y_at_x_left_imp(const X_monotone_curve_2& xcv1,
                             const X_monotone_curve_2& xcv2,
                             const Point_2& CGAL_assertion_code(p),
                             Tag_false) const
    {
      Parameter_space_in_x_2  ps_x = m_self->parameter_space_in_x_2_object();
      Parameter_space_in_y_2  ps_y = m_self->parameter_space_in_y_2_object();
      Construct_min_vertex_2  min_vertex =
        m_self->construct_min_vertex_2_object();
      Equal_2                 equal = m_self->equal_2_object();

      // Check if the left ends of the curves are bounded endpoints.
      const Arr_parameter_space  ps_x1 = ps_x(xcv1, ARR_MIN_END);
      const Arr_parameter_space  ps_y1 =
        (ps_x1 != ARR_INTERIOR ? ARR_INTERIOR : ps_y(xcv1, ARR_MIN_END));
      const bool has_left1 = (ps_x1 == ARR_INTERIOR && ps_y1 == ARR_INTERIOR);

      const Arr_parameter_space  ps_x2 = ps_x(xcv2, ARR_MIN_END);
      const Arr_parameter_space  ps_y2 =
        (ps_x2 != ARR_INTERIOR ? ARR_INTERIOR : ps_y(xcv2, ARR_MIN_END));
      const bool has_left2 = (ps_x2 == ARR_INTERIOR && ps_y2 == ARR_INTERIOR);

      CGAL_assertion(ps_x1 != ARR_RIGHT_BOUNDARY &&
                     ps_x2 != ARR_RIGHT_BOUNDARY);

      // Make sure that p lies on both curves, and that both are defined to its
      // right (so their right endpoint is lexicographically larger than p).
      CGAL_precondition_code(
        Compare_xy_2       compare_xy = m_self->compare_xy_2_object();
        Compare_y_at_x_2   compare_y_at_x = m_self->compare_y_at_x_2_object();
      );

      CGAL_precondition(compare_y_at_x(p, xcv1) == EQUAL &&
                         compare_y_at_x(p, xcv2) == EQUAL);

      CGAL_precondition((! has_left1 ||
                          compare_xy(min_vertex(xcv1), p) == SMALLER) &&
                        (! has_left2 ||
                          compare_xy(min_vertex(xcv2), p) == SMALLER));

      // If one of the curves is vertical, it is below the other one.
      Is_vertical_2       is_vertical = m_self->is_vertical_2_object();

      if (is_vertical(xcv1)) return (is_vertical(xcv2)) ? EQUAL : SMALLER;
      else if (is_vertical(xcv2)) return (LARGER);

      // Perform the comparison based on the existence of bounded left
      // endpoints.
      if (has_left1 && has_left2) {
        // Obtain the left endpoints of xcv1 and xcv2.
        Point_2        left1 = min_vertex(xcv1);
        Point_2        left2 = min_vertex(xcv2);

        if (equal(left1, left2))
          // The two curves have a common left endpoint:
          // Compare them to the right of this point.
          return (m_self->compare_y_at_x_right_2_object()(xcv1, xcv2, left1));
      }

      // We now that the curves do not share a common endpoint, and we can
      // compare their relative y-position (which does not change to the left
      // of the given point p).
      return (m_self->compare_y_position_2_object()(xcv1, xcv2));
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 function object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(this); }

  /*! A functor that determines whether two x-monotone curves intersect.
   */
  class Do_intersect_2 {
  public:
    /*!
     * Determine whether two x-monotone curves intersect.
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \return true if xcv1 and xcv2 intersect false otherwise.
     */
    bool operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2) const
    {
      // The function is implemented based on the Has_do_intersect category.
      // If the category indicates that "do_intersect" is available, it calls
      // the function with same name defined in the base class. Otherwise, it
      // uses the intersection construction to implement this predicate.
      return _do_intersect_imp(xcv1, xcv2, Has_do_intersect_category());
    }

  protected:
    //! The base traits.
    const Self* m_self;

    /*! Constructor.
     * \param self The traits adaptor class. It must be passed, to handle
     *             non stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Do_intersect_2(const Self* self) : m_self(self) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

    /*!
     * Implementation of the operator() in case the HasDoIntersect tag is true.
     */
    bool _do_intersect_imp(const X_monotone_curve_2& xcv1,
                           const X_monotone_curve_2& xcv2,
                           Tag_true) const
    {
      const Base* base = m_self;
      return (base->do_intersect_2_object()(xcv1, xcv2));
    }

    /*!
     * Implementation of the operator() in case the HasDoIntersect tag is false.
     */
    bool _do_intersect_imp(const X_monotone_curve_2& xcv1,
                           const X_monotone_curve_2& xcv2,
                           Tag_false) const
    {
      typedef std::pair<Point_2, Multiplicity>          Intersection_point;
      typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                        Intersection_result;
      std::list<Intersection_result> intersections;
      m_self->intersect_2_object()(xcv1, xcv2, back_inserter(intersections));
      return ! intersections.empty();
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 function object. */
  Do_intersect_2 do_intersect_2_object() const { return Do_intersect_2(this); }

  //@}


  /// \name Overriden functors for boundaries.
  //@{

  // left-right

  /*! A functor that determines the location of a geometric object
   * with respect to the parameter space along the x axis.
   */
  class Parameter_space_in_x_2 {
  public:

    /*!
     * Obtain the location of the given curve end in x.
     * \param xcv The curve.
     * \param ind ARR_MIN_END if we refer to xcv's minimal end,
     *            ARR_MAX_END if we refer to its maximal end.
     * \return The location of the curve end.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                    Arr_curve_end ind) const
    {
      // The function is implemented based on the tag dispatching
      // If the traits class does not support special boundaries, we just
      // return ARR_INTERIOR.
      return parameter_space_in_x(xcv, ind, Psx_2_curve_end_tag());
    }

    /*!
     * Obtain the location of the given curve end in x.
     * \param xcv The curve.
     * \return The location of the curve end in x direction.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv) const
    {
      // The function is implemented based on the tag dispatching.
      // If the traits class does not support special boundaries, we just
      // return ARR_INTERIOR.
      return parameter_space_in_x(xcv, Psx_2_curve_tag());
    }

    /*!
     * Obtain the location of the given point end in x.
     * \param p The point.
     * \return The location of the point end in x direction.
     */
    Arr_parameter_space operator()(const Point_2& p) const
    {
      // The function is implemented based on the tag dispatching
      // If the traits class does not support special boundaries, we just
      // return ARR_INTERIOR.
      return parameter_space_in_x(p, Psx_2_point_tag());
    }


  protected:
    //! The base traits.
    const Base* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Parameter_space_in_x_2(const Base* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

    /*!
     * Implementation of the operator() in case the base should be used.
     */
    Arr_parameter_space parameter_space_in_x(const X_monotone_curve_2& xcv,
                                             Arr_curve_end ind,
                                             Arr_use_traits_tag) const
    { return (m_base->parameter_space_in_x_2_object()(xcv, ind)); }

    /*!
     * Implementation of the operator() in case the dummy should be used.
     */
    Arr_parameter_space parameter_space_in_x(const X_monotone_curve_2&,
                                             Arr_curve_end,
                                             Arr_use_dummy_tag) const
    { return ARR_INTERIOR; }

    /*!
     * Implementation of the operator() in case the base should be used.
     */
    Arr_parameter_space parameter_space_in_x(const X_monotone_curve_2& xcv,
                                              Arr_use_traits_tag) const
    { return (m_base->parameter_space_in_x_2_object()(xcv)); }

    /*!
     * Implementation of the operator() in case the dummy should be used.
     */
    Arr_parameter_space parameter_space_in_x(const X_monotone_curve_2&,
                                              Arr_use_dummy_tag) const
    { return ARR_INTERIOR; }

    /*!
     * Implementation of the operator() in case the base should be used.
     */
    Arr_parameter_space parameter_space_in_x(const Point_2& p,
                                             Arr_use_traits_tag) const
    { return m_base->parameter_space_in_x_2_object()(p); }

    /*!
     * Implementation of the operator() in case the dummy should be used.
     */
    Arr_parameter_space parameter_space_in_x(const Point_2&,
                                             Arr_use_dummy_tag) const
    { return ARR_INTERIOR; }
  };

  /*! Obtain an Parameter_space_in_x_2 function object. */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  {
    return Parameter_space_in_x_2(this);
  }

  /*! A function object that determines whether an x-monotone curve or a
   * point coincide with the vertical identification curve.
   */
  class Is_on_y_identification_2 {
  protected:
    //! The base traits.
    const Base* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Is_on_y_identification_2(const Base* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

  public:
    /*! Determones whether a point lies on the vertical identification curve
     * \param p the point.
     * \return true if p lies on the vertical identification curve, and
     * false otherwise.
     */
    bool operator()(const Point_2& p) const
    { return is_on_y_idn(p, Ioyi_2_point_tag()); }

    /*! Determones whether an x-monotone curve coicide with the vertical
     * identification curve
     * \param xcv the point.
     * \return true if xcv coincides with an identification curve,
     * and false otherwise.
     */
    bool operator()(const X_monotone_curve_2& xcv) const
    { return is_on_y_idn(xcv,  Ioyi_2_curve_tag()); }

  private:
    bool is_on_y_idn(const Point_2& p, Arr_use_traits_tag) const
    { return m_base->is_on_y_identification_2_object()(p); }

    bool is_on_y_idn(const Point_2&, Arr_use_dummy_tag) const
    {
      CGAL_error();
      return SMALLER;
    }

    bool is_on_y_idn(const X_monotone_curve_2& xcv, Arr_use_traits_tag) const
    { return m_base->is_on_y_identification_2_object()(xcv); }

    bool is_on_y_idn(const X_monotone_curve_2&, Arr_use_dummy_tag) const
    {
      CGAL_error();
      return SMALLER;
    }
  };

  /*! Obtain a Is_on_y_identification_2 function object. */
  Is_on_y_identification_2 is_on_y_identification_2_object() const
  {
    return Is_on_y_identification_2(this);
  }

  /*! A function object that compares the y-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  protected:
    //! The base traits.
    const Base* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_near_boundary_2(const Base* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

    /*!
     * Implementation of the operator() in case the base should used.
     */
    Comparison_result comp_y_near_bnd(const X_monotone_curve_2& xcv1,
                                      const X_monotone_curve_2& xcv2,
                                      Arr_curve_end ce,
                                      Arr_use_traits_tag) const
    { return m_base->compare_y_near_boundary_2_object()(xcv1, xcv2, ce); }

    /*!
     * Implementation of the operator() in case the dummy should be used.
     */
    Comparison_result comp_y_near_bnd(const X_monotone_curve_2&,
                                      const X_monotone_curve_2&,
                                      Arr_curve_end, Arr_use_dummy_tag) const
    {
      CGAL_error();
      return EQUAL;
    }

  public:
    /*!
     * Compare the relative y-positions of two curve ends.
     * \param xcv1 The first curve.
     * \param xcv2 The second curve.
     * \param ce The relevant end of xcv1 and xcv2.
     * \pre Both curve ends have a special boundary in x.
     * \return SMALLER if xcv1 lies below xcv2;
     *         LARGER if xcv1 lies above xcv2;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const
    {
      // The function is implemented based on the tag disatching
      // If the traits class does not support open curves, we just
      // return EQUAL, as this comparison will not be invoked anyway.
      return comp_y_near_bnd(xcv1, xcv2, ce, Cmp_y_nb_2_curve_ends_tag());
    }
  };

  /*! Obtain a Compare_y_near_boundary_2 functor. */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(this); }

  /*! A function object that compares the y-coordinate of two given points
   * that lie on vertical boundaries.
   */
  class Compare_y_on_boundary_2 {
  protected:
    //! The base traits.
    const Base* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_on_boundary_2(const Base* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

    /*!
     * Implementation of the operator() in case the base should be used.
     */
    Comparison_result comp_y_on_bnd(const Point_2& p1, const Point_2& p2,
                                    Arr_use_traits_tag) const
    { return m_base->compare_y_on_boundary_2_object()(p1, p2); }

    /*!
     * Implementation of the operator() in case the dummy should be used.
     */
    Comparison_result comp_y_on_bnd(const Point_2&, const Point_2&,
                                    Arr_use_dummy_tag) const
    {
      CGAL_error();
      return SMALLER;
    }

  public:
    /*! Compare the relative y-positions of two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \pre Both points lie on vertical boundaries.
     * \return SMALLER if xcv1 lies below xcv2;
     *         LARGER if xcv1 lies above xcv2;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      // The function is implemented based on the tag dispatching.
      // If the traits class does not support open curves, we just
      // return EQUAL, as this comparison will not be invoked anyway.
      return comp_y_on_bnd(p1, p2, Cmp_y_ob_2_points_tag());
    }
  };

  /*! Obtain a Compare_y_on_boundary_2 function object. */
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const
  { return Compare_y_on_boundary_2(this); }

  // bottom-top

  /*! A functor that determines the location of a geometric object
   * with respect to the parameter space along the y axis.
   */
  class Parameter_space_in_y_2 {
  public:
    /*!
     * Obtain the location of the given curve end in y.
     * \param xcv The curve.
     * \param ind ARR_MIN_END if we refer to xcv's minimal end,
     *            ARR_MAX_END if we refer to its maximal end.
     * \return The location of the curve end.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                   Arr_curve_end ind) const
    {
      // The function is implemented based on the tag dispatching.
      // If the traits class does not support special boundaries, we just
      // return ARR_INTERIOR.
      return parameter_space_in_y(xcv, ind, Psy_2_curve_end_tag());
    }

    /*!
     * Obtain the location of the given curve end in y.
     * \param xcv The curve.
     * \return The location of the curve end in y direction.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv) const
    {
      // The function is implemented based on the tag dispatching.
      // If the traits class does not support special boundaries, we just
      // return ARR_INTERIOR.
      return parameter_space_in_y(xcv, Psy_2_curve_tag());
    }

    /*!
     * Obtain the location of the given point end in y.
     * \param p The point.
     * \return The location of the point end in y direction.
     */
    Arr_parameter_space operator()(const Point_2& p) const
    { return parameter_space_in_y(p, Psy_2_point_tag()); }

  protected:
    //! The base traits.
    const Base* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Parameter_space_in_y_2(const Base* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

    /*!
     * Implementation of the operator() in case the base should be used.
     */
    Arr_parameter_space parameter_space_in_y(const X_monotone_curve_2& xcv,
                                             Arr_curve_end ind,
                                             Arr_use_traits_tag) const
    { return m_base->parameter_space_in_y_2_object()(xcv, ind); }

    /*!
     * Implementation of the operator() in case the dummy should be used.
     */
    Arr_parameter_space parameter_space_in_y(const X_monotone_curve_2&,
                                             Arr_curve_end,
                                             Arr_use_dummy_tag) const
    { return ARR_INTERIOR; }

    /*!
     * Implementation of the operator() in case the base should be used.
     */
    Arr_parameter_space parameter_space_in_y(const X_monotone_curve_2& xcv,
                                             Arr_use_traits_tag) const
    { return (m_base->parameter_space_in_x_2_object()(xcv)); }

    /*!
     * Implementation of the operator() in case the dummy should be used.
     */
    Arr_parameter_space parameter_space_in_y(const X_monotone_curve_2&,
                                             Arr_use_dummy_tag) const
    { return ARR_INTERIOR; }

    /*!
     * Implementation of the operator() in case the base should be used.
     */
    Arr_parameter_space parameter_space_in_y(const Point_2& p,
                                             Arr_use_traits_tag) const
    { return m_base->parameter_space_in_y_2_object()(p); }

    /*!
     * Implementation of the operator() in case the dummy should be used.
     */
    Arr_parameter_space parameter_space_in_y(const Point_2&,
                                             Arr_use_dummy_tag) const
    { return ARR_INTERIOR; }
  };

  /*! Obtain an Parameter_space_in_y_2 function object. */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(this); }

  /*! A function object that determines whether an x-monotone curve or a
   * point coincide with the horizontal identification curve.
   */
  class Is_on_x_identification_2 {
  protected:
    //! The base traits.
    const Base* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Is_on_x_identification_2(const Base* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

  public:
    /*! Determones whether a point lies on the horizontal identification curve
     * \param p the point.
     * \return true if p lies on the vertical identification curve, and
     * false otherwise.
     */
    bool operator()(const Point_2& p) const
    { return is_on_idn(p, Ioxi_2_point_tag()); }

    /*! Determones whether an x-monotone curve coicide with the horizontal
     * identification curve
     * \param xcv the point.
     * \return true if xcv coincides with an identification curve,
     * and false otherwise.
     */
    bool operator()(const X_monotone_curve_2& xcv) const
    { return is_on_x_idn(xcv,  Ioxi_2_curve_tag()); }

  private:
    bool is_on_x_idn(const Point_2& p, Arr_use_traits_tag) const
    { return m_base->is_on_x_identification_2_object()(p); }

    bool is_on_x_idn(const Point_2&, Arr_use_dummy_tag) const
    {
      CGAL_error();
      return SMALLER;
    }

    bool is_on_x_idn(const X_monotone_curve_2& xcv, Arr_use_traits_tag) const
    { return m_base->is_on_x_identification_2_object()(xcv); }

    bool is_on_x_idn(const X_monotone_curve_2&, Arr_use_dummy_tag) const
    {
      CGAL_error();
      return SMALLER;
    }
  };

  /*! Obtain a Is_on_x_identification_2 function object. */
  Is_on_x_identification_2 is_on_x_identification_2_object() const
  { return Is_on_x_identification_2(this); }

  /*! A functor that compares the x-limits of curve ends near the
   * boundary of the parameter space
   */
  class Compare_x_at_limit_2 {
  protected:
    //! The base traits.
    const Base * m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_at_limit_2(const Base * base) : m_base(base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

    /*!
     * Implementation of the operator() in case the base should be used.
     */
    Comparison_result _compare_point_curve(const Point_2& p,
                                           const X_monotone_curve_2& xcv,
                                           Arr_curve_end ce,
                                           Arr_use_traits_tag) const
    { return (m_base->compare_x_at_limit_2_object()(p, xcv, ce)); }

    /*!
     * Implementation of the operator() in case the dummy should be used.
     */
    Comparison_result _compare_point_curve(const Point_2&,
                                           const X_monotone_curve_2&,
                                           Arr_curve_end,
                                           Arr_use_dummy_tag) const
    {
      CGAL_error();
      return EQUAL;
    }

    /*!
     * Implementation of the operator() in case the base should be used.
     */
    Comparison_result _compare_curves(const X_monotone_curve_2& xcv1,
                                      Arr_curve_end ce1,
                                      const X_monotone_curve_2& xcv2,
                                      Arr_curve_end ce2,
                                      Arr_use_traits_tag) const
    { return m_base->compare_x_at_limit_2_object()(xcv1, ce1, xcv2, ce2); }

    /*!
     * Implementation of the operator() in case the dummy should be used.
     */
    Comparison_result _compare_curves(const X_monotone_curve_2&,
                                      Arr_curve_end,
                                      const X_monotone_curve_2&,
                                      Arr_curve_end,
                                      Arr_use_dummy_tag) const
    {
      CGAL_error();
      return EQUAL;
    }

  public:
    /*! Compare the x-limits of a vertical curve and another given
     * curve end.
     * \param p A reference point; we refer to a vertical line incident to p.
     * \param xcv The compared curve.
     * \param ind ARR_MIN_END if we refer to xcv's minimal end;
     *            ARR_MAX_END if we refer to its maximal end.
     * \pre xcv's relevant end has a special boundary in y.
     * \return SMALLER if p lies to the left of xcv;
     *         LARGER if p lies to the right xcv;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce) const
    {
      return _compare_point_curve(p, xcv, ce, Cmp_x_al_2_point_curve_end_tag());
    }

    /*! Compare the x-limits of two curve ends on the boundary.
     * \param xcv1 The first curve.
     * \param ind1 ARR_MIN_END if we refer to xcv1's minimal end;
     *             ARR_MAX_END if we refer to its maximal end.
     * \param xcv2 The second curve.
     * \param ind2 ARR_MIN_END if we refer to xcv2's minimal end;
     *             ARR_MAX_END if we refer to its maximal end.
     * \pre Both curve ends have a special boundary in y.
     * \return SMALLER if xcv1 lies to the left of xcv2;
     *         LARGER if xcv1 lies to the right xcv2;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce2) const
    {
      return _compare_curves(xcv1, ce1, xcv2, ce2, Cmp_x_al_2_curve_ends_tag());
    }
  };

  /*! Obtain a Compare_x_at_limit_2 function object. */
  Compare_x_at_limit_2 compare_x_at_limit_2_object() const
  { return Compare_x_at_limit_2(this); }

  /*! A functor that compares the x-coordinates of curve ends near the
   * boundary of the parameter space
   */
  class Compare_x_near_limit_2 {
  protected:
    //! The base traits.
    const Base* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects,(which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_near_limit_2(const Base* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

    /*!
     * Implementation of the operator() in case the base should be used.
     */
    Comparison_result _compare_curves(const X_monotone_curve_2& xcv1,
                                      const X_monotone_curve_2& xcv2,
                                      Arr_curve_end ce,
                                      Arr_use_traits_tag) const
    { return m_base->compare_x_near_limit_2_object()(xcv1, xcv2, ce); }

    /*!
     * Implementation of the operator() in case the dummy should be used.
     */
    Comparison_result _compare_curves(const X_monotone_curve_2&,
                                      const X_monotone_curve_2&,
                                      Arr_curve_end,
                                      Arr_use_dummy_tag) const
    {
      CGAL_error();
      return EQUAL;
    }

  public:
    /*! Compare the relative x-positions of two curve ends.
     * \param xcv1 The first curve.
     * \param xcv2 The second curve.
     * \param ce ARR_MIN_END if we refer to the curves' minimal end;
     *           ARR_MAX_END if we refer to the curves' maximal end.
     * \pre Both curve ends have a special boundary in y.
     * \return SMALLER if xcv1 lies to the left of xcv2;
     *         LARGER if xcv1 lies to the right xcv2;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const
    { return _compare_curves(xcv1, xcv2, ce, Cmp_x_nl_2_curve_ends_tag()); }
  };

  /*! Obtain a Compare_x_near_limit_2 function object. */
  Compare_x_near_limit_2 compare_x_near_limit_2_object() const
  { return Compare_x_near_limit_2(this); }

  /*! A function object that compares the x-coordinate of two given points
   * that lie on horizontal boundaries.
   */
  class Compare_x_on_boundary_2 {
  protected:
    //! The base traits.
    const Base* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_on_boundary_2(const Base* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

  public:
    /*! Compare the x-coordinate of two given points projected onto the
     * horizontal boundaries
     * \param p1 the first point.
     * \param p2 the second point.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { return comp_x_on_bnd(p1, p2, Cmp_x_ob_2_points_tag()); }

    /*! Compare the x-coordinate of a point and a curve-end projected onto the
     * horizontal boundaries
     * \param pt the point.
     * \param xcv the curve
     * \param ce the curve-end
     */
    Comparison_result operator()(const Point_2& pt,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce) const
    { return comp_x_on_bnd(pt, xcv, ce, Cmp_x_ob_2_point_curve_end_tag()); }

    /*! Compare the x-coordinates of two curve-ends projected onto the horizontal
     * boundaries
     * \param xcv1 the curve
     * \param ce1 the curve-end
     * \param xcv2 the curve
     * \param ce2 the curve-end
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce2) const
    { return comp_x_on_bnd(xcv1, ce1, xcv2, ce2, Cmp_x_ob_2_curve_ends_tag()); }

  private:
    Comparison_result comp_x_on_bnd(const Point_2& p1, const Point_2& p2,
                                    Arr_use_traits_tag) const
    { return m_base->compare_x_on_boundary_2_object()(p1, p2); }

    Comparison_result comp_x_on_bnd(const Point_2&, const Point_2&,
                                    Arr_use_dummy_tag) const
    {
      CGAL_error();
      return SMALLER;
    }

    Comparison_result comp_x_on_bnd(const Point_2& pt,
                                    const X_monotone_curve_2& xcv,
                                    Arr_curve_end ce,
                                    Arr_use_traits_tag) const
    { return m_base->compare_x_on_boundary_2_object()(pt, xcv, ce); }

    Comparison_result comp_x_on_bnd(const Point_2&,
                                    const X_monotone_curve_2& /* xcv */,
                                    Arr_curve_end /* ce */,
                                    Arr_use_dummy_tag) const
    {
      CGAL_error();
      return SMALLER;
    }


    Comparison_result comp_x_on_bnd(const X_monotone_curve_2& xcv1,
                                    Arr_curve_end ce1,
                                    const X_monotone_curve_2& xcv2,
                                    Arr_curve_end ce2,
                                     Arr_use_traits_tag) const
    { return m_base->compare_x_on_boundary_2_object()(xcv1, ce1, xcv2, ce2); }

    Comparison_result comp_x_on_bnd(const X_monotone_curve_2& /* xcv1 */,
                                    Arr_curve_end /* ce1 */,
                                    const X_monotone_curve_2& /* xcv2 */,
                                    Arr_curve_end /* ce2 */,
                                    Arr_use_dummy_tag) const
    {
      CGAL_error();
      return SMALLER;
    }
  };

  /*! Obtain a Compare_x_on_boundary_2 function object. */
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const
  { return Compare_x_on_boundary_2(this); }

  /*! A functor that compares the x-coordinates of curve ends near the
   * boundary of the parameter space
   */
  class Compare_x_near_boundary_2 {
  protected:
    //! The base traits.
    const Base* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_near_boundary_2(const Base* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

    /*!
     * Implementation of the operator() in case the base should be used.
     */
    Comparison_result _compare_curves(const X_monotone_curve_2& xcv1,
                                      const X_monotone_curve_2& xcv2,
                                      Arr_curve_end ce,
                                      Arr_use_traits_tag) const
    { return m_base->compare_x_near_boundary_2_object()(xcv1, xcv2, ce); }

    /*!
     * Implementation of the operator() in case the dummy should be used.
     */
    Comparison_result _compare_curves(const X_monotone_curve_2&,
                                      const X_monotone_curve_2&,
                                      Arr_curve_end,
                                      Arr_use_dummy_tag) const
    {
      CGAL_error();
      return EQUAL;
    }

  public:

    /*! Compare the relative x-positions of two curve ends.
     * \param xcv1 The first curve.
     * \param xcv2 The second curve.
     * \param ce ARR_MIN_END if we refer to the curves' minimal end;
     *           ARR_MAX_END if we refer to the curves' maximal end.
     * \pre Both curve ends have a special boundary in y.
     * \return SMALLER if xcv1 lies to the left of xcv2;
     *         LARGER if xcv1 lies to the right xcv2;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const
    { return _compare_curves(xcv1, xcv2, ce, Cmp_x_nb_2_curve_ends_tag()); }
  };

  /*! Obtain a Compare_x_near_boundary_2 function object. */
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const
  { return Compare_x_near_boundary_2(this); }


  // special non-public comparison functors


  /*! A functor that compares the y-coordinates of curve ends near the
   * boundary of the parameter space
   */
  class Compare_y_curve_ends_2 {
  protected:
    //! The base traits.
    const Self* m_self;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_curve_ends_2(const Self* self) : m_self(self) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

  protected:
    // for open
    Comparison_result _compare_curve_ends(const X_monotone_curve_2& xcv1,
                                          const X_monotone_curve_2& xcv2,
                                          Arr_curve_end ce,
                                          Arr_open_side_tag) const
    {
      Comparison_result res =
        m_self->compare_y_near_boundary_2_object()(xcv1, xcv2, ce);
      return res;
    }

    // for non-open
    Comparison_result _compare_curve_ends(const X_monotone_curve_2& xcv1,
                                          const X_monotone_curve_2& xcv2,
                                          Arr_curve_end ce,
                                          Arr_oblivious_side_tag) const
    {
      Comparison_result res =
        m_self->compare_y_on_boundary_2_object()(
          m_self->construct_vertex_at_curve_end_2_object()(xcv1, ce),
          m_self->construct_vertex_at_curve_end_2_object()(xcv2, ce)
      );
      return res;
    }

  public:
    /*! Compare the relative y-positions of two curve ends.
     * \param xcv1 The first curve.
     * \param ind1 ARR_MIN_END if we refer to xcv1's minimal end;
     *             ARR_MAX_END if we refer to its maximal end.
     * \param xcv2 The second curve.
     * \param ind2 ARR_MIN_END if we refer to xcv2's minimal end;
     *             ARR_MAX_END if we refer to its maximal end.
     * \pre Both curve ends have a special boundary in y.
     * \return SMALLER if xcv1 lies to the left of xcv2;
     *         LARGER if xcv1 lies to the right xcv2;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const
    {

      Arr_parameter_space ps1 =
        m_self->parameter_space_in_x_2_object()(xcv1, ce);

      CGAL_precondition(ps1 != ARR_INTERIOR);

      CGAL_assertion_code(
        Arr_parameter_space ps2 =
        m_self->parameter_space_in_x_2_object()(xcv2, ce);
      );
      CGAL_assertion(ps1 == ps2);

      switch (ps1) {
       case ARR_LEFT_BOUNDARY:
        return _compare_curve_ends(xcv1, xcv2, ce, Left_side_category());

       case ARR_RIGHT_BOUNDARY:
        return _compare_curve_ends(xcv1, xcv2, ce, Right_side_category());

       case ARR_INTERIOR: // fall-through
       default:
        CGAL_error(); // never reached
        return CGAL::EQUAL;
      }
    }
  };

  /*! Obtain a Compare_y_curve_ends_2 function object. */
  Compare_y_curve_ends_2 compare_y_curve_ends_2_object() const
  { return Compare_y_curve_ends_2(this); }

  /*! A functor that compares the x-coordinates of curve ends near the
   * boundary of the parameter space
   */
  class Compare_x_point_curve_end_2 {
  protected:
    //! The base traits.
    const Self* m_self;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_point_curve_end_2(const Self* self) : m_self(self) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

  protected:
    // for open
    Comparison_result _compare_point_curve_end(const Point_2& pt,
                                               const X_monotone_curve_2& xcv,
                                               Arr_curve_end ce,
                                               Arr_open_side_tag) const
    {
      Comparison_result res =
        m_self->compare_x_at_limit_2_object()(pt, xcv, ce);
      if ((res != EQUAL) || m_self->is_vertical_2_object()(xcv)) return res;

      // look at the side from which the
      // vertical asymptote is approached
      res = ((ce == CGAL::ARR_MIN_END) ? CGAL::SMALLER : CGAL::LARGER);
      return res;
    }

    // for contraction
    Comparison_result _compare_point_curve_end(const Point_2& pt,
                                               const X_monotone_curve_2& xcv,
                                               Arr_curve_end ce,
                                               Arr_contracted_side_tag) const
    {
      Comparison_result res =
        m_self->compare_x_on_boundary_2_object()(pt, xcv, ce);
      if ((res != EQUAL) || m_self->is_vertical_2_object()(xcv)) return res;

      // look at the side from which the
      // vertical asymptote is approached
      res = ((ce == CGAL::ARR_MIN_END) ? CGAL::SMALLER : CGAL::LARGER);
      return res;
    }

    // for others
    Comparison_result _compare_point_curve_end(const Point_2& pt,
                                               const X_monotone_curve_2& xcv,
                                               Arr_curve_end ce,
                                               Arr_oblivious_side_tag) const
    {
      Comparison_result res =
        m_self->compare_x_on_boundary_2_object()
          (pt, m_self->construct_vertex_at_curve_end_2_object()(xcv, ce));
      return res;
    }

  public:
    /*! Compare the relative x-positions of a point and a curve end.
     * \param pt The point
     * \param xcv The curve.
     * \param ce ARR_MIN_END if we refer to xcv's minimal end;
     *           ARR_MAX_END if we refer to its maximal end.
     * \pre The curve end has a special boundary in y.
     * \return SMALLER if pt lies to the left of xcv;
     *         LARGER if pt lies to the right xcv;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator()(const Point_2& pt,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce) const
    {
      Arr_parameter_space ps = m_self->parameter_space_in_y_2_object()(xcv, ce);

      CGAL_precondition(ps != ARR_INTERIOR);

      switch (ps) {

      case ARR_BOTTOM_BOUNDARY:
        return _compare_point_curve_end(pt, xcv, ce, Bottom_side_category());

      case ARR_TOP_BOUNDARY:
        return _compare_point_curve_end(pt, xcv, ce, Top_side_category());

      case ARR_INTERIOR:
        // fall-through
      default:
        CGAL_error(); // never reached
        return CGAL::EQUAL;
      }
    }
  };

  /*! Obtain a Compare_x_point_curve_end_2 function object. */
  Compare_x_point_curve_end_2 compare_x_point_curve_end_2_object() const
  { return Compare_x_point_curve_end_2(this); }

  /*! A functor that compares the x-coordinates of curve ends near the
   * boundary of the parameter space
   */
  class Compare_x_curve_ends_2 {
  protected:
    //! The base traits.
    const Self* m_self;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_curve_ends_2(const Self* self) : m_self(self) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

  protected:
    Comparison_result _compare_curve_ends_same_x(const X_monotone_curve_2& xcv1,
                                                 Arr_curve_end ce1,
                                                 const X_monotone_curve_2& xcv2,
                                                 Arr_curve_end ce2) const
    {
      CGAL::Comparison_result res = CGAL::EQUAL;

      CGAL::Arr_parameter_space loc1 =
        m_self->parameter_space_in_y_2_object()(xcv1, ce1);
      CGAL::Arr_parameter_space loc2 =
        m_self->parameter_space_in_y_2_object()(xcv2, ce2);
      bool vert1 = m_self->is_vertical_2_object()(xcv1);
      bool vert2 = m_self->is_vertical_2_object()(xcv2);

      // now we are in the open case: ARR_MIN_END > vertical > ARR_MAX_END
      if (vert1) {
        if (!vert2) {
          res = ((ce2 == CGAL::ARR_MIN_END) ? CGAL::SMALLER : CGAL::LARGER);
          return res;
        }
        // both are vertical
        if (loc1 == loc2) { // both ends converge to the same infinity
          res = CGAL::EQUAL;
          return res;
        }
        res = (loc1 == CGAL::ARR_BOTTOM_BOUNDARY ?
               CGAL::SMALLER : CGAL::LARGER);
        return res;
      }

      if (vert2) {
        res = ((ce1 == CGAL::ARR_MIN_END) ? CGAL::LARGER : CGAL::SMALLER);
        return res;
      }

      // otherwise: both ends have asymptotic behaviour
      if (ce1 == ce2) { // both ends approach asymptote from one side
        if (loc1 == loc2) { // need special y-comparison
          res = m_self->compare_x_near_limit_2_object()(xcv1, xcv2, ce2);
          return res;
        }
        // else: order can be determined without y-comparison
        res = CGAL::EQUAL;
        return res;
      }
      // curve ends approach vertical asymptote (or singularity) from
      // different sides => no comparisons required
      res = ((ce1 == CGAL::ARR_MIN_END) ? CGAL::LARGER : CGAL::SMALLER);
      return res;
    }

    // for open
    Comparison_result _compare_curve_ends(const X_monotone_curve_2& xcv1,
                                          Arr_curve_end ce1,
                                          const X_monotone_curve_2& xcv2,
                                          Arr_curve_end ce2,
                                          Arr_open_side_tag) const {

      Comparison_result res =
        m_self->compare_x_at_limit_2_object()(xcv1, ce1, xcv2, ce2);
      if (res == EQUAL) res = _compare_curve_ends_same_x(xcv1, ce1, xcv2, ce2);
      return res;
    }

    // for contracted
    Comparison_result _compare_curve_ends(const X_monotone_curve_2& xcv1,
                                          Arr_curve_end ce1,
                                          const X_monotone_curve_2& xcv2,
                                          Arr_curve_end ce2,
                                          Arr_contracted_side_tag) const {

      Comparison_result res =
        m_self->compare_x_on_boundary_2_object()(xcv1, ce1, xcv2, ce2);
      if (res == EQUAL) res = _compare_curve_ends_same_x(xcv1, ce1, xcv2, ce2);
      return res;
    }

    // for others
    Comparison_result _compare_curve_ends(const X_monotone_curve_2& xcv1,
                                          Arr_curve_end ce1,
                                          const X_monotone_curve_2& xcv2,
                                          Arr_curve_end ce2,
                                          Arr_oblivious_side_tag) const
    {
      Comparison_result res =
        m_self->compare_x_on_boundary_2_object()
        (m_self->construct_vertex_at_curve_end_2_object()(xcv1, ce1),
         m_self->construct_vertex_at_curve_end_2_object()(xcv2, ce2));
      return res;
    }

  public:
    /*! Compare the relative x-positions of two curve ends.
     * \param xcv1 The first curve.
     * \param ind1 ARR_MIN_END if we refer to xcv1's minimal end;
     *             ARR_MAX_END if we refer to its maximal end.
     * \param xcv2 The second curve.
     * \param ind2 ARR_MIN_END if we refer to xcv2's minimal end;
     *             ARR_MAX_END if we refer to its maximal end.
     * \pre Both curve ends have a special boundary in y.
     * \return SMALLER if xcv1 lies to the left of xcv2;
     *         LARGER if xcv1 lies to the right xcv2;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce2) const
    {
      bool first_open = !m_self->is_closed_2_object()(xcv1, ce1);
      bool second_open = !m_self->is_closed_2_object()(xcv2, ce2);

      if (first_open) {
        if (second_open)
          return _compare_curve_ends(xcv1, ce1, xcv2, ce2,
                                     // both sides are open, so pick one
                                     Bottom_side_category());

        return CGAL::opposite
          (m_self->compare_x_point_curve_end_2_object()
           (m_self->construct_vertex_at_curve_end_2_object()(xcv2, ce2),
            xcv1, ce1));
      }

      if (second_open)
        return m_self->compare_x_point_curve_end_2_object()
          (m_self->construct_vertex_at_curve_end_2_object()(xcv1, ce1),
           xcv2, ce2);

      return _compare_curve_ends(xcv1, ce1, xcv2, ce2,
                                 // both sides are non-open, so pick one
                                 Bottom_side_category());
    }
  };

  /*! Obtain a Compare_x_curve_ends_2 function object. */
  Compare_x_curve_ends_2 compare_x_curve_ends_2_object() const
  { return Compare_x_curve_ends_2(this); }

  class Construct_vertex_at_curve_end_2 {
  protected:
    //! The base traits.
    const Self* m_self;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Construct_vertex_at_curve_end_2(const Self* self) : m_self(self) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

  public:
    Point_2 operator()(const X_monotone_curve_2& xcv, Arr_curve_end ce) const
    {
      return (ce == ARR_MIN_END) ?
        m_self->construct_min_vertex_2_object()(xcv) :
        m_self->construct_max_vertex_2_object()(xcv);
    }
  };

  /*! Obtain a Construct_vertex_at_curve_end_2 function object. */
  Construct_vertex_at_curve_end_2 construct_vertex_at_curve_end_2_object() const
  { return Construct_vertex_at_curve_end_2(this); }

  /*! A function object that determines whether a curve end is closed. */
  class Is_closed_2 {
  protected:
    //! The self traits.
    const Self* m_self;

    /*! Constructor.
     * \param self The traits class itself. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Is_closed_2(const Self* self) : m_self(self) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;

    inline bool _is_closed(Arr_oblivious_side_tag) const { return true; }

    inline bool _is_closed(Arr_open_side_tag) const { return false; }

    inline
    bool _is_closed(const X_monotone_curve_2& xcv, Arr_curve_end ce) const
    {
      Arr_parameter_space ps = m_self->parameter_space_in_x_2_object()(xcv, ce);
      if (ARR_INTERIOR == ps)
        ps = m_self->parameter_space_in_y_2_object()(xcv, ce);

      switch (ps) {
       case ARR_LEFT_BOUNDARY: return _is_closed(Left_side_category());
       case ARR_BOTTOM_BOUNDARY: return _is_closed(Bottom_side_category());
       case ARR_TOP_BOUNDARY: return _is_closed(Top_side_category());
       case ARR_RIGHT_BOUNDARY: return _is_closed(Right_side_category());
       case ARR_INTERIOR: // fall-through
       default: return true;
      }
    }

  public:
    /*! Is the end of an x-monotone curve bounded?
     * \param xcv The x-monotone curve.
     * \param ce The end of xcv identifier.
     * \return true is the curve end is bounded, and false otherwise
     */
    bool operator()(const X_monotone_curve_2& xcv, Arr_curve_end ce) const
    { return _is_closed(xcv, ce); }
  };

  /*! Obtain a Is_closed_2 function object. */
  Is_closed_2 is_closed_2_object() const
  { return Is_closed_2(this); }

  //@}

  /// \name Additional auxiliary functors.
  //@{
  class Is_in_x_range_2 {
  public:
    /*!
     * Check whether the given point is in the x-range of the given x-monotone
     * curve.
     * \param xcv The x-monotone curve.
     * \param p The point.
     * \return true if x(xcv_left) <= x(p) <= x(xcv_right), false otherwise.
     */
    bool operator()(const X_monotone_curve_2& xcv, const Point_2& p) const
    {
      Parameter_space_in_x_2 ps_x = m_self->parameter_space_in_x_2_object();
      Parameter_space_in_y_2 ps_y = m_self->parameter_space_in_y_2_object();
      Compare_x_2            compare_x =  m_self->compare_x_2_object();
      Compare_x_point_curve_end_2
        compare_x_point_curve_end = m_self->compare_x_point_curve_end_2_object();

      // Compare p to the position of the left end of the curve.
      // Note that if the left end of xcv lies at x boundary, p is obviously to
      // its right.
      Arr_parameter_space bx = ps_x(xcv, ARR_MIN_END);
      if (bx == ARR_INTERIOR) {
        Arr_parameter_space by = ps_y(xcv, ARR_MIN_END);

        Comparison_result res = (by == ARR_INTERIOR) ?
          // The left endpoint of xcv is a normal point.
          compare_x(p, m_self->construct_min_vertex_2_object()(xcv)) :
          // The left end of xcv lies at y boundary.
          compare_x_point_curve_end(p, xcv, ARR_MIN_END);

        // If p is to the left of the x-range return false.
        // If p is equal, return true.
        if (res == SMALLER) return false;
        else if (res == EQUAL) return true;
      }

      // If necessary, compare p to the right end of the curve.
      // Note that if this end lies at x boundary, p is obviously to its left.
      bx = ps_x(xcv, ARR_MAX_END);
      if (bx != ARR_INTERIOR) return true;

      Arr_parameter_space by = ps_y(xcv, ARR_MAX_END);

      Comparison_result res = (by == ARR_INTERIOR) ?
        // The right endpoint of xcv is a normal point.
        compare_x(p, m_self->construct_max_vertex_2_object()(xcv)) :
        // The right end of xcv lies at y boundary:
        compare_x_point_curve_end(p, xcv, ARR_MAX_END);

      return (res != LARGER);
    }

    /*!
     * Check whether the x-ranges of the given x-monotone curves overlap.
     * \param xcv1 The first x-monotone curve.
     * \param xcv2 The second x-monotone curve.
     * \return (true) if there is an overlap in the x-ranges of the given
     *         curves; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2) const
    {
      Parameter_space_in_x_2  ps_x = m_self->parameter_space_in_x_2_object();
      Parameter_space_in_y_2  ps_y = m_self->parameter_space_in_y_2_object();
      Compare_x_2             compare_x = m_self->compare_x_2_object();
      Construct_min_vertex_2  min_vertex =
        m_self->construct_min_vertex_2_object();
      Construct_max_vertex_2  max_vertex =
        m_self->construct_max_vertex_2_object();
      Compare_x_point_curve_end_2 compare_x_point_curve_end =
        m_self->compare_x_point_curve_end_2_object();
      Compare_x_curve_ends_2 compare_x_curve_ends =
        m_self->compare_x_curve_ends_2_object();

      // Locate the rightmost of the two left endpoints of the two curves.
      // Note that we guard for curve ends with special boundary.
      Arr_parameter_space       ps_x1, ps_y1;
      Arr_parameter_space       ps_x2, ps_y2;
      const X_monotone_curve_2* xcv_left;
      Arr_parameter_space       by_left;
      Comparison_result         res;

      ps_x1 = ps_x(xcv1, ARR_MIN_END);
      ps_x2 = ps_x(xcv2, ARR_MIN_END);

      if (ps_x1 != ARR_INTERIOR) {
        // If both curves are defined at x boundary, they obviously overlap in
        // their x-ranges.
        if (ps_x2 != ARR_INTERIOR) return true;

        // As xcv2 is not defined at x boundary, take its left end as the
        // rightmost of the two left curve ends.
        xcv_left = &xcv2;
        by_left = ps_y(xcv2, ARR_MIN_END);
      }
      else if (ps_x2 != ARR_INTERIOR) {
        // As xcv1 is not defined at x boundary, take its left end as the
        // rightmost of the two left curve ends.
        xcv_left = &xcv1;
        by_left = ps_y(xcv1, ARR_MIN_END);
      }
      else {
        // Compare the (finite) x-coordinates of the two left ends.
        // We take special care of the case of boundaries in y.
        ps_y1 = ps_y(xcv1, ARR_MIN_END);
        ps_y2 = ps_y(xcv2, ARR_MIN_END);

        if (ps_y1 == ARR_INTERIOR) {
          res = (ps_y2 == ARR_INTERIOR) ?
            compare_x(min_vertex(xcv1), min_vertex(xcv2)) :
            compare_x_point_curve_end(min_vertex(xcv1), xcv2, ARR_MIN_END);
        }
        else {
          res = (ps_y2 == ARR_INTERIOR) ?
            opposite(compare_x_point_curve_end(min_vertex(xcv2), xcv1,
                                               ARR_MIN_END)):
            compare_x_curve_ends(xcv1, ARR_MIN_END, xcv2, ARR_MIN_END);
        }

        if (res == LARGER) {
          xcv_left = &xcv1;
          by_left = ps_y1;
        }
        else {
          xcv_left = &xcv2;
          by_left = ps_y2;
        }
      }

      // Locate the leftmost of the two right endpoints of the two curves.
      // Note that we guard for curve ends with special boundary.
      const X_monotone_curve_2* xcv_right;
      Arr_parameter_space       by_right;

      ps_x1 = ps_x(xcv1, ARR_MAX_END);
      ps_x2 = ps_x(xcv2, ARR_MAX_END);

      if (ps_x1 != ARR_INTERIOR) {
        // If both curves are defined at x boundary, they obviously overlap in
        // their x-ranges.
        if (ps_x2 != ARR_INTERIOR) return true;

        // As xcv2 is not defined at x boundary, take its right end as the
        // leftmost of the two right curve ends.
        xcv_right = &xcv2;
        by_right = ps_y(xcv2, ARR_MAX_END);
      }
      else if (ps_x2 != ARR_INTERIOR) {
        // As xcv1 is not defined at x boundary, take its right end as the
        // leftmost of the two right curve ends.
        xcv_right = &xcv1;
        by_right = ps_y(xcv1, ARR_MAX_END);
      }
      else {
        // Compare the (finite) x-coordinates of the two right ends.
        // We take special care of the case of boundaries in y.
        ps_y1 = ps_y(xcv1, ARR_MAX_END);
        ps_y2 = ps_y(xcv2, ARR_MAX_END);

        if (ps_y1 == ARR_INTERIOR) {
          res = (ps_y2 == ARR_INTERIOR) ?
            compare_x(max_vertex(xcv1), max_vertex(xcv2)) :
            compare_x_point_curve_end(max_vertex(xcv1), xcv2, ARR_MAX_END);
        }
        else {
          res = (ps_y2 == ARR_INTERIOR) ?
            opposite(compare_x_point_curve_end(max_vertex(xcv2), xcv1,
                                               ARR_MAX_END)):
            compare_x_curve_ends(xcv1, ARR_MAX_END, xcv2, ARR_MAX_END);
        }

        if (res == SMALLER) {
          xcv_right = &xcv1;
          by_right = ps_y1;
        }
        else {
          xcv_right = &xcv2;
          by_right = ps_y2;
        }
      }

      // Now compare the (finite) x-coordiates of the left end of xcv_left and
      // the right end of xcv_right.
      if (by_left == ARR_INTERIOR) {
        res = (by_right == ARR_INTERIOR) ?
          compare_x(min_vertex(*xcv_left), max_vertex(*xcv_right)) :
          compare_x_point_curve_end(min_vertex(*xcv_left), *xcv_right,
                                    ARR_MAX_END);
      }
      else {
        res =  (by_right == ARR_INTERIOR) ?
          opposite(compare_x_point_curve_end(max_vertex(*xcv_right),
                                             *xcv_left, ARR_MIN_END)) :
          compare_x_curve_ends(*xcv_left, ARR_MIN_END, *xcv_right, ARR_MAX_END);
      }

      // The two curves overlap in their x-range if and only if the left end
      // of xcv_left is not to the right if the right end of xcv_right.
      return (res != LARGER);
    }

  protected:
    //! The traits itself.
    const Self* m_self;

    /*! Constructor.
     * \param self The traits class itself. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Is_in_x_range_2(const Self* self) : m_self(self) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;
  };

  /*! Obtain an Is_in_x_range_2 function object. */
  Is_in_x_range_2 is_in_x_range_2_object() const
  { return Is_in_x_range_2(this); }

  class Compare_y_position_2 {
  public:
    /*!
     * Obtain the relative of two x-monotone curves with overlapping x-ranges
     * that are disjoint in their interiors.
     * \param xcv1 The first x-monotone curve.
     * \param xcv2 The second x-monotone curve.
     * \pre The x-ranges of the two curves overlap.
     * \return SMALLER if xcv1 lies below xcv2;
     *         LARGER if xcv1 lies above xcv2;
     *         EQUAL in case the common x-range is a single point.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2) const
    {
      CGAL_precondition_code(
        Is_in_x_range_2  is_in_x_range = m_self->is_in_x_range_2_object();
      );
      CGAL_precondition(is_in_x_range(xcv1, xcv2));

      /* The traits class which the basic traits adaptor accepts as a template
       * parameter is a model of the ArrangementBasicTraits_2 concept so it
       * needs not to support intersections at all, therefor it is complicated
       * to check if the x-curves are disjoint in their interiors. Moreover,
       * compare_y_position functor is called only from the arrangement class
       * itself (and some related point-location algorithms), and used only
       * for two curves associated with two arrangement halfedges. These curves
       * are guaranteed to be interior-disjoint. So, it seems that there is no
       * gain in checking the precondition, and it is left unimplemented.
       */

      Parameter_space_in_x_2 ps_x = m_self->parameter_space_in_x_2_object();
      Parameter_space_in_y_2 ps_y = m_self->parameter_space_in_y_2_object();
      Compare_y_at_x_2       compare_y_at_x = m_self->compare_y_at_x_2_object();
      Construct_min_vertex_2 min_vertex =
        m_self->construct_min_vertex_2_object();
      Compare_x_point_curve_end_2
        compare_x_point_curve_end = m_self->compare_x_point_curve_end_2_object();
      Compare_x_curve_ends_2
        compare_x_curve_ends = m_self->compare_x_curve_ends_2_object();
      Compare_y_near_boundary_2
        compare_y_near_bnd = m_self->compare_y_near_boundary_2_object();

      // First check whether any of the curves is defined at x boundary.
      const Arr_parameter_space ps_x1 = ps_x(xcv1, ARR_MIN_END);
      const Arr_parameter_space ps_x2 = ps_x(xcv2, ARR_MIN_END);
      Comparison_result         res;

      CGAL_assertion((ps_x1 != ARR_RIGHT_BOUNDARY) &&
                     (ps_x2 != ARR_RIGHT_BOUNDARY));

      if (ps_x1 != ARR_INTERIOR) {
        if (ps_x2 != ARR_INTERIOR)
          // Compare the relative position of the curves at x boundary.
          return (compare_y_near_bnd(xcv1, xcv2, ARR_MIN_END));

        // Check if the left end of xcv2 lies at y boundary.
        const Arr_parameter_space ps_y2 = ps_y(xcv2, ARR_MIN_END);

        // if xcv2 is below xcv1, return LARGER.
        // if xcv2 is above xcv1, return SMALLER.
        if (ps_y2 == ARR_BOTTOM_BOUNDARY) return (LARGER);
        else if (ps_y2 == ARR_TOP_BOUNDARY) return (SMALLER);

        // Compare the position of the left end of xcv2 (which is a normal
        // point) to xcv1.
        res = compare_y_at_x(min_vertex(xcv2), xcv1);

        // Swap the result.
        if (res == EQUAL) return (EQUAL);
        return ((res == SMALLER) ? LARGER : SMALLER);
      }

      if (ps_x2 != ARR_INTERIOR) {
        // Check if the left end of xcv1 lies at y boundary.
        const Arr_parameter_space ps_y1 = ps_y(xcv1, ARR_MIN_END);

        // If xcv1 is below xcv2 return SMALLER.
        // If xcv1 is above xcv2 return LARGER.
        if (ps_y1 == ARR_BOTTOM_BOUNDARY) return (SMALLER);
        else if (ps_y1 == ARR_TOP_BOUNDARY) return (LARGER);

        // Compare the position of the left end of xcv1 (which is a normal
        // point) to xcv2.
        res = compare_y_at_x(min_vertex(xcv1), xcv2);
        return (res);
      }

      // Check if the left curve end lies at y = +/- oo.
      const Arr_parameter_space ps_y1 = ps_y(xcv1, ARR_MIN_END);
      const Arr_parameter_space ps_y2 = ps_y(xcv2, ARR_MIN_END);
      Comparison_result l_res;

      if (ps_y1 != ARR_INTERIOR) {
        if (ps_y2 != ARR_INTERIOR) {
          // The curve ends have special boundary with oposite signs in y,
          // we readily know their relative position (recall that they do not
          // instersect).
          if ((ps_y1 == ARR_BOTTOM_BOUNDARY) && (ps_y2 == ARR_TOP_BOUNDARY))
            return (SMALLER);
          else if ((ps_y1 == ARR_TOP_BOUNDARY) && (ps_y2 == ARR_BOTTOM_BOUNDARY))
            return (LARGER);

          // Both curves have vertical asymptotes with the same sign in y.
          // Check which asymptote is the rightmost. Note that in this case
          // the vertical asymptotes cannot be equal.
          l_res = compare_x_curve_ends(xcv1, ARR_MIN_END, xcv2, ARR_MIN_END);
          CGAL_assertion(l_res != EQUAL);

          if (ps_y1 == ARR_TOP_BOUNDARY) return (l_res);
          else return ((l_res == SMALLER) ? LARGER : SMALLER);
        }

        // xcv1 has a vertical asymptote and xcv2 has a normal left endpoint.
        // Compare the x-positions of this endpoint and the asymptote.
        const Point_2& left2 = min_vertex(xcv2);

        l_res = compare_x_point_curve_end(left2, xcv1, ARR_MIN_END);
        if (l_res == LARGER) {
          // left2 lies in the x-range of xcv1, so it is safe to compare:
          res = compare_y_at_x(left2, xcv1);

          // Swap the result.
          if (res == EQUAL) return (EQUAL);
          return ((res == SMALLER) ? LARGER : SMALLER);
        }
        return ((ps_y1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER);
      }

      if (ps_y2 != ARR_INTERIOR) {
        // xcv2 has a vertical asymptote and xcv1 has a normal left endpoint.
        // Compare the x-positions of this endpoint and the asymptote.
        const Point_2& left1 = min_vertex(xcv1);

        l_res = compare_x_point_curve_end(left1, xcv2, ARR_MIN_END);
        // left1 lies in the x-range of xcv2, so it is safe to compare:
        if (l_res == LARGER) return (compare_y_at_x(left1, xcv2));
        else return (ps_y2 == ARR_BOTTOM_BOUNDARY) ? LARGER : SMALLER;
      }

      // In this case we compare two normal points.
      Compare_xy_2            compare_xy = m_self->compare_xy_2_object();
      Compare_y_at_x_right_2  compare_y_at_x_right =
        m_self->compare_y_at_x_right_2_object();

      // Obtain the left endpoints of xcv1 and xcv2.
      const Point_2& left1 = min_vertex(xcv1);
      const Point_2& left2 = min_vertex(xcv2);

      // Locate the rightmost point of left1 and left2 and compare its position
      // to the other curve.
      l_res = compare_xy(left1, left2);

      if (l_res != SMALLER) {
        // left1 is in the x-range of xcv2:
        res = compare_y_at_x(left1, xcv2);
        if (res == EQUAL)
          // The two curves intersect at left1. If both curves are defined to
          // the right of the reference point, we can compare them to its
          // right. Otherwise, their share a common endpoint (which is the only
          // overlap in their x-ranges) and are really equal.
          if (l_res == EQUAL) res = compare_y_at_x_right(xcv1, xcv2, left1);
        return (res);
      }
      // left2 is in the x-range of xcv1:
      res = compare_y_at_x(left2, xcv1);

      // The two curves share a common endpoint (which is the only overlap
      // in their x-ranges) and are really equal.
      if (res == EQUAL) return (EQUAL);

      // Swap the result:
      return ((res == SMALLER) ? LARGER : SMALLER);
    }

  protected:
    //! The base traits.
    const Self* m_self;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_position_2(const Self* self) : m_self(self) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;
  };

  /*! Obtain a Compare_y_position_2 function object. */
  Compare_y_position_2 compare_y_position_2_object() const
  { return Compare_y_position_2(this); }

  class Is_between_cw_2 {
  public:
    /*!
     * Check whether the given query curve is encountered when rotating the
     * first curve in a clockwise direction around a given point until reaching
     * the second curve.
     * \param xcv The query curve.
     * \param xcv_to_right Is xcv directed from left to right (that is, the
     *                    common vertex is xcv's left endpoint).
     * \param xcv1 The first curve.
     * \param xcv1_to_right Is xcv1 directed from left to right.
     * \param xcv2 The second curve.
     * \param xcv2_to_right Is xcv2 directed from left to right.
     * \param p The point around which we rotate xcv1.
     * \param xcv_equal_xcv1 Output: does xcv equal xcv1.
     * \param xcv_equal_xcv2 Output: does xcv equal xcv2.
     * \pre p is an end-point of all three curves.
     * \return (true) if xcv is between xcv1 and xcv2; (false) otherwise.
     *         If xcv overlaps xcv1 or xcv2 the result is always (false).
     *         If xcv1 and xcv2 overlap, the result is (true), unless xcv
     *         also overlaps them.
     */
    bool operator()(const X_monotone_curve_2& xcv, bool xcv_to_right,
                     const X_monotone_curve_2& xcv1, bool xcv1_to_right,
                     const X_monotone_curve_2& xcv2, bool xcv2_to_right,
                     const Point_2& p,
                     bool& xcv_equal_xcv1,
                     bool& xcv_equal_xcv2) const
    {
      Compare_y_at_x_left_2   compare_y_at_x_left =
        m_self->compare_y_at_x_left_2_object();
      Compare_y_at_x_right_2  compare_y_at_x_right =
        m_self->compare_y_at_x_right_2_object();

      // Initialize output flags.
      xcv_equal_xcv1 = false;
      xcv_equal_xcv2 = false;

      // Take care of the general 4 cases:
      Comparison_result  l_res, r_res;
      Comparison_result  res1, res2;

      if (!xcv1_to_right && !xcv2_to_right) {
        // Case 1: Both xcv1 and xcv2 are defined to the left of p.
        l_res = compare_y_at_x_left(xcv1, xcv2, p);

        if (l_res == LARGER) {
          // Case 1(a) : xcv1 is above xcv2.
          if (!xcv_to_right) {
            res1 = compare_y_at_x_left(xcv1, xcv, p);
            res2 = compare_y_at_x_left(xcv2, xcv, p);
            if (res1 == EQUAL) xcv_equal_xcv1 = true;
            if (res2 == EQUAL) xcv_equal_xcv2 = true;
            return (res1 == SMALLER || res2 == LARGER);
          }
          return true;
        }

        if (l_res == SMALLER) {
          // Case 1(b): xcv1 is below xcv2.
          if (!xcv_to_right) {
            res1 = compare_y_at_x_left(xcv1, xcv, p);
            res2 = compare_y_at_x_left(xcv2, xcv, p);
            if (res1 == EQUAL) xcv_equal_xcv1 = true;
            if (res2 == EQUAL) xcv_equal_xcv2 = true;
            return (res1 == SMALLER && res2  == LARGER);
          }
          return false;
        }

        // Overlapping segments.
        if (!xcv_to_right) {
          res1 = compare_y_at_x_left(xcv1, xcv, p);
          if (res1 == EQUAL) {
            xcv_equal_xcv1 = true;
            xcv_equal_xcv2 = true;
            return false;
          }
          return true;
        }
        return true;
      }

      if (xcv1_to_right && xcv2_to_right) {
        // Case 2: Both xcv1 and xcv2 are defined to the right of p.
        r_res = compare_y_at_x_right(xcv1, xcv2, p);

        if (r_res == LARGER) {
          // Case 2(a) : xcv1 is above xcv2.
          if (xcv_to_right) {
            res1 = compare_y_at_x_right(xcv1, xcv, p);
            res2 = compare_y_at_x_right(xcv2, xcv, p);
            if (res1 == EQUAL) xcv_equal_xcv1 = true;
            if (res2 == EQUAL) xcv_equal_xcv2 = true;
            return (res1 == LARGER && res2 == SMALLER);
          }
          return false;
        }

        if (r_res == SMALLER) {
          // Case 2(b): xcv1 is below xcv2.
          if (xcv_to_right) {
            res1 = compare_y_at_x_right(xcv1, xcv, p);
            res2 = compare_y_at_x_right(xcv2, xcv, p);
            if (res1 == EQUAL) xcv_equal_xcv1 = true;
            if (res2 == EQUAL) xcv_equal_xcv2 = true;
            return ((res1 == LARGER) || (res2 == SMALLER));
          }
          return true;
        }

        // Overlapping segments.
        if (xcv_to_right) {
          res1 = compare_y_at_x_right(xcv1, xcv, p);
          if (res1 == EQUAL) {
            xcv_equal_xcv1 = true;
            xcv_equal_xcv2 = true;
            return false;
          }
          return true;
        }
        return true;
      }

      if (!xcv1_to_right && xcv2_to_right) {
        // Case 3: xcv1 is defined to the left of p, and xcv2 to its right.
        if (!xcv_to_right) {
          res1 = compare_y_at_x_left(xcv1, xcv, p);
          if (res1 == EQUAL) xcv_equal_xcv1 = true;
          return (res1 == SMALLER);
        }

        res2 = compare_y_at_x_right(xcv2, xcv, p);
        if (res2 == EQUAL) xcv_equal_xcv2 = true;
        return (res2 == SMALLER);
      }

      CGAL_assertion(xcv1_to_right && !xcv2_to_right);

      // Case 4: xcv1 is defined to the right of p, and xcv2 to its left.
      if (xcv_to_right) {
        res1 = compare_y_at_x_right(xcv1, xcv, p);
        if (res1 == EQUAL) xcv_equal_xcv1 = true;
        return (res1  == LARGER);
      }

      res2 = compare_y_at_x_left(xcv2, xcv, p);
      if (res2 == EQUAL) xcv_equal_xcv2 = true;
      return (res2 == LARGER);
    }

  protected:
    //! The base traits.
    const Self* m_self;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Is_between_cw_2(const Self* self) : m_self(self) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;
  };

  /*! Obtain an Is_between_cw_2 function object. */
  Is_between_cw_2 is_between_cw_2_object() const
  { return Is_between_cw_2(this); }

  class Compare_cw_around_point_2 {
  public:
    /*!
     * Compare the two interior disjoint x-monotone curves in a clockwise
     * order around their common endpoint.
     * \param xcv1 The first curve.
     * \param xcv1_to_right Is xcv1 directed from left to right.
     * \param xcv2 The second curve.
     * \param xcv2_to_right Is xcv2 directed from left to right.
     * \param p The common endpoint.
     * \param from_top (true) if we start from 12 o'clock,
     *                 (false) if we start from 6 o'clock.
     * \pre The point p is an endpoint of both curves.
     * \return SMALLER if we encounter xcv1 before xcv2;
     *         LARGER if we encounter xcv2 before xcv1;
     *         EQUAL otherwise.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 bool xcv1_to_right,
                                 const X_monotone_curve_2& xcv2,
                                 bool xcv2_to_right,
                                 const Point_2& p,
                                 bool from_top = true) const
    {
      // Act according to where xcv1 and xcv2 lie.
      if (!xcv1_to_right && !xcv2_to_right)
        // Both are defined to the left of p, and we encounter xcv1 before
        // xcv2 if it is below xcv2:
        return (m_self->compare_y_at_x_left_2_object()(xcv1, xcv2, p));

      if (xcv1_to_right && xcv2_to_right)
        // Both are defined to the right of p, and we encounter xcv1 before
        // xcv2 if it is above xcv2. We therefore reverse the order of the
        // curves when we invoke compare_y_at_x_right:
        return (m_self->compare_y_at_x_right_2_object()(xcv2, xcv1, p));

      if (!xcv1_to_right && xcv2_to_right)
        // If we start from the top, we encounter the right curve (which
        // is xcv2) first. If we start from the bottom, we encounter xcv1 first.
        return (from_top ? LARGER : SMALLER);

      CGAL_assertion(xcv1_to_right && !xcv2_to_right);

      // If we start from the top, we encounter the right curve (which
      // is xcv1) first. If we start from the bottom, we encounter xcv2 first.
      return (from_top ? SMALLER : LARGER);
    }

  protected:
    //! The base traits.
    const Self* m_self;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_cw_around_point_2(const Self* self) : m_self(self) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_basic_adaptor_2<Base>;
  };

  /*! Obtain a Compare_cw_around_point_2 function object. */
  Compare_cw_around_point_2 compare_cw_around_point_2_object() const
  { return Compare_cw_around_point_2(this); }
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
  typedef ArrangementTraits_                           Base_traits_2;
  typedef Arr_traits_basic_adaptor_2<Base_traits_2>    Base;
  typedef Arr_traits_adaptor_2<Base_traits_2>          Self;

  typedef typename Base_traits_2::Curve_2              Curve_2;
  typedef typename Base::X_monotone_curve_2            X_monotone_curve_2;
  typedef typename Base::Point_2                       Point_2;
  typedef typename Base::Multiplicity                  Multiplicity;

  // Categories.
  typedef typename Base::Has_left_category             Has_left_category;
  typedef typename Base::Has_merge_category            Has_merge_category;
  typedef typename Base::Has_do_intersect_category     Has_do_intersect_category;

  typedef typename Base::Left_side_category            Left_side_category;
  typedef typename Base::Bottom_side_category          Bottom_side_category;
  typedef typename Base::Top_side_category             Top_side_category;
  typedef typename Base::Right_side_category           Right_side_category;

  typedef typename Base::Are_all_sides_oblivious_category
    Are_all_sides_oblivious_category;

  /// \name Construction.
  //@{
  /*! Default constructor. */
  Arr_traits_adaptor_2() : Base() {}

  /*! Constructor from a base-traits class. */
  Arr_traits_adaptor_2(const Base_traits_2& traits) : Base(traits) {}
  //@}

  // Inherited functors:
  typedef typename Base::Compare_x_2            Compare_x_2;
  // typedef typename Base::Compare_xy_2           Compare_xy_2;
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

  /*! A functor that compares two points or two x-monotone curves
   * lexigoraphically. Twp points are compared firest by their x-coordinates,
   * then by their y-coordinates. Two curves are compared first their left-most
   * endpoint, then by the graphs, and finally by their right-most endpoint.
   */
  class Compare_xy_2 {
      typedef std::pair<Point_2, Multiplicity>          Intersection_point;
      typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                        Intersection_result;

  public:
    /*! Compare two points lexigoraphically: by x, then by y.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return SMALLER - x(p1) < x(p2);
     *         SMALLER - x(p1) = x(p2) and y(p1) < y(p2);
     *         EQUAL   - x(p1) = x(p2) and y(p1) = y(p2);
     *         LARGER  - x(p1) = x(p2) and y(p1) > y(p2);
     *         LARGER  - x(p1) > x(p2).
     * \pre p1 does not lie on the boundary.
     * \pre p2 does not lie on the boundary.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      Base base(m_self);
      return base.compare_xy_2_object()(p1, p2);
    }

    /*! Compare two x-monotone curves lexigoraphically.
     * Input: C1, C2, intersections = empty
     * compare(C1, C2, intersections)
     *   Compare the left-most points of C1 and C2.
     *   If not equal return the comparison result.
     *   Otherwise (the left-most points are equal)
     *     Compare C1 and C2 to the right of the point.
     *     If not equal return the comparison result.
     *     Otherwise (they overlap)
     *       If intersection is empty, compute the intersections,
     *       Remove the first overlapping section from c1, c2, and intersections.
     *       If intersections is empty
     *         Compare the right-most point and return the result.
     *       return compare(C1, C2, intersections).
     */
    Comparison_result operator()(const X_monotone_curve_2& c1,
                                 const X_monotone_curve_2& c2) const
    {
      std::list<Intersection_result> intersections;
      return operator()(c1, c2, intersections,
                        Are_all_sides_oblivious_category());
    }

  protected:
    //! The base traits.
    const Self& m_self;

    /*! Constructor.
     * \param trait The base traits class. It must be passed, to handle non
     *              stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_xy_2(const Self& self) : m_self(self) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_adaptor_2<Base_traits_2>;

    /// Point-curve
    //@{
    /*! Compare a point and a curve end.
     * Dispatch calls to traits that handle open and close boundaries, resp.
     * The only reason for this dispatcher is the poor choice of different
     * names for the Traits functors that handle close and open boundaries:
     * Open boundary traits use: Compare_x_at_limit_2
     * Close boundary traits use: Compare_x_on_boundary_2
     */
    Comparison_result cmp_x_on_bnd(const Point_2& p,
                                   const X_monotone_curve_2& c,
                                   Arr_curve_end ce) const
    {
      // The complete code would need to check whether ce is bottom or top and
      // use the corresponding dispatching criterion, but
      //  (i) in all out traits if bottom is open, so is top, and
      // (ii) this is going to change soon....
      return cmp_x_on_bnd(p, c, ce, Bottom_side_category());
    }

    // Open
    Comparison_result cmp_x_on_bnd(const Point_2& p,
                                   const X_monotone_curve_2& c,
                                   Arr_curve_end ce,
                                   Arr_open_side_tag) const
    { return m_self.compare_x_at_limit_2_object()(p, c, ce); }

    // Close
    Comparison_result cmp_x_on_bnd(const Point_2& p,
                                   const X_monotone_curve_2& c,
                                   Arr_curve_end ce,
                                   Arr_oblivious_side_tag)  const
    { return m_self.compare_x_on_boundary_2_object()(p, c, ce); }
    //@}

    /// curve-curve
    //@{
    /*! Compare a curve end and a curve end.
     * Dispatch calls to traits that handle open and close boundaries, resp.
     * The only reason for this dispatcher is the poor choice of different
     * names for the Traits functors that handle close and open boundaries:
     * Open boundary traits use: Compare_x_at_limit_2
     * Close boundary traits use: Compare_x_on_boundary_2
     */
    Comparison_result cmp_x_on_bnd(const X_monotone_curve_2& c1,
                                   Arr_curve_end ce1,
                                   const X_monotone_curve_2& c2,
                                   Arr_curve_end ce2) const
    {
      // The complete code would need to check the combination of ce1 and ce2,
      // whther they are (booton,bottom), (bottom,top), (top,bottom) or
      // (top,top) and use the corresponding dispatching criteria, but we don't
      // even have the interface to handle mixed (e.g., open bottom and
      // closed top. Anyway, this is going to change soon....
      return cmp_x_on_bnd(c1, ce1, c2, ce2, Bottom_side_category());
    }

    // Open
    Comparison_result cmp_x_on_bnd(const X_monotone_curve_2& c1,
                                   Arr_curve_end ce1,
                                   const X_monotone_curve_2& c2,
                                   Arr_curve_end ce2,
                                   Arr_open_side_tag) const
    { return m_self.compare_x_at_limit_2_object()(c1, ce1, c2, ce2); }

    // Close
    Comparison_result cmp_x_on_bnd(const X_monotone_curve_2& c1,
                                   Arr_curve_end ce1,
                                   const X_monotone_curve_2& c2,
                                   Arr_curve_end ce2,
                                   Arr_oblivious_side_tag)  const
    { return m_self.compare_x_on_boundary_2_object()(c1, ce1, c2, ce2); }
    //@}

    /*! Compare the max end of two x-monotone curves lexigoraphically.
     */
    Comparison_result compare_max_end(const X_monotone_curve_2& c1,
                                      const X_monotone_curve_2& c2,
                                      Arr_all_sides_oblivious_tag) const
    {
      typedef typename Self::Construct_max_vertex_2 Construct_max_vertex_2;
      Construct_max_vertex_2 ctr_max =
        m_self.construct_max_vertex_2_object();
      const Point_2& p1 = ctr_max(c1);
      const Point_2& p2 = ctr_max(c2);
      return operator()(p1, p2);
    }

    /*! Compare the max (right) end of two x-monotone curves lexigoraphically.
     * \pre the curve overlap
     */
    Comparison_result compare_max_end(const X_monotone_curve_2& c1,
                                      const X_monotone_curve_2& c2,
                                      Arr_not_all_sides_oblivious_tag) const
    {
      typedef typename Base::Parameter_space_in_x_2     Parameter_space_in_x_2;
      typedef typename Base::Parameter_space_in_y_2     Parameter_space_in_y_2;
      Parameter_space_in_x_2 psx = m_self.parameter_space_in_x_2_object();
      Parameter_space_in_y_2 psy = m_self.parameter_space_in_y_2_object();

      Arr_parameter_space px1 = psx(c1, ARR_MAX_END);
      Arr_parameter_space py1 = psy(c1, ARR_MAX_END);

      Arr_parameter_space px2 = psx(c2, ARR_MAX_END);
      Arr_parameter_space py2 = psy(c2, ARR_MAX_END);

      // Handle the trivial cases:
      if ((px1 == ARR_LEFT_BOUNDARY) && (px2 != ARR_LEFT_BOUNDARY))
        return SMALLER;

      if ((px2 == ARR_LEFT_BOUNDARY) && (px1 != ARR_LEFT_BOUNDARY))
        return LARGER;

      if ((px1 == ARR_RIGHT_BOUNDARY) && (px2 != ARR_RIGHT_BOUNDARY))
        return LARGER;

      if ((px2 == ARR_RIGHT_BOUNDARY) && (px1 != ARR_RIGHT_BOUNDARY))
        return SMALLER;

      // The only casese left are px1,px2 = (I,I), (L,L), and (R,R)

      if (px1 == ARR_INTERIOR) {
        CGAL_assertion(px2 == ARR_INTERIOR);

        if ((py1 == ARR_INTERIOR) && (py2 == ARR_INTERIOR))
          return compare_max_end(c1, c2, Arr_all_sides_oblivious_tag());

        // px1, px2, py1 are interior
        if (py1 == ARR_INTERIOR) {
          CGAL_assertion(py2 != ARR_INTERIOR);
#if 0
          // This code is retained (commented out) because in the future, when
          // Compare_x_at_limit_2 and Compare_x_on_boundary will be consolidated,
          // say, to Compare_x_on_boundary, it will replace the call below.
          typedef typename Self::Compare_x_on_boundary_2 Compare_x_on_boundary_2;
          Compare_x_on_boundary_2 cmp_x_on_bnd =
            m_self.compare_x_on_boundary_2_object();
#endif
          const Point_2& c1_max = m_self.construct_max_vertex_2_object()(c1);
          Comparison_result res = cmp_x_on_bnd(c1_max, c2, ARR_MAX_END);
          if (res != EQUAL) return res;

          return (py2 == ARR_TOP_BOUNDARY) ? SMALLER : LARGER;
        }

        // px1, px2, py2 are interior
        if (py2 == ARR_INTERIOR) {
          CGAL_assertion(py1 != ARR_INTERIOR);
#if 0
          // This code is retained (commented out) because in the future, when
          // Compare_x_at_limit_2 and Compare_x_on_boundary will be consolidated,
          // say, to Compare_x_on_boundary, it will replace the call below.
          typedef typename Self::Compare_x_on_boundary_2 Compare_x_on_boundary_2;
          Compare_x_on_boundary_2 cmp_x_on_bnd =
            m_self.compare_x_on_boundary_2_object();
#endif
          const Point_2& c2_max = m_self.construct_max_vertex_2_object()(c2);
          Comparison_result res = cmp_x_on_bnd(c2_max, c1, ARR_MAX_END);
          if (res != EQUAL) return CGAL::opposite(res);

          return (py1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;
        }

        // Both py1 and py2 not interior
#if 0
        // This code is retained (commented out) because in the future, when
        // Compare_x_at_limit_2 and Compare_x_on_boundary will be consolidated,
        // say, to Compare_x_on_boundary, it will replace the call below.
        typedef typename Self::Compare_x_on_boundary_2 Compare_x_on_boundary_2;
        Compare_x_on_boundary_2 cmp_x_on_bnd =
            m_self.compare_x_on_boundary_2_object();
#endif
        Comparison_result res = cmp_x_on_bnd(c1, ARR_MAX_END, c2, ARR_MAX_END);
        if (res != EQUAL) return res;

        if ((py1 == ARR_BOTTOM_BOUNDARY) && (py2 != ARR_BOTTOM_BOUNDARY))
          return SMALLER;
        if ((py1 == ARR_TOP_BOUNDARY) && (py2 != ARR_TOP_BOUNDARY))
          return LARGER;
        return EQUAL;
      }

      // Both endpoints lie either on the left boundary or on the right
      // boundary, which means that their x-coordinates are equal.
      // Handle the trivial cases:
      if ((py1 == ARR_BOTTOM_BOUNDARY) && (py2 != ARR_BOTTOM_BOUNDARY))
        return SMALLER;
      if ((py1 == ARR_TOP_BOUNDARY) && (py2 != ARR_TOP_BOUNDARY))
        return LARGER;

      typedef typename Self::Compare_y_on_boundary_2  Compare_y_on_boundary_2;
      Compare_y_on_boundary_2 cmp_y_on_bnd =
        m_self.compare_y_on_boundary_2_object();
      const Point_2& c1_max = m_self.construct_max_vertex_2_object()(c1);
      const Point_2& c2_max = m_self.construct_max_vertex_2_object()(c2);
      Comparison_result res = cmp_y_on_bnd(c1_max, c2_max);
      return res;
    }

    /*! Split 2 given curves that overlap and have a common sub-curve on their
     * right. Then compare the remaining portions of the curves, respectively.
     */
    Comparison_result
    compare_remainder(const X_monotone_curve_2& c1,
                      const X_monotone_curve_2& c2,
                      std::list<Intersection_result>& intersections) const
    {
      // Right-most sections are equal.
      // Advance to the next respective sections:
      if (intersections.empty()) {
        typedef typename Self::Intersect_2          Intersect_2;
        Intersect_2 intersect = m_self.intersect_2_object();
        intersect(c1, c2, std::back_inserter(intersections));
      }
      // Verify the first intersection is an overlap, remove it, and
      // recursively call.
      const X_monotone_curve_2* xcv =
        boost::get<X_monotone_curve_2>(&(intersections.front()));
      if (! xcv) {
        CGAL_error_msg("The first intersection is not an overlap!");
        return SMALLER;
      }
      intersections.pop_front();
      if (intersections.empty())
        return compare_max_end(c1, c2, Are_all_sides_oblivious_category());

      // Otherwise, split and continue.
      typedef typename Self::Split_2        Split_2;
      Split_2 split = m_self.split_2_object();
      X_monotone_curve_2 c11, c12, c21, c22;
      Construct_max_vertex_2 ctr_max = m_self.construct_max_vertex_2_object();
      const Point_2& p1 = ctr_max(*xcv);
      const Point_2& p2 = ctr_max(*xcv);
      split(c1, p1, c11, c12);
      split(c2, p2, c21, c22);
      return operator()(c12, c22, intersections,
                        Are_all_sides_oblivious_category());
    }

    /*! Compare two x-monotone curves lexigoraphically.
     */
    Comparison_result operator()(const X_monotone_curve_2& c1,
                                 const X_monotone_curve_2& c2,
                                 std::list<Intersection_result>& intersections,
                                 Arr_all_sides_oblivious_tag) const
    {
      const Point_2& c1_min = m_self.construct_min_vertex_2_object()(c1);
      const Point_2& c2_min = m_self.construct_min_vertex_2_object()(c2);

      Comparison_result res = operator()(c1_min, c2_min);
      if (res != EQUAL) return res;

      // Left-most points are equal.
      // Compare their y-coordinates to their right:
      res = m_self.compare_y_at_x_right_2_object()(c1, c2, c1_min);
      if (res != EQUAL) return res;

      return compare_remainder(c1, c2, intersections);
    }

    /*! Compare two x-monotone curves lexigoraphically.
     */
    Comparison_result operator()(const X_monotone_curve_2& c1,
                                 const X_monotone_curve_2& c2,
                                 std::list<Intersection_result>& intersections,
                                 Arr_not_all_sides_oblivious_tag) const
    {
      typedef typename Base::Parameter_space_in_x_2     Parameter_space_in_x_2;
      typedef typename Base::Parameter_space_in_y_2     Parameter_space_in_y_2;
      Parameter_space_in_x_2 psx = m_self.parameter_space_in_x_2_object();
      Parameter_space_in_y_2 psy = m_self.parameter_space_in_y_2_object();

      Arr_parameter_space min_px1 = psx(c1, ARR_MIN_END);
      Arr_parameter_space min_py1 = psy(c1, ARR_MIN_END);

      Arr_parameter_space min_px2 = psx(c2, ARR_MIN_END);
      Arr_parameter_space min_py2 = psy(c2, ARR_MIN_END);

      // Handle the trivial cases:
      if ((min_px1 == ARR_LEFT_BOUNDARY) && (min_px2 != ARR_LEFT_BOUNDARY))
        return SMALLER;

      if ((min_px2 == ARR_LEFT_BOUNDARY) && (min_px1 != ARR_LEFT_BOUNDARY))
        return LARGER;

      if ((min_px1 == ARR_RIGHT_BOUNDARY) && (min_px2 != ARR_RIGHT_BOUNDARY))
        return LARGER;

      if ((min_px2 == ARR_RIGHT_BOUNDARY) && (min_px1 != ARR_RIGHT_BOUNDARY))
        return SMALLER;

      // The only casese left are px1,px2 = (I,I), (L,L), and (R,R)

      if (min_px1 == ARR_INTERIOR) {
        CGAL_assertion(min_px2 == ARR_INTERIOR);

        if ((min_py1 == ARR_INTERIOR) && (min_py2 == ARR_INTERIOR)) {
          return operator()(c1, c2, intersections,
                            Arr_all_sides_oblivious_tag());
        }

        //
        if (min_py1 == ARR_INTERIOR) {
          CGAL_assertion(min_py2 != ARR_INTERIOR);
#if 0
          // This code is retained (commented out) because in the future, when
          // Compare_x_at_limit_2 and Compare_x_on_boundary will be consolidated,
          // say, to Compare_x_on_boundary, it will replace the call below.
          typedef typename Self::Compare_x_on_boundary_2 Compare_x_on_boundary_2;
          Compare_x_on_boundary_2 cmp_x_on_bnd =
            m_self.compare_x_on_boundary_2_object();
#endif
          const Point_2& c1_min = m_self.construct_min_vertex_2_object()(c1);
          Comparison_result res = cmp_x_on_bnd(c1_min, c2, ARR_MIN_END);
          if (res != EQUAL) return res;

          return (min_py2 == ARR_TOP_BOUNDARY) ? SMALLER : LARGER;
        }

        if (min_py2 == ARR_INTERIOR) {
          CGAL_assertion(min_py1 != ARR_INTERIOR);
#if 0
          // This code is retained (commented out) because in the future, when
          // Compare_x_at_limit_2 and Compare_x_on_boundary will be consolidated,
          // say, to Compare_x_on_boundary, it will replace the call below.
          typedef typename Self::Compare_x_on_boundary_2 Compare_x_on_boundary_2;
          Compare_x_on_boundary_2 cmp_x_on_bnd =
            m_self.compare_x_on_boundary_2_object();
#endif
          const Point_2& c2_min = m_self.construct_min_vertex_2_object()(c2);

          Comparison_result res = cmp_x_on_bnd(c2_min, c1, ARR_MIN_END);
          if (res != EQUAL) return CGAL::opposite(res);

          return (min_py1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;
        }

        // Both min_py1 and min_py2 not interior
#if 0
        // This code is retained (commented out) because in the future, when
        // Compare_x_at_limit_2 and Compare_x_on_boundary will be consolidated,
        // say, to Compare_x_on_boundary, it will replace the call below.
        typedef typename Self::Compare_x_on_boundary_2 Compare_x_on_boundary_2;
        Compare_x_on_boundary_2 cmp_x_on_bnd =
            m_self.compare_x_on_boundary_2_object();
#endif
        Comparison_result res = cmp_x_on_bnd(c1, ARR_MIN_END, c2, ARR_MIN_END);
        if (res != EQUAL) return res;

        if ((min_py1 == ARR_BOTTOM_BOUNDARY) && (min_py2 == ARR_TOP_BOUNDARY))
          return SMALLER;
        if ((min_py1 == ARR_TOP_BOUNDARY) && (min_py2 == ARR_BOTTOM_BOUNDARY))
          return LARGER;

        // Left-most points are equal.
        // Compare their x-coordinates near the common left endpoint.
        res = m_self.compare_x_near_boundary_2_object()(c1, c2, ARR_MIN_END);
        if (res != EQUAL) return res;

        return compare_remainder(c1, c2, intersections);
      }


      if ((min_py1 == ARR_BOTTOM_BOUNDARY) && (min_py2 != ARR_BOTTOM_BOUNDARY))
        return SMALLER;
      if ((min_py1 == ARR_TOP_BOUNDARY) && (min_py2 != ARR_TOP_BOUNDARY))
        return LARGER;

      if (min_px1 == ARR_LEFT_BOUNDARY) {
        // The min points of the two curves lie on the left boundary.
        CGAL_assertion(min_px2 == ARR_LEFT_BOUNDARY);
        typedef typename Self::Compare_y_on_boundary_2  Compare_y_on_boundary_2;
        Compare_y_on_boundary_2 cmp_y_on_bnd =
          m_self.compare_y_on_boundary_2_object();
        const Point_2& c1_min = m_self.construct_min_vertex_2_object()(c1);
        const Point_2& c2_min = m_self.construct_min_vertex_2_object()(c2);
        Comparison_result res = cmp_y_on_bnd(c1_min, c2_min);
        if (res != EQUAL) return res;

        typedef typename Self::Is_vertical_2    Is_vertical_2;
        Is_vertical_2 is_vert = m_self.is_vertical_2_object();
        bool vert1 = is_vert(c1);
        bool vert2 = is_vert(c2);
        if (vert1 && ! vert2) return SMALLER;
        if (vert2 && ! vert1) return LARGER;
        if (vert1 && vert2) {
          const Point_2& c1_max = m_self.construct_max_vertex_2_object()(c1);
          const Point_2& c2_max = m_self.construct_max_vertex_2_object()(c2);
          res = cmp_y_on_bnd(c1_max, c2_max);
          if (res == SMALLER) return SMALLER;
        }
        // Compare slightly to the right near the boundary.
        typedef typename Self::Compare_y_near_boundary_2
          Compare_y_near_boundary_2;
        Compare_y_near_boundary_2 cmp_y_near_bnd =
          m_self.compare_y_near_boundary_2_object();
        res = cmp_y_near_bnd(c1, c2, CGAL::ARR_MIN_END);
        if (res != EQUAL) return res;

        // The two curves overlap on their right.
        // Intersect them, and recursively compare the remaining portions,
        // respectively.
        return compare_remainder(c1, c2, intersections);
      }

      CGAL_assertion(min_px1 == ARR_RIGHT_BOUNDARY);
      CGAL_assertion(min_px2 == ARR_RIGHT_BOUNDARY);
      // The min points of the two curves lie on the right boundary.
      // It implies that the entire curves lie on the right boundary, and
      // thus both are vertical.

      typedef typename Self::Compare_y_on_boundary_2  Compare_y_on_boundary_2;
      Compare_y_on_boundary_2 cmp_y_on_bnd =
        m_self.compare_y_on_boundary_2_object();
      const Point_2& c1_min = m_self.construct_min_vertex_2_object()(c1);
      const Point_2& c2_min = m_self.construct_min_vertex_2_object()(c2);
      Comparison_result res = cmp_y_on_bnd(c1_min, c2_min);
      if (res != EQUAL) return res;

      const Point_2& c1_max = m_self.construct_max_vertex_2_object()(c1);
      const Point_2& c2_max = m_self.construct_max_vertex_2_object()(c2);
      res = cmp_y_on_bnd(c1_max, c2_max);
      return res;
    }
  };

  /*! Obtain a Compare_xy_2 function object */
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(*this); }

  /*! A functor that tests whether two x-monotone curves can be merged. */
  class Are_mergeable_2 {
  public:
    /*!
     * Check whether it is possible to merge two given x-monotone curves.
     * \param xcv1 The first curve.
     * \param xcv2 The second curve.
     * \return (true) if the two curves are mergeable - if they are supported
     *         by the same line and share a common endpoint; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2) const
    {
      // The function is implemented based on the Has_merge category.
      return (_are_mergeable_imp(xcv1, xcv2, Has_merge_category()));
    }

  protected:
    //! The base traits.
    const Base* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Are_mergeable_2(const Base* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_adaptor_2<Base_traits_2>;

    /*!
     * Implementation of the operator() in case the Has_merge tag is true.
     */
    bool _are_mergeable_imp(const X_monotone_curve_2& xcv1,
                             const X_monotone_curve_2& xcv2, Tag_true) const
    { return (m_base->are_mergeable_2_object()(xcv1, xcv2)); }

    /*!
     * Implementation of the operator() in case the Has_merge tag is false.
     */
    bool _are_mergeable_imp(const X_monotone_curve_2&,
                            const X_monotone_curve_2&, Tag_false) const
    {
      // Curve merging is not supported:
      return false;
    }
  };

  /*! Obtain an Are_mergeable_2 function object. */
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(this); }

  /*! A functor that merges two x-monotone curves into one. */
  class Merge_2 {
  public:
    /*!
     * Merge two given x-monotone curves into a single curve.
     * \param xcv1 The first curve.
     * \param xcv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable, that is they are supported by the
     *      curve line and share a common endpoint.
     */
    void operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2,
                    X_monotone_curve_2& c) const
    {
      // The function is implemented based on the Has_merge category.
      _merge_imp(xcv1, xcv2, c, Has_merge_category());
    }

  protected:
    //! The base traits.
    const Base* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle non
     *             stateless traits objects, (which stores data).
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Merge_2(const Base* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_traits_adaptor_2<Base_traits_2>;

    /*!
     * Implementation of the operator() in case the HasMerge tag is true.
     */
    void _merge_imp(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2,
                    X_monotone_curve_2& c, Tag_true) const
    { return (m_base->merge_2_object()(xcv1, xcv2, c)); }

    /*!
     * Implementation of the operator() in case the HasMerge tag is false.
     */
    void _merge_imp(const X_monotone_curve_2&, const X_monotone_curve_2&,
                    X_monotone_curve_2&, Tag_false) const
    {
      // This function should never be called!
      CGAL_error_msg( "Merging curves is not supported.");
      return;
    }
  };

  /*! Obtain a Merge_2 function object. */
  Merge_2 merge_2_object() const { return Merge_2(this); }
  //@}
};

} //namespace CGAL

#endif
