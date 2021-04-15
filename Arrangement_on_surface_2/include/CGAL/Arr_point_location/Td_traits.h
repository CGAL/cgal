// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)         : Oren Nechushtan <theoren@math.tau.ac.il>
#ifndef CGAL_TD_TRAITS_H
#define CGAL_TD_TRAITS_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <CGAL/Arr_point_location/Td_active_trapezoid.h>
#include <CGAL/Arr_point_location/Td_inactive_trapezoid.h>
#include <CGAL/Arr_point_location/Td_active_edge.h>
#include <CGAL/Arr_point_location/Td_inactive_edge.h>
#include <CGAL/Arr_point_location/Td_active_vertex.h>
#include <CGAL/Arr_point_location/Td_active_fictitious_vertex.h>
#include <CGAL/Arr_point_location/Td_inactive_vertex.h>
#include <CGAL/Arr_point_location/Td_inactive_fictitious_vertex.h>

namespace CGAL {

template <class Pm_traits_,class Arrangement_>
class Td_traits : public Pm_traits_
{
public:

  //type of td map items type
  enum Type
  {
      NIL = 0,
      TD_ACTIVE_TRAPEZOID,
      TD_INACTIVE_TRAPEZOID,
      TD_ACTIVE_EDGE,
      TD_INACTIVE_EDGE,
      TD_ACTIVE_VERTEX,
      TD_ACTIVE_FICTITIOUS_VERTEX,
      TD_INACTIVE_VERTEX,
      TD_INACTIVE_FICTITIOUS_VERTEX
  };

  //! type of base class
  typedef Pm_traits_                      Traits_base;

  //! type of X_monotone_curve_2
  typedef typename Traits_base::X_monotone_curve_2
                                          X_monotone_curve_2;

  //! type of Arrangement_on_surface_2
  typedef Arrangement_                    Arrangement_on_surface_2;

  //!type of Halfedge_handle
  typedef typename Arrangement_on_surface_2::Halfedge_handle
                                          Halfedge_handle;
  //!type of Halfedge_const_handle
  typedef typename Arrangement_on_surface_2::Halfedge_const_handle
                                          Halfedge_const_handle;
  //!type of Vertex_const_handle
  typedef typename Arrangement_on_surface_2::Vertex_const_handle
                                          Vertex_const_handle;
  //!type of Halfedge_around_vertex_const_circulator
  typedef typename
    Arrangement_on_surface_2::Halfedge_around_vertex_const_circulator
      Halfedge_around_vertex_const_circulator;
  //!type of side tags
  typedef typename Arrangement_on_surface_2::Left_side_category
                                          Left_side_category;
  typedef typename Arrangement_on_surface_2::Bottom_side_category
                                          Bottom_side_category;
  typedef typename Arrangement_on_surface_2::Top_side_category
                                          Top_side_category;
  typedef typename Arrangement_on_surface_2::Right_side_category
                                          Right_side_category;

  //! myself
  typedef Td_traits<Traits_base,Arrangement_on_surface_2>
                                          Self;
  //! type of point
  typedef typename Traits_base::Point_2   Point;

  //! type of Td_active_trapezoid
  typedef CGAL::Td_active_trapezoid<Self> Td_active_trapezoid;

  //! type of Td_inactive_trapezoid
  typedef CGAL::Td_inactive_trapezoid     Td_inactive_trapezoid;

  typedef int                             Td_nothing;

  //! type of Td_active_edge
  typedef CGAL::Td_active_edge<Self>      Td_active_edge;

  //! type of Td_inactive_edge
  typedef CGAL::Td_inactive_edge<Self>    Td_inactive_edge;

  //! type of Td_active_vertex
  typedef CGAL::Td_active_vertex<Self>    Td_active_vertex;

  //! type of Td_active_fictitious_vertex
  typedef CGAL::Td_active_fictitious_vertex<Self>
                                          Td_active_fictitious_vertex;

  //! type of Td_inactive_vertex
  typedef CGAL::Td_inactive_vertex<Self>  Td_inactive_vertex;

  //! type of Td_inactive_fictitious_vertex
  typedef CGAL::Td_inactive_fictitious_vertex<Self>
                                          Td_inactive_fictitious_vertex;

  //! type of td map item (Td_halfedge, Td_vertex or Td_trapezoid)
  typedef boost::variant< Td_nothing,
                          Td_active_trapezoid, Td_inactive_trapezoid,
                          Td_active_edge, Td_inactive_edge,
                          Td_active_vertex, Td_active_fictitious_vertex,
                          Td_inactive_vertex, Td_inactive_fictitious_vertex >  Td_map_item;

    //! type of Curve end pair
  typedef std::pair<const X_monotone_curve_2*, Arr_curve_end>
                                          Curve_end_pair;


  //!Curve_end class represents an X_monotone_curve_2 end
  //  (could be a point or an unbounded curve end)
  //  holds a pointer to the X_monotone_curve_2 and an indicator for the end (min/max)
  class Curve_end
  {
  protected:

    //! pair of pointer to the X_monotone_curve_2 and an indicator
    //    for ARR_MIN_END or ARR_MAX_END
    Curve_end_pair m_pair;

  public:

    //Constructor based on a Curve_end_pair
    Curve_end(Curve_end_pair pr) : m_pair(pr)
    { }

    //Constructor based on a Curve & a Curve-end
    //Curve_end(const X_monotone_curve_2& cv, Arr_curve_end ce) : m_cv(cv), m_ce(ce)
    Curve_end(const X_monotone_curve_2& cv, Arr_curve_end ce)
      : m_pair(std::make_pair(&cv,ce))
    { }

    //Constructor based on a Halfedge & a Curve-end
    //Curve_end(Halfedge_const_handle he, Arr_curve_end ce) : m_cv(he->curve()), m_ce(ce)
    Curve_end(Halfedge_const_handle he, Arr_curve_end ce)
      : m_pair(std::make_pair(&(he->curve()),ce))
    { }

    //access the X-monotone curve
    const X_monotone_curve_2&  cv() const  { return *(m_pair.first);  }

    //access the Curve-end
    Arr_curve_end   ce() const  { return m_pair.second;  }

  };


   /*! A functor that compares the x-coordinates of two edge ends
   */
  class Compare_curve_end_x_2 {
  protected:
    typedef Traits_base Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_curve_end_x_2(const Traits* traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Td_traits<Traits_base, Arrangement_on_surface_2>;

  public:

    /*!
     * Compare the x-coordinates of two edge ends.
     * \param ee1 The first edge end
     * \param ee2 The second edge end
     * \return LARGER if x(ee1) > x(ee2);
     *         SMALLER if x(ee1) < x(ee2);
     *         EQUAL if x(ee1) = x(ee2).
     */
    Comparison_result operator() (const Curve_end& ce1,
                                  const Curve_end& ce2) const
    {
      Arr_parameter_space ce1_x_prm_spc =
        m_traits->parameter_space_in_x_2_object()(ce1.cv(),ce1.ce());
      Arr_parameter_space ce1_y_prm_spc =
        m_traits->parameter_space_in_y_2_object()(ce1.cv(),ce1.ce());
      Arr_parameter_space ce2_x_prm_spc =
        m_traits->parameter_space_in_x_2_object()(ce2.cv(),ce2.ce());
      Arr_parameter_space ce2_y_prm_spc =
        m_traits->parameter_space_in_y_2_object()(ce2.cv(),ce2.ce());

      bool is_ce1_interior = (( ce1_x_prm_spc == ARR_INTERIOR) &&
                              ( ce1_y_prm_spc == ARR_INTERIOR));
      bool is_ce2_interior = (( ce2_x_prm_spc == ARR_INTERIOR) &&
                              ( ce2_y_prm_spc == ARR_INTERIOR));

      //if both are interior
      if (is_ce1_interior && is_ce2_interior)
      {
        return m_traits->compare_x_2_object()
                   ( ((ce1.ce() == ARR_MIN_END) ?
                      m_traits->construct_min_vertex_2_object()(ce1.cv()) :
                      m_traits->construct_max_vertex_2_object()(ce1.cv())  ),
                     ((ce2.ce() == ARR_MIN_END) ?
                      m_traits->construct_min_vertex_2_object()(ce2.cv()) :
                      m_traits->construct_max_vertex_2_object()(ce2.cv())  ));
      }

      //if only ce1 is interior
      if (is_ce1_interior)
        return operator()
                  (((ce1.ce() == ARR_MIN_END) ?
                      m_traits->construct_min_vertex_2_object()(ce1.cv()) :
                      m_traits->construct_max_vertex_2_object()(ce1.cv())  ),
                    ce2);

      //if only ce2 is interior
      if (is_ce2_interior)
        return operator()
                  (ce1,
                   ((ce2.ce() == ARR_MIN_END) ?
                      m_traits->construct_min_vertex_2_object()(ce2.cv()) :
                      m_traits->construct_max_vertex_2_object()(ce2.cv()) ));

      //both are not interior:

      //if both are x-interior
      //   (then both are NOT y-interior, since both are not interior)
      if ( (ce1_x_prm_spc == ARR_INTERIOR) &&
           (ce2_x_prm_spc == ARR_INTERIOR)  )
      {
        //both ce1 and ce2 are not y prm spc interior
        Comparison_result res = m_traits->compare_x_at_limit_2_object()
                                            (ce1.cv(),ce1.ce(),
                                             ce2.cv(),ce2.ce());
        if (res != EQUAL)
          return res;
        //if equal need to compare near limit

        //if param space in y is not the same (one is top and one is bottom)
        //  the bottom is smaller than the top
        if (ce1_y_prm_spc != ce2_y_prm_spc)
          return (ce1_y_prm_spc == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;

        //if the Curve end is not the same the one with the MAX is smaller
        if (ce1.ce() != ce2.ce())
          return (ce1.ce() == ARR_MIN_END) ? LARGER : SMALLER;

        //both have the same Curve end
        return (m_traits->compare_x_near_limit_2_object()
                                              (ce1.cv(),ce2.cv(),ce2.ce()));
      }

      //not both are x-interior

      //set ind value according to the location :
      //  left-bndry  = -1
      //  interior    =  0
      //  right-bndry =  1
      // if (ind1 - ind2) ==0 ->EQUAL
      // if (ind1 - ind2) < 0 ->SMALLER
      // if (ind1 - ind2) > 0 ->LARGER

      int ind1 = (ce1_x_prm_spc == ARR_INTERIOR)? 0
                  : ((ce1_x_prm_spc == ARR_LEFT_BOUNDARY) ? -1 : 1 );

      int ind2 = (ce2_x_prm_spc == ARR_INTERIOR)? 0
                  : ((ce2_x_prm_spc == ARR_LEFT_BOUNDARY) ? -1 : 1 );

      int res = ind1 - ind2;
      if (res == 0)
        return EQUAL;
      return ((res < 0) ? SMALLER : LARGER);
    }

    /*!
     * Compare the x-coordinates of a point and a curve end.
     * \param p The point
     * \param ce The curve end
     * \return LARGER if x(p) > x(ce);
     *         SMALLER if x(p) < x(ce);
     *         EQUAL if x(p) = x(ce).
     */
    Comparison_result operator() (const Point& p,
                                  const Curve_end& ce) const
    {
      bool is_ce_interior =
              ((m_traits->parameter_space_in_x_2_object()
                           (ce.cv(),ce.ce()) == ARR_INTERIOR)   &&
               (m_traits->parameter_space_in_y_2_object()
                           (ce.cv(),ce.ce()) == ARR_INTERIOR));

      //if curve end is interior
      if (is_ce_interior)
      {
        return m_traits->compare_x_2_object()
                   ( p,
                     ((ce.ce() == ARR_MIN_END) ?
                      m_traits->construct_min_vertex_2_object()(ce.cv()) :
                      m_traits->construct_max_vertex_2_object()(ce.cv())  ));
      }

      //if curve end is x-interior but not y-interior
      if (m_traits->parameter_space_in_x_2_object()
                              (ce.cv(),ce.ce()) == ARR_INTERIOR)
      {
          //if curve end y prm space is not interior
          return (m_traits->compare_x_at_limit_2_object()
                                (p,ce.cv(),ce.ce()));
      }

      //if curve end is on the left or right boundaries
      if (m_traits->parameter_space_in_x_2_object()
                            (ce.cv(),ce.ce()) == ARR_LEFT_BOUNDARY)
      {
        return LARGER;
      }
      return SMALLER;
    }

    /*!
     * Compare the x-coordinates of an curve end and a point
     * \param ce The curve end
     * \param p The point
     * \return LARGER if x(ce) > x(p);
     *         SMALLER if x(ce) < x(p);
     *         EQUAL if x(ce) = x(p).
     */
    Comparison_result operator() (const Curve_end& ce,
                                  const Point& p) const
    {
      Comparison_result res = operator()(p, ce);
      if (res == EQUAL)
        return res;
      return (res == LARGER) ? SMALLER : LARGER;
    }

  };

  /*! Obtain a Compare_y_at_x_2 functor object. */
  Compare_curve_end_x_2 compare_curve_end_x_2_object () const
  {
    return Compare_curve_end_x_2(this);
  }




 /*! A functor that compares the y-coordinates of an edge end and a curve at
   * the point x-coordinate
   */
  class Compare_curve_end_y_at_x_2
  {

  protected:

    typedef Td_traits<Traits_base, Arrangement_on_surface_2> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_curve_end_y_at_x_2(const Traits * traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Td_traits<Traits_base, Arrangement_on_surface_2>;

  public:

    /*!
     * Return the location of the given curve end with respect to the input
     *  Halfedge_const_handle.
     * \param ce1 The curve end.
     * \param he  The Halfedge_const_handle.
     * \pre ce1 is in the x-range of he.
     * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    Comparison_result operator() (const Curve_end& ce1,
                                  Halfedge_const_handle  he) const
    {
      CGAL_precondition_code(Halfedge_const_handle invalid_he);
      CGAL_precondition(he != invalid_he);
      return operator()(ce1,he->curve());
    }

    /*!
     * Return the location of the given curve end with respect to the input cv.
     * \param ce1 The curve end.
     * \param cv  The X_monotone_curve_2.
     * \pre ce1 is in the x-range of cv.
     * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    Comparison_result operator() (const Curve_end& ce1,
                                  const X_monotone_curve_2&  cv2) const
    {
      //precondition: ce1 is in the x-range of cv2
      CGAL_precondition (
        (m_traits->compare_curve_end_x_2_object()
                    (ce1, Curve_end(cv2, ARR_MIN_END)) != SMALLER)
           &&
        (m_traits->compare_curve_end_x_2_object()
                    (ce1, Curve_end(cv2, ARR_MAX_END)) != LARGER));

      //get the curve end parameter space in x & y
      Arr_parameter_space ce1_x_prm_spc =
               m_traits->parameter_space_in_x_2_object()(ce1.cv(), ce1.ce());

      Arr_parameter_space ce1_y_prm_spc =
               m_traits->parameter_space_in_y_2_object()(ce1.cv(), ce1.ce());


      if (ce1_x_prm_spc != ARR_INTERIOR)
      {
        //assuming that the edge end is on the same boundary according to
        // the precondition.
        //comparing the curve that contains the given
        //  edge end  and the curve cv2
        return m_traits->compare_y_near_boundary_2_object()
                          (ce1.cv(), cv2, ce1.ce());
      }

      //if ce1_x_prm_spc == ARR_INTERIOR
      if (ce1_y_prm_spc == ARR_INTERIOR)
      {
        //ce1 is interior
        return m_traits->compare_y_at_x_2_object()
                         (((ce1.ce() == ARR_MIN_END) ?
                           m_traits->construct_min_vertex_2_object()(ce1.cv()) :
                           m_traits->construct_max_vertex_2_object()(ce1.cv()) ),
                           cv2);
      }

      //ce1 is an end point of a curve with a vertical asymptote.

      //if the other curve is also vertical or has a vertical asymptote
      //  at the x value of ep

      if ( ((m_traits->parameter_space_in_x_2_object()
                         (cv2, ARR_MIN_END) == ARR_INTERIOR) &&
            (m_traits->parameter_space_in_y_2_object()
                         (cv2, ARR_MIN_END) == ce1_y_prm_spc) &&
            (m_traits->compare_curve_end_x_2_object()
                         (ce1, Curve_end(cv2, ARR_MIN_END)) == EQUAL)) ||
           ((m_traits->parameter_space_in_x_2_object()
                         (cv2, ARR_MAX_END) == ARR_INTERIOR) &&
            (m_traits->parameter_space_in_y_2_object()
                         (cv2, ARR_MAX_END) == ce1_y_prm_spc) &&
            (m_traits->compare_curve_end_x_2_object()
                         (ce1, Curve_end(cv2, ARR_MAX_END)) == EQUAL))  )
      {
        return EQUAL;
      }
      if (ce1_y_prm_spc == ARR_TOP_BOUNDARY)
        return LARGER;
      else //if ce1_y_prm_spc == ARR_BOTTOM_BOUNDARY
        return SMALLER;


    }

  };

  /*! Obtain a Compare_y_at_x_2 functor object. */
  Compare_curve_end_y_at_x_2 compare_curve_end_y_at_x_2_object () const
  {
    return Compare_curve_end_y_at_x_2(this);
  }



  class Equal_curve_end_2
  {
  protected:
    typedef Td_traits<Traits_base, Arrangement_on_surface_2>  Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;
    const Traits_base* m_traits_base; //MICHAL: rational-upd

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    //Equal_curve_end_2(const Traits* traits) : m_traits(traits) {}//MICHAL: rational-upd
    Equal_curve_end_2(const Traits* traits) : m_traits(traits), m_traits_base(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Td_traits<Traits_base, Arrangement_on_surface_2>;

  public:

    bool operator() (const Curve_end& ce1,
                     const Curve_end& ce2) const
    {
      //Kernel kernel; //MICHAL: rational-upd

      bool is_ce1_interior =
              ((m_traits->parameter_space_in_x_2_object()(ce1.cv(),ce1.ce())
                                              == ARR_INTERIOR)      &&
               (m_traits->parameter_space_in_y_2_object()(ce1.cv(),ce1.ce())
                                              == ARR_INTERIOR));
      bool is_ce2_interior =
              ((m_traits->parameter_space_in_x_2_object()(ce2.cv(),ce2.ce())
                                              == ARR_INTERIOR)      &&
               (m_traits->parameter_space_in_y_2_object()(ce2.cv(),ce2.ce())
                                              == ARR_INTERIOR));

      if (is_ce1_interior && is_ce2_interior) //both edge-ends are interior
  {
        return m_traits_base->equal_2_object()
                  ( ((ce1.ce() == ARR_MIN_END) ?
                       m_traits->construct_min_vertex_2_object()(ce1.cv()) :
                       m_traits->construct_max_vertex_2_object()(ce1.cv())  ),
                    ((ce2.ce() == ARR_MIN_END) ?
                       m_traits->construct_min_vertex_2_object()(ce2.cv()) :
                       m_traits->construct_max_vertex_2_object()(ce2.cv())  ));
    }


      //at least one of the edge ends is on the parameter space boundaries

      //if not both are on the boundaries return false
      if (is_ce1_interior || is_ce2_interior)
        return false;

      //both are on the boundaries - so compare the edge ends
      return ( m_traits->compare_curve_end_xy_2_object()(ce1,ce2) == EQUAL);
    }

    bool operator() (const Point& p1,
                     const Point& p2) const
    {
      return m_traits_base->equal_2_object()(p1, p2);
    }

    bool operator() (const Curve_end& ce,
                     const Point& p) const
    {
      return operator()(p, ce);
    }

    bool operator() (const Point& p,
                     const Curve_end& ce) const
    {
      bool is_ce_interior =
              ((m_traits->parameter_space_in_x_2_object()(ce.cv(),ce.ce())
                                              == ARR_INTERIOR)      &&
               (m_traits->parameter_space_in_y_2_object()(ce.cv(),ce.ce())
                                              == ARR_INTERIOR));

      //if ce is on the parameter space boundaries - return false
      // since p is interior
      if (!is_ce_interior)
        return false;

      //else - if ce is interior
      return m_traits_base->equal_2_object()
             ( p,
               ((ce.ce() == ARR_MIN_END) ?
                  m_traits->construct_min_vertex_2_object()(ce.cv()) :
                  m_traits->construct_max_vertex_2_object()(ce.cv())  ));
    }

  };

  /*! Obtain an Equal_curve_end_2 functor object. */
  Equal_curve_end_2 equal_curve_end_2_object () const
  {
    return Equal_curve_end_2(this);
  }

  /*! A functor that compares the coordinates of two edge ends */
  class Compare_curve_end_xy_2
  {
  protected:
    typedef Td_traits<Traits_base, Arrangement_on_surface_2>  Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_curve_end_xy_2(const Traits* traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Td_traits<Traits_base, Arrangement_on_surface_2>;

  public:
    /*!
     * Compare two edge ends lexigoraphically: by x, then by y.
     * \param cv1, cv1_end The first cv end.
     * \param cv2, cv2_end The second cv end.
     * \return LARGER if x(cv1-end) > x(cv2-end),
     *             or if x(cv1-end) = x(cv2-end) and y(cv1-end) > y(cv2-end);
     *         SMALLER if x(cv1-end) < x(cv2-end),
     *             or if x(cv1-end) = x(cv2-end) and y(cv1-end) < y(cv2-end);
     *         EQUAL if the two cv ends are equal.
     */
    Comparison_result operator() (const Curve_end& ce1,
                                  const Curve_end& ce2) const
    {
      Comparison_result res;

      bool is_ce1_interior =
              ((m_traits->parameter_space_in_x_2_object()(ce1.cv(),ce1.ce())
                                              == ARR_INTERIOR)       &&
               (m_traits->parameter_space_in_y_2_object()(ce1.cv(),ce1.ce())
                                              == ARR_INTERIOR));
      bool is_ce2_interior =
              ((m_traits->parameter_space_in_x_2_object()(ce2.cv(),ce2.ce())
                                              == ARR_INTERIOR)       &&
               (m_traits->parameter_space_in_y_2_object()(ce2.cv(),ce2.ce())
                                              == ARR_INTERIOR));

      bool is_ce1_vertical =
              ((m_traits->parameter_space_in_x_2_object()(ce1.cv(),ce1.ce())
                                              == ARR_INTERIOR)       &&
               (m_traits->parameter_space_in_y_2_object()(ce1.cv(),ce1.ce())
                                              != ARR_INTERIOR));
      bool is_ce2_vertical =
              ((m_traits->parameter_space_in_x_2_object()(ce2.cv(),ce2.ce())
                                              == ARR_INTERIOR)       &&
               (m_traits->parameter_space_in_y_2_object()(ce2.cv(),ce2.ce())
                                              != ARR_INTERIOR));

      //if the edge ends are parameter space interior on both x & y
      if ( is_ce1_interior && is_ce2_interior )
      {
        return m_traits->compare_xy_2_object()
                  ( ((ce1.ce() == ARR_MIN_END) ?
                      m_traits->construct_min_vertex_2_object()(ce1.cv()) :
                      m_traits->construct_max_vertex_2_object()(ce1.cv())  ),
                    ((ce2.ce() == ARR_MIN_END) ?
                      m_traits->construct_min_vertex_2_object()(ce2.cv()) :
                      m_traits->construct_max_vertex_2_object()(ce2.cv())  ));
      }

      //at least one curve end is on the parameter space boundaries:

      //if the first curve end is interior
      if ( is_ce1_interior )
      {

        //if the second curve end is of a curve with a vertical asymptote:
        //  x prm spc is interior, y prm spc is not interior
        if ( is_ce2_vertical)
        {
          //res = m_traits->compare_x_near_boundary_2_object()
          res = m_traits->compare_x_at_limit_2_object()
                     (((ce1.ce() == ARR_MIN_END) ?
                       m_traits->construct_min_vertex_2_object()(ce1.cv()) :
                       m_traits->construct_max_vertex_2_object()(ce1.cv())  ),
                      ce2.cv(), ce2.ce());

          if (res != EQUAL)
            return res;
          else
            return (m_traits->parameter_space_in_y_2_object()
                   (ce2.cv(),ce2.ce()) == ARR_TOP_BOUNDARY) ? SMALLER : LARGER;
        }
        else //the second curve end is of an unbounded cv which is not vertical
        {
          //we compare an interior curve end to a curve end on the left/right
          // boundaries. the comparison is simply by the x-coord
          return (ce2.ce() == ARR_MIN_END) ? LARGER : SMALLER;
        }

      }

      //if the second curve end is interior
      if ( is_ce2_interior )
      {
         //if the first curve end is of a vertical line:
        //  x prm spc is interior, y prm spc is not interior
        if ( is_ce1_vertical )
        {
          //res = m_traits->compare_x_near_boundary_2_object()
          res = m_traits->compare_x_at_limit_2_object()
                     ((( ce2.ce() == ARR_MIN_END) ?
                       m_traits->construct_min_vertex_2_object()( ce2.cv()) :
                       m_traits->construct_max_vertex_2_object()( ce2.cv())  ),
                       ce1.cv(),  ce1.ce());
          //need to return the opposite because the function received
          // the curve ends in a reverse order
          if (res != EQUAL)
            return (res == SMALLER) ? LARGER : SMALLER;
          else
            return (m_traits->parameter_space_in_y_2_object()
                        (ce1.cv(),ce1.ce()) ==  ARR_TOP_BOUNDARY) ? LARGER : SMALLER;
        }
        else //the first curve end is of an unbounded cv which is not vertical
        {
          //we compare an interior curve end to an curve end on the left/right
          // boundaries. the comparison is simply by the x-coord
          return (ce1.ce() == ARR_MIN_END) ? SMALLER : LARGER;
        }

      }

      //both curve ends are not interior

      //if both curve ends are of unbounded curves with a vertical asymptote
      if ( is_ce1_vertical && is_ce2_vertical )
      {
        Comparison_result res = m_traits->compare_x_at_limit_2_object()
                                            (ce1.cv(),ce1.ce(),
                                             ce2.cv(),ce2.ce());

        if (res != EQUAL)
          return res;

        //res == EQUAL
        //if equal - need to compare near limit

        Arr_parameter_space ce1_y_prm_spc =
          m_traits->parameter_space_in_y_2_object()(ce1.cv(),ce1.ce()) ;
        Arr_parameter_space ce2_y_prm_spc =
          m_traits->parameter_space_in_y_2_object()(ce2.cv(),ce2.ce()) ;

        //if param space in y is not the same (one is top and one is bottom)
        //  the bottom is smaller than the top
        if (ce1_y_prm_spc != ce2_y_prm_spc)
          return (ce1_y_prm_spc == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;

        //if the Curve end is not the same, the one with the MAX is smaller
        if (ce1.ce() != ce2.ce())
          return (ce1.ce() == ARR_MIN_END) ? LARGER : SMALLER;

        //both have the same Curve end
        return (m_traits->compare_x_near_limit_2_object()
                                          (ce1.cv(),ce2.cv(),ce2.ce()));
      }

      //if only the first curve end is of a curve with a vertical asymptote
      if ( is_ce1_vertical )
      {
        return (ce2.ce() == ARR_MIN_END) ? LARGER : SMALLER;
      }
      //if only the second curve end is of a curve with a vertical asymptote
      if ( is_ce2_vertical )
      {
        return (ce1.ce() == ARR_MIN_END) ? SMALLER : LARGER;
      }

      //both curve ends are not of curves with a vertical asymptote:

      //if not both on left or both on right boundaries
      if (ce1.ce() != ce2.ce())
      {
        return (ce1.ce() == ARR_MIN_END) ? SMALLER : LARGER;
      }

      //both on the same boundary, need to compare the y near the boundary
      //  (ce1.ce() == ce2.ce())
      return m_traits->compare_y_near_boundary_2_object()
                            (ce1.cv(), ce2.cv(), ce1.ce());
    }

    /*!
     * Compare a point and a curve end lexigoraphically: by x, then by y.
     * \param cv1, cv1_end The first cv end.
     * \param cv2, cv2_end The second cv end.
     * \return LARGER if x(p) > x(cv2-end),
     *             or if x(p) = x(cv2-end) and y(p) > y(cv2-end);
     *         SMALLER if x(p) < x(cv2-end),
     *             or if x(p) = x(cv2-end) and y(p) < y(cv2-end);
     *         EQUAL if the point and the cv end are equal.
     */
    Comparison_result operator() (const Point& p,
                                  const Curve_end& ce) const
    {
      Comparison_result res;

      bool is_ce_interior =
              ((m_traits->parameter_space_in_x_2_object()(ce.cv(),ce.ce())
                                              == ARR_INTERIOR)       &&
               (m_traits->parameter_space_in_y_2_object()(ce.cv(),ce.ce())
                                              == ARR_INTERIOR));

      bool is_ce_vertical =
              ((m_traits->parameter_space_in_x_2_object()(ce.cv(),ce.ce())
                                              == ARR_INTERIOR)       &&
               (m_traits->parameter_space_in_y_2_object()(ce.cv(),ce.ce())
                                              != ARR_INTERIOR));

      //if the edge end is parameter space interior on both x & y
      if ( is_ce_interior)
      {
        return m_traits->compare_xy_2_object()
                  ( p,
                    ((ce.ce() == ARR_MIN_END) ?
                      m_traits->construct_min_vertex_2_object()(ce.cv()) :
                      m_traits->construct_max_vertex_2_object()(ce.cv())  ));
      }

      // edge end is on the parameter space boundaries:

      //if edge end is of a vertical line:
      //  x prm spc is interior, y prm spc is not interior
      if ( is_ce_vertical)
      {
        res = m_traits->compare_x_at_limit_2_object()
                           (p, ce.cv(), ce.ce());

        if (res != EQUAL)
          return res;
        else
          return (m_traits->parameter_space_in_y_2_object()
                    (ce.cv(), ce.ce()) == ARR_TOP_BOUNDARY) ? SMALLER : LARGER;
      }
      else //edge end is of an unbounded cv which is not vertical
      {
          //we compare an interior point to an edge end on the left/right
          // boundaries. the comparison is simply by the x-coord
          return (ce.ce() == ARR_MIN_END) ? LARGER : SMALLER;
      }
    }

    Comparison_result operator() (const Curve_end& ce,
                                  const Point& p) const
    {
      Comparison_result res = operator()(p,ce);
      if (res == EQUAL)
        return res;
      return (res == SMALLER) ? LARGER : SMALLER;
    }


    Comparison_result operator() (const Point& p1,
                                  const Point& p2) const
    {
      return m_traits->compare_xy_2_object()(p1,p2);
    }
  };

  /*! Obtain a Compare_curve_end_xy_2 functor object. */
  Compare_curve_end_xy_2 compare_curve_end_xy_2_object () const
  {
    return Compare_curve_end_xy_2(this);
  }




  // Td_traits class ctors and dtor

  Td_traits(const Traits_base& t) : Traits_base(t)
  { }

  Td_traits()
  { }

  ~Td_traits(void)
  {
  }

public:
  /*
    note:
    The traits assume that the trapezoid is active,non empty,
    and planar, that is no two curves intersect in non degenerate curve.
  */

  /* returns true if bottom halfedges of input are the same */
  inline bool is_trpz_bottom_equal(Td_map_item& left_item,
                                                                 Td_map_item& right_item) const
  {
    CGAL_precondition(is_active(left_item) && is_active(right_item));
    CGAL_precondition(is_td_trapezoid(left_item) && is_td_trapezoid(right_item));

    Td_active_trapezoid left (boost::get<Td_active_trapezoid>(left_item));
    Td_active_trapezoid right(boost::get<Td_active_trapezoid>(right_item));

    if (left.is_on_bottom_boundary())
      return (right.is_on_bottom_boundary());

    if (right.is_on_bottom_boundary())
      return (false);

    return (left.bottom() == right.bottom() ||
            left.bottom()->twin() == right.bottom());
  }

  /* returns true if top halfedges of input are the same */
  inline bool is_trpz_top_equal(Td_map_item& left_item,
                                                              Td_map_item& right_item) const
  {
    CGAL_precondition(is_active(left_item) && is_active(right_item));
    CGAL_precondition(is_td_trapezoid(left_item) && is_td_trapezoid(right_item));

    Td_active_trapezoid left (boost::get<Td_active_trapezoid>(left_item));
    Td_active_trapezoid right(boost::get<Td_active_trapezoid>(right_item));

    if (left.is_on_top_boundary())
      return (right.is_on_top_boundary());

    if (right.is_on_top_boundary())
      return (false);

    return (left.top() == right.top() || left.top()->twin() == right.top());
  }

  /* returns true if bottom halfedges of input are the same */
  inline bool is_trapezoids_bottom_equal(const Td_active_trapezoid& left,
                                                                       const Td_active_trapezoid& right) const
  {
    if (left.is_on_bottom_boundary())
      return (right.is_on_bottom_boundary());

    if (right.is_on_bottom_boundary())
      return (false);

    return (left.bottom() == right.bottom() ||
            left.bottom()->twin() == right.bottom());
  }

  /* returns true if top halfedges of input are the same */
  inline bool is_trapezoids_top_equal(const Td_active_trapezoid& left,
                                                                    const Td_active_trapezoid& right) const
  {
    if (left.is_on_top_boundary())
      return (right.is_on_top_boundary());

    if (right.is_on_top_boundary())
      return (false);

    return (left.top() == right.top() || left.top()->twin() == right.top());
  }

  //returns true if the trapezoid is a curve
  bool is_empty_item(const Td_map_item& tr) const
  {
    return (tr.which() == 0);
  }

  //returns true if the trapezoid is a point or a curve
  bool is_trapezoid(const Td_map_item& tr) const
  {
    switch (tr.which())
    {
    case TD_ACTIVE_TRAPEZOID:
    case TD_INACTIVE_TRAPEZOID:
      return true;
    default:
      return false;
    }
  }

  //returns true if the map item is a vertex
  bool is_td_vertex(const Td_map_item& tr) const
  {
    switch (tr.which())
    {
    case TD_ACTIVE_VERTEX:
    case TD_ACTIVE_FICTITIOUS_VERTEX:
    case TD_INACTIVE_VERTEX:
    case TD_INACTIVE_FICTITIOUS_VERTEX:
      return true;
    default:
      return false;
    }
  }

  //returns true if the map item is an edge
  bool is_td_edge(const Td_map_item& tr) const
  {
    switch (tr.which())
    {
    case TD_ACTIVE_EDGE:
    case TD_INACTIVE_EDGE:
      return true;
    default:
      return false;
    }
  }

  //returns true if the map item is an edge
  bool is_td_trapezoid(const Td_map_item& tr) const
  {
    switch (tr.which())
    {
    case TD_ACTIVE_TRAPEZOID:
    case TD_INACTIVE_TRAPEZOID:
      return true;
    default:
      return false;
    }
  }

  //returns true if the trapezoid is a curve
  bool is_fictitious_vertex(const Td_map_item& tr) const
  {
    switch (tr.which())
    {
    case TD_ACTIVE_FICTITIOUS_VERTEX:
    case TD_INACTIVE_FICTITIOUS_VERTEX:
      return true;
    default:
      return false;
    }
  }

  //returns true if the trapezoid is a curve
  bool is_active(const Td_map_item& tr) const
  {
    switch (tr.which())
    {
      case TD_ACTIVE_TRAPEZOID:
      case TD_ACTIVE_EDGE:
      case TD_ACTIVE_VERTEX:
      case TD_ACTIVE_FICTITIOUS_VERTEX:
        return true;
    default:
      return false;
    }
  }

  //returns true if the trapezoid is vertical
  bool is_vertical(Td_map_item& item) const
  {
    CGAL_precondition(is_td_edge(item));
    CGAL_precondition(is_active(item));
    //MICHAL: assumes item is of active edge item - also fails in case of a vertical asymptote
    //MICHAL: check when this is used exactly
    Td_active_edge& e (boost::get<Td_active_edge>(item));
    Halfedge_const_handle he = e.halfedge();
    return (this->compare_curve_end_x_2_object()
               (Curve_end(he,ARR_MIN_END), Curve_end(he,ARR_MAX_END))== EQUAL);
  }

  /* returns whether given edge end is inside the given trapezoid using
    lexicographic order */
  bool is_inside  (Td_map_item& item, const Curve_end& ce) const
  {
    CGAL_precondition( is_active(item) );
    CGAL_precondition( is_td_trapezoid(item) );
    Td_active_trapezoid tr (boost::get<Td_active_trapezoid>(item));

    return
      ( tr.is_on_left_boundary() ||
          (compare_curve_end_xy_2_object()
                    (vtx_to_ce(tr.left()),ce) == SMALLER) )  &&
       ( tr.is_on_right_boundary() ||
          (compare_curve_end_xy_2_object()
                    (vtx_to_ce(tr.right()),ce) == LARGER) )  &&
       ( tr.is_on_bottom_boundary() ||
          (compare_curve_end_y_at_x_2_object()(ce, tr.bottom()) == LARGER) ) &&
       ( tr.is_on_top_boundary() ||
          (compare_curve_end_y_at_x_2_object()(ce, tr.top()) == SMALLER) );
  }

  /*! returns true if the end point is inside the closure of the trapezoid
      (inlcude all boundaries) */
  bool is_in_closure  (const Td_active_trapezoid& tr, const Curve_end& ce ) const
  {
    // test left and right sides
    if ((tr.is_on_left_boundary()   ||
         (compare_curve_end_xy_2_object()
                   (ce,vtx_to_ce(tr.left())) != SMALLER))  &&
        (tr.is_on_right_boundary()  ||
         (compare_curve_end_xy_2_object()
                   (ce,vtx_to_ce(tr.right())) != LARGER))  )
    {
      // test bottom side
      if (!tr.is_on_bottom_boundary() &&
          compare_curve_end_y_at_x_2_object()(ce,tr.bottom()) == SMALLER )
      {
        return false;
      }

      // test top side
      if (!tr.is_on_top_boundary() &&
                compare_curve_end_y_at_x_2_object()(ce,tr.top()) == LARGER)
      {
        return false;
      }

      return true;
    }
    return false;
  }
  /*! returns true if the end point is inside the closure of the trapezoid
      (inlcude all boundaries) */
  bool is_in_closure (const Td_active_edge& e, const Curve_end& ce ) const
  {
    // test left and right sides
    if (compare_curve_end_xy_2_object()(ce,Curve_end(e.halfedge(),ARR_MIN_END)) == SMALLER)
      return false;
    if (compare_curve_end_xy_2_object()(ce,Curve_end(e.halfedge(),ARR_MAX_END)) == LARGER)
      return false;
    return true;
  }


  Curve_end vtx_to_ce(Vertex_const_handle v) const
  {
    //the circulator is of incoming halfedges
    Halfedge_around_vertex_const_circulator he = v->incident_halfedges();
    //if the vertex is associated with a point on the bounded coords,
    // we can take any incident halfedge. o/w if the vertex lies at infinity,
    //  it has 2 fictitious incident halfedges
    if (v->is_at_open_boundary() && he->source()->is_at_open_boundary()) ++he;
    if (v->is_at_open_boundary() && he->source()->is_at_open_boundary()) ++he;

    return Curve_end(he->curve(),
                     (he->direction() == ARR_RIGHT_TO_LEFT)?
                      ARR_MIN_END : ARR_MAX_END);
  }

public:

  static Vertex_const_handle empty_vtx_handle() {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(Vertex_const_handle, m_empty_vtx_handle);
    return m_empty_vtx_handle;
  }

  static Halfedge_const_handle empty_he_handle() {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(Halfedge_const_handle, m_empty_he_handle);
    return m_empty_he_handle;
  }

};

} //namespace CGAL


#endif
