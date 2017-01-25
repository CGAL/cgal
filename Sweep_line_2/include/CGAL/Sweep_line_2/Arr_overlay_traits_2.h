// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_OVERLAY_TRAITS_2_H
#define CGAL_ARR_OVERLAY_TRAITS_2_H

#include <CGAL/license/Sweep_line_2.h>


/*!
 * Defintion of the Arr_overlay_traits_2 class-template.
 */

#include <boost/variant.hpp>
#include <boost/optional.hpp>

#include <CGAL/Arr_tags.h>
#include <CGAL/Object.h>

namespace CGAL {

/*! \class
 * A meta-traits class that stores a red or a blue halfedge handle with every
 * x-monotone curve, and a red or blue vertex handle with each point. This
 * information is used to speed up the overlay of a red arrangement and a blue
 * arrangement one on top of the other.
 */
template <typename GeometryTraits_, typename ArrangementRed_,
          typename ArrangementBlue_>
class Arr_overlay_traits_2 {
public:
  typedef GeometryTraits_                           Traits_2;
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
  typedef typename Traits_2::Multiplicity           Multiplicity; 

  typedef typename Traits_2::Compare_x_2            Base_compare_x_2;
  typedef typename Traits_2::Compare_xy_2           Base_compare_xy_2;
  typedef typename Traits_2::Construct_min_vertex_2 Base_construct_min_vertex_2;
  typedef typename Traits_2::Construct_max_vertex_2 Base_construct_max_vertex_2;
  typedef typename Traits_2::Is_vertical_2          Base_is_vertical_2;
  typedef typename Traits_2::Compare_y_at_x_2       Base_compare_y_at_x_2;
  typedef typename Traits_2::Compare_y_at_x_right_2 Base_compare_y_at_x_right_2;
  typedef typename Traits_2::Intersect_2            Base_intersect_2;
  typedef typename Traits_2::Split_2                Base_split_2;
  typedef typename Traits_2::Equal_2                Base_equal_2;

  typedef typename Traits_2::Has_do_intersect_category
                                                    Has_do_intersect_category;

  typedef typename internal::Arr_complete_left_side_category<Traits_2>::Category
                                                    Left_side_category;
  typedef typename internal::Arr_complete_bottom_side_category<Traits_2>::Category
                                                    Bottom_side_category;
  typedef typename internal::Arr_complete_top_side_category<Traits_2>::Category
                                                    Top_side_category;
  typedef typename internal::Arr_complete_right_side_category<Traits_2>::Category
                                                    Right_side_category;

  /* Overlay is implemented as sweep-line visitor. The sweep-line algorithm
   * never uses Compare_y_at_x_left_2, and it never performs merging of curves.
   * Thus, AreMergeable_2 and Merge_2 are not needed either.
   */
  typedef Tag_false                                 Has_left_category;
  typedef Tag_false                                 Has_merge_category;
  
  // The color of a feature.
  enum Color {
    RED,          // From the "red" arrangement.
    BLUE,         // From the "blue" arrangement.
    RB_OVERLAP    // Red-blue overlap.
  };

  typedef boost::variant<Halfedge_handle_red, Vertex_handle_red,
                         Face_handle_red>           Cell_handle_red;
  typedef boost::optional<Cell_handle_red>          Optional_cell_red;

  typedef boost::variant<Halfedge_handle_blue, Vertex_handle_blue,
                         Face_handle_blue>          Cell_handle_blue;
  typedef boost::optional<Cell_handle_blue>         Optional_cell_blue;
      
  template <typename Handle_red>
  Optional_cell_red make_optional_cell_red(Handle_red handle_red)
  { return boost::make_optional(Cell_handle_red(handle_red)); }

  template <typename Handle_blue>
  Optional_cell_red make_optional_cell_blue(Handle_blue handle_blue)
  { return boost::make_optional(Cell_handle_blue(handle_blue)); }
  
private:
  const Traits_2* m_base_traits;        // The base traits object.

public:
  /*! Default constructor. */
  Arr_overlay_traits_2() {}

  /*! Constructor from a base traits class. */
  Arr_overlay_traits_2(const Traits_2& base_tr) : m_base_traits(&base_tr) {}

  const Traits_2* base_traits() const { return m_base_traits; }

  /*! \class
   * Nested extension of the x-monotone curve type.
   */
  class Ex_x_monotone_curve_2 {
  public:
    typedef Base_x_monotone_curve_2     Base;

  protected:
    Base                  m_base_xcv;               // The base curve.
    Halfedge_handle_red   m_red_halfedge_handle;    // The red halfedge.
    Halfedge_handle_blue  m_blue_halfedge_handle;   // The blue halfedge.

  public:
    /*! Default constructor. */
    Ex_x_monotone_curve_2() : 
      m_base_xcv(),
      m_red_halfedge_handle(),
      m_blue_halfedge_handle()
    {}

    /*! Constructor from a curve. */
    Ex_x_monotone_curve_2(const Base& xcv) :
      m_base_xcv(xcv),
      m_red_halfedge_handle(),
      m_blue_halfedge_handle()
    {}

    /*! Constructor from a curve and halfedge handles. */
    Ex_x_monotone_curve_2(const Base& xcv,
                          Halfedge_handle_red  he_r,
                          Halfedge_handle_blue he_b) :
      m_base_xcv(xcv),
      m_red_halfedge_handle(he_r),
      m_blue_halfedge_handle(he_b)
    {
      CGAL_precondition((he_r == Halfedge_handle_red()) ||
                        (he_r->direction() == ARR_RIGHT_TO_LEFT));
      CGAL_precondition((he_b == Halfedge_handle_blue()) ||
                        (he_b->direction() == ARR_RIGHT_TO_LEFT));
    }

    /*! Get the base curve (const version). */
    const Base& base() const { return m_base_xcv; }

    /*! Get the base curve (non-const version). */
    Base& base() { return m_base_xcv; }

    /*! Casting to a base curve (const version). */
    operator const Base&() const { return m_base_xcv; }

    /*! Casting to a base curve (const version). */
    operator Base&() { return m_base_xcv; }

    /*! Get the red halfedge handle. */
    Halfedge_handle_red red_halfedge_handle() const
    { return m_red_halfedge_handle; }
    
    /*! Get the blue halfedge handle. */
    Halfedge_handle_blue blue_halfedge_handle() const
    { return m_blue_halfedge_handle; }

    /*! Set the red halfedge handle. */
    void set_red_halfedge_handle(Halfedge_handle_red he_r)
    {
      CGAL_precondition((he_r == Halfedge_handle_red()) ||
                        (he_r->direction() == ARR_RIGHT_TO_LEFT));

      m_red_halfedge_handle = he_r;
    }

    /*! Set the blue halfedge handle. */
    void set_blue_halfedge_handle(Halfedge_handle_blue he_b)
    {
      CGAL_precondition((he_b == Halfedge_handle_blue()) ||
                        (he_b->direction() == ARR_RIGHT_TO_LEFT));

      m_blue_halfedge_handle = he_b;
    }

    /*! Get the color of the subcurve. */
    Color color() const
    {
      Halfedge_handle_red     null_red_he;
      Halfedge_handle_blue    null_blue_he;

      if ((m_red_halfedge_handle != null_red_he) &&
          (m_blue_halfedge_handle == null_blue_he))
        return RED;

      if ((m_blue_halfedge_handle != null_blue_he) &&
          (m_red_halfedge_handle == null_red_he))
        return BLUE;
      
      // Overlap, return the RB_OVERLAP color:
      CGAL_assertion((m_red_halfedge_handle != null_red_he) && 
                     (m_blue_halfedge_handle != null_blue_he));
      return RB_OVERLAP;
    }
  }; // nested class Ex_x_monotone_curve_2 - END
  
  typedef Ex_x_monotone_curve_2                     X_monotone_curve_2;

  // For debugging purposes:
  friend std::ostream& operator<<(std::ostream& os,
                                  const X_monotone_curve_2& xcv)
  {
    os << xcv.base();
    return os;
  }
 
  /*! \class
   * Nested extension of the point type.
   */
  class Ex_point_2 {
  public:
    typedef Base_point_2    Base;

  protected:
    Base m_base_pt;                     // The base point.
    Optional_cell_red m_red_cell;       // The "red" object.
    Optional_cell_blue m_blue_cell;     // The "blue" object.

  public:
    /*! Default constructor. */
    Ex_point_2() :
      m_base_pt(),
      m_red_cell(),
      m_blue_cell()
    {}

    /*! Constructor from a base point. */
    Ex_point_2(const Base& pt) :
      m_base_pt(pt),
      m_red_cell(),
      m_blue_cell()
    {}

    /*! Constructor from a base point with red and blue objects. */
    Ex_point_2(const Base& pt, const Optional_cell_red& cell_red,
               const Optional_cell_blue& cell_blue) :
      m_base_pt(pt),
      m_red_cell(cell_red),
      m_blue_cell(cell_blue)
    {}

    /*! Get the base point (const version). */
    const Base& base() const { return m_base_pt; }

    /*! Get the base point (non-const version). */
    Base& base() { return m_base_pt; }

    /*! Casting to a base point (const version). */
    operator const Base&() const { return m_base_pt; }

    /*! Casting to a base point (non-const version). */
    operator Base&() { return m_base_pt; }

    /*! Get the red object. */
    const Optional_cell_red& red_cell() const { return m_red_cell; }

    /*! Get the blue object. */
    const Optional_cell_blue& blue_cell() const { return m_blue_cell; }

    /*! Check if the red object is empty. */
    bool is_red_cell_empty() const { return !m_red_cell; }

    /*! Check if the blue object is empty. */
    bool is_blue_cell_empty() const { return !m_blue_cell; }

    /*! Set the red object. */
    void set_red_cell(const Optional_cell_red& cell_red)
    { m_red_cell = cell_red; }

    /*! Set the blue object. */
    void set_blue_cell(const Optional_cell_blue& cell_blue)
    { m_blue_cell = cell_blue; }

    /*! Obtain the red cell handle or NULL if it doesn't exist. */
    const Cell_handle_red* red_cell_handle() const
    { return m_red_cell ? &(*m_red_cell) : NULL; }

    /*! Obtain the blue cell handle or NULL if it doesn't exist. */
    const Cell_handle_blue* blue_cell_handle() const
    { return m_blue_cell ? &(*m_blue_cell) : NULL; }
    
    /*! Obtain the red vertex handle or NULL if it doesn't exist. */
    const Vertex_handle_red* red_vertex_handle() const
    {
      return m_red_cell ? boost::get<Vertex_handle_red>(&(*m_red_cell)) : NULL;
    }

    /*! Obtain the blue vertex handle or NULL if it doesn't exist. */
    const Vertex_handle_blue* blue_vertex_handle() const
    {
      return
        m_blue_cell ? boost::get<Vertex_handle_blue>(&(*m_blue_cell)) : NULL;
    }
  };

  typedef Ex_point_2                                Point_2;

  // For debugging purposes:
  friend std::ostream& operator<<(std::ostream& os, const Point_2& pt)
  {
    os << pt.base();
    // os << ", red? " << pt.is_red_cell_empty() << ", blue? "
    //    << pt.is_blue_cell_empty();
    return os;
  }

  /*! A functor that computes intersections between x-monotone curves. */
  class Intersect_2 {
  protected:
    //! The base traits.
    const Arr_overlay_traits_2* m_traits;

    /*! Constructor.
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Intersect_2(const Arr_overlay_traits_2* traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    template<class OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xcv1,
                              const X_monotone_curve_2& xcv2,
                              OutputIterator oi)
    {
      // In case the curves originate from the same arrangement, they are
      // obviously interior-disjoint.
      if (xcv1.color() == xcv2.color()) return oi;
      
      if ((xcv1.color() == RB_OVERLAP) || (xcv2.color() == RB_OVERLAP))
        return oi;

      // Compute the intersection points between the curves. Note that if
      // xcv1 and xcv2 are subcruves of x-monotone curves that had intersected
      // before the current point on the status line, we may get a filter
      // failure if we send the subcurve whose left endpoint is to the left
      // of the other curve - this is because their previously computed
      // intersection point p may be equal to the this left endpoint. As many
      // traits classes start by computing the intersection between the
      // supporting curves and then check whether the result is in the x-range
      // of both subcurves, this will result in a filter failure. However, if
      // we send xcv1 first, then p is obviusly not in its x-range and there is
      // no filter failure.
      //
      //              / xcv1
      //             /
      //            /
      //       ----+--
      //          /
      //         /
      //      p +------------- xcv2
      //              ^
      //              |
      //              status line
      //
      // Note that we do not bother with curves whose left ends are open,
      // since such curved did not intersect before.

      const std::pair<Base_point_2, unsigned int>* base_ipt;
      const Base_x_monotone_curve_2* overlap_xcv;
      bool send_xcv1_first = true;
      OutputIterator oi_end;

      Parameter_space_in_x_2 ps_x_op = m_traits->parameter_space_in_x_2_object();
      Parameter_space_in_y_2 ps_y_op = m_traits->parameter_space_in_y_2_object();
      const Arr_parameter_space bx1 = ps_x_op(xcv1, ARR_MIN_END);
      const Arr_parameter_space by1 = ps_y_op(xcv1, ARR_MIN_END);
      const Arr_parameter_space bx2 = ps_x_op(xcv2, ARR_MIN_END);
      const Arr_parameter_space by2 = ps_y_op(xcv2, ARR_MIN_END);

      const Traits_2* m_base_tr = m_traits->base_traits();

      if ((bx1 == ARR_INTERIOR) && (by1 == ARR_INTERIOR) &&
          (bx2 == ARR_INTERIOR) && (by2 == ARR_INTERIOR))
      {
        send_xcv1_first =
          (m_base_tr->compare_xy_2_object()
           (m_base_tr->construct_min_vertex_2_object()(xcv1.base()),
            m_base_tr->construct_min_vertex_2_object()(xcv2.base())) == LARGER);
      }

      oi_end = (send_xcv1_first) ?
        m_base_tr->intersect_2_object()(xcv1.base(), xcv2.base(), oi) :
        m_base_tr->intersect_2_object()(xcv2.base(), xcv1.base(), oi);

      // Convert objects that are associated with Base_x_monotone_curve_2 to
      // the exteneded X_monotone_curve_2. 
      while (oi != oi_end) {
        base_ipt = object_cast<std::pair<Base_point_2, unsigned int> >(&(*oi));

        if (base_ipt != NULL) {
          // We have a red-blue intersection point, so we attach the
          // intersecting red and blue halfedges to it.
          Optional_cell_red red_cell;
          Optional_cell_blue blue_cell;

          if (xcv1.color() == RED) {
            CGAL_assertion(xcv2.color() == BLUE);
            red_cell =
              boost::make_optional(Cell_handle_red(xcv1.red_halfedge_handle()));
            blue_cell =
              boost::make_optional(Cell_handle_blue(xcv2.blue_halfedge_handle()));
          }
          else {
            CGAL_assertion((xcv2.color() == RED) && (xcv1.color() == BLUE));
            red_cell =
              boost::make_optional(Cell_handle_red(xcv2.red_halfedge_handle()));
            blue_cell =
              boost::make_optional(Cell_handle_blue(xcv1.blue_halfedge_handle()));
          }

          // Create the extended point and add the multiplicity.
          Point_2 ex_point(base_ipt->first, red_cell, blue_cell);
          *oi++ = CGAL::make_object(std::make_pair(ex_point, base_ipt->second));
        }
        else {
          overlap_xcv = object_cast<Base_x_monotone_curve_2>(&(*oi));
          CGAL_assertion(overlap_xcv != NULL);

          // We have a red-blue overlap, so we mark the curve accordingly.
          Halfedge_handle_red        red_he;
          Halfedge_handle_blue       blue_he;
          
          if (xcv1.color() == RED) {
            red_he = xcv1.red_halfedge_handle();
            
            // Overlap can occur only between curves from a different color.
            CGAL_assertion(xcv2.color() == BLUE);
            blue_he = xcv2.blue_halfedge_handle();
          }
          else {
            CGAL_assertion((xcv1.color() == BLUE) && (xcv2.color() == RED));

            red_he = xcv2.red_halfedge_handle();
            blue_he = xcv1.blue_halfedge_handle();
          }
          
          *oi++ = CGAL::make_object(X_monotone_curve_2(*overlap_xcv,
                                                       red_he, blue_he));
        }
      }

      // Return the past-the-end iterator.
      return oi_end;
    }
  };

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object() const { return Intersect_2(this); }

  /*! A functor that splits an arc at a point. */
  class Split_2 {
  protected:
    //! The base operator.
    Base_split_2    m_base_split;

    /*! Constructor.
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Split_2(const Base_split_2& base) : m_base_split(base) {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    void operator()(const X_monotone_curve_2& xcv, const Point_2& p,
                    X_monotone_curve_2& c1, X_monotone_curve_2& c2)
    {
      // Split the base curve.
      m_base_split(xcv.base(), p.base(), c1.base(), c2.base());

      // Duplicate the halfedge data to the resulting curves.
      c1.set_red_halfedge_handle(xcv.red_halfedge_handle());
      c1.set_blue_halfedge_handle(xcv.blue_halfedge_handle());

      c2.set_red_halfedge_handle(xcv.red_halfedge_handle());
      c2.set_blue_halfedge_handle(xcv.blue_halfedge_handle());
    }
  };

  /*! Obtain a Split_2 functor object. */
  Split_2 split_2_object() const
  { return Split_2(m_base_traits->split_2_object()); }

  /*! A functor that obtains the left endpoint of an x-monotone curve. */
  class Construct_min_vertex_2 {
  protected:
    //! The base operators.
    Base_construct_min_vertex_2  m_base_min_v;
    Base_equal_2                 m_base_equal;

    /*! Constructor.
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Construct_min_vertex_2(const Base_construct_min_vertex_2& base_min_v,
                           const Base_equal_2& base_equal) :
      m_base_min_v(base_min_v),
      m_base_equal(base_equal)
    {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Point_2 operator()(const X_monotone_curve_2& xcv) 
    {
      // Create the objects that wrap the arrangement vertex.
      // Note that the halfedges associated with the curves are always
      // directed from right to left, so their target is the smaller end.
      const Base_point_2& base_p = m_base_min_v(xcv.base());
      Optional_cell_red red_cell;
      Optional_cell_blue blue_cell;

      if ((xcv.color() == RED) || (xcv.color() == RB_OVERLAP))
        red_cell =
          (! xcv.red_halfedge_handle()->target()->is_at_open_boundary() &&
           m_base_equal(base_p, xcv.red_halfedge_handle()->target()->point())) ?
          boost::make_optional(Cell_handle_red(xcv.red_halfedge_handle()->target())) :
          boost::make_optional(Cell_handle_red(xcv.red_halfedge_handle()));

      if ((xcv.color() == BLUE) || (xcv.color() == RB_OVERLAP))
        blue_cell =
          (! xcv.blue_halfedge_handle()->target()->is_at_open_boundary() &&
           m_base_equal(base_p, xcv.blue_halfedge_handle()->target()->point())) ?
          boost::make_optional(Cell_handle_blue(xcv.blue_halfedge_handle()->target())) :
          boost::make_optional(Cell_handle_blue(xcv.blue_halfedge_handle()));

      return Point_2(base_p, red_cell, blue_cell);
    }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  {
    return 
      Construct_min_vertex_2(m_base_traits->construct_min_vertex_2_object(),
                             m_base_traits->equal_2_object());
  }
 
  /*! A functor that obtains the right endpoint of an x-monotone curve. */
  class Construct_max_vertex_2 {
  protected:
    //! The base operators.
    Base_construct_max_vertex_2  m_base_max_v;
    Base_equal_2                 m_base_equal;

    /*! Constructor.
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Construct_max_vertex_2(const Base_construct_max_vertex_2& base_max_v,
                           const Base_equal_2& base_equal) :
      m_base_max_v(base_max_v),
      m_base_equal(base_equal)
    {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Point_2 operator()(const X_monotone_curve_2& xcv) const
    {
      // Create the objects that wrap the arrangement vertex.
      // Note that the halfedges associated with the curves are always
      // directed from right to left, so their target is the smaller end.
      const Base_point_2& base_p = m_base_max_v(xcv.base());
      Optional_cell_red red_cell;
      Optional_cell_blue blue_cell;

      if ((xcv.color() == RED) || (xcv.color() == RB_OVERLAP))
        red_cell = 
          (! xcv.red_halfedge_handle()->source()->is_at_open_boundary() &&
           m_base_equal(base_p, xcv.red_halfedge_handle()->source()->point())) ?
          boost::make_optional(Cell_handle_red(xcv.red_halfedge_handle()->source())) :
          boost::make_optional(Cell_handle_red(xcv.red_halfedge_handle()));
      
      if ((xcv.color() == BLUE) || (xcv.color() == RB_OVERLAP))
        blue_cell = 
          (! xcv.blue_halfedge_handle()->source()->is_at_open_boundary() &&
           m_base_equal(base_p, xcv.blue_halfedge_handle()->source()->point())) ?
          boost::make_optional(Cell_handle_blue(xcv.blue_halfedge_handle()->source())) :
          boost::make_optional(Cell_handle_blue(xcv.blue_halfedge_handle()));

      return (Point_2(base_p, red_cell, blue_cell));
    }
  };

  /*! Obtain a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  {
    return
      Construct_max_vertex_2(m_base_traits->construct_max_vertex_2_object(),
                             m_base_traits->equal_2_object());
  }

  /*! A functor that checks whether a given x-monotone curve is vertical. */
  class Is_vertical_2 {
  protected:
    //! The base operator.
    Base_is_vertical_2 m_base_is_vert;

    /*! Constructor.
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Is_vertical_2(const Base_is_vertical_2& base) : m_base_is_vert(base) {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    bool operator()(const X_monotone_curve_2& xcv) const
    { return m_base_is_vert(xcv.base()); }
  };

  /*! Obtain a Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const
  { return Is_vertical_2(m_base_traits->is_vertical_2_object()); }

  /*! A functor that checks whether two points and two x-monotone curves are
   * identical.
   */
  class Equal_2 {
  protected:
    //! The base operator.
    Base_equal_2 m_base_equal;

    /*! Constructor.
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Equal_2(const Base_equal_2& base) : m_base_equal(base) {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    bool operator()(const Point_2& p1, const Point_2& p2) const
    { return m_base_equal(p1.base(), p2.base()); }

    bool operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2) const
    { return m_base_equal(xcv1.base(), xcv2.base()); }
  };

  /*! Obtain a Equal_2 functor object. */
  Equal_2 equal_2_object() const
  { return Equal_2(m_base_traits->equal_2_object()); }
  
  /*! A functor that compares the x-coordinates of two points */
  class Compare_x_2 {
  protected:
    //! The base operator.
    Base_compare_x_2 m_base_cmp_x;

    /*! Constructor.
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_2(const Base_compare_x_2& base) : m_base_cmp_x(base) {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { return m_base_cmp_x(p1.base(), p2.base()); }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Compare_x_2 compare_x_2_object() const
  { return Compare_x_2(m_base_traits->compare_x_2_object()); }

  /*! A functor that compares two points lexigoraphically: by x, then by y. */
  class Compare_xy_2 {
  protected:
    //! The base operator.
    Base_compare_xy_2 m_base_cmp_xy;

    /*! Constructor.
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_xy_2(const Base_compare_xy_2& base) : m_base_cmp_xy(base) {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      // Check if there wither points represent red or blue vertices.
      const Vertex_handle_red* vr1 = p1.red_vertex_handle();
      const Vertex_handle_red* vr2 = p2.red_vertex_handle();
      const Vertex_handle_blue* vb1 = p1.blue_vertex_handle();
      const Vertex_handle_blue* vb2 = p2.blue_vertex_handle();

      // EFEF: The code below looks strange. I think it should be:
      // if (assign_v1_red && assign_v2_red) return (vr1 == vr2);
      // if (assign_v1_blue && assign_v2_blue) return (vb1 == vb2);
      // return (m_base_cmp_xy (p1.base(), p2.base()));
      // If both points, p1 and p2, belong to the same arrangement, then:
      // - the points are not equal if the corresponding handles are not equal.
      // - if the corresponding handles are equal, the points must be equal
      //   but this holds for any two points.
      
      if ((vr1 && vb1) || (vr2 && vb2))
        // In case of an overlapping vertex, just perform the comparison.
        return (m_base_cmp_xy(p1.base(), p2.base()));

      // If both points are vertices of the same color, avoid the geometric
      // comparison if they refer to the same vertex.
      if ((vr1 && vr2 && (*vr1 == *vr2)) || (vb1 && vb2 && (*vb1 == *vb2)))
        return (EQUAL);

      return (m_base_cmp_xy(p1.base(), p2.base()));
    }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const
  { return (Compare_xy_2(m_base_traits->compare_xy_2_object())); }

  /*! A functor that compares the y-coordinates of a point and an
   * x-monotone curve at the point x-coordinate.
   */
  class Compare_y_at_x_2 {
  protected:
    //! The base operator.
    Base_compare_y_at_x_2 m_base_cmp_y_at_x;

    /*! Constructor.
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_at_x_2(const Base_compare_y_at_x_2& base) :
        m_base_cmp_y_at_x(base)
    {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv) const
    { return m_base_cmp_y_at_x(p.base(), xcv.base()); }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return (Compare_y_at_x_2(m_base_traits->compare_y_at_x_2_object())); }

  /*! A functor that compares compares the y-coordinates of two x-monotone
   * curves immediately to the right of their intersection point.
   */
  class Compare_y_at_x_right_2 {
  protected:
    //! The base operator.
    Base_compare_y_at_x_right_2    m_base_cmp_y_at_x_right;

    /*! Constructor.
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_at_x_right_2(const Base_compare_y_at_x_right_2& base):
      m_base_cmp_y_at_x_right(base)
    {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 const Point_2& p) const
    { return m_base_cmp_y_at_x_right(xcv1.base(), xcv2.base(), p.base()); }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  {
    return 
      Compare_y_at_x_right_2(m_base_traits->compare_y_at_x_right_2_object());
  }

  // left-right

  /*! A functor that determines whether an endpoint of an x-monotone curve lies
   * on a boundary of the parameter space along the x axis.
   */
  class Parameter_space_in_x_2 {
  protected:
    //! The base traits.
    const Traits_2* m_base;

    /*! Constructor.
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Parameter_space_in_x_2(const Traits_2* tr) : m_base(tr) {}
 
    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                   Arr_curve_end ce) const
    { return m_base->parameter_space_in_x_2_object()(xcv.base(), ce); }

    Arr_parameter_space operator()(const Point_2& p) const
    { return m_base->parameter_space_in_x_2_object()(p.base()); }

    Arr_parameter_space operator()(const X_monotone_curve_2& xcv) const
    { return m_base->parameter_space_in_x_2_object()(xcv.base()); }
  };

  /*! Obtain an Parameter_space_in_x_2 functor object. */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { return Parameter_space_in_x_2(m_base_traits); } 
 
  /*! A function object that determines whether an x-monotone curve or a
   * point coincide with the vertical identification curve.
   */
  class Is_on_x_identification_2 {
  protected:
    //! The base traits.
    const Traits_2* m_base;

    /*! Constructor.
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Is_on_x_identification_2(const Traits_2* tr) : m_base(tr) {}
 
    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Arr_parameter_space operator()(const Point_2& p) const
    { return m_base->is_on_x_identification_2_object()(p.base()); }

    Arr_parameter_space operator()(const X_monotone_curve_2& xcv) const
    { return m_base->is_on_x_identification_2_object()(xcv.base()); }
  };

  /*! Obtain an Is_on_x_identification_2 functor object. */
  Is_on_x_identification_2 is_on_x_identification_2_object() const
  { return Is_on_x_identification_2(m_base_traits); } 

  /*! A functor that compares the y-values of pointss on the
   * boundary of the parameter space.
   */
  class Compare_y_on_boundary_2 {
  protected:
    //! The base traits.
    const Traits_2* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle the
     *             case it is not stateless (e.g., it stores data).
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_on_boundary_2(const Traits_2* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Comparison_result operator()(const Point_2& pt1, const Point_2& pt2) const
    { return m_base->compare_y_on_boundary_2_object()(pt1.base(), pt2.base()); }
  };
  
  /*! Obtain a Compare_y_on_boundary_2 functor. */
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const
  { return Compare_y_on_boundary_2(m_base_traits); }

  /*! A functor that compares the y-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  protected:
    //! The base traits.
    const Traits_2* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle the
     *             case it is not stateless (e.g., it stores data).
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_near_boundary_2(const Traits_2* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2, 
                                 Arr_curve_end ce) const
    {
      // If the traits class does not support open curves, we just
      // return EQUAL, as this comparison will not be invoked anyway.
      return m_base->compare_y_near_boundary_2_object()(xcv1.base(),
                                                        xcv2.base(), ce);
    }
  };
  
  /*! Obtain a Compare_y_near_boundary_2 functor. */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(m_base_traits); }
  

  // bottom-top
  
  /*! A functor that determines whether an endpoint of an x-monotone arc lies
   * on a boundary of the parameter space along the y axis.
   */
  class Parameter_space_in_y_2 {
  protected:
    //! The base traits.
    const Traits_2* m_base;

    /*! Constructor.
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Parameter_space_in_y_2(const Traits_2* tr) : m_base(tr) {}
   
    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                   Arr_curve_end ce) const
    { return m_base->parameter_space_in_y_2_object()(xcv.base(), ce); }

    Arr_parameter_space operator()(const Point_2& p) const
    { return m_base->parameter_space_in_y_2_object()(p.base()); }

    Arr_parameter_space operator()(const X_monotone_curve_2& xcv) const
    { return m_base->parameter_space_in_y_2_object()(xcv.base()); }
  };

  /*! Obtain an Parameter_space_in_y_2 functor object. */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(m_base_traits); } 

  /*! A function object that determines whether an x-monotone curve or a
   * point coincide with the vertical identification curve.
   */
  class Is_on_y_identification_2 {
  protected:
    //! The base traits.
    const Traits_2* m_base;

    /*! Constructor.
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Is_on_y_identification_2(const Traits_2* tr) : m_base(tr) {}
 
    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Arr_parameter_space operator()(const Point_2& p) const
    { return m_base->is_on_y_identification_2_object()(p.base()); }

    Arr_parameter_space operator()(const X_monotone_curve_2& xcv) const
    { return m_base->is_on_y_identification_2_object()(xcv.base()); }
  };

  /*! Obtain an Is_on_y_identification_2 functor object. */
  Is_on_y_identification_2 is_on_y_identification_2_object() const
  { return Is_on_y_identification_2(m_base_traits); } 

  /*! A functor that compares the x-limits of curve ends on the
   * boundary of the parameter space.
   */
  class Compare_x_at_limit_2 {
  protected:
    //! The base traits.
    const Traits_2* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle the
     *             case it is not stateless (e.g., it stores data).
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_at_limit_2(const Traits_2* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce) const
    { return m_base->compare_x_at_limit_2_object()(p.base(), xcv.base(), ce); }

    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce2) const
    {
      return m_base->compare_x_at_limit_2_object()(xcv1.base(), ce1,
                                                   xcv2.base(), ce2);
    }
  };

  /*! Obtain a Compare_x_at_limit_2 functor. */
  Compare_x_at_limit_2 compare_x_at_limit_2_object() const
  { return Compare_x_at_limit_2(m_base_traits); }

  /*! A functor that compares the x-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_x_near_limit_2 {
  protected:
    //! The base traits.
    const Traits_2* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle the
     *             case it is not stateless (e.g., it stores data).
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_near_limit_2(const Traits_2* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const
    {
      return m_base->compare_x_near_limit_2_object()(xcv1.base(), xcv2.base(),
                                                     ce);
    }
  };

  /*! Obtain a Compare_x_near_limit_2 functor. */
  Compare_x_near_limit_2 compare_x_near_limit_2_object() const
  { return Compare_x_near_limit_2(m_base_traits); }

  /*! A functor that compares the y-values of pointss on the
   * boundary of the parameter space.
   */
  class Compare_x_on_boundary_2 {
  protected:
    //! The base traits.
    const Traits_2* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle the
     *             case it is not stateless (e.g., it stores data).
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_on_boundary_2(const Traits_2* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Comparison_result operator()(const Point_2& pt1, const Point_2& pt2) const
    { return m_base->compare_x_on_boundary_2_object()(pt1.base(), pt2.base()); }

    Comparison_result operator()(const Point_2& pt,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce) const
    {
      return m_base->compare_x_on_boundary_2_object()(pt.base(), xcv.base(),
                                                      ce);
    }

    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce2) const
    {
      return m_base->compare_x_on_boundary_2_object()(xcv1.base(), ce1,
                                                      xcv2.base(), ce2);
    }
  };
  
  /*! Obtain a Compare_x_on_boundary_2 functor. */
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const
  { return Compare_x_on_boundary_2(m_base_traits); }

  /*! A functor that compares the x-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_x_near_boundary_2 {
  protected:
    //! The base traits.
    const Traits_2* m_base;

    /*! Constructor.
     * \param base The base traits class. It must be passed, to handle the
     *             case it is not stateless (e.g., it stores data).
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_near_boundary_2(const Traits_2* base) : m_base(base) {}

    //! Allow its functor obtaining function calling the protected constructor.
    friend class Arr_overlay_traits_2<Traits_2,
                                      Arrangement_red_2, Arrangement_blue_2>;
    
  public:
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const
    {
      return m_base->compare_x_near_boundary_2_object()(xcv1.base(), xcv2.base(),
							ce);
    }
  };

  /*! Obtain a Compare_x_near_boundary_2 functor. */
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const
  { return Compare_x_near_boundary_2(m_base_traits); }
};

} //namespace CGAL

#endif
