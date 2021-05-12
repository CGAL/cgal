// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_INSERTION_TRAITS_2_H
#define CGAL_ARR_INSERTION_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Defintion of the Arr_insertion_traits_2<Traits,Arrangement> class.
 */

#include <CGAL/Surface_sweep_2/Arr_basic_insertion_traits_2.h>

namespace CGAL {

/*! \class Arr_insertion_traits_2
 *
 * A meta-traits class that stores a halfedge handle with every x-monotone
 * curve, and a vertex handle with each point. This information is used to
 * speed up the aggregated insertion process.
 */
template <typename GeometryTraits_2, typename Arrangement_2_>
class Arr_insertion_traits_2 :
  public Arr_basic_insertion_traits_2<GeometryTraits_2, Arrangement_2_>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Arrangement_2_                                Arrangement_2;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Arrangement_2                                 Arr2;
  typedef Arr_basic_insertion_traits_2<Gt2, Arr2>       Base;

public:
  typedef typename Gt2::Intersect_2                     Base_intersect_2;
  typedef typename Gt2::Split_2                         Base_split_2;
  typedef typename Base::Base_x_monotone_curve_2        Base_x_monotone_curve_2;
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;
  typedef typename Base::Halfedge_handle                Halfedge_handle;
  typedef typename Base::Base_point_2                   Base_point_2;
  typedef typename Base::Point_2                        Point_2;
  typedef typename Base::Multiplicity                   Multiplicity;

  typedef typename Base::Has_left_category            Has_left_category;
  typedef typename Base::Has_do_intersect_category    Has_do_intersect_category;

  // should be ok, as basic_insertion (=Base) completes incomplete tags
  typedef typename Base::Left_side_category             Left_side_category;
  typedef typename Base::Bottom_side_category           Bottom_side_category;
  typedef typename Base::Top_side_category              Top_side_category;
  typedef typename Base::Right_side_category            Right_side_category;

  /* Insertion is implemented as sweep-line visitor. The sweep-line algorithm
   * never performs merging of curves. Therefore, AreMergeable_2 and
   * Merge_2 are not needed either.
   */
  typedef Tag_false                                   Has_merge_category;

public:
  /*! Constructor with a traits class. */
  Arr_insertion_traits_2(const Gt2& tr) : Base(tr) {}

  /*! A functor that compares compares the y-coordinates of two x-monotone
   * curves immediately to the right of their intersection point.
   */
  class Intersect_2 {
  protected:
    //! The base operators.
    Base_intersect_2 m_base_intersect;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Intersect_2(const Base_intersect_2& base) :  m_base_intersect (base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_insertion_traits_2<Gt2, Arrangement_2>;

  public:
    template<typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi)
    {
      typedef std::pair<Point_2, Multiplicity>          Intersection_point;
      typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                        Intersection_result;
      typedef boost::variant<Intersection_point, Base_x_monotone_curve_2>
                                                        Intersection_base_result;

      Halfedge_handle invalid_he;

      if ((cv1.halfedge_handle() != invalid_he) &&
          (cv2.halfedge_handle() != invalid_he) &&
          (cv1.halfedge_handle() != cv2.halfedge_handle()))
      {
        // The curves are interior-disjoint as both of them are already in
        //  the arrangement.
        return oi;
      }

      std::vector<Intersection_base_result> xections;
      m_base_intersect(cv1.base(), cv2.base(), std::back_inserter(xections));
      // convert objects that are associated with Base_x_monotone_curve_2 to
      // X_monotone_curve_2
      for (const auto& xection : xections) {
        const Intersection_point*
          p_p = boost::get<Intersection_point>(&xection);
        if (p_p != nullptr) {
          *oi++ = Intersection_result(xection);
          continue;
        }
        const Base_x_monotone_curve_2* base_cv_p =
          boost::get<Base_x_monotone_curve_2>(&xection);
        CGAL_assertion(base_cv_p);

        // Add halfedge handles to the resulting curve.
        Halfedge_handle he;
        if (cv1.halfedge_handle() != invalid_he) he = cv1.halfedge_handle();
        else if (cv2.halfedge_handle() != invalid_he)
          he = cv2.halfedge_handle();
        X_monotone_curve_2 cv(*base_cv_p, he);
        cv.set_overlapping();
        *oi++ = Intersection_result(cv);
      }
      xections.clear();
      return oi;
    }
  };

  /*! Obtain a Intersect_2 function object */
  Intersect_2 intersect_2_object () const
  { return (Intersect_2(this->m_base_traits->intersect_2_object())); }

  /*! A functor that splits an arc at a point. */
  class Split_2 {
  protected:
    //! The base operator.
    Base_split_2 m_base_split;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Split_2(const Base_split_2& base) : m_base_split (base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_insertion_traits_2<Gt2, Arrangement_2>;

  public:
    void operator()(const X_monotone_curve_2& cv, const Point_2 & p,
                    X_monotone_curve_2& c1, X_monotone_curve_2& c2)
    {
      m_base_split(cv.base(), p.base(), c1.base(), c2.base());
      c1.set_halfedge_handle(cv.halfedge_handle());
      c2.set_halfedge_handle(cv.halfedge_handle());
    }
  };

  /*! Obtain a Split_2 function object */
  Split_2 split_2_object() const
  { return (Split_2(this->m_base_traits->split_2_object())); }
};

} // namespace CGAL

#endif
