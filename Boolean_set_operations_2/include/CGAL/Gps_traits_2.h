// Copyright (c) 2005  Tel-Aviv University (Israel).
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

#ifndef CGAL_GPS_TRAITS_2_H
#define CGAL_GPS_TRAITS_2_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/General_polygon_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>
#include <CGAL/Boolean_set_operations_2/Gps_traits_adaptor.h>

namespace CGAL {

template <typename ArrTraits_2,
          typename GeneralPolygon_2 = General_polygon_2<ArrTraits_2> >
class Gps_traits_2 : public ArrTraits_2 {
  typedef ArrTraits_2                                   Base;
  typedef Gps_traits_2<ArrTraits_2, GeneralPolygon_2>   Self;

public:
  typedef typename Base::Point_2                        Point_2;
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;
  typedef typename Base::Multiplicity                   Multiplicity;

  // Polygon_2 type is required by GeneralPolygonSetTraits Concept
  typedef GeneralPolygon_2                              Polygon_2;
  // Backward compatibility
  typedef Polygon_2                                     General_polygon_2;

  // Polygon_with_holes_2 type required by GeneralPolygonSetTraits Concept.
  typedef CGAL::General_polygon_with_holes_2<Polygon_2> Polygon_with_holes_2;
  // Backward compatibility
  typedef Polygon_with_holes_2
    General_polygon_with_holes_2;

  typedef typename Polygon_2::Curve_const_iterator      Curve_const_iterator;

  typedef typename Polygon_with_holes_2::Hole_const_iterator
                                                        Hole_const_iterator;

  typedef typename Base::Compare_endpoints_xy_2         Compare_endpoints_xy_2;
  typedef typename Base::Construct_min_vertex_2         Construct_min_vertex_2;
  typedef typename Base::Construct_max_vertex_2         Construct_max_vertex_2;

  /*!
   * A functor for constructing a polygon from a range of x-monotone curves.
   */
  class Construct_polygon_2 {
  public:
    template <class XCurveIterator>
    void
    operator()(XCurveIterator begin, XCurveIterator end, Polygon_2& pgn) const
    { pgn.init(begin, end); }
  };

  Construct_polygon_2 construct_polygon_2_object() const
  { return Construct_polygon_2(); }

  /*!
   * A functor for scanning all x-monotone curves that form a polygon boundary.
   */
  class Construct_curves_2 {
  public:
    std::pair<Curve_const_iterator, Curve_const_iterator>
    operator()(const Polygon_2& pgn) const
    { return std::make_pair(pgn.curves_begin(), pgn.curves_end()); }
  };

  Construct_curves_2 construct_curves_2_object() const
  { return Construct_curves_2(); }

  /*!
   * An auxiliary functor used for validity checks.
   */
  typedef Gps_traits_adaptor<Base>                      Traits_adaptor;

  /*typedef CGAL::Is_valid_2<Self, Traits_adaptor>           Is_valid_2;
    Is_valid_2 is_valid_2_object() const
    {
    Traits_adaptor   tr_adp;

    return (Is_valid_2 (*this, tr_adp));
    }*/

  //Added Functionality from GeneralPolygonWithHoles Concept to the traits.

  /*A functor for constructing the outer boundary of a polygon with holes*/
  class Construct_outer_boundary {
  public:
    Polygon_2 operator()(const Polygon_with_holes_2& pol_wh) const
    { return pol_wh.outer_boundary(); }
  };

  Construct_outer_boundary construct_outer_boundary_object() const
  { return Construct_outer_boundary(); }

  /* A functor for constructing the container of holes of a polygon with holes
   */
  class Construct_holes {
  public:
    std::pair<Hole_const_iterator, Hole_const_iterator>
    operator()(const Polygon_with_holes_2& pol_wh) const
    { return std::make_pair(pol_wh.holes_begin(), pol_wh.holes_end()); }
  };

  Construct_holes construct_holes_object() const { return Construct_holes(); }

  /* A functor for constructing a Polygon_with_holes_2 from a
   * Polygon_2 (and possibly a range of holes).
   *
   * constructs a general polygon with holes using a given general polygon
   * outer as the outer boundary and a given range of holes. If outer is an
   * empty general polygon, then an unbounded polygon with holes will be
   * created. The holes must be contained inside the outer boundary, and the
   * polygons representing the holes must be strictly simple and pairwise
   * disjoint, except perhaps at the vertices.
   */
  class Construct_polygon_with_holes_2 {
  public:
    Polygon_with_holes_2 operator()(const Polygon_2& pgn_boundary) const
    { return Polygon_with_holes_2(pgn_boundary); }

    template <class HolesInputIterator>
    Polygon_with_holes_2 operator()(const Polygon_2& pgn_boundary,
                                    HolesInputIterator h_begin,
                                    HolesInputIterator h_end) const
    { return Polygon_with_holes_2(pgn_boundary, h_begin,h_end); }
  };

  Construct_polygon_with_holes_2 construct_polygon_with_holes_2_object() const
  { return Construct_polygon_with_holes_2(); }

  // Return true if the outer boundary is empty, and false otherwise.
  class Is_unbounded {
  public:
    bool operator()(const Polygon_with_holes_2& pol_wh) const
    { return pol_wh.is_unbounded(); }
  };

  Is_unbounded is_unbounded_object() const { return Is_unbounded(); }

  // Equality operator
  class Equal_2 {
  protected:

    /*! The base traits (in case it has state) */
    const Base& m_traits;

    /*! Construct from base traits.
     * \param traits the traits (in case it has state)
     * The constructor is declared protected to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Equal_2(const Base& traits) : m_traits(traits) {}

    friend class Gps_traits_2;

  public:
    /*! Determine whether the two points are the same.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const
    { return m_traits.equal_2_object()(p1, p2); }

    /*! Check whether the two x-monotone curves are the same (have the same
     * graph).
     * \param cv1 the first curve.
     * \param cv2 the second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    { return m_traits.equal_2_object()(cv1, cv2); }

    //! Compare two general polygons
    bool operator()(const Polygon_2& pgn1, const Polygon_2& pgn2) const {
      if (pgn1.size() != pgn2.size()) return false;
      if (pgn1.is_empty() && ! pgn2.is_empty()) return false;
      if (pgn2.is_empty()) return false;

      auto eql = m_traits.equal_2_object();
      auto it1 = pgn1.curves_begin();
      auto it2 = pgn2.curves_begin();
      for (; it2 != pgn2.curves_end(); ++it2) if (eql(*it1, *it2)) break;
      if (it2 == pgn2.curves_end()) return false;
      ++it1;
      ++it2;
      while (it1 != pgn1.curves_end()) {
        if (! eql(*it1++, *it2++)) return false;
        if (it2 == pgn2.curves_end()) it2 = pgn2.curves_begin();
      }
      return true;
    }

    //! Compare two general polygons
    bool operator()(const Polygon_with_holes_2& pgn1,
                    const Polygon_with_holes_2& pgn2) const {
      if (! operator()(pgn1.outer_boundary(), pgn2.outer_boundary()))
        return false;
      if (pgn1.number_of_holes(), pgn2.number_of_holes()) return false;
      auto it1 = pgn1.holes_begin();
      auto it2 = pgn2.holes_begin();
      while (it1 != pgn1.holes_end())
        if (! operator()(*it1++, *it2++)) return false;
      return true;
    }
  };

  Equal_2 equal_2_object() const { return Equal_2(*this); }

};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
