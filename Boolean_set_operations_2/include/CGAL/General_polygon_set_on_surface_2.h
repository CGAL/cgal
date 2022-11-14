// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>
//                 Ophir Setter    <ophir.setter@cs.tau.ac.il>

#ifndef CGAL_GENERAL_POLYGON_SET_ON_SURFACE_2_H
#define CGAL_GENERAL_POLYGON_SET_ON_SURFACE_2_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Boolean_set_operations_2/Gps_on_surface_base_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>

namespace CGAL {

namespace Boolean_set_operation_2_internal
{
  struct PreconditionValidationPolicy
  {
   /*! is_valid - Checks if a Traits::Polygon_2 OR
       Traits::Polygon_with_holes_2 are valid.
       This validation policy checks that polygons are valid in a
       CGAL_precondition macro. We inherit from Gps_on_surface_base_2
       and use preconditions to validate the input polygons.
   */
    template <class Polygon, class Traits>
    inline static void is_valid(const Polygon& p, const Traits& t)
    {
      CGAL_precondition(is_valid_unknown_polygon(p, t));
      CGAL_USE(p); CGAL_USE(t);
    }
  };
}

// General_polygon_set_on_surface_2
/*
  This class is derived from Gps_on_surface_base_2.
  It enforces the validation conditions for general polygons, and is therefore
  the basic implementation that should be used by the user
*/
template <class Traits_, class TopTraits_>
  class General_polygon_set_on_surface_2 :
  public Gps_on_surface_base_2<Traits_, TopTraits_,
                 Boolean_set_operation_2_internal::PreconditionValidationPolicy>
{
protected:
  typedef Traits_                                                   Traits_2;
  typedef General_polygon_set_on_surface_2<Traits_2, TopTraits_>    Self;
  typedef Gps_on_surface_base_2<Traits_2, TopTraits_,
   Boolean_set_operation_2_internal::PreconditionValidationPolicy>  Base;

public:
  typedef typename Base::Polygon_2                                  Polygon_2;
  typedef typename Base::Polygon_with_holes_2
    Polygon_with_holes_2;
  typedef typename Base::Arrangement_on_surface_2
    Arrangement_on_surface_2;

public:

  // default costructor
  General_polygon_set_on_surface_2() : Base()
  {}

  // constructor from a traits object
  General_polygon_set_on_surface_2(const Traits_2& traits) : Base(traits) {}

  // copy constructor
  General_polygon_set_on_surface_2(const Self& ps) : Base(ps) {}

  // assignment operator
  General_polygon_set_on_surface_2& operator=(const Self& ps)
  {
    Base::operator=(ps);
    return (*this);
  }

  // constructor from a polygon
  explicit General_polygon_set_on_surface_2(const Polygon_2& pgn) : Base(pgn) {}

  // constructor from a polygon with holes
  explicit
  General_polygon_set_on_surface_2(const Polygon_with_holes_2& pwh) :
    Base(pwh)
  {}

  // constructor from a polygon and a traits object
  explicit General_polygon_set_on_surface_2(const Polygon_2& pgn,
                                            const Traits_2& traits) :
    Base(pgn, traits) {}

  // constructor from a polygon with holes and a traits object
  explicit
  General_polygon_set_on_surface_2(const Polygon_with_holes_2& pwh,
                                   const Traits_2& traits) :
    Base(pwh, traits)
  {}

protected:
  General_polygon_set_on_surface_2(Arrangement_on_surface_2* arr) : Base(arr) {}

public:
  // destructor
  virtual ~General_polygon_set_on_surface_2() {}

  void intersection(const Self& gps1, const Self& gps2)
  { Base::intersection(gps1.base(), gps2.base()); }

  void join(const Self& gps1, const Self& gps2)
  { Base::join(gps1.base(), gps2.base()); }

  void symmetric_difference(const Self& gps1, const Self& gps2)
  { Base::symmetric_difference(gps1.base(), gps2.base()); }

  // For some reason the below functions (the ones that we call "using" for)
  // are hidden by the function in this class and are not found in the parent's
  // class (General_polygon_set_on_surface_2) when they are called on an
  // object of type General_polygon_set_2.
  // Check in the Vandervoorde / Stroustrup books what is the exact reason.
  // (There may be a better and more correct solution.)
  using Base::intersection;
  using Base::join;
  using Base::symmetric_difference;

private:
  const Base& base() const { return static_cast<const Base&> (*this); }

  Base& base() { return static_cast<Base&> (*this); }
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_GENERAL_POLYGON_SET_ON_SURFACE_2_H
