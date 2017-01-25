// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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

#ifndef CGAL_OVERLAY_SUBCURVE_H
#define CGAL_OVERLAY_SUBCURVE_H

#include <CGAL/license/Sweep_line_2.h>


/*! \file
 * Definition of the Arr_overlay_subcurve class-template.
 */

#include <CGAL/Sweep_line_2/Arr_construction_subcurve.h>

namespace CGAL {

/*! \class
 * Representation of a subcurve in the overlay process. A subcurve can
 * originate from a "red" arrangement or from a "blue" arrangement that are
 * overlaid one on top of the other. It stores information of the halfedge
 * from the opposite color it "sees" from below, so it can be easily located
 * in the proper place in the resulting arrangement.
 */
template<class Traits_>
class Arr_overlay_subcurve : public Arr_construction_subcurve<Traits_>
{
public:

  typedef Traits_                                        Traits_2;
  typedef typename Traits_2::Point_2                     Point_2;
  typedef typename Traits_2::X_monotone_curve_2          X_monotone_curve_2;

  typedef Arr_construction_subcurve<Traits_2>            Base;
  typedef Arr_overlay_subcurve<Traits_2>                 Self;

  typedef typename Base::Status_line_iterator            Status_line_iterator;

  typedef typename Traits_2::Color                       Color;

  typedef typename Traits_2::Halfedge_handle_red         Halfedge_handle_red;
  typedef typename Traits_2::Face_handle_red             Face_handle_red;
  typedef typename Face_handle_red::value_type           Face_red;

  typedef typename Traits_2::Halfedge_handle_blue        Halfedge_handle_blue;
  typedef typename Traits_2::Face_handle_blue            Face_handle_blue;
  typedef typename Face_handle_blue::value_type          Face_blue;

  typedef Sweep_line_event<Traits_2, Self>               Event;

protected:

  // Data members:
  Self* m_above;        // A subcurve of an opposite color that lies above.

  union {
    const Face_red* red;
    const Face_blue* blue;
  } m_top_face;         // If m_above is NULL, points the top face in
                        // the arrangement of the opposite color that
                        // contains the subcurve.
  
public:
  /*! Constructor. */
  Arr_overlay_subcurve() :
    Base(),
    m_above(NULL)
  {
    m_top_face.red = NULL;
  }

  /*! constructor given a curve. */
  Arr_overlay_subcurve(const X_monotone_curve_2& curve) :
    Base(curve),
    m_above(NULL)
  {
    m_top_face.red = NULL;
  }

  /*! Get the subcurve lying above above this subcurve in the status line. */
  Self* subcurve_above() const { return (m_above); }

  /*! Set the subcurve above. */
  void set_subcurve_above(Self* sc) { m_above = sc; }

  /*! Get the color of the associated curve. */
  Color color() const { return (this->m_lastCurve.color()); }

  /*! Check if two subcurves have the same color. */
  bool has_same_color(const Self* sc) const
  { return (this->m_lastCurve.color() == sc->color()); }

  /*! Get the red halfedge that represents the subcurve. */
  Halfedge_handle_red red_halfedge_handle() const
  { return (this->m_lastCurve.red_halfedge_handle()); }

  /*! Get the blue halfedge that represents the subcurve. */
  Halfedge_handle_blue blue_halfedge_handle() const
  { return (this->m_lastCurve.blue_halfedge_handle()); }

  /*! Get the red top face that contains the subcurve. */
  const Face_handle_red red_top_face() const
  { return (Face_handle_red(m_top_face.red)); }

  /*! Get the blue top face that contains the subcurve. */
  const Face_handle_blue blue_top_face() const
  { return (Face_handle_blue(m_top_face.blue)); }

  /*! Set the red top face. */
  void set_red_top_face(Face_handle_red fh) { m_top_face.red = &(*fh); }

  /*! Set the blue top face. */
  void set_blue_top_face(Face_handle_blue fh) { m_top_face.blue = &(*fh); }

  /*! Copy the top face from the given subcurve. */
  void set_top_face(const Self* sc)
  {
    CGAL_precondition(sc->m_above == NULL);
   
    // Mark there is no curve above and copy the face pointer.
    m_above = NULL;
    m_top_face.red = sc->m_top_face.red;
  }
};
  
} //namespace CGAL

#endif
