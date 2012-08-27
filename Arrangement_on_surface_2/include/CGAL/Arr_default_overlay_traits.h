// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_DEFAULT_OVERLAY_TRAITS_H
#define CGAL_ARR_DEFAULT_OVERLAY_TRAITS_H

/*! \file
 * Definition of default overlay-traits classes.
 */

#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Sweep_line_2/Arr_default_overlay_traits_base.h>

namespace CGAL {

/*!
 * \class
 * An overlay-traits class for computing the overlay of two arrangement that
 * are templated with the default DCEL classes, namely they store no extra
 * data with their DCEL features. The resulting arrangement is also assumed
 * to be templated with the default DCEL.
 */
template <class Arrangement_>
class Arr_default_overlay_traits :
  public _Arr_default_overlay_traits_base<Arrangement_,
                                          Arrangement_,
                                          Arrangement_>
{};

/*!
 * \class
 * An overlay-traits class for computing the overlay of two arrangement whose
 * face records are extended with auxiliary data fields, of type Data1 and
 * Data2, respectively. The resulting arrangement is also assumed to be
 * templated with the face-extended DCEL, where each face stores an auxiliart
 * Res_data field.
 * The resulting data object that corresponds to the overlay of two data
 * object of type Data1 and Data2 is computed using the functor
 * Overlay_face_data.
 */
template <class ArrangementA, class ArrangementB, class ArrangementR,
	  class OverlayFaceData_>
class Arr_face_overlay_traits :
  public _Arr_default_overlay_traits_base<ArrangementA,
                                          ArrangementB,
                                          ArrangementR> 
{
public:

  typedef typename ArrangementA::Face_const_handle    Face_handle_A;
  typedef typename ArrangementB::Face_const_handle    Face_handle_B;
  typedef typename ArrangementR::Face_handle          Face_handle_R;

  typedef OverlayFaceData_                            Overlay_face_data;

private:

  Overlay_face_data         overlay_face_data;

public:

  /*!
   * Create a face f that matches the overlapping region between f1 and f2.
   */
  virtual void create_face (Face_handle_A f1,
			    Face_handle_B f2,
			    Face_handle_R f) const
  {
    // Overlay the data objects associated with f1 and f2 and store the result
    // with f.
    f->set_data (overlay_face_data (f1->data(), f2->data()));
    return;
  }

};

} //namespace CGAL

#endif
