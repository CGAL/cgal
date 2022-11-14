// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_DEFAULT_OVERLAY_TRAITS_H
#define CGAL_ARR_DEFAULT_OVERLAY_TRAITS_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 *
 * Definition of default overlay-traits classes.
 */

#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Surface_sweep_2/Arr_default_overlay_traits_base.h>

namespace CGAL {

/*! \class
 *
 * An overlay-traits class for computing the overlay of two arrangement that
 * are templated with the default DCEL classes, namely they store no extra
 * data with their DCEL features. The resulting arrangement is also assumed
 * to be templated with the default DCEL.
 */
template <typename Arrangement_>
class Arr_default_overlay_traits :
  public _Arr_default_overlay_traits_base<Arrangement_, Arrangement_,
                                          Arrangement_>
{};

/*! \class
 *
 * An overlay-traits class for computing the overlay of two arrangement whose
 * face records are extended with auxiliary data fields, of type Data1 and
 * Data2, respectively. The resulting arrangement is also assumed to be
 * templated with the face-extended DCEL, where each face stores an auxiliary
 * Res_data field.
 * The resulting data object that corresponds to the overlay of two data
 * object of type Data1 and Data2 is computed using the functor
 * Overlay_face_data.
 */
template <typename ArrangementA, typename ArrangementB, typename ArrangementR,
          typename OverlayFaceData_>
class Arr_face_overlay_traits :
  public _Arr_default_overlay_traits_base<ArrangementA, ArrangementB,
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
  /*! Create a face f that matches the overlapping region between f1 and f2.
   */
  virtual void create_face(Face_handle_A f1, Face_handle_B f2,
                           Face_handle_R f) const
  {
    // Overlay the data objects associated with f1 and f2 and store the result
    // with f.
    f->set_data (overlay_face_data (f1->data(), f2->data()));
    return;
  }
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
