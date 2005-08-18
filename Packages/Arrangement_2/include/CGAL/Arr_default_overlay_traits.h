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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Ron Wein <baruchzu@post.tau.ac.il>

#ifndef ARR_DEFAULT_OVERLAY_TRAITS_H
#define ARR_DEFAULT_OVERLAY_TRAITS_H

#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2/Arr_overlay_traits.h>

CGAL_BEGIN_NAMESPACE

/*!
 * \class
 * An overlay-traits class for computing the overlay of two arrangement that
 * are templated with the default DCEL classes, namely they store no extra
 * data with their DCEL features. The resulting arrangement is also assumed
 * to be templated with the default DCEL.
 */
template <class Arrangement_>
class Arr_default_overlay_traits :
  public _Arr_default_overlay_traits<Arrangement_, Arrangement_, Arrangement_>
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
template <class Arrangement1_, class Arrangement2_, class ResArr_,
	  class OverlayFaceData_>
class Arr_face_overlay_traits :
  public _Arr_default_overlay_traits<Arrangement1_, Arrangement2_, ResArr_> 
{
public:

  typedef typename Arrangement1_::Face_const_handle    Face_handle1;
  typedef typename Arrangement2_::Face_const_handle    Face_handle2;
  typedef typename ResArr_::Face_handle                Res_face_handle;

  typedef OverlayFaceData_                             Overlay_face_data;

private:

  Overlay_face_data         overlay_face_data;

public:

  /*!
   * Create a face f that matches the overlapping region between f1 and f2.
   */
  virtual void create_face (Face_handle1 f1,
			    Face_handle2 f2,
			    Res_face_handle f) const
  {
    // Overlay the data objects associated with f1 and f2 and store the result
    // with f.
    f->set_data (overlay_face_data (f1->data(), f2->data()));
    return;
  }

};

CGAL_END_NAMESPACE

#endif
