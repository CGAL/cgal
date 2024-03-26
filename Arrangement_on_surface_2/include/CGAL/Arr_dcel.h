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
// Author(s): Ron Wein        <wein@post.tau.ac.il>
//            Baruch Zukerman <baruchzu@post.tau.ac.il>
//            Efi Fogel       <efifogel@gmail.com>

#ifndef CGAL_ARR_DCEL_H
#define CGAL_ARR_DCEL_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The definition of the Arr_dcel<Traits> class.
 */

#include <CGAL/Arr_dcel_base.h>

namespace CGAL {

/*! \class
 * The arrangement DCEL class.
 * The Traits parameters corresponds to a geometric traits class, which
 * defines the Point_2 and X_monotone_curve_2 types.
 */
template <typename Traits,
          typename V = Arr_vertex_base<typename Traits::Point_2>,
          typename H = Arr_halfedge_base<typename Traits::X_monotone_curve_2>,
          typename F = Arr_face_base>
class Arr_dcel : public Arr_dcel_base<V, H, F> {
public:
  /*! \struct
   * An auxiliary structure for rebinding the DCEL with a new traits class.
   */
  template <typename T>
  struct rebind {
  private:
    using Pnt = typename T::Point_2;
    using Xcv = typename T::X_monotone_curve_2;
    using Rebind_v = typename V::template rebind<Pnt>;
    using V_other = typename Rebind_v::other;
    using Rebind_h = typename H::template rebind<Xcv>;
    using H_other = typename Rebind_h::other;

  public:
    using other = Arr_dcel<T, V_other, H_other, F>;
  };

  /*! Default constructor. */
  Arr_dcel() {}

  /*! Destructor. */
  virtual ~Arr_dcel() {}
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
