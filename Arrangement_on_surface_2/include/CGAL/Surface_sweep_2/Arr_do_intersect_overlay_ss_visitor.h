// Copyright (c) 2006,2007,2009,2010,2011,2025 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s):      Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_DO_INTERSECT_ARR_OVERLAY_SS_VISITOR_H
#define CGAL_DO_INTERSECT_ARR_OVERLAY_SS_VISITOR_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_do_intersect_overlay_ss_visitor class-template.
 */

#include <CGAL/Default.h>
#include <CGAL/Surface_sweep_2/Arr_overlay_ss_visitor.h>

namespace CGAL {

/*! \class Arr_do_intersect_overlay_ss_visitor
 *
 * A sweep-line visitor for overlaying a "red" arrangement and a "blue"
 * arrangement as long as the edges do not intersect in their interiors. If
 * there are no intersections, the overlay arrangement is constructed. All three
 * arrangements are embedded on the same type of surface and use the same
 * geometry traits. Otherwise, the process is terminated without any delay (that
 * is, once an intersection is detected).
 */
template <typename OverlayHelper, typename OverlayTraits, typename Visitor_ = Default>
class Arr_do_intersect_overlay_ss_visitor :
    public Arr_overlay_ss_visitor<
      OverlayHelper, OverlayTraits,
      typename Default::Get<Visitor_,
                            Arr_do_intersect_overlay_ss_visitor<OverlayHelper, OverlayTraits, Visitor_> >::type> {

private:
  using Overlay_helper = OverlayHelper;
  using Overlay_traits = OverlayTraits;

  using Self = Arr_do_intersect_overlay_ss_visitor<Overlay_helper, Overlay_traits, Visitor_>;
  using Visitor = typename Default::Get<Visitor_, Self>::type;
  using Base = Arr_overlay_ss_visitor<Overlay_helper, Overlay_traits, Visitor>;

public:
  using Arrangement_red_2 = typename Base::Arrangement_red_2;
  using Arrangement_blue_2 = typename Base::Arrangement_blue_2;
  using Arrangement_2 = typename Base::Arrangement_2;

  /*! Constructor */
  Arr_do_intersect_overlay_ss_visitor(const Arrangement_red_2* red_arr,
                                      const Arrangement_blue_2* blue_arr,
                                      Arrangement_2* res_arr,
                                      Overlay_traits* overlay_traits) :
    Base(red_arr, blue_arr, res_arr, overlay_traits)
  {}

  /*! Destructor */
  virtual ~Arr_do_intersect_overlay_ss_visitor() {}
};

} // namespace CGAL

#endif
