// Copyright (c) 1997-2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Radu Ursu

#ifndef CGAL_QT_WIDGET_OPTIMISATION_CIRCLE_2_H
#define CGAL_QT_WIDGET_OPTIMISATION_CIRCLE_2_H

#include <CGAL/IO/Qt_widget.h>

namespace CGAL{
template<class Traits>
Qt_widget&
operator << (Qt_widget &ws,
	     const CGAL::Optimisation_circle_2<Traits>& oc){
  typedef typename Traits::Point_2  Point_2;
  typedef typename Traits::Circle_2 Circle_2;

  double cx( CGAL::to_double( oc.center().x()));
  double cy( CGAL::to_double( oc.center().y()));
  double sr( CGAL::to_double( oc.squared_radius()));

  if( ! CGAL_NTS is_negative(sr))
    ws << Circle_2( Point_2(cx, cy), sr);
  return ws;
}

}//end namespace

#endif
