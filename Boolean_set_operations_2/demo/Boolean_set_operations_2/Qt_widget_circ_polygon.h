// Copyright (c) 2005  Tel-Aviv University (Israel).
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

#ifndef CGAL_QT_WIDGET_CIRC_POLYGON_H
#define CGAL_QT_WIDGET_CIRC_POLYGON_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <utility>

namespace CGAL
{
  template <class Tr>
  Qt_widget& operator<<(Qt_widget& w, const General_polygon_2<Tr>& pgn)
  {
    typedef typename General_polygon_2<Tr>::Curve_const_iterator  CI;
    typedef Simple_cartesian<double>    DK;
    typedef DK::Point_2                 DPT;
    typedef CGAL::Polygon_2<DK>         DPGN;

    std::list<std::pair<double, double> > pair_list;
    for(CI citr = pgn.curves_begin();citr != pgn.curves_end(); ++citr)
    {
      if(citr->is_linear())
      {
        // when the curve is linear approximate will allways return
        // two pairs (for each endpoint) regardless the parameter of number
        // of points
        citr->approximate(std::back_inserter(pair_list), 0);
        continue;
      }

      // circular arc
      double sx = CGAL::to_double(citr->source().x());
      double tx = CGAL::to_double(citr->target().x());
      int x_min;
      int x_max;
      if(citr->is_directed_right())
      {
        x_min =  w.x_pixel(sx);
        x_max =  w.x_pixel(tx);
      }
      else
      {
        x_min = w.x_pixel(tx);
        x_max = w.x_pixel(sx);
      }
      const int     n = x_max - x_min + 1;
      if (n <= 0)
        continue;

      citr->approximate(std::back_inserter(pair_list), n);
    }

    DPGN app_pgn;
    for(std::list<std::pair<double, double> >::iterator it = pair_list.begin();
        it != pair_list.end();
        ++it)
    {
      DPT pt(it->first, it->second);
      app_pgn.push_back(pt);
    }

    w<<app_pgn;
    return w;
  }


}//end namespace CGAL


#endif
