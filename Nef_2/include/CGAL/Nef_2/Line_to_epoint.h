// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_LINE_TO_EPOINT_H
#define CGAL_LINE_TO_EPOINT_H

#include <CGAL/license/Nef_2.h>


namespace CGAL {

template <class Kernel_>
struct Line_to_epoint {
  typedef Kernel_ Kernel;
  typedef typename Kernel::RT RT;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Line_2 Line_2;
  enum Point_type { SWCORNER=1, LEFTFRAME, NWCORNER, 
                    BOTTOMFRAME, STANDARD, TOPFRAME,
                    SECORNER, RIGHTFRAME, NECORNER };


  static RT dx(const Line_2& l) { return l.b(); }
  static RT dy(const Line_2& l) { return -l.a(); }

  static FT ordinate_distance(const Line_2& l) { 
    return FT(-l.c()) / l.b();
  }

  static Point_type determine_type(const Line_2& l)
  {
    RT adx = CGAL_NTS abs(dx(l)), ady = CGAL_NTS abs(dy(l));
    int sdx = CGAL_NTS sign(dx(l)), sdy = CGAL_NTS sign(dy(l));
    int cmp_dx_dy = CGAL_NTS compare(adx,ady), s(1);
    if ( (sdx < 0) && ( (cmp_dx_dy > 0) || ( (cmp_dx_dy == 0) && 
                                             (sdy != (s = CGAL_NTS sign(ordinate_distance(l))))))) {
      if (0 == s) return ( sdy < 0 ? SWCORNER : NWCORNER );
      else        return LEFTFRAME;
    } else if ( (sdx > 0) && ( (cmp_dx_dy > 0) || ( (cmp_dx_dy == 0) && 
                                                    (sdy != (s = CGAL_NTS sign(ordinate_distance(l))))))) { 
      if (0 == s) return ( sdy < 0 ? SECORNER : NECORNER );
      else        return RIGHTFRAME;
    } else if ( (sdy < 0) && ( (cmp_dx_dy < 0) || ( (cmp_dx_dy == 0) && 
                                                    (ordinate_distance(l) < FT(0))))) {
      return BOTTOMFRAME;
    } else if ( (sdy > 0) && ( (cmp_dx_dy < 0) || ( (cmp_dx_dy == 0) && 
                                                    (ordinate_distance(l) > FT(0))))) {
      return TOPFRAME;
    }
    CGAL_error_msg(" determine_type: degenerate line.");
    return (Point_type)-1; // never come here
  }


};

} //namespace CGAL
#endif //CGAL_LINE_TO_EPOINT_H
