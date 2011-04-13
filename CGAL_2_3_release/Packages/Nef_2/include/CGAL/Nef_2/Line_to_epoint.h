// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_2/Line_to_epoint.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// source        : nef_2d/Simple_extended_kernel.lw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Simple converter
// ============================================================================

#ifndef CGAL_LINE_TO_EPOINT_H
#define CGAL_LINE_TO_EPOINT_H

CGAL_BEGIN_NAMESPACE

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

  static FT ordinate_distance(const Line_2& l)
  { return Kernel::make_FT(-l.c(),l.b()); }

  static Point_type determine_type(const Line_2& l)
  {
    RT adx = CGAL_NTS abs(dx(l)), ady = CGAL_NTS abs(dy(l));
    int sdx = CGAL_NTS sign(dx(l)), sdy = CGAL_NTS sign(dy(l));
    int cmp_dx_dy = CGAL_NTS compare(adx,ady), s(1);
    if (sdx < 0 && ( cmp_dx_dy > 0 || cmp_dx_dy == 0 && 
        sdy != (s = CGAL_NTS sign(ordinate_distance(l))))) {
      if (0 == s) return ( sdy < 0 ? SWCORNER : NWCORNER );
      else        return LEFTFRAME;
    } else if (sdx > 0 && ( cmp_dx_dy > 0 || cmp_dx_dy == 0 && 
               sdy != (s = CGAL_NTS sign(ordinate_distance(l))))) { 
      if (0 == s) return ( sdy < 0 ? SECORNER : NECORNER );
      else        return RIGHTFRAME;
    } else if (sdy < 0 && ( cmp_dx_dy < 0 || cmp_dx_dy == 0 && 
               ordinate_distance(l) < FT(0))) {
      return BOTTOMFRAME;
    } else if (sdy > 0 && ( cmp_dx_dy < 0 || cmp_dx_dy == 0 && 
               ordinate_distance(l) > FT(0))) {
      return TOPFRAME;
    }
    CGAL_assertion_msg(false," determine_type: degenerate line.");
    return (Point_type)-1; // never come here
  }


};

CGAL_END_NAMESPACE
#endif //CGAL_LINE_TO_EPOINT_H

