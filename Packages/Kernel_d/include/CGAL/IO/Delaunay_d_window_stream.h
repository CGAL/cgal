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
// file          : include/CGAL/Delaunay_d_window_stream.h
// package       : Kernel_d
// chapter       : Basic
//
// source        : ddgeo/Delaunay_d.lw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Susan Hert <hert@mpi-sb.mpg.de>
//
// implementation: Higher dimensional geometry
// ============================================================================
#ifndef CGAL_DELAUNAY_D_WINDOW_STREAM_H
#define CGAL_DELAUNAY_D_WINDOW_STREAM_H

#include <CGAL/Delaunay_d.h>
#include <CGAL/IO/Convex_hull_d_window_stream.h>
#include <CGAL/IO/Window_stream.h>

CGAL_BEGIN_NAMESPACE

/*{\Mtext \headerline{Low Dimensional Output Routines}
include |<CGAL/IO/Delaunay_d_window_stream.h>|
\setopdims{2cm}{1cm}}*/

template <typename R, typename Lifted_R>
void d2_show(const Delaunay_d<R,Lifted_R>& D,
             CGAL::Window_stream& W, 
             typename Delaunay_d<R,Lifted_R>::Delaunay_voronoi_kind k = 
             Delaunay_d<R,Lifted_R>::NEAREST)
/*{\Mfunc draws the underlying simplicial complex |D| into window |W|.\\
\precond |dim == 2|. }*/
{ 
  CGAL_assertion_msg(D.dimension() == 2, "d2_map: dim != 2.");
  Regular_complex_d<R> RC(2);
  D.project(RC, (k == Delaunay_d<R,Lifted_R>::NEAREST ? -1 : +1));
  CGAL::d2_show(RC,W);
}

template <typename R, typename Lifted_R>
void d2_map(const Delaunay_d<R,Lifted_R>& D, 
            GRAPH< typename Delaunay_d<R,Lifted_R>::Point_d, int >& DTG, 
            typename Delaunay_d<R,Lifted_R>::Delaunay_voronoi_kind k = 
            Delaunay_d<R,Lifted_R>::NEAREST)
/*{\Mfunc constructs a LEDA graph representation of the nearest 
(|kind = NEAREST| or the furthest (|kind = FURTHEST|) site
Delaunay triangulation.\\ \precond |dim() == 2|. }*/
{ 
  CGAL_assertion_msg(D.dimension() == 2, "d2_map: dim != 2.");
  Regular_complex_d<R> RC(2);
  D.project(RC, (k == NEAREST ? -1 : +1));
  d2_map(RC,DTG);
}


CGAL_END_NAMESPACE
#endif //CGAL_DELAUNAY_D_WINDOW_STREAM_H


