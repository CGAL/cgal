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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef CGAL_DELAUNAY_D_WINDOW_STREAM_H
#define CGAL_DELAUNAY_D_WINDOW_STREAM_H

#include <CGAL/LEDA_basic.h>
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
            CGAL_LEDA_SCOPE::GRAPH< typename Delaunay_d<R,Lifted_R>::Point_d, 
                                    int >& DTG, 
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


