// Copyright (c) 2001  Martin-Luther-University Halle-Wittenberg (Germany).
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
// Author(s)     : Matthias Baesken, Algorithmic Solutions

#include <CGAL/LEDA/bitmaps/button32/a_annulus.xbm>
#include <CGAL/LEDA/bitmaps/button32/w_annulus.xbm>
#include <CGAL/LEDA/bitmaps/button32/circle.xbm>
#include <CGAL/LEDA/bitmaps/button32/dice.xbm>
#include <CGAL/LEDA/bitmaps/button32/empty_circle.xbm>
#include <CGAL/LEDA/bitmaps/button32/encl_circle.xbm>
#include <CGAL/LEDA/bitmaps/button32/exit.xbm>
#include <CGAL/LEDA/bitmaps/button32/f_triang.xbm>
#include <CGAL/LEDA/bitmaps/button32/f_voro.xbm>
#include <CGAL/LEDA/bitmaps/button32/grid.xbm>
#include <CGAL/LEDA/bitmaps/button32/help.xbm>
#include <CGAL/LEDA/bitmaps/button32/hull.xbm>
#include <CGAL/LEDA/bitmaps/button32/inside.xbm>
#include <CGAL/LEDA/bitmaps/button32/intersect.xbm>
#include <CGAL/LEDA/bitmaps/button32/line.xbm>
#include <CGAL/LEDA/bitmaps/button32/point.xbm>
#include <CGAL/LEDA/bitmaps/button32/poly.xbm>
#include <CGAL/LEDA/bitmaps/button32/rect.xbm>
#include <CGAL/LEDA/bitmaps/button32/tree.xbm>
#include <CGAL/LEDA/bitmaps/button32/triang.xbm>
#include <CGAL/LEDA/bitmaps/button32/triangle.xbm>
#include <CGAL/LEDA/bitmaps/button32/union.xbm>
#include <CGAL/LEDA/bitmaps/button32/voro.xbm>

#define num_xbm_button32 23

static unsigned char* xbm_button32[] = {
a_annulus_bits,
w_annulus_bits,
circle_bits,
dice_bits,
empty_circle_bits,
encl_circle_bits,
exit_bits,
f_triang_bits,
f_voro_bits,
grid_bits,
help_bits,
hull_bits,
inside_bits,
intersect_bits,
line_bits,
point_bits,
poly_bits,
rect_bits,
tree_bits,
triang_bits,
triangle_bits,
union_bits,
voro_bits
};


static const char* name_xbm_button32[] = {
"a_annulus",
"w_annulus",
"circle",
"dice",
"empty_circle",
"encl_circle",
"exit",
"f_triang",
"f_voro",
"grid",
"help",
"hull",
"inside",
"intersect",
"line",
"point",
"poly",
"rect",
"tree",
"triang",
"triangle",
"union",
"voro"
};


#if defined(__GNUC__)
inline char xbm_button32_unused_warning()
{ return xbm_button32[0][0] + name_xbm_button32[0][0]; }
#endif
