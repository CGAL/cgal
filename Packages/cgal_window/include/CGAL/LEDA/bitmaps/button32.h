// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 2001, May 23
//
// file          : ?/include/CGAL/LEDA/bitmaps/button32.h
// package       : cgal_window (0.9.7)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// author(s)     : ?
// coordinator   : ?
//
// ======================================================================

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
