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
// release       : $CGAL_Revision: CGAL-2.3-I-75 $
// release_date  : $CGAL_Date: 2001/06/21 $
//
// file          : include/CGAL/LEDA/pixmaps/texture.h
// package       : cgal_window (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.9.8
// revision_date : 12 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>)
// ======================================================================
#include <CGAL/LEDA/pixmaps/texture/blue.xpm>
#include <CGAL/LEDA/pixmaps/texture/clouds.xpm>
#include <CGAL/LEDA/pixmaps/texture/dirt.xpm>
#include <CGAL/LEDA/pixmaps/texture/pool.xpm>
#include <CGAL/LEDA/pixmaps/texture/sharks.xpm>
#include <CGAL/LEDA/pixmaps/texture/slate.xpm>
#include <CGAL/LEDA/pixmaps/texture/space.xpm>
#include <CGAL/LEDA/pixmaps/texture/wall.xpm>
#include <CGAL/LEDA/pixmaps/texture/waves.xpm>
#include <CGAL/LEDA/pixmaps/texture/wood.xpm>

#include <CGAL/LEDA/pixmaps/leda.xpm>
#include <CGAL/LEDA/pixmaps/algosol.xpm>

#define num_texture 12

static const char** xpm_texture[] = {
algosol_xpm,
leda_xpm,
blue_xpm,
clouds_xpm,
dirt_xpm,
pool_xpm,
sharks_xpm,
slate_xpm,
space_xpm,
wall_xpm,
waves_xpm,
wood_xpm
};

static const char* name_texture[] = {
"algosol",
"leda",
"blue",
"clouds",
"dirt",
"pool",
"sharks",
"slate",
"space",
"wall",
"waves",
"wood"
};

#if defined(__GNUC__)
inline char texture_unused_warning()
{ return xpm_texture[0][0][0] + name_texture[0][0]; }
#endif


