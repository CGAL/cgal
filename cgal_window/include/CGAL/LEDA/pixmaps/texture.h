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


