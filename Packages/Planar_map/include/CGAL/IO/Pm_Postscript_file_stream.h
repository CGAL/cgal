// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// Author(s)     : Eti Ezra <estere@post.tau.ac.il>

#ifndef CGAL_IO_PM_POSTSCRIPT_FILE_STREAM_H
#define CGAL_IO_PM_POSTSCRIPT_FILE_STREAM_H

#include <CGAL/IO/leda_window.h>
#include <CGAL/IO/Postscript_file_stream.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/IO/Pm_drawer.h>
#include  <CGAL/IO/draw_pm.h>

CGAL_BEGIN_NAMESPACE

template <class Dcel,class Traits>
Postscript_file_stream & operator << (Postscript_file_stream & ps,
                                      const Planar_map_2<Dcel,Traits> & pm)
{
  Pm_drawer< Planar_map_2<Dcel,Traits> , Postscript_file_stream>  drawer(ps);
  draw_pm(pm, drawer, ps);
  return ps;
}  

CGAL_END_NAMESPACE

#endif
