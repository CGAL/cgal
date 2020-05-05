// Copyright (c) 2015  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BOOST_GRAPH_IO_H
#define CGAL_BOOST_GRAPH_IO_H

#include <CGAL/boost/graph/IO/3MF.h>
#include <CGAL/boost/graph/IO/GOCAD.h>
#include <CGAL/boost/graph/IO/INP.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/boost/graph/IO/PLY.h>
#include <CGAL/boost/graph/IO/STL.h>
#include <CGAL/boost/graph/IO/VTK.h>
#include <CGAL/boost/graph/IO/WRL.h>

namespace CGAL{
template <class FaceGraph, typename NamedParameters>
bool read_polygon_mesh(std::istream& is,
                       FaceGraph& g,
                       const NamedParameters& np)
{
  bool ok = false;
  ok = read_OFF(is, g, np);
  if(ok)
    return true;
  g.clear();
  is.clear();//reset the error state
  is.seekg (0, is.beg);
  ok = read_OBJ(is, g, np);
  if(ok)
    return true;
  g.clear();
  is.clear();
  is.seekg (0, is.beg);
  ok = read_PLY(is, g, np);
  if(ok)
    return true;
  g.clear();
  is.clear();
  is.seekg (0, is.beg);
  ok = read_STL(is, g, np);
  if(ok)
    return true;
  g.clear();
  is.clear();
  is.seekg (0, is.beg);
  ok = read_GOCAD(is, g, np);
  return ok;
}
template <class FaceGraph>
bool read_polygon_mesh(std::istream& is,
                       FaceGraph& g)
{
  return read_polygon_mesh(is, g, parameters::all_default());
}
}//end CGAL
#endif // CGAL_BOOST_GRAPH_IO_H
