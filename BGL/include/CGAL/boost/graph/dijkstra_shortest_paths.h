// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot


#ifndef CGAL_BOOST_GRAPH_DIJKSTRA_SHORTEST_PATHS_H
#define CGAL_BOOST_GRAPH_DIJKSTRA_SHORTEST_PATHS_H

// This will push/pop a VC15 warning
#include <CGAL/boost/graph/Named_function_parameters.h>

#include <boost/version.hpp>
#include <climits>

#if BOOST_VERSION == 105400
  #ifdef BOOST_GRAPH_DIJKSTRA_HPP
  #    pragma message \
      "Warning: the header file boost/graph/dijkstra_shortest_paths.hpp "       \
      "of boost 1.54 contains a bug that may impact some functions in CGAL. "   \
      "Please consider including CGAL/boost/graph/dijkstra_shortest_paths.hpp "  \
      "before boost header"
  #endif
  #include <CGAL/boost/graph/dijkstra_shortest_paths.hpp>
#else
  #include <boost/graph/dijkstra_shortest_paths.hpp>
#endif

#endif // CGAL_BOOST_GRAPH_DIJKSTRA_SHORTEST_PATHS_H
