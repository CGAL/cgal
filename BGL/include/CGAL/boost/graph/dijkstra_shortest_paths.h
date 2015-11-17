// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Sebastien Loriot


#ifndef CGAL_BOOST_GRAPH_DIJKSTRA_SHORTEST_PATHS_H
#define CGAL_BOOST_GRAPH_DIJKSTRA_SHORTEST_PATHS_H

// This will push/pop a VC15 warning
#include <CGAL/boost/graph/named_function_params.h>

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
