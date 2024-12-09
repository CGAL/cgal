// Copyright (c) 2021  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_BOOST_GRAPH_BREADTH_FIRST_SEARCH_H
#define CGAL_BOOST_GRAPH_BREADTH_FIRST_SEARCH_H

// This will push/pop a VC++ warning
#include <CGAL/Named_function_parameters.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4172) // Address warning inside boost named parameters
#endif

#include <boost/graph/breadth_first_search.hpp>

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_BOOST_GRAPH_BREADTH_FIRST_SEARCH_H
