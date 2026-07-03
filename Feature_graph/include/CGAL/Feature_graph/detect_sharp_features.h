// Copyright (c) 2026 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ange Clement

#ifndef CGAL_FEATURE_GRAPH_DETECT_SHARP_FEATURES_H
#define CGAL_FEATURE_GRAPH_DETECT_SHARP_FEATURES_H

#include <boost/graph/adjacency_list.hpp>

#include <CGAL/license/Feature_graph.h>

namespace CGAL
{

namespace Feature_graph
{

template <typename K>
struct Detect_sharp_features_on_labeled_image
{
  typedef typename K::Point_3 Point_3;
  typedef std::size_t Size;
  typedef typename K::FT FT;

  typedef Point_3 VertexProperty;
  typedef boost::adjacency_list<boost::vecS,boost::listS,boost::undirectedS,VertexProperty> Feature_graph;

  template <typename Image, typename CGAL_NP_TEMPLATE_PARAMETERS>
  Feature_graph operator()(const Image& image, const CGAL_NP_CLASS& np = parameters::default_values()) const
  {
    // do nothing for now
    return Feature_graph();
  }
};

} // end namespace Feature_graph

} // end namespace CGAL

#endif // CGAL_FEATURE_GRAPH_DETECT_SHARP_FEATURES_H