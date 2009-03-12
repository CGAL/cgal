// Copyright (c) 2007-2008  INRIA (France).
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
// $URL: svn+ssh://ggael@scm.gforge.inria.fr/svn/cgal/trunk/Surface_reconstruction_3/include/CGAL/Neighborhood/knn_point_neighbor_search.h $
// $Id: knn_point_neighbor_search.h 47904 2009-02-04 18:46:54Z ggael $
//
//
// Author(s)     : Gael Guennebaud

#ifndef CGAL_KNN_POINT_NEIGHBOR_SEARCH_H
#define CGAL_KNN_POINT_NEIGHBOR_SEARCH_H

#include <vector>
#include "point_kdtree_3.h"

namespace CGAL {

template <class SearchTraits>
class knn_point_neighbor_search {

public:

  typedef PointKdtree_3<SearchTraits>  Tree;
  typedef typename SearchTraits::Point_d Point_d;
  typedef typename SearchTraits::FT FT;
  typedef std::pair<Point_d,FT> Point_with_transformed_distance;
  typedef std::pair<Point_d,FT> Neighbor;
  typedef std::vector<Neighbor> NN_list;

public:

  typedef typename NN_list::const_iterator iterator;

  inline iterator begin() const { return m_list.begin(); }

  inline iterator end() const { return m_list.end(); }

  // constructor
  knn_point_neighbor_search(Tree& tree, const Point_d& q, unsigned int k=1) : m_list(k)
  {
    KnnNeighborhood<SearchTraits> knn(k);
    tree.doQueryK(q,knn);
    int r = knn.size()-1;
    for (int i=0; i<knn.size(); ++i)
    {
      m_list[r-i] = Neighbor(tree.object(knn.index(i)), knn.value(i));
    }
  }

  ~knn_point_neighbor_search() {}

private:

  NN_list m_list;

}; // class

} // namespace CGAL

#endif  // CGAL_KNN_POINT_NEIGHBOR_SEARCH_H
