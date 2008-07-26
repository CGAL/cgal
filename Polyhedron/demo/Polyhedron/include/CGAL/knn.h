/*

USAGE

typedef CNeighbor_search<kernel> Neighbor_search;
Neighbor_search search;

std::list<Point> points;
search.init(points);

// query
Point p = nearest_point(query);

*/

#ifndef _NEIGHBOR_SEARCH_
#define _NEIGHBOR_SEARCH_
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <list>

template <class kernel>
class CNeighbor_search
{
public:
	typedef typename kernel::FT FT;
	typedef typename kernel::Point_3 Point;
  typedef typename CGAL::Search_traits_3<kernel> TreeTraits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;

private:
  Tree m_tree;

public:
  CNeighbor_search() {}
  ~CNeighbor_search() {}

  void init(const std::list<Point>& points)
  {
    m_tree = Tree(points.begin(),points.end());
  }

  FT distance_nearest_point(const Point& query)
  {
    Neighbor_search search(m_tree,query,1); // only first nearest neighbor
    return (FT)std::sqrt(CGAL_NTS to_double(search.begin()->second));
  }

  Point nearest_point(const Point& query)
  {
    Neighbor_search search(m_tree,query,1); // only first nearest neighbor
    return search.begin()->first;
  }
};
#endif // _NEIGHBOR_SEARCH_

