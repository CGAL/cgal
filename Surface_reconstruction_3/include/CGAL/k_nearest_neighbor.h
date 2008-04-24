#ifndef CGAL_K_NEIGHBOR_NEIGHBOR_H
#define CGAL_K_NEIGHBOR_NEIGHBOR_H

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_vertex_handle_3.h>
#include <list>

CGAL_BEGIN_NAMESPACE


/// Wrapper around Orthogonal_k_neighbor_search for Triangulation_3 Vertex_handles.
template <class Gt, class Vertex_handle>
class K_nearest_neighbor
{
public:
  typedef Gt  Geom_traits;

  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point;
  typedef Point_vertex_handle_3<Vertex_handle> Point_vertex_handle_3;
  typedef Search_traits_vertex_handle_3<Vertex_handle> Traits;
  typedef Euclidean_distance_vertex_handle_3<Vertex_handle> KDistance;
  typedef Orthogonal_k_neighbor_search<Traits,KDistance> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;

private:
  Tree m_tree;

public:
  K_nearest_neighbor() {}

  template <class InputIterator> ///< InputIterator value_type is Vertex_handle.
  K_nearest_neighbor(InputIterator first, InputIterator beyond)
  {
    m_tree = Tree(first, beyond);
  }

  /// Default copy constructor, operator =() and destructor are fine.

  bool get_k_nearest_neighbors(const Point_vertex_handle_3& query,
                               const unsigned int nb,
                               std::list<Vertex_handle>& kvertices)
  {
    Neighbor_search search(m_tree,query,nb); // only nb nearest neighbors
    Search_iterator it = search.begin();
    for(unsigned int i=0;i<nb;i++,it++)
    {
      if(it == search.end())
        return false;
      kvertices.push_back((Vertex_handle)it->first);
    }
    return true;
  }

}; // K_nearest_neighbor


CGAL_END_NAMESPACE

#endif // CGAL_K_NEIGHBOR_NEIGHBOR_H

