
#ifndef CGAL_TRIANGULATION_HASH_FUNCTIONS_H
#define CGAL_TRIANGULATION_HASH_FUNCTIONS_H

#include <CGAL/config.h>
#include <boost/functional/hash.hpp>
#include <CGAL/Mesh_3/comparison_operators.h>
#include <CGAL/array.h>
#include <CGAL/Triangulation_utils_3.h>


namespace CGAL {
  namespace Mesh_3 {

  template<typename Tr>
  struct Hash_facet
    : public std::unary_function<typename Tr::Facet, std::size_t>
  {
    typedef typename Tr::Point Point;
    typedef typename Tr::Facet Facet;
    typedef typename Tr::Vertex_handle Vertex_handle;

    std::size_t operator()(const Facet& f) const
    {
      CGAL_PROFILER("Hash facet");
      CGAL::cpp11::array<Vertex_handle,3> v;
      for(int i = 0; i < 3; ++i)
        v[i] = f.first->vertex(
            Triangulation_utils_3::vertex_triple_index(f.second,i));
      
      Vertex_handle_comparator<Vertex_handle> vcomp;
      sort3(v[0], v[1], v[2], vcomp);
      
      std::size_t seed = 0;
      for(std::size_t i = 0; i < 3; ++i)
      {
        const Point& p = v[i]->point();
        boost::hash_combine(seed, p.x());
        boost::hash_combine(seed, p.y());
        boost::hash_combine(seed, p.z());
      }
      return seed;
    }
  };

  template<typename Tr>
  struct Hash_cell
    : public std::unary_function<typename Tr::Cell_handle, std::size_t>
  {
    typedef typename Tr::Point Point;
    typedef typename Tr::Cell_handle Cell_handle;
    typedef typename Tr::Vertex_handle Vertex_handle;

    std::size_t operator()(const Cell_handle& c) const
    {
      CGAL_PROFILER("Hash cell");
      CGAL::cpp11::array<Vertex_handle,4> v;
      for(int i = 0; i < 4; ++i)
        v[i] = c->vertex(i);
        
      Vertex_handle_comparator<Vertex_handle> vcomp;
      sort4(v[0], v[1], v[2], v[3], vcomp);

      std::size_t seed = 0;
      for(std::size_t i = 0; i < 4; ++i)
      { 
        const Point& p = v[i]->point();
        boost::hash_combine(seed, p.x());
        boost::hash_combine(seed, p.y());
        boost::hash_combine(seed, p.z());
      }
      return seed;
    }
  };


  }
}

#endif //CGAL_TRIANGULATION_HASH_FUNCTIONS_H
