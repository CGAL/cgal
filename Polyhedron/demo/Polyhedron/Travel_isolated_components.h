#ifndef TRAVEL_ISOLATED_COMPONENTS_H
#define TRAVEL_ISOLATED_COMPONENTS_H

#include <boost/optional.hpp>
#include <vector>
#include "One_ring_iterators.h"
#include <CGAL/Default.h>


template<typename Mesh>
struct Id_getter{
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor mesh_vd;
  typedef typename boost::graph_traits<Mesh>::face_descriptor mesh_fd;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor mesh_hd;
  typedef typename boost::graph_traits<Mesh>::edge_descriptor mesh_ed;
  const Mesh& mesh;
  typename boost::property_map<Mesh, boost::vertex_index_t >::type vmap;
  typename boost::property_map<Mesh, boost::face_index_t >::type fmap;
  typename boost::property_map<Mesh, boost::halfedge_index_t >::type hmap;
  typename boost::property_map<Mesh, boost::edge_index_t >::type emap;

  Id_getter(const Mesh& mesh)
    :mesh(mesh)
  {
    vmap = get(boost::vertex_index, mesh);
    fmap = get(boost::face_index, mesh);
    hmap = get(boost::halfedge_index, mesh);
    emap = get(boost::edge_index, mesh);
  }

  std::size_t id(mesh_vd vd)
  {
    return get(vmap, vd);
  }
  std::size_t id(mesh_fd fd)
  {
    return get(fmap, fd);
  }
  std::size_t id(mesh_hd hd)
  {
    return get(hmap, hd);
  }
  std::size_t id(mesh_ed ed)
  {
    return get(emap, ed);
  }
};

template<typename Mesh>
class Travel_isolated_components {
  const Mesh& mesh;
  Id_getter<Mesh> id_getter;
public:
  Travel_isolated_components(const Mesh& mesh)
    :mesh(mesh), id_getter(mesh){}
  // for transform iterator
  template<class Descriptor>
  struct Get_handle {
    typedef Descriptor result_type;
    template<class Iterator>
    result_type operator()(Iterator it) const
    { return it; }
  };

  // to be used in get_minimum_isolated_component function
  struct Minimum_visitor
  {
    Minimum_visitor()
      : minimum(-1)
    {}
    template<class Descriptor>
    void operator()(const std::vector<Descriptor>& C)
    {
      minimum = (std::min)(minimum, C.size());
    }

    std::size_t minimum;
  };

  // to be used in select_isolated_components function
  template<class OutputIterator>
  struct Selection_visitor
  {
    Selection_visitor(std::size_t threshold_size, OutputIterator out)
      : threshold_size(threshold_size), out(out), any_inserted(false) { }

    template<class Descriptor>
    void operator()(const std::vector<Descriptor>& C) {
      if(C.size() <= threshold_size) {
        any_inserted = true;
        out = std::copy(C.begin(), C.end(), out);
      }
      else {
        minimum_visitor(C);
      }
    }

    std::size_t     threshold_size;
    OutputIterator  out;
    bool            any_inserted;
    Minimum_visitor minimum_visitor; // hold minimum of NOT inserted components
  };

  template <class Handle>
  std::size_t id(Handle h) { return id_getter.id(h); }
  /* AF shouldn't the edge property map do that right?
  std::size_t id(boost::graph_traits<Polyhedron>::edge_descriptor ed)
  {
    return id_getter.id(halfedge(ed, id_getter.mesh)) /2;
  }
  */
  // NOTE: prior to call this function, id fields should be updated
  template<class Descriptor, class InputIterator, class IsSelected, class Visitor>
  void travel(InputIterator begin,
              InputIterator end,
              std::size_t size,
              const IsSelected& selection,
              Visitor& visitor)
  {
    std::vector<bool> mark(size, false);

    for(; begin != end; ++begin)
    {
      Descriptor h = *begin;

      if(mark[id(h)] || selection.count(h)) { continue; }

      std::vector<Descriptor> C;
      C.push_back(h);
      mark[id(h)] = true;
      std::size_t current_index = 0;

      bool neigh_to_selection = false;
      while(current_index < C.size()) {
        Descriptor current = C[current_index++];

        for(One_ring_iterator<Mesh, Descriptor> circ(current, mesh); circ; ++circ)
        {
          Descriptor nv = circ;
          neigh_to_selection |= (selection.count(nv)!=0);
          if(!mark[id(nv)] && !selection.count(nv)) {
            mark[id(nv)] = true;
            C.push_back(nv);
          }
        }
      }
      if(neigh_to_selection) { visitor(C); }
    }
  }
};
#endif
