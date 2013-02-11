#ifndef PROPERTY_MAPS_FOR_EDIT_PLUGIN
#define PROPERTY_MAPS_FOR_EDIT_PLUGIN

#include "Polyhedron_type.h"

template<class P>
class Polyhedron_vertex_deformation_index_map
{
private:

  typedef P Polyhedron ;


public:

  typedef boost::read_write_property_map_tag                                  category;
  typedef std::size_t                                                       value_type;
  typedef std::size_t                                                       reference;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor key_type;

  Polyhedron_vertex_deformation_index_map()
  {    
  }
 

};

template<class P>
std::size_t
get(  Polyhedron_vertex_deformation_index_map<P>, typename P::Vertex_handle vh)
{
  return vh->id();
}

template<class P>
void
put(  Polyhedron_vertex_deformation_index_map<P>&, typename P::Vertex_handle vh, std::size_t s)
{
  vh->id() = s;
}



template<class P>
class Polyhedron_edge_deformation_index_map
{
private:

  typedef P Polyhedron ;


public:

  typedef boost::read_write_property_map_tag                                  category;
  typedef std::size_t                                                       value_type;
  typedef std::size_t                                                       reference;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor key_type;

  Polyhedron_edge_deformation_index_map()
  {}
  
};

namespace boost {

  template<class P>
  std::size_t
    get(  Polyhedron_edge_deformation_index_map<P> /* pmap */, typename P::Halfedge_handle eh)
  {
    return eh->id();
  }

  template<class P>
  void
    put(  Polyhedron_edge_deformation_index_map<P>& /* pmap */, typename P::Halfedge_handle eh, std::size_t s)
  {
    eh->id() = s;
  }
}

template<class P>
class Polyhedron_edge_deformation_length_map
{
private:

  typedef P Polyhedron ;
  

public:
  std::vector<double>& m_lengths;

  typedef boost::read_write_property_map_tag                                  category;
  typedef double                                                       value_type;
  typedef double&                                                       reference;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor key_type;

  Polyhedron_edge_deformation_length_map(const P&,std::vector<double>& lengths)
    : m_lengths(lengths) {}
};

namespace boost {

  template<class P>
  double
    get(  Polyhedron_edge_deformation_length_map<P> pmap, typename P::Halfedge_handle eh)
  {
    CGAL_assertion( pmap.m_lengths.size() > eh->id() );
    return pmap.m_lengths[eh->id()];
  }

  template<class P>
  void
    put(  Polyhedron_edge_deformation_length_map<P>& pmap, typename P::Halfedge_handle eh, double s)
  {
    CGAL_assertion( pmap.m_lengths.size() > eh->id() );
    pmap.m_lengths[eh->id()] = s;
  }
}

#endif //PROPERTY_MAPS_FOR_EDIT_PLUGIN
