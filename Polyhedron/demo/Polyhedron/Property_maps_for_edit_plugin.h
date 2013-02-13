#ifndef PROPERTY_MAPS_FOR_EDIT_PLUGIN
#define PROPERTY_MAPS_FOR_EDIT_PLUGIN

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
namespace boost {

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

/**
 * A custom property map for deformation which behaves like zero initialized pmap.
 * It is designed for deforming very large polyhedrons with small ROS by not requiring any 
 * actual initialization. 
 * It just stores ids for ROS vertices/edges and when any other vertex/edge is queried, 0 will be returned.
 */
template<class Key>
class Polyhedron_zero_default_index_map
{
public:
  typedef boost::read_write_property_map_tag  category;
  typedef std::size_t                         value_type;
  typedef std::size_t                         reference;
  typedef Key                                 key_type;

  Polyhedron_zero_default_index_map(std::map<key_type, size_t>& internal_map)
  : internal_map(&internal_map) { }

  std::map<key_type, size_t>* internal_map;
};

namespace boost {

  template<class Key> std::size_t get( Polyhedron_zero_default_index_map<Key>& pmap
                                   , typename Key k)
  {
    std::map<Key, size_t>::iterator found = pmap.internal_map->find(k);
    // if the key doesn't exist in the map, then retun 0 to simulate zero initialization
    if(found == pmap.internal_map->end()) { return 0; }
    return found->second;
  }

  template<class Key> void put( Polyhedron_zero_default_index_map<Key>& pmap
                            , typename Key k, std::size_t s)
  {
    // also provide cleaning facility (it will be useful when ROS is cleaned and another ROS is added to the system)
    if(s == 0) { pmap.internal_map->erase(k); }
    else       { pmap.internal_map->operator[](k) = s; }
  }
}

#endif //PROPERTY_MAPS_FOR_EDIT_PLUGIN
