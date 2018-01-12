#ifndef CGAL_POLYGON_MESH_PROCESSING_CONSTRAINTS_MAP_H
#define CGAL_POLYGON_MESH_PROCESSING_CONSTRAINTS_MAP_H

#include <CGAL/property_map.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

  template<typename Descriptor>
  struct Constrained_vertices_map
  {
    typedef Descriptor                          key_type;
    typedef bool                                value_type;
    typedef value_type&                         reference;
    typedef boost::read_write_property_map_tag  category;

    // to change this to boost::shared_ptr
    std::shared_ptr<std::set<Descriptor>> const_things;

  public:
    Constrained_vertices_map() : const_things(new std::set<Descriptor>) {}

    friend bool get(const Constrained_vertices_map& map, const key_type& d)
    {
      typename std::set<Descriptor>::iterator it = map.const_things->find(d);
      return it != map.const_things->end() ? true : false;
    }

    friend void put(Constrained_vertices_map& map, const key_type& d)
    {
      map.const_things->insert(d);
    }
  };

}
}
}

#endif //CGAL_POLYGON_MESH_PROCESSING_CONSTRAINTS_MAP_H
