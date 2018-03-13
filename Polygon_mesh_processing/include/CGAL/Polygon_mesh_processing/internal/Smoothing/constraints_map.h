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


  template<typename PolygonMesh, typename VertexPointMap,
           typename CotangentValue = CGAL::internal::Cotangent_value_Meyer<PolygonMesh, VertexPointMap>>
  struct Edge_cotangent_weight : CotangentValue
  {
      Edge_cotangent_weight(PolygonMesh& pmesh_, VertexPointMap vpmap_)
        : CotangentValue(pmesh_, vpmap_)
      {}

      PolygonMesh& pmesh()
      {
        return CotangentValue::pmesh();
      }

      typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;
      typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor     vertex_descriptor;

      double operator()(halfedge_descriptor he)
      {
        if(is_border_edge(he, pmesh()))
        {
          halfedge_descriptor h1 = next(he, pmesh());
          vertex_descriptor vs = source(he, pmesh());
          vertex_descriptor vt = target(he, pmesh());
          vertex_descriptor v1 = target(h1, pmesh());
          return (CotangentValue::operator ()(vs, v1, vt));
        }
        else
        {
          halfedge_descriptor h1 = next(he, pmesh());
          halfedge_descriptor h2 = prev(opposite(he, pmesh()), pmesh());
          vertex_descriptor vs = source(he, pmesh());
          vertex_descriptor vt = target(he, pmesh());
          vertex_descriptor v1 = target(h1, pmesh());
          vertex_descriptor v2 = source(h2, pmesh());
          return ( CotangentValue::operator()(vs, v1, vt) + CotangentValue::operator()(vs, v2, vt) ) / 2.0;
        }
      }
  };

  template<typename PolygonMesh>
  struct Incident_area
  {
    Incident_area(PolygonMesh& mesh) : pmesh(mesh){}

    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;

    double operator()(halfedge_descriptor he)
    {
      halfedge_descriptor hopp = opposite(he, pmesh);
      face_descriptor f1 = face(he, pmesh);
      face_descriptor f2 = face(hopp, pmesh);

      double A1 = f1 == boost::graph_traits<PolygonMesh>::null_face() ? 0 : face_area(f1, pmesh);
      double A2 = f2 == boost::graph_traits<PolygonMesh>::null_face() ? 0 : face_area(f2, pmesh);
      return A1 + A2;
    }
    PolygonMesh& pmesh;
  };






}
}
}

#endif //CGAL_POLYGON_MESH_PROCESSING_CONSTRAINTS_MAP_H
