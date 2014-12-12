#include <fstream>

#include <boost/graph/graph_traits.hpp>
#include <CGAL/property_map.h>

#include <CGAL/basic.h>
#include <CGAL/Kernel/global_functions.h>

// Polyhedron
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

// Graph traits adaptors
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

typedef CGAL::Cartesian<double>                                      Kernel;
typedef Kernel::Point_3                                              Point;
typedef Kernel::Vector_3                                             Vector;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

template<typename HalfedgeGraph, 
         typename PointMap, 
         typename NormalMap>
void calculate_face_normals(const HalfedgeGraph& g, 
                            PointMap pm, 
                            NormalMap nm) 
{
  typedef boost::graph_traits<HalfedgeGraph> GraphTraits;
  typedef typename GraphTraits::face_iterator face_iterator;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::property_traits<PointMap>::value_type point;
  typedef typename boost::property_traits<NormalMap>::value_type normal;

  face_iterator fb, fe;
  for(boost::tie(fb, fe) = faces(g); fb != fe; ++fb)
  {
    halfedge_descriptor edg = halfedge(*fb, g);
    halfedge_descriptor edgb = edg;

    point p0 = pm[target(edg, g)];
    edg = next(edg, g);
    point p1 = pm[target(edg, g)];
    edg = next(edg, g);
    point p2 = pm[target(edg, g)];
    edg = next(edg, g);
      
    if(edg == edgb) {
      // triangle
      nm[*fb] = CGAL::unit_normal(p1, p2, p0);
    } else {
      // not a triangle
      normal n(CGAL::NULL_VECTOR);

      do {
        n = n + CGAL::normal(p1, p2, p0);
        p0 = p1;
        p1 = p2;
        
        edg = next(edg, g);
        p2 = pm[target(edg, g)];
      } while(edg != edgb);
      
      nm[*fb] = n / CGAL::sqrt(n.squared_length());
    }
  }
}

int main(int, char** argv)
{
  typedef boost::property_map< 
    Polyhedron,
    CGAL::face_index_t 
    >::const_type Face_index_map;

  std::ifstream in(argv[1]);
  Polyhedron P;
  in >> P ;
  
  // initialize facet indices
  std::size_t i = 0;
  for(Polyhedron::Facet_iterator it = P.facets_begin(); it != P.facets_end(); ++it, ++i) 
  {
    it->id() = i;
  }

  // Ad hoc property_map to store normals. Face_index_map is used to
  // map face_descriptors to a contiguous range of indices. See
  // http://www.boost.org/libs/property_map/doc/vector_property_map.html
  // for details.
  boost::vector_property_map<Vector, Face_index_map> 
    normals(get(CGAL::face_index, P));

  calculate_face_normals(
    P // Graph
    , get(CGAL::vertex_point, P) // map from vertex_descriptor to point
    , normals // map from face_descriptor to Vector_3
    );

  std::cout << "Normals" << std::endl;
  for(Polyhedron::Facet_iterator it = P.facets_begin(); it != P.facets_end(); ++it) { 
    // Facet_iterator is a face_descriptor, so we can use it as the
    // key here
    std::cout << normals[it] << std::endl;
  }

  return 0;
}
