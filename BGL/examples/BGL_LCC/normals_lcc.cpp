#include <fstream>

#include <boost/graph/graph_traits.hpp>
#include <CGAL/property_map.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3                Point;
typedef Kernel::Vector_3               Vector;

typedef CGAL::Linear_cell_complex_traits<3, Kernel> LCC_traits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
          <2, 3, LCC_traits>::type LCC;

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

int main(int argc, char** argv)
{
  typedef boost::property_map<LCC, CGAL::face_index_t>::const_type
                 Face_index_map;

  LCC lcc;
  CGAL::read_off((argc>1)?argv[1]:"cube.off", lcc);

  // Ad hoc property_map to store normals. Face_index_map is used to
  // map face_descriptors to a contiguous range of indices. See
  // http://www.boost.org/libs/property_map/doc/vector_property_map.html
  // for details.
  boost::vector_property_map<Vector, Face_index_map>
    normals(static_cast<unsigned>(num_faces(lcc)), get(CGAL::face_index, lcc));

  calculate_face_normals(
    lcc // Graph
    , get(CGAL::vertex_point, lcc) // map from vertex_descriptor to point
    , normals // map from face_descriptor to Vector_3
    );

  std::cout << "Normals" << std::endl;
  for(LCC::Attribute_range<2>::type::iterator it=lcc.attributes<2>().begin();
      it!=lcc.attributes<2>().end(); ++it)
  {
    // Facet_iterator is a face_descriptor, so we can use it as the
    // key here
    std::cout << normals[it] << std::endl;
  }

  return 0;
}
