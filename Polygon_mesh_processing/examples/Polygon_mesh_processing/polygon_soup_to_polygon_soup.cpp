#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/boost/graph/helpers.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::FT                                                   FT;
typedef K::Point_3                                              Point_3;

typedef CGAL::Surface_mesh<Point_3>                             Mesh;

typedef std::array<FT, 3>                                       Point;
typedef std::array<std::size_t,3>                               Polygon;


struct Soup {
   std::vector<Point> points;
   std::vector<Polygon> polygons;
};



namespace  boost {
  template <>
  struct graph_traits<Soup> {
    typedef std::size_t vertex_descriptor;
    typedef std::size_t face_descriptor;
    typedef std::size_t edge_descriptor;
    typedef std::size_t halfedge_descriptor;

    typedef std::size_t vertices_size_type;
    typedef std::size_t faces_size_type;
    typedef std::size_t  edges_size_type;

    static vertex_descriptor null_vertex() { return static_cast<vertex_descriptor>(-1); }
    static face_descriptor null_face() { return static_cast<face_descriptor>(-1); }
    static edge_descriptor null_edge() { return static_cast<edge_descriptor>(-1); }
 };

} // namespace boost

namespace CGAL {
  namespace BGL {
template <>
bool is_valid_face_descriptor(typename boost::graph_traits<::Soup>::face_descriptor ,
                              const Soup & ,
                              const bool verb)
{
  return true;
}
  } // namespace BGL
} // namespace CGAL

template<typename Mesh>
class Soup_point_pmap
{
public:
  typedef boost::read_write_property_map_tag category;
  typedef std::array<FT, 3> value_type;
  typedef std::array<FT, 3> reference;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor key_type;

  Soup_point_pmap()
    : sm_(nullptr)
  {}

  Soup_point_pmap(const Mesh& sm)
    : sm_(const_cast<Mesh*>(&sm))
    {}

  Soup_point_pmap(const Soup_point_pmap& pm)
    : sm_(pm.sm_)
    {}

  reference operator[](key_type v) const
  {
    CGAL_assertion(sm_!=nullptr);
    return sm_->points[v];
  }

  inline friend reference get(const Soup_point_pmap<Mesh>& pm, key_type v)
  {
    CGAL_precondition(pm.sm_!=nullptr);
    return pm.sm_->points[v];
  }

  inline friend void put(const Soup_point_pmap<Mesh>& pm, key_type v, const value_type& p)
  {
    CGAL_precondition(pm.sm_!=nullptr);
    pm.sm_->points[v] = p;
  }

  private:
  mutable Mesh* sm_;
};


void reserve(Soup& soup, std::size_t num_points, std::size_t /* num_edges */ , std::size_t num_faces)
{
  soup.points.reserve(num_points);
  soup.polygons.reserve(num_faces);
}

std::size_t add_vertex(Soup& soup)
{
  std::array<FT, 3> p;
  soup.points.push_back(p);
  return soup.points.size() - 1; // return the index of the newly added vertex
}

namespace boost {

template <>
struct property_map<Soup, boost::vertex_point_t >
{
  typedef Soup_point_pmap<Soup> type;
  typedef Soup_point_pmap<Soup> const_type;
};

template <>
struct property_traits<Soup_point_pmap<Soup> >
{
  typedef std::array<FT, 3> value_type;
  typedef std::array<FT, 3> reference;
  typedef std::size_t key_type;

  typedef boost::read_write_property_map_tag category;

  static value_type get(const Soup_point_pmap<Soup>& pm, key_type v)
  {
    return pm[v];
  }

  static void put(Soup_point_pmap<Soup>& pm, key_type v, const value_type& p)
  {
    pm[v] = p;
  }
};

Soup_point_pmap<Soup> get(CGAL::vertex_point_t, const Soup& soup)
{
  return Soup_point_pmap<Soup>(soup);
}

} // namespace boost


namespace CGAL {
  namespace Euler {

    template <typename VertexRange>
    typename boost::graph_traits<Soup>::face_descriptor
    add_face(const VertexRange& vr, Soup& soup)
    {
      Polygon polygon = { vr[0], vr[1], vr[2] }; // assuming triangular faces
      soup.polygons.push_back(polygon);
      return soup.polygons.size() - 1; // return the index of the newly added
    }
  }
}


namespace PMP = CGAL::Polygon_mesh_processing;

int main(int, char**)
{
  std::vector<Point> points;
  std::vector<Polygon> polygons;

  points.push_back(CGAL::make_array<FT>( 0,  0, 0)); // 0
  points.push_back(CGAL::make_array<FT>( 1,  0, 0)); // 1
  points.push_back(CGAL::make_array<FT>( 0,  1, 0)); // 2
  points.push_back(CGAL::make_array<FT>(-1,  0, 0)); // 3


  polygons.push_back({0,1,2});
  polygons.push_back({1,0,3});

  PMP::orient_polygon_soup(points, polygons);

  Soup soup;
  PMP::polygon_soup_to_polygon_mesh(points, polygons, soup); // , CGAL::parameters::default_values(), CGAL::parameters::vertex_point_map(Soup_point_pmap<Soup>(soup)));

  std::cout << "Soup has " << soup.points.size() << " points and " << soup.polygons.size() << " polygons." << std::endl;
  for (std::size_t i = 0; i < soup.points.size(); ++i) {
    const auto& p = soup.points[i];
    std::cout << "Point " << i << ": (" << p[0] << ", " << p[1] << ", " << p[2] << ")" << std::endl;
  }
  for (std::size_t i = 0; i < soup.polygons.size(); ++i) {
    const auto& poly = soup.polygons[i];
    std::cout << "Polygon " << i << ": ";
    for (const auto& v : poly) {
      std::cout << v << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}
