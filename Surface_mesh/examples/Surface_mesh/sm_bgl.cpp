
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::vertex_iterator vertex_iterator;

typedef boost::property_map<Mesh, boost::vertex_index_t >::type          Vertex_index_pmap;
typedef boost::property_map<Mesh, boost::vertex_point_t >::type          Vertex_point_pmap; 

typedef boost::property_map<Mesh, boost::vertex_property_t<bool> >::type Vertex_bool_pmap;
// boost::property_traits<Vertex_bool_pmap>::key_type == vertex_descriptor
// boost::property_traits<Vertex_bool_pmap>::value_type == bool


int main()
{
  Mesh m;
  m.add_vertex(Point_3(1,2,3));

  Vertex_index_pmap vi_pm = get(boost::vertex_index_t(), m);

  Vertex_point_pmap vp_pm = get(boost::vertex_point_t(), m);

  Vertex_bool_pmap vs_pm  = add(boost::vertex_property_t<bool>("visited",false), m);

  vertex_iterator b,e;
  for(boost::tie(b,e) = vertices(m); b!=e; ++b){
    vertex_descriptor v = *b;
    put(vs_pm, v, true);
    std::cout << "vertex " << v << " with index " << vi_pm[v] << " at " << get(vp_pm,v) << std::endl;
  }

  
  Vertex_bool_pmap vs_pm2  = add(boost::vertex_property_t<bool>("visited",false), m);
  if(! vs_pm2){
    std::cerr << "Error: vertex property map named 'visited' already exists\n";
  }

  remove(vs_pm, m);
  
  return 0;
}
               
