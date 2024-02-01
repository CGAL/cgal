#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/Bsurf/locally_shortest_path.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/boost/graph/generators.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>


namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;
typedef PMP::Face_location<Mesh, double>                      Face_location;
typedef PMP::Edge_location<Mesh, double>                      Edge_location;


//TODO: move to header and add doc
template <class TriangleMesh, class FT>
void simplify_sequence(std::vector<PMP::Edge_location<TriangleMesh, FT>>& edge_locations, const TriangleMesh& tm)
{
  std::cout << "  simplify_sequence:\n";
  std::cout << "  input.size() = " << edge_locations.size() << "\n";
  typedef boost::graph_traits<TriangleMesh> BGT;
  typedef typename BGT::vertex_descriptor vertex_descriptor;

  auto print_seq = [&]()
  {
    std::cout << "  ";
    typedef typename BGT::halfedge_descriptor halfedge_descriptor;
    for (auto el : edge_locations)
    {
      auto var = PMP::get_descriptor_from_location(el, tm);
      if (const vertex_descriptor* vd_ptr = std::get_if<vertex_descriptor>(&var))
        std::cout << "  " << *vd_ptr;
      else
      {
        if (const halfedge_descriptor* ed_ptr = std::get_if<halfedge_descriptor>(&var))
          std::cout << "  " << *ed_ptr;
        else
          std::cout << "  ERROR";
      }
    }
    std::cout << "\n";
  };
  print_seq();


  std::vector<PMP::Edge_location<TriangleMesh, FT>> filtered_locations;
  std::size_t end=edge_locations.size();
  filtered_locations.reserve(edge_locations.size());
  for(std::size_t curr=0; curr<end;)
  {
    filtered_locations.push_back(edge_locations[curr]);
    auto var = PMP::get_descriptor_from_location(edge_locations[curr], tm);
    if (const vertex_descriptor* vd_ptr = std::get_if<vertex_descriptor>(&var))
    {
      do{
        vertex_descriptor vd = *vd_ptr;
        ++curr;
        if (curr!=end)
        {
          var = PMP::get_descriptor_from_location(edge_locations[curr], tm);
          vd_ptr = std::get_if<vertex_descriptor>(&var);
          if (vd_ptr==nullptr || *vd_ptr!=vd)
            break;
        }
        else
          break;
      }
      while(true);
    }
    else
      ++curr;
  }

  if (filtered_locations.size()!=edge_locations.size())
    filtered_locations.swap(edge_locations);

  std::cout << "  output.size() = " << edge_locations.size() << "\n";
  print_seq();
}

void test_vertex_pair(unsigned int vid1, unsigned int vid2, const Mesh& mesh)
{
  Mesh::Vertex_index v1(vid1), v2(vid2);
  std::cout << "running " << v1 << " " << v2 << "\n";
  static int rid=-1;
  std::string fname="grid_shortest_path"+std::to_string(++rid)+".polylines.txt";
  std::cout << "  writing in " << fname << "\n";

  std::ofstream out(fname);
  out << std::setprecision(17);

  Mesh::Face_index f1=face(halfedge(v1, mesh), mesh), f2 = face(halfedge(v2, mesh), mesh);

  Face_location src(f1, CGAL::make_array(0.,0.,0.));
  Face_location tgt(f2, CGAL::make_array(0.,0.,0.));

  Mesh::Halfedge_index hf1 = prev(halfedge(f1, mesh), mesh);
  int i1=0;
  while (target(hf1, mesh)!= v1)
  {
    hf1=next(hf1,mesh);
    ++i1;
  }
  src.second[i1]=1;
  Mesh::Halfedge_index hf2 = prev(halfedge(f2, mesh), mesh);
  int i2=0;
  while (target(hf2, mesh)!= v2)
  {
    hf2=next(hf2,mesh);
    ++i2;
  }
  tgt.second[i2]=1;

  std::vector<Edge_location> edge_locations;
  auto src_bk=src, tgt_bk=tgt;
  PMP::locally_shortest_path<double>(src, tgt, mesh, edge_locations);
  assert(get<Mesh::Vertex_index>(PMP::get_descriptor_from_location(src,mesh))==v1);
  assert(get<Mesh::Vertex_index>(PMP::get_descriptor_from_location(tgt,mesh))==v2);

  simplify_sequence(edge_locations, mesh);

  out << edge_locations.size()+2;
  out << " " << PMP::construct_point(src, mesh);
  for (auto el : edge_locations)
    out << " " << PMP::construct_point(el, mesh);
  out << " " << PMP::construct_point(tgt, mesh) << "\n";
  out << std::flush;
  std::size_t expected_size = edge_locations.size();

  /// other direction

  src=src_bk;
  tgt=tgt_bk;
  edge_locations.clear();
  PMP::locally_shortest_path<double>(tgt, src, mesh, edge_locations);
  assert(get<Mesh::Vertex_index>(PMP::get_descriptor_from_location(src,mesh))==v1);
  assert(get<Mesh::Vertex_index>(PMP::get_descriptor_from_location(tgt,mesh))==v2);

  out << edge_locations.size()+2;
  out << " " << PMP::construct_point(tgt, mesh);
  for (auto el : edge_locations)
    out << " " << PMP::construct_point(el, mesh);
  out << " " << PMP::construct_point(src, mesh) << "\n";
  out << std::flush;
  CGAL_assertion(edge_locations.size() == expected_size);
}

int main()
{
  Mesh tmp;
  CGAL::make_grid(10, 10, tmp, [](int i, int j){return K::Point_3(i,j,0);}, true);
  Mesh mesh;
  CGAL::Polygon_mesh_processing::extrude_mesh(tmp, mesh, K::Vector_3(0,0,-1));
  std::ofstream("grid.off") << mesh;

  test_vertex_pair(84,36,mesh);
  test_vertex_pair(40,36,mesh);
  test_vertex_pair(78,28,mesh);

  return 0;
}
