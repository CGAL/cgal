#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <cstring>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_3                                       Point;

typedef CGAL::Surface_mesh<Point>                             SMesh;
typedef boost::graph_traits<SMesh>::vertex_descriptor         vertex_descriptor;
typedef boost::graph_traits<SMesh>::face_descriptor           face_descriptor;

int main()
{
  std::ifstream in(CGAL::data_file_path("meshes/colored_tetra.ply"));
  SMesh mesh;
  CGAL::IO::read_PLY(in, mesh);

  std::cerr << "Read mesh with " << mesh.number_of_vertices() << " vertices and "
            << mesh.number_of_faces() << " faces" << std::endl;

  std::cerr << "Properties associated with vertices:" << std::endl;
  std::vector<std::string> properties = mesh.properties<SMesh::Vertex_index>();
  for(std::size_t i = 0; i < properties.size(); ++ i)
    std::cerr << " * " << properties[i] << std::endl;

  std::cerr << "Properties associated with faces:" << std::endl;
  properties = mesh.properties<SMesh::Face_index>();
  for(std::size_t i = 0; i < properties.size(); ++ i)
    std::cerr << " * " << properties[i] << std::endl;

  mesh.add_property_map<SMesh::Edge_index, short>("id", 42);
  mesh.add_property_map<SMesh::Halfedge_index, float>("u", 13.f);
  mesh.add_property_map<SMesh::Halfedge_index, float>("v", 37.f);

  // Append second mesh
  std::ifstream in2("tetra.ply");
  CGAL::IO::read_PLY(in2, mesh);

  std::ofstream out("out.ply");
//  CGAL::IO::set_binary_mode(out);
  CGAL::IO::write_PLY(out, mesh);

  // extra test for read/write of properties
  mesh = SMesh();
  in.close();
  in.open(CGAL::data_file_path("meshes/colored_tetra.ply"));
  CGAL::IO::read_PLY(in, mesh);
  float fvalue=1001;
  auto v_uvmap = mesh.add_property_map<SMesh::Vertex_index, std::vector<float>>("v:uv").first;
  auto v_umap = mesh.add_property_map<SMesh::Vertex_index, float>("v:u").first;
  auto v_vmap = mesh.add_property_map<SMesh::Vertex_index, float>("v:v").first;
  for (SMesh::Vertex_index v : vertices(mesh))
  {
    v_uvmap[v]={fvalue, -fvalue};
    v_umap[v]=fvalue;
    v_vmap[v]=-fvalue;
    ++fvalue;
  }

  double dvalue=2001;
  auto f_uvmap = mesh.add_property_map<SMesh::Face_index, std::vector<double>>("f:uv").first;
  auto f_umap = mesh.add_property_map<SMesh::Face_index, double>("f:u").first;
  auto f_vmap = mesh.add_property_map<SMesh::Face_index, double>("f:v").first;
  for (SMesh::Face_index f : faces(mesh))
  {
    f_uvmap[f]={dvalue, -dvalue};
    f_umap[f]=dvalue;
    f_vmap[f]=-dvalue;
    ++dvalue;
  }

  out.close();
  out.open("out_ascii.ply");
  CGAL::IO::write_PLY(out, mesh);
  out.close();
  out.open("out_binary.ply");
  CGAL::IO::set_binary_mode(out);
  CGAL::IO::write_PLY(out, mesh);
  out.close();
  mesh.clear();

  const std::array<std::string,2> fnames = {"out_ascii.ply", "out_binary.ply"};
  for (std::string fn : fnames)
  {
    std::cout << "Reading " << fn << "\n";
    in.close();
    in.open(fn);
    SMesh mesh_bis;
    CGAL::IO::read_PLY(in, mesh_bis);

    assert((mesh_bis.property_map<SMesh::Vertex_index, float>("v:u").second));
    assert((mesh_bis.property_map<SMesh::Vertex_index, float>("v:v").second));
    assert((mesh_bis.property_map<SMesh::Vertex_index, std::vector<float>>("v:uv").second));

    v_uvmap = mesh_bis.property_map<SMesh::Vertex_index, std::vector<float>>("v:uv").first;
    v_umap = mesh_bis.property_map<SMesh::Vertex_index, float>("v:u").first;
    v_vmap = mesh_bis.property_map<SMesh::Vertex_index, float>("v:v").first;

    fvalue=1001;
    for (SMesh::Vertex_index v : vertices(mesh_bis))
    {
      assert(v_uvmap[v].size()==2);
      assert(v_uvmap[v][0]==fvalue);
      assert(v_uvmap[v][1]==-fvalue);
      assert(v_umap[v]==fvalue);
      assert(v_vmap[v]==-fvalue);
      ++fvalue;
    }

    assert((mesh_bis.property_map<SMesh::Face_index, double>("f:u").second));
    assert((mesh_bis.property_map<SMesh::Face_index, double>("f:v").second));
    assert((mesh_bis.property_map<SMesh::Face_index, std::vector<double>>("f:uv").second));

    f_uvmap = mesh_bis.property_map<SMesh::Face_index, std::vector<double>>("f:uv").first;
    f_umap = mesh_bis.property_map<SMesh::Face_index, double>("f:u").first;
    f_vmap = mesh_bis.property_map<SMesh::Face_index, double>("f:v").first;

    dvalue=2001;
    for (SMesh::Face_index f : faces(mesh_bis))
    {
      assert(f_uvmap[f].size()==2);
      assert(f_uvmap[f][0]==dvalue);
      assert(f_uvmap[f][1]==-dvalue);
      assert(f_umap[f]==dvalue);
      assert(f_vmap[f]==-dvalue);
      ++dvalue;
    }
  }

  return 0;
}
