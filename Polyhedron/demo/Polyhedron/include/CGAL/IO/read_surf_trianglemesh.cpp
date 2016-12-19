#include <iostream>
#include <string>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

/*!
 * read_surf_file reads a file which extension is .surf and fills `output` with one `Mesh` per patch.
 * Mesh is a model of FaceListGraph.
 */
template<class Mesh, class NamedParameters>
bool read_surf_file(std::istream& input, std::vector<Mesh>& output, const NamedParameters& np)
{
  typedef typename GetGeomTraits<Mesh,
      NamedParameters>::type Kernel;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  std::vector<Point_3> vertices;
  std::string line;
  std::istringstream iss;
  int nb_vertices(0);
  //ignore header
  while(std::getline(input, line))
    if(line.compare(0, 8, "Vertices") != 0)
      continue;
    else
    {
      iss.clear();
      iss.str(line);
      std::string dump;
      iss >> dump >> nb_vertices;
      vertices.reserve(nb_vertices);
      break;
    }
  //get vertices
  while(std::getline(input, line))
  {
    double x(0),y(0),z(0);
    if(line.compare(0, 16, "NBranchingPoints") == 0)
    {
      break;
    }
    iss.clear();
    iss.str(line);
    iss >> x >> y >> z;
    vertices.push_back(Point_3(x,y,z));
  }
  std::cout<<nb_vertices<<" vertices"<<std::endl;
  int nb_patches = 0;
  //get number of patches
  while(std::getline(input, line))
    if(line.compare(0, 7, "Patches") != 0)
      continue;
    else
    {
      //reserve size for points
      iss.clear();
      iss.str(line);
      std::string dump;
      iss >> dump >> nb_patches;
      break;
    }
  std::cout<<nb_patches<<" patches"<<std::endl;
  for(int i=0; i < nb_patches; ++i)
  {
    Mesh mesh;
    while(std::getline(input, line))
      if(line.compare(0, 9, "Triangles") != 0)
        continue;
      else
      {
        break;
      }
    std::getline(input, line);
    //connect triangles
    while(std::getline(input, line))
    {
      std::size_t a(0),b(0),c(0);
      if(line.compare(0, 1, "}") == 0)
      {
        break;
      }
      iss.clear();
      iss.str(line);
      iss >> a >> b >> c;
      std::vector<vertex_descriptor> face;
      face.push_back(add_vertex(vertices[a-1], mesh));
      face.push_back(add_vertex(vertices[b-1], mesh));
      face.push_back(add_vertex(vertices[c-1], mesh));
      CGAL::Euler::add_face(face, mesh);
    }
    output.push_back(mesh);
  }
}

