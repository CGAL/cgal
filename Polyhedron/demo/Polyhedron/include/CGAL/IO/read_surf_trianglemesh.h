#include <iostream>
#include <string>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

//a material is composed of an material Id and a material name.
typedef std::pair<int, std::string> material;

struct MaterialData
{
  material innerRegion;
  material outerRegion;
};

/*!
 * read_surf reads a file which extension is .surf and fills `output` with one `Mesh` per patch.
 * Mesh is a model of FaceListGraph.
 */
template<class Mesh, class NamedParameters>
void read_surf(std::istream& input, std::vector<Mesh>& output, std::vector<MaterialData>& metadata,  const NamedParameters&)
{
  typedef typename GetGeomTraits<Mesh,
      NamedParameters>::type Kernel;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  std::vector<Point_3> points;
  std::string line;
  std::istringstream iss;
  int nb_vertices(0);
  std::vector<material> materials;
  //ignore header
  while(std::getline(input, line))
    if(line.compare(0, 13, "    Materials") != 0)
      continue;
    else
    {
      while(std::getline(input, line))
        if(line.compare(0, 5, "    }") == 0)
          break;
        else
        {
          //get materials metadata
          material _material;
          iss.clear();
          iss.str(line);
          std::string dump;
          iss >> _material.second;
          std::getline(input, line);
          std::getline(input, line);
          iss.clear();
          iss.str(line);
          iss >> dump >> _material.first;
          std::getline(input, line);

          materials.push_back(_material);
        }
      break;
    }

  //get number of vertices
  while(std::getline(input, line))
    if(line.compare(0, 8, "Vertices") != 0)
      continue;
    else
    {
      iss.clear();
      iss.str(line);
      std::string dump;
      iss >> dump >> nb_vertices;
      points.reserve(nb_vertices);
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
    points.push_back(Point_3(x,y,z));
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
  std::cout<<nb_patches<<" patch(es)"<<std::endl;
  metadata.resize(nb_patches);
  for(int i=0; i < nb_patches; ++i)
  {
    std::map<std::size_t, vertex_descriptor> vertices;
    Mesh mesh;
    //get metada
    while(std::getline(input, line))
      if(line.compare(0, 11, "InnerRegion") != 0)
        continue;
      else
      {
        std::string name;
        iss.clear();
        iss.str(line);
        std::string dump;
        iss >> dump >> name;
        for(std::size_t j=0; j<materials.size(); ++j)
        {
          if(materials[j].second.compare(name) == 0)
          {
            metadata[i].innerRegion=materials[j];
            break;
          }
        }
        std::getline(input, line);
        iss.clear();
        iss.str(line);
        iss >> dump >> name;
        for(std::size_t j=0; j<materials.size(); ++j)
        {
          if(materials[j].second.compare(name) == 0)
          {
            metadata[i].outerRegion=materials[j];
            break;
          }
        }
        break;
      }

    while(std::getline(input, line))
      if(line.compare(0, 9, "Triangles") != 0)
        continue;
      else
      {
        break;
      }
    //std::getline(input, line);
    //connect triangles
    while(std::getline(input, line))
    {
      std::size_t index[3];
      if(line.compare(0, 1, "}") == 0)
      {
        break;
      }
      iss.clear();
      iss.str(line);
      iss >> index[0] >> index[1] >> index[2];
      for(int id=0; id<3; ++id)
      {
        if(vertices.find(index[id]) == vertices.end())
        {
          vertices.insert(std::make_pair(index[id], add_vertex(points[index[id]-1], mesh)));
        }
      }
      std::vector<vertex_descriptor> face;
      for(int id=0; id<3; ++id)
        face.push_back(vertices[index[id]]);
      CGAL::Euler::add_face(face, mesh);
    }
    output.push_back(mesh);
  }
}

template<class Mesh>
void read_surf(std::istream& input, std::vector<Mesh>& output, std::vector<MaterialData>& metadata)
{
  read_surf(input, output, metadata, CGAL::Polygon_mesh_processing::parameters::all_default());
}

