#include <iostream>
#include <string>
#include <CGAL/array.h>
#include <CGAL/Bbox_3.h>
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

bool get_material_metadata(std::istream& input,
                           std::string& line,
                           material& _material)
{
  std::istringstream iss;
  iss.str(line);

  iss >> _material.second;//name

  while (std::getline(input, line))
  {
    std::string prop; //property
    iss.clear();
    iss.str(line);
    iss >> prop;

    if (prop.compare("Id") == 0)
      iss >> _material.first;
    else if (prop.compare("}") == 0)
      return true; //end of this material
    else
      CGAL_assertion(prop.compare("Color") == 0);
  }
  return false;
}

bool line_starts_with(const std::string& line, const char* cstr)
{
  std::string line_copy = line;
  std::size_t fnws = line_copy.find_first_not_of(" \t");
  if (fnws != std::string::npos)
  {
    line_copy.erase(0, fnws);
    return (line_copy.find(cstr) == 0);
  }
  return false;
}

/*!
 * read_surf reads a file which extension is .surf and fills `output` with one `Mesh` per patch.
 * Mesh is a model of FaceListGraph.
 */
template<class Mesh, class NamedParameters>
bool read_surf(std::istream& input, std::vector<Mesh>& output,
    std::vector<MaterialData>& metadata,
    CGAL::Bbox_3& grid_box,
    CGAL::cpp11::array<unsigned int, 3>& grid_size,
    const NamedParameters&)
{
  typedef typename CGAL::GetGeomTraits<Mesh,
      NamedParameters>::type Kernel;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  std::vector<Point_3> points;
  std::string line;
  std::istringstream iss;
  std::size_t nb_vertices(0);
  std::vector<material> materials;
  //ignore header
  while(std::getline(input, line))
  {
    std::size_t fnws=line.find_first_not_of(" \t");
    if(fnws != std::string::npos)
      line.erase(0, fnws);

    if (line_starts_with(line, "Materials"))
    {
      while(std::getline(input, line))
      {
        std::size_t fnws=line.find_first_not_of(" \t");
        if(fnws != std::string::npos)
          line.erase(0, fnws);
        if(line.compare(0, 1, "}") == 0)
          break;
        else
        {
          material _material;
          if (!get_material_metadata(input, line, _material))
            return false;
          materials.push_back(_material);
        }
      }
      break;
    }
  }

  //get grid box
  while (std::getline(input, line))
  {
    line.erase(0, line.find_first_not_of(" \t"));
    if (line.compare(0, 7, "GridBox") != 0)
      continue;
    else
    {
      iss.clear();
      line.erase(0, 7);
      iss.str(line);
      double xmin, xmax, ymin, ymax, zmin, zmax;
      iss >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;
      grid_box = CGAL::Bbox_3(xmin, ymin, zmin, xmax, ymax, zmax);
      break;
    }
  }

  //get grid size
  while (std::getline(input, line))
  {
    line.erase(0, line.find_first_not_of(" \t"));
    if (line.compare(0, 8, "GridSize") != 0)
      continue;
    else
    {
      iss.clear();
      line.erase(0, 8);
      iss.str(line);
      iss >> grid_size[0] >> grid_size[1] >> grid_size[2];
      break;
    }
  }

  //get number of vertices
  while(std::getline(input, line))
  {
    std::size_t fnws=line.find_first_not_of(" \t");
    if(fnws != std::string::npos)
      line.erase(0, fnws);
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
  }
  //get vertices
  while(std::getline(input, line))
  {
    std::size_t fnws=line.find_first_not_of(" \t");
    if(fnws != std::string::npos)
      line.erase(0, fnws);
    double x(0),y(0),z(0);
    if(line.compare(0, 16, "NBranchingPoints") == 0)
    {
      break;
    }
    iss.clear();
    iss.str(line);
    iss >> CGAL::iformat(x) >> CGAL::iformat(y) >> CGAL::iformat(z);
    points.push_back(Point_3(x,y,z));
  }
  std::cout<<nb_vertices<<" vertices"<<std::endl;
  int nb_patches = 0;
  //get number of patches
  while(std::getline(input, line))
  {
    std::size_t fnws=line.find_first_not_of(" \t");
    if(fnws != std::string::npos)
      line.erase(0, fnws);
    if(line.compare(0, 7, "Patches") != 0)
      continue;
    else
    {
      iss.clear();
      iss.str(line);
      std::string dump;
      iss >> dump >> nb_patches;
      break;
    }
  }
  std::cout<<nb_patches<<" patch(es)"<<std::endl;
  metadata.resize(nb_patches);
  output.resize(nb_patches);
  const vertex_descriptor null_vertex = boost::graph_traits<Mesh>::null_vertex();
  std::vector<vertex_descriptor> vertices(nb_vertices, null_vertex);
  for(int i=0; i < nb_patches; /* i is incremented in the body */)
  {
    Mesh& mesh = output[i];
    //get metada
    while(std::getline(input, line))
    {
      std::size_t fnws=line.find_first_not_of(" \t");
      if(fnws != std::string::npos)
        line.erase(0, fnws);
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
    }

    while(std::getline(input, line))
    {
      std::size_t fnws=line.find_first_not_of(" \t");
      if(fnws != std::string::npos)
        line.erase(0, fnws);
      if(line.compare(0, 9, "Triangles") != 0)
        continue;
      else
      {
        break;
      }
    }
    //std::getline(input, line);
    //connect triangles
    while(std::getline(input, line))
    {
      std::size_t fnws=line.find_first_not_of(" \t");
      if(fnws != std::string::npos)
        line.erase(0, fnws);
      std::size_t index[3];
      if(line.compare(0, 1, "}") == 0)
      {
        break;
      }
      iss.clear();
      iss.str(line);
      iss >> index[0] >> index[1] >> index[2];
      CGAL::cpp11::array<vertex_descriptor, 3> face;
      for(int id=0; id<3; ++id)
      {
        if(vertices[index[id]-1] == null_vertex)
        {
          vertices[index[id]-1] = add_vertex(points[index[id]-1], mesh);
        }
        face[id] = vertices[index[id]-1];
      }
      CGAL::Euler::add_face(face, mesh);
    }
    // reset the `vertices` vector, unless that was the last iteration of
    // the loop
    if(++i < nb_patches) std::fill(vertices.begin(),
                                   vertices.end(),
                                   null_vertex);
  } // end loop on patches

  return true;
}

template<class Mesh>
bool read_surf(std::istream& input, std::vector<Mesh>& output,
  std::vector<MaterialData>& metadata,
  CGAL::Bbox_3& grid_box,
  CGAL::cpp11::array<unsigned int, 3>& grid_size)
{
  return read_surf(input, output, metadata, grid_box, grid_size,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

