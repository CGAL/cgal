#include <iostream>
#include <string>
#include <algorithm>

#include <CGAL/array.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

//a material is composed of an material Id and a material name.
typedef std::pair<int, std::string> material;

struct MaterialData
{
  material innerRegion;
  material outerRegion;
};

std::string to_lower_case(const std::string& str)
{
  std::string r = str;
  std::transform(r.begin(), r.end(), r.begin(), ::tolower);
  return r;
}

bool get_material_metadata(std::istream& input,
                           std::string& line,
                           material& _material,
                           int material_id)
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
    {
      int tmp_id;
      iss >> tmp_id;
      _material.first = material_id;

      if ((0 == to_lower_case(_material.second).compare("exterior"))
          && _material.first != 0)
      {
        std::cerr << "Exterior should have index 0. ";
        std::cerr << "In this file it has index " << _material.first << "." << std::endl;
        std::cerr << "Reader failed, because Meshing will fail to terminate." << std::endl;
        return false;
      }
    }
    else if (prop.compare("}") == 0)
      return true; //end of this material
  }
  return false;
}

bool line_starts_with(const std::string& line, const char* cstr)
{
  const std::size_t fnws = line.find_first_not_of(" \t");
  if (fnws != std::string::npos)
    return (line.compare(fnws, strlen(cstr), cstr) == 0);
  return false;
}
namespace IO{
namespace internal{
bool treat_surf_materials(std::istream& input,
                          std::vector<material>& materials,
                          int& material_id)
{
  std::string line;
  while(std::getline(input, line))
  {
    if(line_starts_with(line, "}"))
      break;
    else
    {
      material _material;
      if (!get_material_metadata(input, line, _material, material_id++))
        return false;
      materials.push_back(_material);
    }
  }
  return true;
}

void treat_surf_grid_box(const std::string& line,
                        CGAL::Bbox_3& grid_box)
{
  std::istringstream iss;
  iss.str(line);
  std::string dump;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  iss >> dump >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;
  grid_box = CGAL::Bbox_3(xmin, ymin, zmin, xmax, ymax, zmax);
}

void treat_surf_grid_size(const std::string& line,
                         std::array<unsigned int, 3>& grid_size)
{
  std::istringstream iss;
  iss.str(line);
  std::string dump;
  iss >> dump >> grid_size[0] >> grid_size[1] >> grid_size[2];
}

template<typename Point_3>
void treat_surf_vertices(std::istream& input,
                         const std::string& input_line,
                         std::size_t& nb_vertices,
                         std::vector<Point_3>& points)
{
  std::istringstream iss;
  iss.str(input_line);
  std::string dump;
  iss >> dump >> nb_vertices;
  points.reserve(nb_vertices);
  std::cout<<nb_vertices<<" vertices"<<std::endl;

  std::string line;
  //get vertices
  while(std::getline(input, line))
  {
    if (line_starts_with(line, "NBranchingPoints"))
      break;
    double x(0),y(0),z(0);
    iss.clear();
    iss.str(line);
    iss >> CGAL::iformat(x) >> CGAL::iformat(y) >> CGAL::iformat(z);
    points.push_back(Point_3(x,y,z));
  }
}

template<typename MaterialData, typename Mesh>
void get_surf_patches(const std::string& line,
                      int& nb_patches,
                      std::vector<MaterialData>& metadata,
                      std::vector<Mesh>& output)
{
  std::istringstream iss;
  iss.str(line);
  std::string dump;
  iss >> dump >> nb_patches;
  std::cout<<nb_patches<<" patch(es)"<<std::endl;
  metadata.resize(nb_patches);
  output.resize(nb_patches);
}

void treat_surf_inner_region(std::istream& input,
                             const std::string& input_line,
                             std::vector<material>& materials,
                             std::vector<MaterialData>& metadata,
                             const int& i)
{
  std::istringstream iss;
  std::string name;
  iss.clear();
  iss.str(input_line);
  std::string dump;
  iss >> dump >> name;
  name.erase(std::remove(name.begin(), name.end(), '"'), name.end());
  name.erase(std::remove(name.begin(), name.end(), ','), name.end());
  for(std::size_t j=0; j<materials.size(); ++j)
  {
    if(materials[j].second.compare(name) == 0)
    {
      metadata[i].innerRegion=materials[j];
      break;
    }
  }
  std::string line;
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
}

void treat_surf_triangles(const std::string& line,
                          std::size_t& nb_triangles)
{
  std::istringstream iss;
  iss.str(line);
  std::string dump;
  iss >> dump >> nb_triangles;
}

void connect_surf_triangles(std::istream& input,
                            const std::size_t& nb_triangles,
                            std::vector<std::array<std::size_t, 3> >& polygons)
{
  typedef std::array<std::size_t, 3> Triangle_ind;

  polygons.reserve(nb_triangles);
  std::string line;
  std::istringstream iss;
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

    Triangle_ind polygon;
    polygon[0] = index[0] - 1;
    polygon[1] = index[1] - 1;
    polygon[2] = index[2] - 1;
    polygons.push_back(polygon);
  }
  CGAL_assertion(nb_triangles == polygons.size());
}

template<typename Point_3, typename DuplicatedPointsOutIterator, typename Mesh>
void build_surf_patch(DuplicatedPointsOutIterator out,
                      std::vector<Mesh>& output,
                      std::vector<Point_3> points,
                      std::vector<std::array<std::size_t, 3> >& polygons,
                      const int& i)
{
  namespace PMP = CGAL::Polygon_mesh_processing;
  typedef std::array<std::size_t, 3> Triangle_ind;
  if (!PMP::is_polygon_soup_a_polygon_mesh(polygons))
  {
    std::cout << "Orientation of patch #" << (i + 1) << "...";
    std::cout.flush();

    const std::size_t nbp_init = points.size();
    bool no_duplicates =
        PMP::orient_polygon_soup(points, polygons);//returns false if some points
    //were duplicated

    std::cout << "\rOrientation of patch #" << (i + 1) << " done";

    if(!no_duplicates) //collect duplicates
    {
      for (std::size_t i = nbp_init; i < points.size(); ++i)
        *out++ = points[i];
      std::cout << " (non manifold -> "
                << (points.size() - nbp_init) << " duplicated vertices)";
    }
    std::cout << "." << std::endl;
  }

  Mesh& mesh = output[i];

  PMP::internal::PS_to_PM_converter<std::vector<Point_3>, std::vector<Triangle_ind> > converter(points, polygons);
  converter(mesh, false/*insert_isolated_vertices*/);

  CGAL_assertion(PMP::remove_isolated_vertices(mesh) == 0);
  CGAL_assertion(is_valid_polygon_mesh(mesh));
}

template<class Point>
bool fill_binary_vertices(std::istream& is, std::vector<Point>& points, const std::size_t& nb_points)
{
  float f[3];
  std::size_t pos = 0;
  while(pos < nb_points)
  {
    is.read(reinterpret_cast<char*>(&f[0]), sizeof(f[0]));
    if(!is.good())
      return false;
    is.read(reinterpret_cast<char*>(&f[1]), sizeof(f[1]));
    if(!is.good())
      return false;
    is.read(reinterpret_cast<char*>(&f[2]), sizeof(f[2]));
    if(!is.good())
      return false;
    points.push_back(Point(f[0], f[1], f[2]));
    ++pos;
  }
  return true;
}

bool fill_binary_triangles(std::istream& is, std::vector<std::array<std::size_t, 3> >& triangles, const std::size_t& nb_trs)
{
  int tr[7];
  std::size_t pos = 0;
  while(pos < nb_trs)
  {
    for(int i = 0; i < 7; ++i)
    {
      is.read(reinterpret_cast<char*>(&tr[i]), sizeof(tr[i]));
      if(!is.good())
        return false;
    }
    std::array<std::size_t, 3> tri;
    for(int i = 0; i < 3; ++i)
      tri[i] = tr[i];
    triangles.push_back(tri);
    ++pos;
  }
  return true;
}

bool fill_binary_patch(std::istream& is, std::vector<std::size_t>& patch, std::size_t size)
{
  int id = 0;
  std::size_t pos = 0;
  while(pos < size)
  {
    is.read(reinterpret_cast<char*>(&id), sizeof(int));
    if(!is.good())
      return false;
    patch.push_back(id);
    ++pos;
  }
  return true;
}

template<typename Point_3, typename DuplicatedPointsOutIterator, typename Mesh>
bool build_binary_surf_patch(DuplicatedPointsOutIterator& out,
                             std::vector<Mesh>& output,
                             std::vector<Point_3>& points,
                             const std::vector<std::array<std::size_t, 3> >& polygons,
                             const std::vector<std::size_t>& patch,
                             const std::size_t i)
{
  namespace PMP = CGAL::Polygon_mesh_processing;
  typedef std::array<std::size_t, 3> Triangle_ind;
  std::vector<Triangle_ind> triangles;
  triangles.reserve(patch.size());
  for(const std::size_t& id : patch)
  {
    triangles.push_back(polygons[id]);
  }
  if (!PMP::is_polygon_soup_a_polygon_mesh(triangles))
  {
    std::cout << "Orientation of patch #" << (i + 1) << "...";
    std::cout.flush();

    const std::size_t nbp_init = points.size();
    bool no_duplicates =
        PMP::orient_polygon_soup(points, triangles);//returns false if some points
    //were duplicated

    std::cout << "\rOrientation of patch #" << (i + 1) << " done";

    if(!no_duplicates) //collect duplicates
    {
      for (std::size_t i = nbp_init; i < points.size(); ++i)
        *out++ = points[i];
      std::cout << " (non manifold -> "
                << (points.size() - nbp_init) << " duplicated vertices)";
    }
    std::cout << "." << std::endl;
  }

  Mesh& mesh = output[i];

  PMP::internal::PS_to_PM_converter<std::vector<Point_3>, std::vector<Triangle_ind> > converter(points, triangles);
  converter(mesh, false/*insert_isolated_vertices*/);

  CGAL_assertion(PMP::remove_isolated_vertices(mesh) == 0);
  CGAL_assertion(is_valid_polygon_mesh(mesh));
  return true;
}
}//end internal
}//end IO

/*!
 * read_surf reads a file which extension is .surf and fills `output` with one `Mesh` per patch.
 * Mesh is a model of FaceListGraph.
 */
template<class Mesh, typename DuplicatedPointsOutIterator, class NamedParameters>
bool read_surf(std::istream& input, std::vector<Mesh>& output,
    std::vector<MaterialData>& metadata,
    CGAL::Bbox_3& grid_box,
    std::array<unsigned int, 3>& grid_size,
    DuplicatedPointsOutIterator out,
    const NamedParameters&)
{
  typedef typename CGAL::GetGeomTraits<Mesh, NamedParameters>::type Kernel;
  typedef typename Kernel::Point_3 Point_3;

  std::vector<Point_3> points;
  std::string line;
  std::size_t nb_vertices(0);
  std::vector<material> materials;
  //ignore header
  int material_id = 0;
  int nb_patches = 0;
  if(!std::getline(input, line))
    return false;
  bool binary = (line.find("BINARY") != std::string::npos);
  if(!binary)
  {
    while(std::getline(input, line))
    {
      if (line_starts_with(line, "Materials"))
      {
        if(!IO::internal::treat_surf_materials(input, materials, material_id))
          return false;
      }

      //get grid box
      if (line_starts_with(line, "GridBox"))
      {
        IO::internal::treat_surf_grid_box(line, grid_box);
      }

      //get grid size
      if (line_starts_with(line, "GridSize"))
      {
        IO::internal::treat_surf_grid_size(line, grid_size);
      }

      //get number of vertices
      if (line_starts_with(line, "Vertices"))
      {
        IO::internal::treat_surf_vertices(input, line, nb_vertices, points);
      }

      //get number of patches
      if (line_starts_with(line, "Patches"))
      {
        IO::internal::get_surf_patches(line, nb_patches, metadata, output);
        break;
      }
    }

    for(int i=0; i < nb_patches; ++i)
    {
      std::size_t nb_triangles(0);
      //get metada
      while(std::getline(input, line))
      {
        if (line_starts_with(line, "InnerRegion"))
        {
          IO::internal::treat_surf_inner_region(input, line, materials, metadata, i);
        }
        if (line_starts_with(line, "Triangles"))
        {
          IO::internal::treat_surf_triangles(line, nb_triangles);
          break;
        }
      }
      //connect triangles
      typedef std::array<std::size_t, 3> Triangle_ind;
      std::vector<Triangle_ind> polygons;
      IO::internal::connect_surf_triangles(input, nb_triangles, polygons);

      //build patch
      IO::internal::build_surf_patch(out, output, points, polygons, i);
    } // end loop on patches
  }
  else
  {
    int nTriangles = 0;
    typedef std::array<std::size_t, 3> Triangle_ind;
    std::vector<Triangle_ind> triangles;
    std::vector<int> tr_per_patches;
    int patch_counter = -1;
    while(std::getline(input, line))
    {
      //get nb vertices
      if (line_starts_with(line, "nVertices"))
      {
        std::istringstream iss;
        iss.str(line);
        std::string dump;
        if(!(iss >> dump >> nb_vertices))
           return false;
      }
    //get nb triangles
      if (line_starts_with(line, "nTriangles"))
      {
        std::istringstream iss;
        iss.str(line);
        std::string dump;
        if(!(iss >> dump >> nTriangles))
          return false;
        std::cout<<nTriangles<<" triangles."<<std::endl;
      }
    //get nb patches and nb triangles per patches
      if (line_starts_with(line, "define Patches"))
      {
        std::istringstream iss;
        iss.str(line);
        std::string dump;
        if(patch_counter <0)
        {
          if(!(iss >> dump >> dump >> nb_patches))
            return false;
          std::cout<<nb_patches<<" patches."<<std::endl;
          ++patch_counter;
          metadata.resize(nb_patches);
          output.resize(nb_patches);
        }
        else if(patch_counter < nb_patches)
        {
          int nb_tr = 0;
          if(!(iss >> dump >> dump >> nb_tr))
            return false;
          tr_per_patches.push_back(nb_tr);
          ++patch_counter;
        }
        else
        {
          std::cerr<<"Error in input file. Incoherent number of materials."<<std::endl;
          return false;
        }
        //reset patch_counter
        if(patch_counter == nb_patches)
          patch_counter = 0;
      }

    //get materials Ids
      if (line_starts_with(line, "Materials"))
      {
        if(!IO::internal::treat_surf_materials(input, materials, material_id))
          return false;
      }
    //get patches  InnerRegion
      if (line_starts_with(line, "InnerRegion"))
      {
        IO::internal::treat_surf_inner_region(input, line, materials,metadata, patch_counter);
      }
      if (line_starts_with(line, "GridBox"))
      {
        IO::internal::treat_surf_grid_box(line, grid_box);
      }

      //get grid size
      if (line_starts_with(line, "GridSize"))
      {
        IO::internal::treat_surf_grid_size(line, grid_size);
      }
      //fill vertices
      if (line_starts_with(line, "@1"))
        IO::internal::fill_binary_vertices(input, points, nb_vertices);
      //fill triangles
      else if (line_starts_with(line, "@2"))
      {
        IO::internal::fill_binary_triangles(input, triangles, nTriangles);
      }
      //fill patches
      else if (line_starts_with(line, "@"))
      {
        std::vector<std::size_t> patch;
        IO::internal::fill_binary_patch(input, patch, tr_per_patches[patch_counter]);
        IO::internal::build_binary_surf_patch(out, output, points, triangles, patch, patch_counter++);
      }
    }
  }
  return true;
}

template<class Mesh, typename DuplicatedPointsOutIterator>
bool read_surf(std::istream& input, std::vector<Mesh>& output,
  std::vector<MaterialData>& metadata,
  CGAL::Bbox_3& grid_box,
  std::array<unsigned int, 3>& grid_size,
  DuplicatedPointsOutIterator out)
{
  return read_surf(input, output, metadata, grid_box, grid_size, out,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

