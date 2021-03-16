#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <fstream>
#include <iostream>

namespace CGAL{
namespace Polygon_mesh_processing{
//todo face_descriptor ça doit etre le facerange::value_type je pense. checker le c++20 de value_type par contre
template<typename TriangleMesh, typename FaceRange, typename FacePairRange, typename NamedParameters1, typename NamedParameters2 >
void compare_meshes(const TriangleMesh& m1, const TriangleMesh& m2,
                    FacePairRange& common, FaceRange& tm1_only, FaceRange& tm2_only,
                    const NamedParameters1& np1,const NamedParameters2& np2)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  typedef typename GetVertexPointMap < TriangleMesh, NamedParameters1>::const_type VPMap1;
  typedef typename GetVertexPointMap < TriangleMesh, NamedParameters2>::const_type VPMap2;
  VPMap1 vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, m1));
  VPMap2 vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, m2));
  typedef typename boost::property_traits<VPMap2>::value_type Point_3;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;



  std::map<Point_3, std::size_t> point_id_map;
  std::vector<std::size_t> m1_vertex_id(num_vertices(m1), -1);
  std::vector<std::size_t> m2_vertex_id(num_vertices(m2), -1);

  //iterate both meshes to set ids to all points, and set vertex/point_id maps.
  std::size_t id =0;
  for(auto v : vertices(m1))
  {
    const Point_3& p = get(vpm1, v);
    auto res = point_id_map.insert(std::make_pair(p, id));
    if(res.second)
      id++;
    m1_vertex_id[(std::size_t)v]=res.first->second;
  }
  for(auto v : vertices(m2))
  {
    const Point_3& p = get(vpm2, v);
    auto res = point_id_map.insert(std::make_pair(p, id));
    if(res.second)
      id++;
    m2_vertex_id[(std::size_t)v]=res.first->second;
  }

  //fill a set with the "faces point-ids" of m1 and then iterate faces of m2 to compare.
  std::set<std::vector<std::size_t> > m1_faces;
  for(auto f : faces(m1))
  {
    std::vector<std::size_t> ids;
    for(auto v : CGAL::vertices_around_face(halfedge(f, m1), m1))
    {
      ids.push_back(m1_vertex_id[(std::size_t)v]);
    }
    std::sort(ids.begin(), ids.end());
    m1_faces.insert(ids);
  }
  std::map<std::vector<std::size_t>, face_descriptor> m2_faces_map;
  for(auto f : faces(m2))
  {
    std::vector<std::size_t> ids;
    for(auto v : CGAL::vertices_around_face(halfedge(f, m2), m2))
    {
      ids.push_back(m2_vertex_id[(std::size_t)v]);
    }
    std::sort(ids.begin(), ids.end());
    m2_faces_map.insert({ids, f});
    if(m1_faces.find(ids) == m1_faces.end())   {
      tm2_only.push_back(f);
    }
  }

  for(auto f : faces(m1))
  {
    std::vector<std::size_t> ids;
    for(auto v : CGAL::vertices_around_face(halfedge(f, m1), m1))
    {
      ids.push_back(m1_vertex_id[(std::size_t)v]);
    }
    std::sort(ids.begin(), ids.end());
    auto m2_face_it = m2_faces_map.find(ids);
    if(m2_face_it == m2_faces_map.end())
    {
      tm1_only.push_back(f);
    }
    else
    {
      common.push_back(std::make_pair(f, m2_face_it->second));
    }
  }
}
}//pmp
}//cgal
typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;

typedef K::Point_3                                                Point;

typedef CGAL::Surface_mesh<Point>                                 Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor      vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor        face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename1 = (argc > 1) ? argv[1] : "data/tet.off";
  const char* filename2 = (argc > 2) ? argv[2] : "data/tet-bis.off";

  Surface_mesh mesh1, mesh2;
  if(!PMP::read_polygon_mesh(filename1, mesh1))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }
  if(!PMP::read_polygon_mesh(filename2, mesh2))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }
  std::vector<std::pair<face_descriptor, face_descriptor> > common;
  std::vector<face_descriptor> m1_only, m2_only;
  PMP::compare_meshes(mesh1, mesh2, common, m1_only, m2_only, CGAL::parameters::all_default(), CGAL::parameters::all_default());
  std::cout<<"Faces only in m1 : "<<std::endl;
  for(const auto& f : m1_only)
  {
    std::cout<<f<<", ";
  }
  std::cout<<"\n Faces only in m2: "<<std::endl;
  for(const auto& f : m2_only)
  {
    std::cout<<f<<", ";
  }
  std::cout<<"\n Faces in both: "<<std::endl;
  for(const auto& f_pair : common)
  {
    std::cout<<f_pair.first<<", "<<f_pair.second<<";;"<<std::endl;
  }
  return 0;
}
