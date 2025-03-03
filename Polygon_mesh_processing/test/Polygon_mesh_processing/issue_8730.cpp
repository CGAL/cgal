#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;

using vertex_index = Surface_mesh::vertex_index;
using halfedge_index = Surface_mesh::halfedge_index;
using face_index = Surface_mesh::face_index;
using edge_index = Surface_mesh::edge_index;


std::vector<Point_3> points = {
    {0,0,0},
    {0,0.5,0},
    {0,1,0},
    {0.5,0,0},
    {0.5,0.5,0},
    {0.5,1,0},
    {1,0,0},
    {1,0.5,0},
    {1,1,0},
    {0,0,0.5},
    {0,0.5,0.5},
    {0,1,0.5},
    {0.5,0,0.5},
    {0.5,1,0.5},
    {1,0,0.5},
    {1,0.5,0.5},
    {1,1,0.5},
    {0,0,1},
    {0,0.5,1},
    {0,1,1},
    {0.5,0,1},
    {0.5,0.5,1},
    {0.5,1,1},
    {1,0,1},
    {1,0.5,1},
    {1,1,1}
};

std::vector<std::array<int, 3>> faces = {
    {7,4,5},
    {7,6,4},
    {7,5,8},
    {4,2,5},
    {2,4,1},
    {11,5,2},
    {11,13,5},
    {11,2,10},
    {2,1,10},
    {11,10,18},
    {11,18,19},
    {10,17,18},
    {10,9,17},
    {20,18,17},
    {20,21,18},
    {12,20,17},
    {12,14,20},
    {12,17,9},
    {3,12,9},
    {3,6,12},
    {3,9,0},
    {3,0,1},
    {9,1,0},
    {1,9,10},
    {3,1,4},
    {6,3,4},
    {12,6,14},
    {6,7,14},
    {15,14,7},
    {14,15,23},
    {15,7,8},
    {15,8,16},
    {24,15,16},
    {16,8,13},
    {22,16,13},
    {13,8,5},
    {22,25,16},
    {22,13,19},
    {13,11,19},
    {22,19,21},
    {21,24,22},
    {19,18,21},
    {21,23,24},
    {22,24,25},
    {24,16,25},
    {23,21,20},
    {24,23,15},
    {14,23,20}
};

int main()
{
  Surface_mesh mesh;
  std::vector<vertex_index> vertices(points.size());
  for (size_t i = 0; i < points.size(); ++i)
  {
    vertices[i] = mesh.add_vertex(points[i]);
  }

  for (size_t i = 0; i < faces.size(); ++i)
  {
    vertex_index v0 = vertices[faces[i][0]];
    vertex_index v1 = vertices[faces[i][1]];
    vertex_index v2 = vertices[faces[i][2]];
    mesh.add_face(v0, v1, v2);
  }

  auto vertexToPatches = mesh.add_property_map<vertex_index, std::set<size_t>>("v:patches").first;
  auto faceToPatch = mesh.add_property_map<face_index, size_t>("f:patch", std::numeric_limits<size_t>::max()).first;
  auto edgeToIsFeature = mesh.add_property_map<edge_index, bool>("e:is_feature", false).first;

  double angle_in_deg = 90;

  CGAL::Polygon_mesh_processing::detect_sharp_edges(mesh, angle_in_deg, edgeToIsFeature);

  auto number_of_patches = CGAL::Polygon_mesh_processing::
    connected_components(mesh, faceToPatch,
      CGAL::parameters::edge_is_constrained_map(edgeToIsFeature));

  assert(number_of_patches==6);

  CGAL::Polygon_mesh_processing::detect_vertex_incident_patches(mesh,
    faceToPatch, vertexToPatches, edgeToIsFeature);

  std::array<int,4> degrees = CGAL::make_array(0,0,0,0);

  for (auto v : mesh.vertices())
  {
    auto d = vertexToPatches[v].size();
    assert(d<4);
    degrees[d]+=1;
  }

  assert(degrees[0]==6);
  assert(degrees[1]==0);
  assert(degrees[2]==12);
  assert(degrees[3]==8);

  return 0;
}