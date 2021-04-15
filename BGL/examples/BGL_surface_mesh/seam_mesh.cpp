
// Enable the operator<< for a seam vertex/halfedge/edge/face
#define CGAL_SEAM_MESH_INSERT_OPERATOR

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <iostream>
#include <fstream>
#include <vector>


typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor         SM_vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor       SM_halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor           SM_edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor           SM_face_descriptor;

typedef Mesh::Property_map<SM_edge_descriptor, bool>            Seam_edge_pmap;
typedef Mesh::Property_map<SM_vertex_descriptor, bool>          Seam_vertex_pmap;
typedef CGAL::Seam_mesh<Mesh, Seam_edge_pmap, Seam_vertex_pmap> Seam_mesh;

typedef boost::graph_traits<Seam_mesh>::vertex_descriptor         vertex_descriptor;
typedef boost::graph_traits<Seam_mesh>::halfedge_descriptor       halfedge_descriptor;
typedef boost::graph_traits<Seam_mesh>::edge_descriptor           edge_descriptor;
typedef boost::graph_traits<Seam_mesh>::face_descriptor           face_descriptor;


int main(int argc, char* argv[])
{
  Mesh sm;
  std::ifstream in((argc>1) ? argv[1] : "data/cube.off");
  in >> sm;

  Seam_edge_pmap seam_edge_pm =
      sm.add_property_map<SM_edge_descriptor, bool>("e:on_seam", false).first;
  Seam_vertex_pmap seam_vertex_pm =
      sm.add_property_map<SM_vertex_descriptor, bool>("v:on_seam",false).first;

  Seam_mesh mesh(sm, seam_edge_pm, seam_vertex_pm);

  // Add seams. Below are two different seam selection files (choose one).
  // 1st: Split the cube in two connected components
  // 2nd: Seams are all edges incident to a vertex
  SM_halfedge_descriptor smhd = mesh.add_seams("data/two_connected_components.selection.txt");
  // SM_halfedge_descriptor smhd = mesh.add_seams("data/flatten.selection.txt");

  assert(smhd != SM_halfedge_descriptor());

  std::cout << "Added: " << mesh.number_of_seam_edges() << " seam edges" << std::endl;

  halfedge_descriptor bhd(smhd);
  std::cout << "First seam halfedge: " << smhd << std::endl;
  std::cout << "source = " << source(smhd, mesh)
            << " (pointing to vertex: " << source(smhd, sm) << ")" << std::endl;
  std::cout << "target = " << target(smhd, mesh)
            << " (pointing to vertex: " << target(smhd, sm) << ")" << std::endl;

  std::cout << "prev of seam halfedge in (base) mesh: " << prev(smhd, sm) << std::endl;
  std::cout << "prev of seam halfedge in seam mesh: " << prev(bhd, mesh) << std::endl;

  std::cout << "next of seam halfedge in (base) mesh: " << next(smhd, sm) << std::endl;
  std::cout << "next of seam halfedge in seam mesh: " << next(bhd, mesh) << std::endl;

  std::cout << "opposite of seam halfedge in (base) mesh: " << opposite(smhd, sm) << std::endl;
  std::cout << "opposite of seam halfedge in seam mesh: " << opposite(bhd, mesh) << std::endl;

  std::cout << "vertices on one of the seams" << std::endl;
  for(halfedge_descriptor hd :
                halfedges_around_face(opposite(bhd, mesh), mesh)){
    std::cout << target(hd.tmhd, sm) << " ";
  }
  std::cout << std::endl;

  std::cout << "vertices around " << target(smhd , sm) << " in (base) mesh" << std::endl;
  for(SM_halfedge_descriptor hd : halfedges_around_target(smhd, sm)){
    std::cout << source(hd, sm) << " ";
  }
  std::cout << std::endl;

  std::cout << "vertices around " << target(bhd , mesh) << " in seam mesh" << std::endl;
  for(halfedge_descriptor hd : halfedges_around_target(bhd, mesh)){
    std::cout << source(hd.tmhd, sm) << " ";
  }
  std::cout << std::endl;

  std::cout << "vertices around " << source(smhd , sm) << " in (base) mesh" << std::endl;
  for(SM_halfedge_descriptor hd :
                halfedges_around_source(source(smhd, sm), sm)){
     std::cout << target(hd, sm) << " ";
  }
  std::cout << std::endl;

  std::cout << "vertices around " << source(bhd , mesh) << " in seam mesh" << std::endl;
  for(halfedge_descriptor hd :
                halfedges_around_source(source(bhd, mesh), mesh)){
    std::cout << target(hd.tmhd, sm) << " ";
  }
  std::cout << std::endl;

  std::cout << "vertices around vertices in seam mesh" << std::endl;
  for(vertex_descriptor vd : vertices(mesh)){
    halfedge_descriptor hd = halfedge(vd, mesh);
    std::cout << " " << vd << " has incident vertices:" << std::endl;
    for(halfedge_descriptor hd2 : halfedges_around_target(hd, mesh)){
      std::cout << "  " << hd2;
    }
    std::cout << std::endl;
  }
  std::cout << "done" << std::endl;

  std::cout << "the (base) mesh has: " << num_halfedges(sm) << " halfedges" << std::endl;
  std::cout << "the seam mesh has: " << num_halfedges(mesh) << " halfedges" << std::endl;
  std::cout << "halfedges in (base) mesh" << std::endl;
  for(SM_halfedge_descriptor hd : halfedges(sm)){
     std::cout << hd << " ";
  }
  std::cout << std::endl;

  std::cout << "halfedges in seam mesh" << std::endl;
  for(halfedge_descriptor hd : halfedges(mesh)){
     std::cout << hd << " ";
  }
  std::cout << std::endl;

  std::cout << "faces of the base and seam meshes" << std::endl;
  for(face_descriptor fd : faces(mesh)){
    std::cout << fd << " ";
  }
  std::cout << std::endl;

  std::vector<face_descriptor> faces;
  CGAL::Polygon_mesh_processing::connected_component(face(bhd, mesh),
                                                     mesh,
                                                     std::back_inserter(faces));
  std::cout << "the connected component (in the seam mesh) given by halfedge: "  << smhd;
  std::cout << " has " << faces.size() << " faces." << std::endl;

  std::cout << "accessing coordinates of the source of halfedge " << smhd << ": ";
  boost::property_map<Seam_mesh, CGAL::vertex_point_t>::type vpm =
                                                    get(CGAL::vertex_point, mesh);
  std::cout << get(vpm, source(smhd, mesh)) << std::endl;

  return 0;
}

