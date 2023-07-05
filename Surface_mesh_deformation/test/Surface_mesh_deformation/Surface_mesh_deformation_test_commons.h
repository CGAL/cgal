#include <boost/graph/graph_traits.hpp>
#include <CGAL/property_map.h>
#include <optional>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>

template<class Polyhedron>
void read_to_polyhedron(const std::string file_name, Polyhedron& mesh)
{
  std::ifstream input(file_name);

  if ( !input || !(input >> mesh) || mesh.empty() ){
    std::cerr << "Error: can not read " << file_name << std::endl;
    assert(false);
  }
}

template<class Polyhedron>
void init_indices(Polyhedron& poly) {
  typedef typename boost::graph_traits<Polyhedron>::vertex_iterator      vertex_iterator;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_iterator    halfedge_iterator;

  vertex_iterator vb, ve;
  std::size_t counter = 0;
  for(boost::tie(vb, ve) = vertices(poly); vb != ve; ++vb, ++counter) {
    (*vb)->id() = counter;
  }

  counter = 0;
  halfedge_iterator heb, hee;
  for(boost::tie(heb, hee) = halfedges(poly); heb != hee; ++heb, ++counter) {
    (*heb)->id() = counter;
  }
}

template<class DeformMesh>
std::vector<typename DeformMesh::vertex_descriptor>
read_rois(DeformMesh& deform_mesh,
  const std::string& roi_file,
  const std::string& handle_file)
{
  std::ifstream roi_stream(roi_file.c_str());
  std::ifstream handle_stream(handle_file.c_str());
  if(!roi_stream || !handle_stream) {
    std::cerr << "Error: can not read roi or handle files" << std::endl;
    assert(false);
  }

  typedef typename boost::graph_traits<typename DeformMesh::Halfedge_graph>::vertex_iterator vertex_iterator;
  typedef typename DeformMesh::vertex_descriptor vertex_descriptor;
  // put all vertices to a vector
  typename DeformMesh::Halfedge_graph const& polyhedron = deform_mesh.halfedge_graph();

  std::vector<vertex_descriptor> vvertices;
  vvertices.reserve(num_vertices(polyhedron));
  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = vertices(polyhedron); vb != ve; ++vb) {
    vvertices.push_back(*vb);
  }
  // load handles and roi from txt

  // put all handles into one handle group
  std::size_t id;
  std::vector<typename DeformMesh::vertex_descriptor> hg;

  while(handle_stream >> id) {
    deform_mesh.insert_control_vertex(vvertices[id]);
    hg.push_back(vvertices[id]);
  }
  while(roi_stream >> id) {
    deform_mesh.insert_roi_vertex(vvertices[id]);
  }

  return hg;
}

template<class DeformMesh>
void preprocess_and_deform(DeformMesh& deform_mesh,
  const std::string& roi_file,
  const std::string& handle_file,
  CGAL::Simple_cartesian<double>::Vector_3 translate,
  int deformation_iteration)
{
  std::vector<typename DeformMesh::vertex_descriptor> hg =
    read_rois(deform_mesh, roi_file, handle_file);

  CGAL::Timer timer; timer.start();

  if(!deform_mesh.preprocess()) {
    std::cerr << "Error: preprocess() failed!" << std::endl;
    assert(false);
  }
  std::cerr << "Preprocess time: " << timer.time() << std::endl;
  timer.reset();

  deform_mesh.translate(hg.begin(), hg.end(), translate);
  deform_mesh.deform(deformation_iteration, 0);

  std::cerr << "Deformation time (one handle translated and deform method iterated " <<
    deformation_iteration << " times): " << timer.time() << std::endl;
}
