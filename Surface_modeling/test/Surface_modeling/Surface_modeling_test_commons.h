#include <boost/property_map/property_map.hpp>
#include <boost/optional.hpp>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

template<class PolyhedronWithId, class KeyType>
struct Polyhedron_with_id_property_map
    : public boost::put_get_helper<std::size_t&,
             Polyhedron_with_id_property_map<PolyhedronWithId, KeyType> >
{
public:
    typedef KeyType      key_type;
    typedef std::size_t  value_type;
    typedef value_type&  reference;
    typedef boost::lvalue_property_map_tag category;
        
    reference operator[](key_type key) const { return key->id(); }
};

template<class Polyhedron>
bool read_to_polyhedron(const char* file_name, Polyhedron& mesh)
{
  std::ifstream input(file_name);

  if ( !input || !(input >> mesh) || mesh.empty() ){
    std::cerr << "Error: can not read " << file_name << std::endl;
    return false;
  }
  return true;
}

template<class DeformMesh>
boost::optional<typename DeformMesh::Handle_group>
read_rois(DeformMesh& deform_mesh, 
  const std::string& roi_file,
  const std::string& handle_file)
{
  std::ifstream roi_stream(roi_file);
  std::ifstream handle_stream(handle_file);
  if(!roi_stream || !handle_stream) {
    std::cerr << "Error: can not read roi or handle files" << std::endl;
    return boost::optional<typename DeformMesh::Handle_group>();
  }

  typedef typename boost::graph_traits<typename DeformMesh::Polyhedron>::vertex_iterator vertex_iterator;
  typedef typename DeformMesh::vertex_descriptor vertex_descriptor;
  // put all vertices to a vector
  typename DeformMesh::Polyhedron const& polyhedron = deform_mesh.halfedge_graph();

  std::vector<vertex_descriptor> vertices;
  vertices.reserve(boost::num_edges(polyhedron));
  vertex_iterator vb, ve;
  for(boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb) {
    vertices.push_back(*vb);
  }
  // load handles and roi from txt

  // put all handles into one handle group
  std::size_t id;
  typename DeformMesh::Handle_group active_handle_group = deform_mesh.create_handle_group();

  while(handle_stream >> id) {
    deform_mesh.insert_handle(active_handle_group, vertices[id]); 
  }
  while(roi_stream >> id) {
    deform_mesh.insert_roi(vertices[id]);
  }

  return boost::make_optional(active_handle_group);
}

template<class DeformMesh>
bool preprocess_and_deform(DeformMesh& deform_mesh, 
  const std::string& roi_file,
  const std::string& handle_file,
  typename DeformMesh::Vector translate, 
  int deformation_iteration) 
{
  boost::optional<typename DeformMesh::Handle_group> active_handle_group = 
    read_rois(deform_mesh, roi_file, handle_file);
  if(!active_handle_group) { return false; }

  CGAL::Timer timer; timer.start();

  if(!deform_mesh.preprocess()) {
    std::cerr << "Error: preprocess() failed!" << std::endl;
    return false;
  }
  std::cerr << "Preprocess time: " << timer.time() << std::endl;
  timer.reset();

  deform_mesh.translate(*active_handle_group, translate);
  deform_mesh.deform(deformation_iteration, 0);

  std::cerr << "Deformation time (one handle translated and deform method iterated " <<
    deformation_iteration << " times): " << timer.time() << std::endl;
  return true;
}
