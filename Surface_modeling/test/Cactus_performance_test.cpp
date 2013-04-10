#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Timer.h>

#define CGAL_EIGEN3_ENABLED
//#define CGAL_SUPERLU_ENABLED
#include <CGAL/Deform_mesh.h>

typedef CGAL::Simple_cartesian<double>   Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>  Polyhedron;

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

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_descriptor      edge_descriptor;

typedef Polyhedron_with_id_property_map<Polyhedron, vertex_descriptor> Vertex_index_map; // use id field of vertices
typedef Polyhedron_with_id_property_map<Polyhedron, edge_descriptor>   Edge_index_map;   // use id field of edges

typedef CGAL::Deform_mesh<Polyhedron, Vertex_index_map, Edge_index_map, CGAL::ORIGINAL_ARAP>  Deform_mesh_original;
typedef CGAL::Deform_mesh<Polyhedron, Vertex_index_map, Edge_index_map, CGAL::SPOKES_AND_RIMS> Deform_mesh_spokes;

typedef boost::graph_traits<Polyhedron>::vertex_iterator     vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_iterator       edge_iterator;

template<class DeformMesh>
typename DeformMesh::Handle_group read_rois(DeformMesh& deform_mesh)
{
   // load handles and roi from txt
  std::ifstream handle_stream("data/cactus_handle.txt"); // there is only one handle in cactus_handle.txt
  std::ifstream roi_stream("data/cactus_roi.txt");
  std::vector<int> handles;
  std::vector<int> rois;
  int id;
  while(handle_stream >> id) { handles.push_back(id); }
  while(roi_stream >> id) { rois.push_back(id); }

  typename DeformMesh::Handle_group active_handle_group = deform_mesh.create_handle_group();

  typename DeformMesh::Polyhedron const& polyhedron = deform_mesh.halfedge_graph();
  id = 0;
  vertex_iterator vb, ve;  
  for(boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb, ++id)
	{
    // not efficient but small poly
    if(std::find(handles.begin(), handles.end(), id) != handles.end()) { 
      deform_mesh.insert_handle(active_handle_group, *vb);
    }

    if(std::find(rois.begin(), rois.end(), id) != rois.end()) {
      deform_mesh.insert_roi(*vb);
    }
  }

  return active_handle_group;
}

int main()
{
  Polyhedron mesh;
  std::ifstream("data/cactus.off") >> mesh;

  Deform_mesh_original deform_mesh(mesh, Vertex_index_map(), Edge_index_map()); 
  // load handles and roi from txt
  Deform_mesh_original::Handle_group active_handle_group = read_rois(deform_mesh);
  deform_mesh.preprocess();

  deform_mesh.translate(active_handle_group, Deform_mesh_original::Vector(-0.55, -0.30, 0.0) );

  CGAL::Timer timer; timer.start();
  deform_mesh.deform(500, 0);
  timer.stop();

  std::cout << "---------------------------------------" << std::endl;
  std::cout << "Deformation performance: " << timer.time() << std::endl;
  std::cout << "---------------------------------------" << std::endl; 

  // std::ofstream("data/deformed_cactus.off") << mesh;

  std::cout << "Press ENTER to continue...";
  std::cin.get();
}

