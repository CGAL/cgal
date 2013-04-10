#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/properties_Polyhedron_3.h>

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

const double squared_threshold = 0.0001; // alert if average difs between precomputed and deformed mesh models is above threshold

#define PERFORMANCE            // define for checking performance
// #define MESH_DIFFERENCE        // define for checking difs with precomputed-saved deformations


template <class T>
std::string toString(const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

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

bool compare_mesh(const Polyhedron& mesh_1, const Polyhedron& mesh_2)
{
  Polyhedron::Vertex_const_iterator it_1 = mesh_1.vertices_begin();
  Polyhedron::Vertex_const_iterator it_2 = mesh_2.vertices_begin();
  Kernel::Vector_3 total_dif(0,0,0);
  for( ; it_1 != mesh_1.vertices_end(); ++it_1 , ++it_2)
  {
    total_dif = total_dif + (it_1->point() - it_2->point());    
  }
  total_dif = total_dif / mesh_1.size_of_vertices();

  std::cout << total_dif << std::endl;

  bool fail = total_dif.squared_length() > squared_threshold;
  return fail;
}

// read deformation session saved as a handle differences
template<class DeformMesh>
bool read_handle_difs_and_deform(DeformMesh& deform_mesh, typename DeformMesh::Handle_group& active_handle_group)
{
  std::ifstream dif_stream("data/cactus_handle_differences.txt");
  std::vector<Kernel::Vector_3> dif_vector;
  double x, y, z;
  while(dif_stream >> x >> y >> z)
  {
    dif_vector.push_back(Kernel::Vector_3(x, y, z));
  }
  CGAL::Timer timer;
  
  for(std::size_t i = 0; i < dif_vector.size(); ++i)
  {
    timer.start();
    deform_mesh.translate(active_handle_group, dif_vector[i]);
    deform_mesh.deform();
    timer.stop();
#ifdef MESH_DIFFERENCE
    // read pre-deformed cactus
    std::string predeformed_cactus_file = "data/cactus_deformed/cactus_deformed_" + toString(i) + ".off";
    Polyhedron predeformed_cactus;
	
    std::ifstream(predeformed_cactus_file) >> predeformed_cactus;
	  bool fail = compare_mesh(predeformed_cactus, deform_mesh.halfedge_graph());
    if( fail ) { return true; }
#endif
	  // for saving deformation
    //std::ofstream(predeformed_cactus_file) << deform_mesh.halfedge_graph();
    //std::cout << predeformed_cactus_file << std::endl;
  }

#ifdef PERFORMANCE
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "Deformation performance: " << timer.time() << std::endl;
  std::cout << "---------------------------------------" << std::endl; 
  std::cout << "Press ENTER to continue...";
  std::cin.get();
#endif

  return false;
}

int main()
{
  Polyhedron mesh;
  std::ifstream("data/cactus.off") >> mesh;

  Deform_mesh_original deform_mesh_1(mesh, Vertex_index_map(), Edge_index_map()); 
  // load handles and roi from txt
  Deform_mesh_original::Handle_group active_handle_group_1 = read_rois(deform_mesh_1);
  deform_mesh_1.preprocess();

  bool fail_1 = read_handle_difs_and_deform(deform_mesh_1, active_handle_group_1);
  if(fail_1) { return EXIT_FAILURE; }
  std::cout << "ORIGINAL ARAP Success!" << std::endl;

  Deform_mesh_spokes deform_mesh_2(mesh, Vertex_index_map(), Edge_index_map()); 
  // load handles and roi from txt
  Deform_mesh_spokes::Handle_group active_handle_group_2 = read_rois(deform_mesh_2);
  deform_mesh_2.preprocess();

  bool fail_2 = read_handle_difs_and_deform(deform_mesh_2, active_handle_group_2);
  if(fail_2) { return EXIT_FAILURE; }
  std::cout << "SPOKES AND RIMS ARAP Success!" << std::endl;
}

