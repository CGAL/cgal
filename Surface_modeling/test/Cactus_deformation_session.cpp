#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Deform_mesh.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include "Property_maps_for_edit_plugin.h"

#include <CGAL/Eigen_solver_traits.h>

#define CGAL_SUPERLU_ENABLED

#ifdef CGAL_SUPERLU_ENABLED
	#include <Eigen/SuperLUSupport>
	typedef CGAL::Eigen_solver_traits<Eigen::SuperLU<CGAL::Eigen_sparse_matrix<double>::EigenType> > DefaultSolver;
#else
	#include <Eigen/SparseLU>
	///////////////////////////////////////////////// 
	namespace CGAL {
	namespace internal {
		template <class FT, class EigenMatrix, class EigenOrdering>
		struct Get_eigen_matrix< ::Eigen::SparseLU<EigenMatrix, EigenOrdering >, FT> {
			typedef Eigen_sparse_matrix<FT, ::Eigen::ColMajor> type;
		};
	} // internal
	} // CGAL
	///////////////////////////////////////////////// 
	typedef CGAL::Eigen_solver_traits<
			Eigen::SparseLU<
			  CGAL::Eigen_sparse_matrix<double, Eigen::ColMajor>::EigenType,
			  Eigen::COLAMDOrdering<int> >  > DefaultSolver;
#endif
  
typedef CGAL::Simple_cartesian<double>   Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>  Polyhedron;

typedef Polyhedron_vertex_deformation_index_map<Polyhedron> Vertex_index_map;
typedef Polyhedron_edge_deformation_index_map<Polyhedron>   Edge_index_map;

typedef CGAL::Deform_mesh<Polyhedron, DefaultSolver, Vertex_index_map, Edge_index_map> Deform_mesh;

template <class T>
std::string toString(const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

Deform_mesh::Handle_group read_rois(Deform_mesh& deform_mesh)
{
   // load handles and roi from txt
  std::ifstream handle_stream("data/cactus_handle.txt"); // there is only one handle in cactus_handle.txt
  std::ifstream roi_stream("data/cactus_roi.txt");
  std::vector<int> handles;
  std::vector<int> rois;
  int id;
  while(handle_stream >> id) { handles.push_back(id); }
  while(roi_stream >> id) { rois.push_back(id); }

  Deform_mesh::Handle_group active_handle_group;

  id = 0;
  for(Polyhedron::Vertex_iterator it = deform_mesh.polyhedron.vertices_begin(); it != deform_mesh.polyhedron.vertices_end();
	  ++it, ++id)
	{
    // not efficient but small poly
    if(std::find(handles.begin(), handles.end(), id) != handles.end()) { 
      active_handle_group = deform_mesh.insert_handle(it);
    }

    if(std::find(rois.begin(), rois.end(), id) != rois.end()) {
      deform_mesh.insert_roi(it);
    }
  }

  return active_handle_group;
}

void compare_mesh(Polyhedron& mesh_1, Polyhedron& mesh_2)
{
  Polyhedron::Vertex_iterator it_1 = mesh_1.vertices_begin();
  Polyhedron::Vertex_iterator it_2 = mesh_2.vertices_begin();
  double total_dif_x, total_dif_y, total_dif_z;
  total_dif_x = total_dif_y = total_dif_z = 0;
  for( ; it_1 != mesh_1.vertices_end(); ++it_1 , ++it_2)
  {
    Kernel::Vector_3 dif = (it_1->point() - it_2->point());
    total_dif_x += std::abs(dif.x());
    total_dif_y += std::abs(dif.y());
    total_dif_z += std::abs(dif.z());
  }

  total_dif_x /= mesh_1.size_of_vertices();
  total_dif_y /= mesh_1.size_of_vertices();
  total_dif_z /= mesh_1.size_of_vertices();

  std::cout << total_dif_x << " " << total_dif_y << " " << total_dif_z << std::endl;
}

// read deformation session saved as a handle differences
void read_handle_difs_and_deform(Deform_mesh& deform_mesh, Deform_mesh::Handle_group& active_handle_group)
{
  std::ifstream dif_stream("data/cactus_handle_differences.txt");
  std::vector<Kernel::Vector_3> dif_vector;
  double x, y, z;
  while(dif_stream >> x >> y >> z)
  {
    dif_vector.push_back(Kernel::Vector_3(x, y, z));
  }

  for(int i = 0; i < dif_vector.size(); ++i)
  {
    deform_mesh.translate(active_handle_group, dif_vector[i]);
    deform_mesh.deform();

    // read pre-deformed cactus
    std::string predeformed_cactus_file = "data/cactus_deformed/cactus_deformed_" + toString(i) + ".off";
    Polyhedron predeformed_cactus;
	
    std::ifstream(predeformed_cactus_file) >> predeformed_cactus;
	compare_mesh(predeformed_cactus, deform_mesh.polyhedron);
	// for saving deformation
    //std::ofstream(predeformed_cactus_file) << deform_mesh.polyhedron;
    //std::cout << predeformed_cactus_file << std::endl;
  }
}

int main()
{
  Polyhedron mesh;
	std::ifstream("data/cactus.off") >> mesh;

  Deform_mesh deform_mesh(mesh, Vertex_index_map(), Edge_index_map()); 
  // load handles and roi from txt
  Deform_mesh::Handle_group active_handle_group = read_rois(deform_mesh);

  deform_mesh.preprocess();

  read_handle_difs_and_deform(deform_mesh, active_handle_group);
}




