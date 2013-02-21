#include "StdAfx.h"

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

int main()
{
  Polyhedron mesh;
	std::ifstream("data/square.off") >> mesh;

  Deform_mesh deform_mesh(mesh, Vertex_index_map(), Edge_index_map()); 
  // load handles and roi from txt
  std::ifstream handle_stream("data/Symmetry_test_handle.txt");
  std::ifstream non_roi_stream("data/Symmetry_test_non_roi.txt");
  std::vector<int> handles;
  std::vector<int> non_rois;
  int id;
  while(handle_stream >> id) { handles.push_back(id); }
  while(non_roi_stream >> id) { non_rois.push_back(id); }
  std::cout << handles.size() << std::endl;
  std::cout << non_rois.size();
  id = 0;
  for(Polyhedron::Vertex_iterator it = mesh.vertices_begin(); it != mesh.vertices_end();
	  ++it, ++id)
	{
    // not efficient but small poly
    if(std::find(handles.begin(), handles.end(), id) != handles.end())
    {
      deform_mesh.insert_handle(it);
    }

    if(std::find(non_rois.begin(), non_rois.end(), id) == non_rois.end())
    {
      deform_mesh.insert_roi(it);
    }
  }

  deform_mesh.preprocess();

  Kernel::Vector_3 dif(-0.45, -0.65, 0);
  deform_mesh(dif);
  for(int i = 0; i < 50; ++i)
  {
    deform_mesh.deform();
  }
  std::ofstream("data/square_deformed.off") << mesh;
}

