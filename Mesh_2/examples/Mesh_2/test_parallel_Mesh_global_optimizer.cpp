#define CGAL_MESH_2_OPTIMIZER_VERBOSE
#define CGAL_MESH_2_OPTIMIZERS_DEBUG
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>

#include <CGAL/lloyd_optimize_mesh_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_selection.h>

#include <CGAL/draw_triangulation_2.h>
#include <limits.h>


#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>


#include <iostream>



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_mesh_vertex_base_2<K>                Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>                  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>  CDT;
typedef CGAL::Delaunay_triangulation_2<K>                   DT;

typedef CDT::Point Point;
typedef CGAL::Creator_uniform_2<double,Point>            Creator;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Vertex_iterator Vertex_iterator;

#define RADIUS 5.5
#define RADIUS_BIRTH 5
#define RADIUS_DEATH 5


void cell_death(CDT &cdt)
{
  CDT::Finite_vertices_iterator it;
  // CDT::Vertex_handle v; // equivalent à l'iterateur

  //  std::vector<CDT::Vertex_handle> vector_vertices;
  std::vector<CDT::Vertex_handle> vector_vertices;
  std::vector<CDT::Vertex_handle> dead_cells;
  std::vector<CDT::Vertex_handle>::iterator it_cells;
  int count=0;
  for (it = cdt.finite_vertices_begin();
  	   it!= cdt.finite_vertices_end();
  	   ++it)
  	{
	  // on ne garde que les points à l'intérieur des contraintes
	  // delete only points without constraints
	  if (!cdt.are_there_incident_constraints(it))
	    {
	      vector_vertices.push_back(it);
	      count++;
	    }
  	}

  int n =10;
 
  random_selection(vector_vertices.begin(), vector_vertices.end(), n,
		   std::back_inserter(dead_cells));


  for (it_cells = dead_cells.begin(); it_cells!= dead_cells.end(); ++it_cells)
    cdt.remove(*it_cells);

}

void cell_birth(CDT &cdt)
{
  // insert random points
  CGAL::Random_points_in_disc_2<Point,Creator> g(RADIUS_BIRTH);
  int new_cells = 10;
  for (int i=0; i<new_cells; ++i)
    {
      cdt.insert(*g++);
    }
  
}

int startSimu(int nb_points, int nb_iter)
{
  
  CDT cdt;
  // insert random points at first
  CGAL::Random_points_in_disc_2<Point,Creator> g(RADIUS);
  for (int i=0; i<nb_points; ++i)
    {
      cdt.insert(*g++);
    }


  // loop on birth death and lloyd optimization
  for (int n=0; n<nb_iter; ++n)
    {
      lloyd_optimize_mesh_2(cdt,
      			    CGAL::parameters::max_iteration_number = nb_iter,
      			    CGAL::parameters::mark = false);

    }

   
  return 0;
}

/// main function
int main(int argc, char**argv)
{
  // arguments
  int nb_points = 20; // number of points
  int nb_iter=5; // nb of iterations


  return startSimu(nb_points, nb_iter);
}
