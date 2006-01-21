//file: examples/Interpolation/nn_coordinates_3.C 
// author(s)     : Raphaelle Chaine <Raphaelle.Chaine@sophia.inria.fr, raphaelle.chaine@liris.cnrs.fr>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h> 
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>

#include <CGAL/natural_neighbor_coordinates_3.h>

#include <fstream>

typedef double NT; //Number Type 

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};

typedef K::Point_3 Point3;
typedef K::Vector_3 Vector3;
typedef K::Sphere_3 Sphere_3;  

typedef CGAL::Triangulation_vertex_base_3<K> Vb3d;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb3d>  Vbh;

typedef CGAL::Triangulation_cell_base_3<K> Cb3d;

typedef CGAL::Triangulation_data_structure_3<Vbh,Cb3d> Tds3d;

typedef CGAL::Delaunay_triangulation_3<K, Tds3d> Delaunay3d;
typedef CGAL::Triangulation_hierarchy_3<Delaunay3d> Dh;

typedef Dh::Facet Facet;  
typedef Dh::Vertex_handle Vertex_handle;
typedef Dh::Cell_handle Cell_handle;
typedef Dh::Finite_vertices_iterator Finite_vertices_iterator;
typedef Dh::Vertex_iterator Vertex_iterator;
typedef Dh::Facet_circulator Facet_circulator;
typedef Dh::Cell_iterator Cell_iterator;

typedef K::Construct_circumcenter_3 Construct_circumcenter_3;

int main()
{  
  Dh triangulation;
  
  std::fstream iFile("data/points3", std::ios::in);
  Point3 p;
  
  while ( iFile >> p )
    triangulation.insert(p);
  
  Point3 pp[3];  
  std::cout << "Consider the natural coordinates of P1, P2 and P3 with regard to the triangulation of data/points3 " << std::endl;
  pp[0]=Point3(106.55,112.57,110.4); //inside data/points3 convex hull
  std::cout << "P1 is inside the convex hull" << std::endl;
  pp[1]=Point3(250,100,140); //on data/points3 convex hull
  std::cout << "P2 is on a vertex of the triangulation" << std::endl;
  pp[2]=Point3(0,0,0); //outside data/points3 convex hull
  std::cout << "P2 is outside the convex hull" << std::endl;

  for(int ii=0;ii<3;ii++)
    {  
      std::vector< std::pair< Vertex_handle,NT> > coor_laplace;
      std::vector< std::pair< Vertex_handle,NT> > coor_sibson;
      NT norm_coeff_laplace, norm_coeff_sibson;

      std::cout << "Point P"<< ii+1 << " : "<<pp[ii].x() << " " 
		<< pp[ii].y() << " " 
		<< pp[ii].z() << std::endl;

      laplace_natural_neighbor_coordinates_3(triangulation,pp[ii],
					     std::back_inserter(coor_laplace),
					     norm_coeff_laplace);
      
      std::cout << "Linear combination of natural neighbors with Laplace natural coordinates";
      std::cout << " + correctness (ensured only with an exact number type supporting sqrt)" << std::endl;  
      std::cout << is_correct_natural_neighborhood(triangulation,pp[ii],
						   coor_laplace.begin(),
						   coor_laplace.end(),
						   norm_coeff_laplace) 
		<< std::endl;

      sibson_natural_neighbor_coordinates_3(triangulation,pp[ii],
					    std::back_inserter(coor_sibson),
					    norm_coeff_sibson);
      std::cout << "Linear combination of natural neighbors with Sibson natural coordinates" << std::endl;  
      std::cout << " + correctness (ensured only with an exact number type)" << std::endl;  
      std::cout << is_correct_natural_neighborhood(triangulation,pp[ii],
						   coor_sibson.begin(),
						   coor_sibson.end(),
						   norm_coeff_sibson) 
		<< std::endl;      
    } 
  return 0;
}

