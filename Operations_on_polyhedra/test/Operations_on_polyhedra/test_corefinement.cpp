#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#define CGAL_TODO_WARNINGS
#define CGAL_COREFINEMENT_DEBUG

#include <CGAL/intersection_of_Polyhedra_3.h>
#include <CGAL/intersection_of_Polyhedra_3_refinement_visitor.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef CGAL::Polyhedron_3<Kernel>                           Polyhedron;

typedef CGAL::Node_visitor_refine_polyhedra<Polyhedron> Split_visitor;
typedef std::vector<Kernel::Point_3> Polyline;
typedef std::vector<std::pair<Polyhedron::Facet_handle,Polyhedron::Facet_handle> > Polyline_info;

int main(int argc,char** argv) {

    if (argc!=3){
      std::cerr << "Usage "<< argv[0] << " file1.off file2.off\n";
      return 1;
    }

    Polyhedron P, Q;    
    std::ifstream file(argv[1]);
    file >> P;
    file.close();
    file.open(argv[2]);
    file >> Q;
    file.close();
    
    CGAL::set_pretty_mode(std::cerr);
    
    if (!P.is_pure_triangle() || !Q.is_pure_triangle() ){
      std::cerr << "Input must be triangulated" << std::endl;
      exit(EXIT_FAILURE);
    }    
    
    CGAL_assertion(P.is_valid());
    CGAL_assertion(Q.is_valid());
    
    std::cout << "Polys " << &P << " " << &Q << std::endl;
    std::cout << "Vertices " <<  P.size_of_vertices() << " " << Q.size_of_vertices() << std::endl;
    std::cout << "Faces " <<  P.size_of_facets() << " " << Q.size_of_facets() << std::endl;
    std::cout << "Hedges " << P.size_of_halfedges() << " " << Q.size_of_halfedges() << std::endl;
  
    std::list<Polyline> polylines;
    std::list<Polyline_info> polyline_infos;
    
    
    CGAL::Intersection_of_Polyhedra_3<Polyhedron,Kernel,Split_visitor> polyline_intersections;

    polyline_intersections( P,Q,std::back_inserter(polylines) );
    
    std::cout << polylines.size() << " polylines: " << std::endl;
    for (std::list<Polyline>::iterator it=polylines.begin();it!=polylines.end();++it)
      std::cout << it->size() << " ";
    std::cout << std::endl;
    
    CGAL_assertion(P.is_valid());
    CGAL_assertion(P.is_pure_triangle());
    CGAL_assertion(Q.is_valid());
    CGAL_assertion(Q.is_pure_triangle());
    
    std::ofstream out1("out1.off");
    std::ofstream out2("out2.off");
    out1 << P;
    out2 << Q;
    
    return 0;
}
