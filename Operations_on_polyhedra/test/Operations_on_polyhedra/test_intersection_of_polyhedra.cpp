#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <vector>
#include <fstream>

#define CGAL_TODO_WARNINGS

#include <CGAL/intersection_of_Polyhedra_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef CGAL::Polyhedron_3<Kernel>                           Polyhedron;


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
    
    std::cout << "Size of polyhedra " << P.size_of_facets() << " " << Q.size_of_facets() << std::endl;
  
    std::list<std::vector<Kernel::Point_3> > polylines;
    
    CGAL::Intersection_of_Polyhedra_3<Polyhedron> polyline_intersections;

    polyline_intersections(P,Q,std::back_inserter(polylines));
    
    std::cout << "Nb polylines found : " << polylines.size() << std::endl;
    
    int k=0;
    for (std::list<std::vector<Kernel::Point_3> >::iterator it_poly=polylines.begin();it_poly!=polylines.end();++it_poly)
    {
      bool is_cycle=it_poly->begin()!=boost::prior(it_poly->end()) && *it_poly->begin()==*boost::prior(it_poly->end());
      std::cout << (is_cycle?std::string("cycle "):std::string("polyline ")) << ++k << " made of " << it_poly->size()-1 << " edges" << std::endl;
      
      //~ std::cout << "polyline " << ++k<< std::endl;
      //~ for (std::vector<Kernel::Point_3>::iterator it=it_poly->begin();it!=it_poly->end();++it)
        //~ std::cout << *it << std::endl;
    }   
    
    return 0;
}
