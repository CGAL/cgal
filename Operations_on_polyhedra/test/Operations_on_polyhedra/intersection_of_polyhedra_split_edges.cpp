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


struct Is_on_polyline{
  bool operator()(Polyhedron::Halfedge_handle) const { return true;}
};

struct Set_vertex_corner{
  typedef CGAL::cpp0x::tuple<CGAL::internal_IOP::Intersection_type,
                             Polyhedron::Halfedge_handle,
                             CGAL::internal_IOP::Intersection_type,
                             Polyhedron::Halfedge_handle>                          Info;  
  
  void operator()(Polyhedron::Vertex_handle,int,Polyhedron*) {}
  void add_info_to_node(int i,Polyhedron*,const Info& info){
    std::cout << i<< " "<< CGAL::cpp0x::get<2>(info) << std::endl;
  }
};


typedef CGAL::Node_visitor_for_polyline_split<Polyhedron,Is_on_polyline,Set_vertex_corner> Split_visitor;
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

    std::cout << P.size_of_facets() << " " << Q.size_of_facets() << std::endl;
  
    std::list<Polyline> polylines;
    std::list<Polyline_info> polyline_infos;
    
    
    CGAL::Intersection_of_Polyhedra_3<Polyhedron,Kernel,Split_visitor> polyline_intersections;

    polyline_intersections( P,Q,
                            CGAL::dispatch_output<Polyline,Polyline_info>
                              (std::back_inserter(polylines),std::back_inserter(polyline_infos))
    );
    
    std::cout << polylines.size() << " polylines and " << polyline_infos.size() << " infos" << std::endl;
    
    CGAL_assertion(P.is_valid());
    CGAL_assertion(Q.is_valid());
    
    std::ofstream out1("out1.off");
    std::ofstream out2("out2.off");
    out1 << P;
    out2 << Q;
    
    return 0;
}
