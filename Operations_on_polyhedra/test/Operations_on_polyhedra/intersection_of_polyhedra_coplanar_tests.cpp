#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/intersection_of_Polyhedra_3.h>
#include <list>
#include <set>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Default_polyhedron_ppmap<Polyhedron> Ppmap;
typedef CGAL::internal_IOP::Intersection_point_with_info<Kernel,Polyhedron::Halfedge_handle, Ppmap> Cpl_inter_pt;



int main(int argc, char** argv)
{
  assert(argc==3);
  
  Polyhedron P,Q;

  std::ifstream input;
  input.open(argv[1]); assert(input);
  input >> P;
  input.close();
  input.open(argv[2]); assert(input);
  input >> Q;
  input.close();
  
  assert(P.size_of_vertices() == 3 && Q.size_of_vertices()==3);
  
  Ppmap ppmap;
  std::multiset<int> reference;
  
  //set up one result as the reference
  std::list<Cpl_inter_pt> inter_pts;
  CGAL::internal_IOP::intersection_coplanar_facets<Kernel>(P.halfedges_begin(),Q.halfedges_begin(), ppmap, inter_pts);
  std::cout << "====> reference: ";
  for (std::list<Cpl_inter_pt>::iterator iti=inter_pts.begin();iti!=inter_pts.end();++iti){
    reference.insert(iti->debug_unique_type_int());
    iti->print_debug();
  }
  std::cout << std::endl;
  
  std::multiset<int> res_to_test;
  for (Polyhedron::Halfedge_iterator hp=P.halfedges_begin();hp!=P.halfedges_end();++hp){
    for (Polyhedron::Halfedge_iterator hq=Q.halfedges_begin();hq!=Q.halfedges_end();++hq){
      //check P against Q
      inter_pts.clear();
      res_to_test.clear();
      CGAL::internal_IOP::intersection_coplanar_facets<Kernel>(hp,hq,ppmap,inter_pts);
      std::cout << "====> P vs Q: ";
      for (std::list<Cpl_inter_pt>::iterator iti=inter_pts.begin();iti!=inter_pts.end();++iti){
        res_to_test.insert(iti->debug_unique_type_int());
        iti->print_debug();
        assert(iti->is_valid());
      }
      std::cout << std::endl;
      assert(reference==res_to_test);
      //check Q against P
      inter_pts.clear();
      res_to_test.clear();
      CGAL::internal_IOP::intersection_coplanar_facets<Kernel>(hq,hp,ppmap,inter_pts);
      std::cout << "====> Q vs P: ";
      for (std::list<Cpl_inter_pt>::iterator iti=inter_pts.begin();iti!=inter_pts.end();++iti){
        res_to_test.insert(iti->debug_unique_type_int());
        iti->print_debug();
        assert(iti->is_valid());
      }
      std::cout << std::endl;
      assert(reference==res_to_test);
    }
  }
  
  
  return 0;
}
