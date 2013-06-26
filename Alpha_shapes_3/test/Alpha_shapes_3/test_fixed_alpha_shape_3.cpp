#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#ifndef CGAL_NO_DEPRECATED_CODE
#  include <CGAL/Weighted_alpha_shape_euclidean_traits_3.h>
#endif
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <iostream>
#include <fstream>
#include <CGAL/Timer.h>

#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>

//#include "Filtered_weighted_alpha_shape_euclidean_traits_3.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel                                     Kernel;
typedef CGAL::Regular_triangulation_euclidean_traits_3<Kernel>                                  WFixed_Gt;
typedef WFixed_Gt                                                                               Gt;
//typedef CGAL::Filtered_weighted_alpha_shape_euclidean_traits_3<Kernel,true>                     Gt;


typedef CGAL::Weighted_point<Kernel::Point_3,Kernel::FT>                                        Weighted_point;
typedef CGAL::Fixed_alpha_shape_vertex_base_3<WFixed_Gt>                                        WFixed_Vb;
typedef CGAL::Fixed_alpha_shape_cell_base_3<WFixed_Gt>                                          WFixed_Cb;
typedef CGAL::Triangulation_data_structure_3<WFixed_Vb,WFixed_Cb>                               WFixed_TDS;
typedef CGAL::Regular_triangulation_3<WFixed_Gt,WFixed_TDS>                                     WFixed_DT;
typedef CGAL::Fixed_alpha_shape_3< WFixed_DT >                                                  WFixed_AS;


typedef CGAL::Alpha_shape_vertex_base_3<Gt>                                                     WVb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>                                                       WCb;
typedef CGAL::Triangulation_data_structure_3<WVb,WCb>                                           WTDS;
typedef CGAL::Regular_triangulation_3<Gt,WTDS>                                                  WDT;
typedef CGAL::Alpha_shape_3< WDT >                                                              WAS;

//Unweighted stuff
typedef CGAL::Fixed_alpha_shape_vertex_base_3<Kernel>                                           Fixed_Vb;
typedef CGAL::Fixed_alpha_shape_cell_base_3<Kernel>                                             Fixed_Cb;
typedef CGAL::Triangulation_data_structure_3<Fixed_Vb,Fixed_Cb>                                 Fixed_TDS;
typedef CGAL::Delaunay_triangulation_3<Kernel,Fixed_TDS>                                        Fixed_DT;
typedef CGAL::Fixed_alpha_shape_3<Fixed_DT>                                                     Fixed_AS;

typedef CGAL::Alpha_shape_vertex_base_3<Kernel>                                                 Vb;
typedef CGAL::Alpha_shape_cell_base_3<Kernel>                                                   Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>                                             TDS;
typedef CGAL::Delaunay_triangulation_3<Kernel,TDS>                                              DT;
typedef CGAL::Alpha_shape_3<DT>                                                                 AS;



template <class Object>
void fill_wp_lists(const char* file_path,std::list<Object>& Ls,double rw=0){
  double x,y,z,r;
  std::ifstream input(file_path);
 
  while(input){
    input >> x;
    if (!input) break;
    input >> y >> z >> r;
    Ls.push_back(Object(Kernel::Point_3(x,y,z),(r+rw)*(r+rw)));
  }
}


template <class Iterator,class Alpha_shape_3>
void print_simplices_classif(Iterator begin,Iterator end, const Alpha_shape_3& As){
  unsigned count[4]={0,0,0,0};
  for (Iterator it=begin;it!=end;++it){
    typename Alpha_shape_3::Classification_type type=As.classify(it);
    switch (type){
      case Alpha_shape_3::REGULAR:
        ++count[0];
      break;
      case Alpha_shape_3::INTERIOR:
        ++count[1];
      break;
      case Alpha_shape_3::SINGULAR:
        ++count[2];
      break;
      case Alpha_shape_3::EXTERIOR:
        ++count[3];
      break;
      default:
        std::cout << type << std::endl;
        assert(false);
    }
  }
  std::cout << "R I S E : "<<  count[0] << " "<< count[1] << " "<< count[2] << " "<< count[3] << std::endl;
}

template <class Iterator,class Alpha_shape_3>
void print_simplices_classif_fe(Iterator begin,Iterator end, const Alpha_shape_3& As){
  unsigned count[4]={0,0,0,0};
  for (Iterator it=begin;it!=end;++it){
    typename Alpha_shape_3::Classification_type type=As.classify(*it);
    switch (type){
      case Alpha_shape_3::REGULAR:
        ++count[0];
      break;
      case Alpha_shape_3::INTERIOR:
        ++count[1];
      break;
      case Alpha_shape_3::SINGULAR:
        ++count[2];
      break;
      case Alpha_shape_3::EXTERIOR:
        ++count[3];
      break;
      default:
        std::cout << type << std::endl;
        assert(false);
    }
  }
  std::cout << "R I S E : "<<  count[0] << " "<< count[1] << " "<< count[2] << " "<< count[3] << std::endl;
}

template <class AS1,class AS2>
void print_comparison(const AS1& as1, const AS2& as2){
  std::cout << "Cells\n";
  print_simplices_classif(as1.finite_cells_begin(),as1.finite_cells_end(),as1); 
  print_simplices_classif(as2.finite_cells_begin(),as2.finite_cells_end(),as2);
  std::cout << "Facets\n";
  print_simplices_classif_fe(as1.finite_facets_begin(),as1.finite_facets_end(),as1); 
  print_simplices_classif_fe(as2.finite_facets_begin(),as2.finite_facets_end(),as2);
  std::cout << "Edges\n";
  print_simplices_classif_fe(as1.finite_edges_begin(),as1.finite_edges_end(),as1); 
  print_simplices_classif_fe(as2.finite_edges_begin(),as2.finite_edges_end(),as2);
  std::cout << "Vertices\n";
  print_simplices_classif(as1.finite_vertices_begin(),as1.finite_vertices_end(),as1); 
  print_simplices_classif(as2.finite_vertices_begin(),as2.finite_vertices_end(),as2);
}

template <class Iterator1,class Iterator2,class AS1,class AS2>
void compare_classif(Iterator1 begin1, Iterator1 end1,const AS1& as1,Iterator2 it2, const AS2& as2,std::string sname,bool debug){
  unsigned nb=0;
  for (Iterator1 it1=begin1;it1!=end1;++it1){
    ++nb;
    if ( static_cast<int>(as1.classify(it1)) != static_cast<int>(as2.classify(it2))){
      std::cerr << nb << " Pb in " << sname << "\n";
      print_comparison(as1,as2);
      exit(EXIT_FAILURE);
    }
    ++it2;
  }
  if (debug)
    std::cout << sname << ": " << nb << " identical\n";
}

template <class Iterator1,class Iterator2,class AS1,class AS2>
void compare_facets_classif(Iterator1 begin1, Iterator1 end1,const AS1& as1,Iterator2 it2, const AS2& as2,std::string sname,bool debug){
  unsigned nb=0;
  for (Iterator1 it1=begin1;it1!=end1;++it1){
    for (int i=0; i<4;++i){
      typename AS1::Facet f1(it1,i);
      typename AS2::Facet f2(it2,i);
      if ( static_cast<int>(as1.classify(f1)) != static_cast<int>(as2.classify(f2))){
        std::cerr << " Pb in " << sname << "\n";
        print_comparison(as1,as2);
        exit(EXIT_FAILURE);
      }
    }
    ++it2;
  }
  if (debug)
    std::cout << sname << ": " << nb << " identical\n";
}

template <class Iterator1,class Iterator2,class AS1,class AS2>
void compare_edges_classif(Iterator1 begin1, Iterator1 end1,const AS1& as1,Iterator2 it2, const AS2& as2,std::string sname,bool debug){
  unsigned nb=0;
  for (Iterator1 it1=begin1;it1!=end1;++it1){
    for (int i=0; i<4;++i){
      for (int j=i+1;j<4;++j){
        typename AS1::Edge e1(it1,i,j);
        typename AS2::Edge e2(it2,i,j);
        if ( static_cast<int>(as1.classify(e1)) != static_cast<int>(as2.classify(e2))){
          std::cerr << " Pb in " << sname << "\n";
          print_comparison(as1,as2);
          exit(EXIT_FAILURE);
        }
      }
    }
    ++it2;
  }
  if (debug)
    std::cout << sname << ": " << nb << " identical\n";
}

template <class AS1,class AS2>
void compare_all(const AS1& as1, const AS2& as2,bool debug=false)
{
  compare_classif(as1.finite_cells_begin(),as1.finite_cells_end(),as1,as2.finite_cells_begin(),as2,"Cells",debug);  
  compare_facets_classif(as1.finite_cells_begin(),as1.finite_cells_end(),as1,as2.finite_cells_begin(),as2,"Facets",debug);
  compare_edges_classif(as1.finite_cells_begin(),as1.finite_cells_end(),as1,as2.finite_cells_begin(),as2,"Edges",debug);
  compare_classif(as1.finite_vertices_begin(),as1.finite_vertices_end(),as1,as2.finite_vertices_begin(),as2,"Vertices",debug); 
}
  
void test_dynamic_insert(const std::list<Weighted_point >& lst)
{
  typedef std::list<Weighted_point >::const_iterator Iterator;
  Iterator min_it=lst.begin();
  WFixed_AS dynamic_as;
  
  while(min_it!=lst.end() && dynamic_as.dimension() != 3)
  {
    dynamic_as.insert(*min_it);
    ++min_it;
  }
//  int k=0;  
  for (Iterator it=min_it ;it!=lst.end(); ++it)
  {
//    std::cout << ++k << " " << std::flush;    
    dynamic_as.insert(*it);
    WFixed_DT tr_copy;    
    tr_copy.set_infinite_vertex( tr_copy.tds().copy_tds( dynamic_as.tds(), dynamic_as.infinite_vertex() ) );
    WFixed_AS static_as (tr_copy);
    compare_all(dynamic_as,static_as);
  }
  std::cout << "done"<< std::endl;
}

void test_dynamic_remove(const std::list<Weighted_point >& lst)
{
   WFixed_AS dynamic_as(lst.begin(),lst.end());
  
//  int k=0;
  while (dynamic_as.dimension() == 3)
  {
//    std::cout << ++k << " " << std::flush;
    WFixed_DT tr_copy;
    tr_copy.set_infinite_vertex( tr_copy.tds().copy_tds( dynamic_as.tds(), dynamic_as.infinite_vertex() ) );
    WFixed_AS static_as( tr_copy );
    compare_all(dynamic_as,static_as);
    dynamic_as.remove(dynamic_as.finite_vertices_begin());
  }
  std::cout << "done"<< std::endl;
  
  while (dynamic_as.number_of_vertices() != 0) dynamic_as.remove(dynamic_as.finite_vertices_begin());
}

void make_one_run(const char* filename){
  std::cout << "== testing with "  << filename << " ==\n";
//read weighted points
  std::list<Weighted_point > lst;
  fill_wp_lists(filename,lst);

  //---Test weighted alpha shape
//build regular triangulation
  WFixed_DT T(lst.begin(),lst.end());
  if (lst.size()!=T.number_of_vertices())
    std::cout << lst.size()-T.number_of_vertices() << " hidden vertices.\n";
  
  std::cout << "Build Fixed weighted alpha complex" << std::endl;
  WFixed_AS wfixed_as(T);
 
  
//copy triangulation for familly alpha-shape
  WDT T1;
  T1.set_infinite_vertex( T1.tds().copy_tds( wfixed_as.tds(),wfixed_as.infinite_vertex() ) );
  std::cout << "Build familly weighted alpha complex" << std::endl;
  WAS w_as(T1,0,WAS::GENERAL);
  

//DEBUG info
//  print_comparison(wfixed_as,w_as);
  
//compare classification of simplices
  std::cout << "Compare both.... ";
  compare_all(wfixed_as,w_as); 
  std::cout << "OK\n";

//---Test alpha shape
  Fixed_DT delaunay0(lst.begin(),lst.end());
  DT delaunay1;
  delaunay1.set_infinite_vertex( delaunay1.tds().copy_tds( delaunay0.tds(),delaunay0.infinite_vertex() ) );

  std::cout << "Build Fixed alpha complex" << std::endl;
  Fixed_AS fixed_as(delaunay0);
  std::cout << "Build familly alpha complex" << std::endl;
  AS as(delaunay1,0,AS::GENERAL); 
  std::cout << "Compare both.... ";
  compare_all(fixed_as,as);
  std::cout << "OK\n";
  
//test dynamic version
  std::cout << "Test dynamic remove \n";
  test_dynamic_remove(lst);

  std::cout << "Test dynamic insert \n";
  test_dynamic_insert(lst);
}

int main(int argc, char** argv){
  if (argc==1){
    std::cerr << "Nothing was tested\n";
    return EXIT_SUCCESS;
  }  
  
  for (int i=1;i<argc;++i) make_one_run(argv[i]);
  return 0;
}
