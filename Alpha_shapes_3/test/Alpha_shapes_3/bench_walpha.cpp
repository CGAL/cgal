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

#include "copy_tds.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel                                     Kernel;
typedef CGAL::Regular_triangulation_euclidean_traits_3<Kernel>                                  EPIC_traits;

typedef CGAL::Weighted_point<Kernel::Point_3,Kernel::FT>                                        Weighted_point;
typedef CGAL::Fixed_alpha_shape_vertex_base_3<EPIC_traits>                                      WFixed_Vb;
typedef CGAL::Fixed_alpha_shape_cell_base_3<EPIC_traits>                                        WFixed_Cb;
typedef CGAL::Triangulation_data_structure_3<WFixed_Vb,WFixed_Cb>                               WFixed_TDS;
typedef CGAL::Regular_triangulation_3<EPIC_traits,WFixed_TDS>                                   WFixed_DT;
typedef CGAL::Fixed_alpha_shape_3< WFixed_DT >                                                  WFixed_AS;


typedef CGAL::Alpha_shape_vertex_base_3<EPIC_traits>                                            WVb;
typedef CGAL::Alpha_shape_cell_base_3<EPIC_traits>                                              WCb;
typedef CGAL::Triangulation_data_structure_3<WVb,WCb>                                           WTDS;
typedef CGAL::Regular_triangulation_3<EPIC_traits,WTDS>                                         WDT;
typedef CGAL::Alpha_shape_3< WDT >                                                              WAS;

typedef CGAL::Alpha_shape_vertex_base_3<EPIC_traits,CGAL::Default,CGAL::Tag_true,CGAL::Tag_true> WVb_f;
typedef CGAL::Alpha_shape_cell_base_3<EPIC_traits,CGAL::Default,CGAL::Tag_true,CGAL::Tag_true>  WCb_f;
typedef CGAL::Triangulation_data_structure_3<WVb_f,WCb_f>                                       WTDS_f;
typedef CGAL::Regular_triangulation_3<EPIC_traits,WTDS_f>                                       WDT_f;
typedef CGAL::Alpha_shape_3< WDT_f,CGAL::Tag_true >                                                            WAS_f;




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



void make_one_run(const char* filename){
  std::cout << "== testing with "  << filename << " ==\n";
//read weighted points
  std::list<Weighted_point > lst;
  fill_wp_lists(filename,lst,1.4);

  //---Test weighted alpha shape
//build regular triangulation
  WFixed_DT T(lst.begin(),lst.end());
  if (lst.size()!=T.number_of_vertices())
    std::cout << lst.size()-T.number_of_vertices() << " hidden vertices.\n";
  
  std::cout << "Build Fixed weighted alpha complex" << std::endl;
  CGAL::Timer time;
  time.start();
  WFixed_AS wfixed_as(T);
  time.stop();
  std::cout << "Fixed "<< time.time() << std::endl;
 
  time.reset();
  
//copy triangulation for familly alpha-shape
  WDT T1;
  copy_tds(wfixed_as.tds(),T1.tds(),wfixed_as.infinite_vertex(),T1.infinite_vertex());
  std::cout << "Build familly weighted alpha complex" << std::endl;
  time.start();
  WAS w_as(T1,0,WAS::GENERAL);
  time.stop();
  std::cout << "Familly "<< time.time() << std::endl;
  
  //copy triangulation for familly alpha-shape
  CGAL::Interval_nt<false>::Protector protector;
  WDT_f T1f;
  copy_tds(wfixed_as.tds(),T1f.tds(),wfixed_as.infinite_vertex(),T1f.infinite_vertex());
  std::cout << "Build familly filtered weighted alpha complex" << std::endl;
  time.start();
  WAS_f w_asf(T1f,0,WAS_f::GENERAL);
  time.stop();
  std::cout << "Familly filtered "<< time.time() << std::endl;
}

int main(int argc, char** argv){
  if (argc==1){
    std::cerr << "Nothing was tested\n";
    return EXIT_SUCCESS;
  }  
  
  for (int i=1;i<argc;++i) make_one_run(argv[i]);
  return 0;
}
