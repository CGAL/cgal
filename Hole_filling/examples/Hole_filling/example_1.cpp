#define CGAL_SUPERLU_ENABLED
#undef NDEBUG
#define DEBUG_TRACE
#include <CGAL/Hole_filling.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <boost/tuple/tuple.hpp>

typedef CGAL::Simple_cartesian<double>   Kernel;
typedef CGAL::Polyhedron_3<Kernel>       Polyhedron;
typedef Polyhedron::Facet_handle         Facet_handle;
typedef Polyhedron::Traits::Point_3      Point_3;

//
//void test_triangulate_polyline(Polyhedron& poly) {
//  typedef std::vector< std::vector<boost::tuple<int, int, int> > > Triangles_list;
//  Triangles_list tris;
//  // construct polyline from border
//  std::vector<Point_3> polyline;
//  for(Polyhedron::Halfedge_iterator it = poly.halfedges_begin(); it != poly.halfedges_end(); ++it){
//    if(!it->is_border()) { continue; }
//    Polyhedron::Halfedge_around_facet_circulator hf_around_facet = it->facet_begin();
//    do {
//      polyline.push_back(hf_around_facet->vertex()->point());
//    } while(++hf_around_facet != it->facet_begin());
//  }
//  std::cout << "tri begin" << std::endl;
//  for(int i = 0; i < 100; ++i) {
//    tris.push_back(std::vector<boost::tuple<int, int, int> >());
//    //CGAL::triangulate_hole_polyline(polyline.begin(), polyline.end(),std::back_inserter(tris.back()));
//  }
//  std::cout << "tri end" << std::endl;
//  bool is_all_equal = true;
//  for(std::size_t i = 0; i + 1 < tris.size() && is_all_equal; ++i)
//  {
//    for(std::size_t j = 0; j < tris[0].size(); ++j) {
//      if(tris[i][j].get<0>() != tris[i+1][j].get<0>()
//        || tris[i][j].get<1>() != tris[i+1][j].get<1>()
//        || tris[i][j].get<2>() != tris[i+1][j].get<2>()) 
//      {
//        is_all_equal = false;
//        std::cout << "---------not equal points---------" << std::endl;
//        std::cout << tris[i][j].get<0>() << tris[i][j].get<1>() << tris[i][j].get<2>() << std::endl;
//        std::cout << tris[i][j+1].get<0>() << tris[i][j+1].get<1>() << tris[i][j+1].get<2>() << std::endl;
//        std::cout << "----------------------------------" << std::endl;
//      }
//    }
//  }
//
//  if(!is_all_equal) {
//    std::cout << "Error: test_triangulate results are different!" << std::endl;
//  }
//  else {
//    std::cout << "OK: test_triangulate results are the same" << std::endl;
//  }
//}
//
//void test_triangulate(Polyhedron& poly) {
//  typedef std::vector< std::vector<Point_3> > Triangles_list; // 3 point for each tri
//  Triangles_list tris; //hold results
//
//  for(int i = 0; i < 100; ++i) {
//    Polyhedron tmp = poly;
//    for(Polyhedron::Halfedge_iterator it = tmp.halfedges_begin(); it != tmp.halfedges_end(); ++it){
//      if(it->is_border()) {
//        std::vector<Facet_handle> facets;
//        CGAL::triangulate_hole(tmp, it, std::back_inserter(facets));
//        tris.push_back(std::vector<Point_3>());
//        for(std::vector<Facet_handle>::iterator it = facets.begin(); it != facets.end(); ++it) {
//          Polyhedron::Halfedge_around_facet_circulator b = (*it)->facet_begin(), e(b);
//          do{
//            tris.back().push_back(b->vertex()->point());
//          }while(++b != e);
//        }
//        break;
//      }
//    }
//  }
//
//  bool is_all_equal = true;
//  for(std::size_t i = 0; i + 1 < tris.size() && is_all_equal; ++i)
//  {
//    for(std::size_t j = 0; j < tris[0].size(); ++j) {
//      if(tris[i][j] != tris[i+1][j]) {
//        is_all_equal = false;
//        std::cout << "---------not equal points---------" << std::endl;
//        std::cout << tris[i][j] << std::endl;
//        std::cout << tris[i+1][j] << std::endl;
//        std::cout << "----------------------------------" << std::endl;
//      }
//    }
//  }
//
//  if(!is_all_equal) {
//    std::cout << "Error: test_triangulate results are different!" << std::endl;
//  }
//  else {
//    std::cout << "OK: test_triangulate results are the same" << std::endl;
//  }
//}

//void test_triangulate_fair(Polyhedron& poly) {
//  typedef std::vector< std::vector<Point_3> > Triangles_list;
//  Triangles_list tris; //hold results
//
//  for(int i = 0; i < 100; ++i) {
//    Polyhedron tmp = poly;
//    for(Polyhedron::Halfedge_iterator it = tmp.halfedges_begin(); it != tmp.halfedges_end(); ++it){
//      if(it->is_border()) {
//        std::vector<Facet_handle> facets;
//        CGAL::triangulate_and_refine_hole(tmp, it, std::back_inserter(facets));
//        tris.push_back(std::vector<Point_3>());
//        for(std::vector<Facet_handle>::iterator it = facets.begin(); it != facets.end(); ++it) {
//          Polyhedron::Halfedge_around_facet_circulator b = (*it)->facet_begin(), e(b);
//          do{
//            tris.back().push_back(b->vertex()->point());
//          }while(++b != e);
//        }
//        break;
//      }
//    }
//  }
//
//  bool is_all_equal = true;
//  for(std::size_t i = 0; i + 1 < tris.size() && is_all_equal; ++i)
//  {
//    if(tris[i] != tris[i+1]) { is_all_equal = false; }
//  }
//
//  if(!is_all_equal) {
//    std::cout << "Error: test_triangulate_fair results are different!" << std::endl;
//  }
//  else {
//    std::cout << "OK: test_triangulate_fair results are the same" << std::endl;
//  }
//}

#include <boost/function_output_iterator.hpp>
struct Nop_functor {
  template<class T>
  void operator()(const T& /* t */) const {}
};
typedef boost::function_output_iterator<Nop_functor> Nop_out;


int main() {
  Polyhedron poly;
  std::ifstream input("data/mech-holes-shark-bug.off");
  if ( !input || !(input >> poly) || poly.empty() ) {
    std::cerr<< "Cannot open file" << std::endl;
    return 1;
  }

  for(Polyhedron::Halfedge_iterator it = poly.halfedges_begin(); it != poly.halfedges_end(); ++it){
    if(it->is_border()) {
      //CGAL::triangulate_and_refine_hole(poly, it, Nop_out(), Nop_out(), 14.21);
      CGAL::triangulate_refine_and_fair_hole(poly, it, Nop_out(), Nop_out(), 10.0);
      break;
    }
  }

  //test_triangulate(poly);
  //test_triangulate_polyline(poly);
  std::ofstream out("data/out.off");
  out << poly;
  out.close();
  std::cout << "Done!" << std::endl;
  std::cin.get();
}

