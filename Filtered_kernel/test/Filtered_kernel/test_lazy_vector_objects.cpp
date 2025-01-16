#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersection_2.h>
#include <CGAL/intersection_3.h>
#include <limits>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Segment_2 Segment_2;
typedef K::Triangle_2 Triangle_2;
typedef K::Triangle_3 Triangle_3;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;


int main()
{
  {
    Point_2 p1(0,0), p2(2,0), p3(1,3), q1(0,2), q2(2,2), q3(1,-1);
    Triangle_2 t1(p1,p2,p3);
    Triangle_2 t2(q1,q2,q3);

    Segment_2 s1(p1,p3), s2(p2, q1);
    CGAL::Object obj = CGAL::intersection(t1,t2);
    const std::vector<Point_2> *V;
    if( !(V = CGAL::object_cast<std::vector<Point_2> > (&obj))  ){
      std::cerr << "ERROR" << std::endl;
      return EXIT_FAILURE;
    }
    else{
      std::cerr << "OK" << std::endl;
      for(std::size_t i = 0; i < V->size(); i++){
        std::cerr << (*V)[i] << std::endl;
        std::cerr << CGAL::exact((*V)[i]) << std::endl;
      }
    }
    obj = CGAL::intersection(s1,s2);

    // check the variant return type
    const auto o_variant = CGAL::intersection(t1,t2);
    if(!o_variant) {
      std::cerr << "ERROR, empty" << std::endl;
      return EXIT_FAILURE;
    }

    V = nullptr;
    if( !(V = std::get_if<std::vector<Point_2> >(&(*o_variant)))  ){
      std::cerr << "ERROR, something other than vector< point_2 >" << std::endl;
      return EXIT_FAILURE;
    }
    else {
      std::cerr << "OK" << std::endl;
      for(std::size_t i = 0; i < V->size(); i++){
        std::cerr << (*V)[i] << std::endl;
        std::cerr << CGAL::exact((*V)[i]) << std::endl;
      }
    }
  }

  {
    Point_3 p1(0,0,1), p2(2,0,1), p3(1,3,1), q1(0,2,1), q2(2,2,1), q3(1,-1,1);
    Triangle_3 t1(p1,p2,p3);
    Triangle_3 t2(q1,q2,q3);

    CGAL::Object obj = CGAL::intersection(t1,t2);
    const std::vector<Point_3> *V;
    if( !(V = CGAL::object_cast<std::vector<Point_3> > (&obj))  ){
      std::cerr << "ERROR" << std::endl;
      return EXIT_FAILURE;
    }
    else{
      std::cerr << "OK" << std::endl;
      for(std::size_t i = 0; i < V->size(); i++){
        std::cerr << (*V)[i] << std::endl;
        std::cerr << CGAL::exact((*V)[i]) << std::endl;
      }
    }

    // check the variant return type
    const auto o_variant = CGAL::intersection(t1,t2);
    if(!o_variant) {
      std::cerr << "ERROR, empty" << std::endl;
      return EXIT_FAILURE;
    }

    V = nullptr;
    if( !(V = std::get_if<std::vector<Point_3> > (&(*o_variant)))  ){
      std::cerr << "ERROR" << std::endl;
      return EXIT_FAILURE;
    }
    else{
      std::cerr << "OK" << std::endl;
      for(std::size_t i = 0; i < V->size(); i++){
        std::cerr << (*V)[i] << std::endl;
        std::cerr << CGAL::exact((*V)[i]) << std::endl;
      }
    }
  }
  //making the interval construction failing
  {
    double eps = std::numeric_limits<double>::epsilon();
    std::cout << "Epsilon is " << eps << std::endl;
    Point_2 p1(0,0), p2(2,0), p3(1,3), q1(0,eps), q2(2,eps), q3(1,-1);
    Triangle_2 t1(p1,p2,p3);
    Triangle_2 t2(q1,q2,q3);

    Segment_2 s1(p1,p3), s2(p2, q1);
    CGAL::Object obj = CGAL::intersection(t2,t1);
    const std::vector<Point_2> *V;
    if( !(V = CGAL::object_cast<std::vector<Point_2> > (&obj))  ){
      std::cerr << "ERROR" << std::endl;
      return EXIT_FAILURE;
    }
    else{
      std::cerr << "OK" << std::endl;
      for(std::size_t i = 0; i < V->size(); i++){
        std::cerr << "A " << (*V)[i] << std::endl;
        std::cerr << "E " << CGAL::exact((*V)[i]) << std::endl;
      }
      assert(V->size()==6);
    }
    obj = CGAL::intersection(s1,s2);
  }


}
