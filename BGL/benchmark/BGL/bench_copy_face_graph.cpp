
#define CGAL_CFG_USE_VECTOR

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Timer.h>

#include <fstream>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Surface_mesh<K::Point_3> SM;



int main(int argc, char* argv[])
{

  CGAL::Timer t;
  t.start();
  Polyhedron A, B, C;
  SM smA, smB, smC;

  {
    std::ifstream in(argv[1]);
    //in >> A;
  }

  {
    std::ifstream in(argv[1]);
    in >> smA;
  }

  t.stop();
  std::cerr << "\nreading " << t.time() << std::endl;
   t.reset(); t.start();
  /*
  copy_face_graph(A,B);
  t.stop();
  std::cerr << "\nP -> P " << t.time();
   t.reset(); t.start();
  */
  copy_face_graph(smA,smB);
  t.stop();
  std::cerr << "\nS -> S " << t.time() << std::endl;
   t.reset(); t.start();
  
 

  copy_face_graph_2(smA,smC);
  t.stop();
  std::cerr << "\nS -> S  V2  " << t.time();

  /*
   t.reset(); t.start();
  
  copy_face_graph(smA,C);
  t.stop();
  std::cerr << "\nS -> P " << t.time() << std::endl;

   */

  return 0;
}
