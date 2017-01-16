#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Timer.h>

#include <fstream>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Surface_mesh<K::Point_3> SM;
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> Polyhedron_with_ids;

template<class Mesh>
void run(const char* fname)
{
  CGAL::Timer t;
  t.start();
  Mesh A, B;

  std::ifstream in(fname);
  in >> A;
  in.close();
  t.stop();
  std::cerr << "  reading " << t.time() << std::endl;

  t.reset(); t.start();
  copy_face_graph_old(A,B);
  t.stop();
  std::cerr << "  S -> S OLD  " << t.time() << std::endl;
  clear(B);

  t.reset(); t.start();
  copy_face_graph(A,B);
  t.stop();
  std::cerr << " S -> S  NEW  " << t.time() << std::endl;
  
}

int main(int argc, char* argv[])
{

  std::cerr << "With Polyhedron\n";
  run<Polyhedron>(argv[1]);
  std::cerr << "With Surface_mesh\n";
  run<SM>(argv[1]);
  std::cerr << "With Polyhedron_with_ids\n";
  run<Polyhedron_with_ids>(argv[1]);

  return 0;
}
