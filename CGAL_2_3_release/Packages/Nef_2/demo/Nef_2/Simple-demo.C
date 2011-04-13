#ifdef CGAL_USE_LEDA
#include <CGAL/basic.h>
#include <CGAL/leda_integer.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/IO/Nef_polyhedron_2_Window_stream.h>

template <>
struct ring_or_field<leda_integer> {
  typedef ring_with_gcd kind;
  typedef leda_integer RT;
  static RT gcd(const RT& r1, const RT& r2) 
  { return ::gcd(r1,r2); }
};

typedef CGAL::Extended_homogeneous<leda_integer> EKernel;
typedef CGAL::Nef_polyhedron_2<EKernel> Nef_polyhedron;
typedef Nef_polyhedron::Point     Point;
typedef Nef_polyhedron::Line      Line;

static CGAL::Window_stream W(600,600);

int main()
{
  SETDTHREAD(71); 
  CGAL::set_pretty_mode ( std::cerr );
  std::cerr << "using " << CGAL::pointlocationversion << std::endl;
  std::cerr << "using " << CGAL::sweepversion << std::endl;

  W.init(-CGAL::frame_default,CGAL::frame_default,-CGAL::frame_default);
  W.set_show_coordinates(true);
  W.set_grid_mode(5);
  W.set_node_width(3);
  W.display(0,0);

  Point p1(-2,-1),p2(2,1);
  Line l1(p1,p2), l2(3,7,23);
  Nef_polyhedron N1(l1), N2(l2);

  W<<N1; W.read_mouse();
  W<<N2; W.read_mouse();
  W<< (N1*N2); W.read_mouse();

  return 0;
}

#else // CGAL_USE_LEDA

int main() { return 0; }

#endif // CGAL_USE_LEDA

