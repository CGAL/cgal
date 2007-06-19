#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#else
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#endif
#include <CGAL/Lazy_kernel.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/Nef_3/convex_decomposition_3.h> 
#include <CGAL/convexity_check_3.h>

#include <fstream>
#include <sstream>

//#define CGAL_WITH_LAZY_KERNEL
#ifdef CGAL_WITH_LAZY_KERNEL
typedef CGAL::Gmpq NT;
//typedef leda_rational NT;
typedef CGAL::Lazy_kernel<CGAL::Simple_cartesian<NT> > Kernel;
#else
#ifdef CGAL_USE_LEDA
typedef leda_integer NT;
#else
typedef CGAL::Gmpz NT;
#endif
typedef CGAL::Homogeneous<NT> Kernel;
#endif
typedef Kernel::RT RT;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Aff_transformation_3 Aff_transformation_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
//typedef Polyhedron_3::Vertex_const_iterator Vertex_const_iterator;
#ifdef CGAL_NEF_INDEXED_ITEMS
typedef CGAL::Nef_polyhedron_3<Kernel,CGAL::SNC_indexed_items>     Nef_polyhedron_3;
#else
typedef CGAL::Nef_polyhedron_3<Kernel>     Nef_polyhedron_3;
#endif
typedef Nef_polyhedron_3::Vertex_const_iterator Vertex_const_iterator;
typedef Nef_polyhedron_3::Vertex_const_handle Vertex_const_handle;
typedef Nef_polyhedron_3::Halfedge_const_handle Halfedge_const_handle;
typedef Nef_polyhedron_3::Halffacet_const_handle Halffacet_const_handle;
typedef Nef_polyhedron_3::Volume_const_iterator Volume_const_iterator;
typedef Nef_polyhedron_3::SHalfedge_const_handle SHalfedge_const_handle;
typedef Nef_polyhedron_3::SHalfloop_const_handle SHalfloop_const_handle;
typedef Nef_polyhedron_3::SFace_const_handle SFace_const_handle;

class Volume_output {
  bool twin;
  std::vector<Plane_3> planes;
public:

  typedef std::vector<Plane_3>::const_iterator CI;

  Volume_output(bool twin_ = true) : twin(twin_){}
  bool is_in(const Plane_3 p) {
    for(CI ci = planes.begin(); ci != planes.end(); ++ci)
      if(*ci == p)
	return true;
    return false;
  }

  void visit(Halffacet_const_handle f) {
    if(twin) f = f->twin();
    Plane_3 p = f->twin()->plane();
    if(!is_in(p))
      planes.push_back(p);
  }
  void visit(SFace_const_handle s) {}
  void visit(Halfedge_const_handle e) {}
  void visit(Vertex_const_handle v) {}
  void visit(SHalfedge_const_handle se) {}
  void visit(SHalfloop_const_handle sl) {}

  void dump() const {
    std::cout << planes.size() << std::endl;
    for(CI ci = planes.begin(); ci != planes.end(); ++ci) {
#ifdef CGAL_WITH_LAZY_KERNEL
      std::cout << CGAL::to_double(ci->a().exact()) << " " 
		<< CGAL::to_double(ci->b().exact()) << " "
		<< CGAL::to_double(ci->c().exact()) << " "
		<< CGAL::to_double(ci->d().exact()) << std::endl;
#else
      std::cout << *ci << std::endl;
#endif
    }
  }
    
  CI planes_begin() const {
    return planes.begin();
  }
    
  CI planes_end() const {
    return planes.end();
  }
};

template <typename CI>
Nef_polyhedron_3 create_from_halfspaces(CI begin, CI end, bool invert) {
  
  Nef_polyhedron_3 cube;
  std::ifstream in("Nef3/centered_cube.nef3");
  in >> cube;

  int mag = 1000000;
  Aff_transformation_3 scale(mag, 0, 0, 
			     0, mag, 0, 
			     0, 0, mag, 1);
  cube.transform(scale);
 
  CI pli;
  for(pli = begin; pli != end; ++pli) {
    if(invert) {
      Plane_3 inverted_plane(-pli->a(), -pli->b(),
			     -pli->c(), -pli->d());
      cube = cube.intersection(inverted_plane);
    } else
      cube = cube.intersection(*pli);
  }  
  return cube;
}

int main(int argc, char* argv[]) {

  Nef_polyhedron_3 CSP;
  std::cin >> CSP;

  if(CSP.number_of_vertices()==0)
    return 0;

  std::list<Point_3> points;
  Vertex_const_iterator v;
  for(v = CSP.vertices_begin(); v != CSP.vertices_end(); ++v)
    points.push_back(v->point());

  Polyhedron_3 CV;
  convex_hull_3( points.begin(), points.end(), CV);

  Nef_polyhedron_3 NCV(CV);
  Nef_polyhedron_3 DIFF = NCV-CSP;

  convex_decomposition_3<Nef_polyhedron_3>(DIFF);

  int nov = 1;
  Volume_const_iterator c;
  for(c=++(DIFF.volumes_begin()); c!=DIFF.volumes_end(); ++c)
    if(c->mark()) ++nov;
  std::cout << nov << std::endl;

  Volume_output vout_ncv;
  NCV.visit_shell_objects(NCV.volumes_begin()->shells_begin(), vout_ncv);
  vout_ncv.dump();

  for(c=++(DIFF.volumes_begin()); c!=DIFF.volumes_end(); ++c) {
    if(c->mark()) {
      Volume_output vout(true);
      DIFF.visit_shell_objects(c->shells_begin(), vout);
      
      Nef_polyhedron_3 obs1 = 
	create_from_halfspaces(vout.planes_begin(), vout.planes_end(), true);

      Polyhedron_3 pobs2;
      DIFF.convert_inner_shell_to_polyhedron(c->shells_begin(), pobs2);
      Nef_polyhedron_3 obs2(pobs2);      

      Nef_polyhedron_3 empty = obs1.symmetric_difference(obs2);
      if(!empty.is_empty())
	std::cerr << "error " << std::endl;
      vout.dump();
    }
  }
}
