//#include <CGAL/leda_integer.h>
//#include <CGAL/leda_rational.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Lazy_kernel.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/Nef_3/convex_decomposition_3.h> 
#include <CGAL/convexity_check_3.h>

#include <CGAL/Nef_3/trunk_offset.h>

#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>

#include <fstream>
#include <sstream>

#ifdef CGAL_WITH_LAZY_KERNEL
typedef CGAL::Gmpq NT;
//typedef leda_rational NT;
typedef CGAL::Lazy_kernel<CGAL::Simple_cartesian<NT> > Kernel;
#else
typedef CGAL::Gmpz NT;
typedef CGAL::Homogeneous<NT> Kernel;
#endif
typedef Kernel::RT RT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 Plane_3;
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
typedef CGAL::Trunk_offset<Nef_polyhedron_3> TO;

class Volume_output {
  bool twin;
  std::vector<Plane_3> planes;
  typedef std::vector<Plane_3>::const_iterator CI;
public:
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

};

void read( const char* name, Polyhedron_3& poly) {
  std::ifstream in( name);
  if ( ! in) { 
    std::cerr << "minkowsky_sum: error: cannot open file '"<< name
	      << "' for reading." << std::endl;
    std::exit( 1);
  }
  in >> poly;
  if ( ! in) { 
    std::cerr << "minkowsky_sum: error: reading from file '"<< name << "'."
	      << std::endl;
    std::exit( 1);
  }
}

int main(int argc, char* argv[]) {

  if ( argc < 3 || argc > 5) {
    std::cerr << "Usage: " << argv[0] << " <infile>" << std::endl;
    std::exit(1);
  }

  std::ifstream trunk( argv[2] );
  char buffer[80];
  int nv,nf,nc,nfl,col,man;
  std::vector<Point_3> points;
  std::list<std::pair<int*, int*> > facets;

  trunk.getline(buffer,80);
  int version = std::atoi(buffer+10);
  trunk.getline(buffer,80);
  trunk >> nv;
  trunk.getline(buffer,80);
  trunk.getline(buffer,80);
  trunk >> nf;
  if(version > 3) {
    trunk.getline(buffer,80);
    trunk.getline(buffer,80);
    trunk >> nfl;
  }
  trunk.getline(buffer,80);
  trunk.getline(buffer,80);
  trunk >> nc;

  trunk.getline(buffer,80);
  for(int i=0;i<nc;++i)
    trunk.getline(buffer,80);

  trunk.getline(buffer,80);

  std::cerr << "number of vertices " << nv << std::endl;
  for(int i=0;i<nv;++i) {
    double a,b,c;
    trunk >> a >> b >> c;
    RT x(round(a)), y(round(b)), z(round(c));
    //    FT x(a), y(b), z(c);
/*
    Point_3 p(x.numerator()   * y.denominator() * z.denominator(),
	      x.denominator() * y.numerator()   * z.denominator(),
	      x.denominator() * y.denominator() * z.numerator(),
	      x.denominator() * y.denominator() * z.denominator() );
*/
    Point_3 p(x,y,z,1);
//    p = normalized(p);
    points.push_back(p);
  }
 
  std::cerr << "number of facets " << nf << std::endl;
  trunk.getline(buffer,80);
  trunk.getline(buffer,80);

  int mod = argc > 3 ? std::atoi(argv[3]) : 256; 
  int suf = argc > 4 ? std::atoi(argv[4]) : 0;
  int off = argc > 5 ? std::atoi(argv[5]) : 0;

  int face[3*nf];
  for(int i=0;i<nf;++i){
    trunk >> face[3*i] >> face[3*i+1] >> face[3*i+2] >> col >> man;
    if(i>=off)
      facets.push_back(std::make_pair(face+(3*i),face+(3*i+3))); 
  }

  Polyhedron_3 P;
  read( argv[1], P);

  CGAL_assertion(is_strongly_convex_3(P));

  CGAL::Timer t;
  t.start();
  TO to(mod, suf);
  Nef_polyhedron_3 CSP = to(points.begin(), points.end(),
			    facets.begin(), facets.end(), P);

  CSP = !CSP;

  t.stop();
  std::cerr << "Runtime CSP: " << t.time() << std::endl;
  std::ofstream outCSP("csp.nef3");
  outCSP << CSP;
  t.start();

  points.clear();
  Vertex_const_iterator v;
  for(v = CSP.vertices_begin(); v != CSP.vertices_end(); ++v)
    points.push_back(v->point());

  Polyhedron_3 CV;
  convex_hull_3( points.begin(), points.end(), CV);

  Nef_polyhedron_3 NCV(CV);
  Nef_polyhedron_3 DIFF = NCV-CSP;

  t.stop();
  std::cerr << "Runtime CSP,CV,(CV-CSP): " << t.time() << std::endl;
  std::ofstream outNCV("ncv.nef3");
  std::ofstream outDIFF("diff.nef3");
  outNCV << NCV;
  outDIFF << DIFF;
  t.start();
  
  convex_decomposition_3<Nef_polyhedron_3>(DIFF);

  t.stop();
  std::cerr << "Runtime CSP,CV,(CV-CSP),Decomposition: " 
	    << t.time() << std::endl;
  t.start();


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
      vout.dump();
    }
  }

  t.stop();
  std::cerr << "Runtime CSP,CV,(CV-CSP),Decomposition,Output: " 
	    << t.time() << std::endl;    

  /*
  std::ofstream out("temp.nef3");
  out << result;

  QApplication a(argc, argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron_3>* w =
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron_3>(result);
  a.setMainWidget(w);
  w->show();
  a.exec();
  */
}
