//#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Lazy_kernel.h>
//#include <CGAL/Fractional_traits.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

#include <CGAL/convexity_check_3.h>

#include <CGAL/Nef_3/trunk_offset.h>

#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>

#include <fstream>
#include <sstream>

//typedef CGAL::Gmpq FT;
typedef leda_rational FT;
//typedef CGAL::Gmpz NT;
//typedef CGAL::Homogeneous<NT> Kernel;
typedef CGAL::Lazy_kernel<CGAL::Simple_cartesian<FT> > Kernel;
//typedef Fractional_traits<Kernel::FT> FracTraits;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef Polyhedron_3::Vertex_const_iterator Vertex_const_iterator;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3;
typedef CGAL::Trunk_offset<Kernel> TO;

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
  int nv,nf,nc,col,man;
  std::vector<Point_3> points;
  std::list<std::pair<int*, int*> > facets;

  trunk.getline(buffer,80);
  trunk.getline(buffer,80);
  trunk >> nv;
  trunk.getline(buffer,80);
  trunk.getline(buffer,80);
  trunk >> nf;
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
    FT x(round(a)), y(round(b)), z(round(c));
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
  Nef_polyhedron_3 result = to(points.begin(), points.end(),
			       facets.begin(), facets.end(), P);

  t.stop();
  std::cerr << "Runtime Boundary Offset: " << t.time() << std::endl;

  std::ofstream out("temp.nef3");
  out << result;

  QApplication a(argc, argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron_3>* w =
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron_3>(result);
  a.setMainWidget(w);
  w->show();
  a.exec();
}
