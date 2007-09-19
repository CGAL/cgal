//#define CGAL_CHECK_EXPENSIVE
//#define CGAL_CHECK_EXACTNESS
#include <CGAL/Arrangement_of_spheres_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Qt_multithreaded_examiner_viewer_2.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_qt_viewer.h>
#include <CGAL/Arrangement_of_spheres_3/read_spheres.h>
#include <vector>

typedef CGAL::Simple_cartesian<CGAL::Gmpq> K;

typedef K::Point_3 P3;
typedef K::Vector_3 V3;
typedef K::Vector_2 V2;
typedef K::FT FT;
typedef K::Sphere_3 S3;
typedef K::Circle_2 C2;
typedef K::Point_2 P2;
typedef K::Line_3 L3;
typedef K::Line_2 L2;
typedef K::Plane_3 Pl3;

K k;
K::Intersect_3 i3= k.intersect_3_object();

Pl3 sweep_plane(FT x) {
  P3 p(x,0,0);
  V3 v(1,0,0);
  return Pl3(p,v);
}

FT squared_depth(P3 p, FT x){
  return CGAL::square(p[0]-x);
}


bool intersect(S3 s, FT x, C2 &c){
  FT r2= s.squared_radius() - squared_depth(s.center(), x);
  if (r2 <0) return false;
  CGAL_assertion(r2>=0);
  c= C2(P2(s.center()[1],
	   s.center()[2]), r2);
  return true;
}

bool intersect(L3 s, FT x, P2 &p2){
  CGAL::Object o= i3(s, sweep_plane(x));
  P3 p;
  if (assign(p, o)) {
    p2=P2(p.y(), p.z());
    return true;
  } else {
    return false;
  }
}

bool intersect(Pl3 s, FT x, L2& l2){
  CGAL::Object o= i3(s, sweep_plane(x));
  L3 p;
  if (assign(p, o)) {
    l2=L2(P2(p.point().y(), p.point().z()), V2(p.to_vector().y(), p.to_vector().z()));
    return true;
  } else {
    return false;
  }
}


struct Show_circles {
  Show_circles(const std::vector<S3> &ss,
	       const std::vector<L3> &l, 
	       const std::vector<Pl3> &p, FT z):
    spheres_(ss), lines_(l), planes_(p), z_(z){}

  void operator()(CGAL::Qt_examiner_viewer_2 *qtv) {
    QMutexLocker lock(qtv->lock());
    for (unsigned int i=0; i< spheres_.size(); ++i) {
      C2 c;
      if (intersect(spheres_[i], z_,c )) {
	*qtv << c;
	std::cout << c << std::endl;
	double ar= std::sqrt(CGAL::to_double(c.squared_radius()));
	P2 pt(c.center().x()+FT(ar), c.center().y());
	*qtv << pt;
	std::ostringstream oss;
	oss << i;
	*qtv << oss.str();
      } else {
	std::cout << "Sphere " << i << " missed" << std::endl;
      }
    }
    for (unsigned int i=0; i< lines_.size(); ++i) {
      P2 p2;
      if (intersect(lines_[i], z_, p2)) {
	*qtv << p2;
	std::ostringstream oss;
	oss << "l" << i;
	*qtv << oss.str();
      } else {
	std::cout << "Line " << i << " missed" << std::endl;
      }
    }

    for (unsigned int i=0; i< planes_.size(); ++i) {
      L2 p2;
      if (intersect(planes_[i], z_, p2)) {
	*qtv << p2;
	V3 offset=(.0001 * planes_[i].orthogonal_vector());
	P3 pt=planes_[i].point() +offset;
	Pl3 np3(pt,
	       planes_[i].orthogonal_vector());
	intersect(np3, z_, p2);
	*qtv << CGAL::GRAY;
	*qtv << p2;
	*qtv << CGAL::BLACK;
      } else {
	std::cout << "Plane " << i << " missed" << std::endl;
      }
    }
    
    //qtv->unlock();
    qtv->show_everything();
  }

  std::vector< S3> spheres_;
  std::vector< L3> lines_;
  std::vector< Pl3> planes_;
  FT z_;
};

 

int main(int argc, char *argv[]) {
#ifdef CGAL_AOS3_USE_TEMPLATES
  //typedef CGAL::Simple_cartesian<double> K;
  typedef CGAL::Arrangement_of_spheres_traits_3<K> Traits;
  typedef CGAL::Arrangement_of_spheres_3<Traits> Arrangement;
#else 
  typedef CGAL::Arrangement_of_spheres_3 Arrangement;
#endif

  /*char buf[]="0/1 0/1 0/1 1/1 1";
  std::istringstream iss(buf);
  S3 s3;
  iss >> s3;
  std::cout << s3 << std::endl;*/


  std::vector<S3> spheres;
  std::vector<L3> lines;
  std::vector<Pl3> planes;

  std::ifstream in(argv[2]);
  
  while (true) {
    char buf[10000];
    in.getline(buf, 10000);
    if (!in) break;
    //std::cout << buf << std::endl;
    if (buf[0]=='#' || buf[0]=='\0') {
      //std::cout << "Line is comment: " << buf << std::endl;
      continue;
    } else if (buf[0]== 'l') {
      //std::cout << "Reading line " << std::endl;
      std::istringstream iss(buf+2);
      L3 l3;
      iss >> l3;
      if (iss) {
	lines.push_back(l3);
	//std::cout << "l " << l3 << std::endl;
      } else {
	std::cout << "Error reading from line " << buf << std::endl;
      }
    } else if (buf[0]== 'p') {
      //std::cout << "Reading plane " << std::endl;
      std::istringstream iss(buf+2);
      Pl3 l3;
      iss >> l3;
      if (iss) {
	planes.push_back(l3);
	//std::cout << "p " << l3 << std::endl;
      } else {
	std::cout << "Error reading from line " << buf << std::endl;
      }
    } else {
      std::istringstream iss(buf);
      S3 l3;
      iss >> l3;
      if (iss) {
	spheres.push_back(l3);
	//std::cout << "s " << l3 << std::endl;
      } else {
	std::cout << "Error reading from line " << buf << std::endl;
      }
    } 
  }
  
  FT z;
  {
    std::istringstream iss(argv[1]);
    iss >> z;
    if (!iss) return EXIT_FAILURE;
  }
  typedef Show_circles CS;
  CS cs(spheres, lines, planes, z);
  
  CGAL::Qt_multithreaded_examiner_viewer_2<CS> qtv(cs, argc, argv);
  

  return qtv();
}
