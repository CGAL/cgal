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


struct Show_circles {
  Show_circles(const std::vector<C2> &ss,
	       const std::vector<L2> &l):
    spheres_(ss), lines_(l){}

  void operator()(CGAL::Qt_examiner_viewer_2 *qtv) {
    
    for (unsigned int i=0; i< spheres_.size(); ++i) {
      *qtv << spheres_[i];
      double ar= std::sqrt(CGAL::to_double(spheres_[i].squared_radius()));
      P2 pt(spheres_[i].center().x()+FT(ar), spheres_[i].center().y());
      *qtv << pt;
      std::ostringstream oss;
      oss << i;
      *qtv << oss.str();
      
    }
    for (unsigned int i=0; i< lines_.size(); ++i) {
      *qtv << lines_[i];
    }
    
    
    
    
    qtv->show_everything();
  }
  
  std::vector< C2> spheres_;
  std::vector< L2> lines_;
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


  std::vector<C2> spheres;
  std::vector<L2> lines;
 
  std::ifstream in(argv[1]);
  
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
      L2 l3;
      iss >> l3;
      if (iss) {
	lines.push_back(l3);
	//std::cout << "l " << l3 << std::endl;
      } else {
	std::cout << "Error reading from line " << buf << std::endl;
      }
    } else {
      std::istringstream iss(buf);
      C2 l3;
      iss >> l3;
      if (iss) {
	spheres.push_back(l3);
	//std::cout << "s " << l3 << std::endl;
      } else {
	std::cout << "Error reading from line " << buf << std::endl;
      }
    } 
  }
  
 
  typedef Show_circles CS;
  CS cs(spheres, lines);
  
  CGAL::Qt_multithreaded_examiner_viewer_2<CS> qtv(cs, argc, argv);
  

  return qtv();
}
