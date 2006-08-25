#include <CGAL/basic.h>
#include <iostream>
#include <string>
#include <list>
#include <sstream>
#include <cassert>
#include <CGAL/Cartesian.h>
#include <boost/variant.hpp>
#include <CGAL/intersections.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Circular_kernel_2.h>

#include <CGAL/point_generators_2.h>

#include <CGAL/MP_Float.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>

using namespace std;
class Dxf_converter {
  //typedef std::pair<double, double> Point_2;
/*
    typedef CGAL::Quotient<CGAL::MP_Float>                       NT1;
  typedef CGAL::Cartesian<NT1>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_2_2<NT1>                      Algebraic_k1;
  typedef CGAL::Curved_kernel<Linear_k1, Algebraic_k1>         CircularKernel;*/
// typedef double                                              type;
// typedef CGAL::Cartesian<type>                       NT1; 
//typedef CGAL::Quotient<CGAL::MP_Float>                       NT1;
typedef CGAL::Gmpq						NT1;
//    typedef CGAL::Quotient<NT1>                       NT1;
  typedef CGAL::Cartesian<NT1>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT1>                      Algebraic_k1;
  typedef CGAL::Circular_kernel_2<Linear_k1, Algebraic_k1>         CK;
 typedef  CK::Circular_arc_point_2 Circular_arc_point_2;
  typedef  CK::Line_2  Line_2;
  typedef CK::Point_2 Point_2;
  typedef  CK::Circle_2 Circle_2;  
  typedef std::list<std::pair<Point_2, double> > Polygon;
  typedef std::list<Polygon> Polygons;
 typedef CK::Intersect_2   Intersect_2;
 typedef CK::FT FT;
private:

  void
  header(std::istream& is)
  {
    int n;
    double xmin, ymin;
    double xmax, ymax;
    is >> n;
    assert(n == 9);
    char c;
    is >> c;
    assert(c == '$');
    std::string str;
    is >> str;
    if(str == std::string("EXTMIN")){
      is >> n;
      assert(n == 10);
    is >> xmin;
    is >> n;
    assert(n == 20);
    is >> ymin;
    }
    is >> n;
    assert(n == 9);
    is >> c;
    assert(c == '$');
    is >> str;
    if(str == "EXTMAX"){
      is >> n;
      assert(n == 10);
      is >> xmax;
      is >> n;
      assert(n == 20);
      is >> ymax;
    }
  }
  
  
  void
  skip_header(std::istream& is)
  {
    int n;
    is >> n;
    assert(n == 0);
    std::string str;
    is >> str;
    assert(str == "SECTION");
    is >> n;
    assert(n == 2);
    is >> str;
    if(str == "HEADER"){
      header(is);
    }
    is >> n;
    assert(n == 0);
    is >> str;
    assert(str == "ENDSEC");
  }
  
  
  
  void 
  read_circle(std::istream& is,std::ostream& os)
  {
    int n;
    double cx, cy, r;
    std::string str;
    is >> n;
    assert(n == 8);
    is >> n;
    assert(n == 0);
  
  is >> n;
  assert(n == 10);
  is >> cx;
  is >> n;
  assert(n == 20);
  is >> cy;
  is >> n;
  assert(n == 40);
  is >> r;
 NT1 x(cx),y(cy),gr(r);
  os<<"Circle_2 ("<<"Point_2("<<"Rational(" <<x.numerator()<< ", "<<x.denominator()<<")"<< ", "<<"Rational(" <<y.numerator()<< ", "<<y.denominator()<<")""),"<<"Rational(" <<gr.numerator()<< ", "<<gr.denominator()<<")"<<")"<<std::endl;
//  circ = typename K::Construct_circle_2()(Point_2(cx,cy), r*r);
}

void
write_polygon(Polygons polygons,std::ostream& os){


  Point_2 first_point;
  Point_2 ps;
  Point_2 pt ;
  Point_2 center;
  FT bulge;
  for( Polygons::iterator it = polygons.begin(); it != polygons.end(); it++){
     Polygon::iterator pit = it->begin();
    first_point = pit->first;    
    while(true){
      ps = pit->first;
      bulge = pit->second;
      pit++;
      if(pit ==it->end()){
	break;
      }
      pt = pit->first;
      if(bulge == FT(0)){
	if(ps != pt){
	  os<<"Line_arc_2("<<"Point_2("<<"Rational(" <<ps.x().numerator()<< ", "<<ps.x().denominator()<<"), "<<"Rational(" <<ps.y().numerator()<< ", "<<ps.y().denominator()<<")),"<<"Point_2("<<"Rational(" <<pt.x().numerator()<< ", "<<pt.x().denominator()<<"), "<<"Rational(" <<pt.y().numerator()<< ", "<<pt.y().denominator()<<")))"<<std::endl;
	}
      } else {
os<<"Circular_arc_2( Point_2( Rational(" <<ps.x().numerator()<< ", "<<ps.x().denominator()<<")"<< ", "<<"Rational(" <<ps.y().numerator()<< ", "<<ps.y().denominator()<<")),"<<"Point_2("<<"Rational(" <<pt.x().numerator()<< ", "<<pt.x().denominator()<<"), "<<"Rational(" <<pt.y().numerator()<< ", "<<pt.y().denominator()<<")),"<<"Rational(" <<bulge.numerator()<< ", "<<bulge.denominator()<<"))"<<std::endl;
//         const FT sqr_bulge = CGAL::square(bulge);
//         const FT common = (FT(1) - sqr_bulge) / (FT(4)*bulge);
//         const FT x_coord = ((ps.x() + pt.x())/FT(2)) + common*(ps.y() - pt.y());
//         const FT y_coord = ((ps.y() + pt.y())/FT(2)) + common*(pt.x() - ps.x());
//         const FT sqr_rad = CGAL::squared_distance(ps, pt) * (FT(1)/sqr_bulge + FT(2) + sqr_bulge) / FT(16); 
//  	Circular_arc_point_2 cps = ps;
//  	Circular_arc_point_2 cpt = pt;
// 	os << "Circular_arc_2("<<"Circle_2 ("<<"Point_2("<< x_coord<<","<<y_coord <<"),"<<sqr_rad <<"),"<< "Point_2("<< cps.x()<<","<<cps.y() <<"),"<<"Point_2("<< cpt.x()<<","<<cpt.y() <<"))"<<std::endl; 

	//Circular_arc_2 arc = Circular_arc_2(Circle_2(Point_2(x_coord, y_coord), sqr_rad), cps, cpt);
// os<<"Circular_arc_2("<<"Circle_2 ("<<"Point_2("<<"Rational(" <<x_coord.numerator()<< ", "<<x_coord.denominator()<<")"<< ", "<<"Rational(" <<y_coord.numerator()<< ", "<<y_coord.denominator()<<")""),"<<"Rational(" <<sqr_rad.numerator()<< ", "<<sqr_rad.denominator()<<")"<<"),"<<"Point_2("<<"Rational(" <<cps.x().numerator()<< ", "<<cps.x().denominator()<<"), "<<"Rational(" <<cps.y().numerator()<< ", "<<cps.y().denominator()<<")),"<<"Point_2("<<"Rational(" <<cpt.x().numerator()<< ", "<<cpt.x().denominator()<<"), "<<"Rational(" <<cpt.y().numerator()<< ", "<<cpt.y().denominator()<<")))"<<std::endl;
      }
    }
    if(bulge == FT(0)){
      if(ps != first_point){
	  os<<"Line_arc_2("<<"Point_2("<<"Rational(" <<ps.x().numerator()<< ", "<<ps.x().denominator()<<"), "<<"Rational(" <<ps.y().numerator()<< ", "<<ps.y().denominator()<<")),"<<"Point_2("<<"Rational(" <<first_point.x().numerator()<< ", "<<first_point.x().denominator()<<"), "<<"Rational(" <<first_point.y().numerator()<< ", "<<first_point.y().denominator()<<")))"<<std::endl;
      }
    } else {
      pt = first_point;
os<<"Circular_arc_2( Point_2( Rational(" <<ps.x().numerator()<< ", "<<ps.x().denominator()<<")"<< ", "<<"Rational(" <<ps.y().numerator()<< ", "<<ps.y().denominator()<<")),"<<"Point_2("<<"Rational(" <<pt.x().numerator()<< ", "<<pt.x().denominator()<<"), "<<"Rational(" <<pt.y().numerator()<< ", "<<pt.y().denominator()<<")),"<<"Rational(" <<bulge.numerator()<< ", "<<bulge.denominator()<<"))"<<std::endl;
//       const FT sqr_bulge = CGAL::square(bulge);
//       const FT common = (FT(1) - sqr_bulge) / (FT(4)*bulge);
//       const FT x_coord = ((ps.x() + pt.x())/FT(2)) + common*(ps.y() - pt.y());
//       const FT y_coord = ((ps.y() + pt.y())/FT(2)) + common*(pt.x() - ps.x());      
//       const FT sqr_rad = CGAL::squared_distance(ps, pt) * (FT(1)/sqr_bulge + FT(2) + sqr_bulge) / FT(16);       
//       Circular_arc_point_2 cps = ps;
//       Circular_arc_point_2 cpt = pt;
/*	os << "Circular_arc_2("<<"Circle_2 ("<<"Point_2("<< x_coord<<","<<y_coord <<"),"<<sqr_rad <<"),"<< "Point_2("<< cps.x()<<","<<cps.y() <<"),"<<"Point_2("<< cpt.x()<<","<<cpt.y() <<"))"<<std::endl;*/ 
// os<<"Circular_arc_2("<<"Circle_2 ("<<"Point_2("<<"Rational(" <<x_coord.numerator()<< ", "<<x_coord.denominator()<<")"<< ", "<<"Rational(" <<y_coord.numerator()<< ", "<<y_coord.denominator()<<")""),"<<"Rational(" <<sqr_rad.numerator()<< ", "<<sqr_rad.denominator()<<")"<<"),"<<"Point_2("<<"Rational(" <<cps.x().numerator()<< ", "<<cps.x().denominator()<<"), "<<"Rational(" <<cps.y().numerator()<< ", "<<cps.y().denominator()<<")),"<<"Point_2("<<"Rational(" <<cpt.x().numerator()<< ", "<<cpt.x().denominator()<<"), "<<"Rational(" <<cpt.y().numerator()<< ", "<<cpt.y().denominator()<<")))"<<std::endl;
    }
  }

}




void
read_polygon(std::istream& is,Polygon& poly)
{
  Polygons polygons;
  int n;
  int i;
  double x, y, len;
  std::string str;
 i=0;
  do {
    is >> n;
    if(n != 0){
      int m;
      is >> m; 
    }
  } while(n != 0);

  do {
    is >> str;
    if(str == "VERTEX"){
      is >> n;
      assert(n == 8);
      is >> n;
      assert(n == 0);
      is >> n;
      assert(n == 10);
      is >> x;
      is >> n;
      assert(n == 20);
      is >> y;
      is >> n;
      len = 0;
      if(n == 42){
	is >> len;
      } else {
	assert(n == 0);
      }
      //std::cout<<"Polygon"<<x<<y<<len<<std::endl;
      poly.push_back(std::make_pair(CK::Construct_point_2()(x,y), len));
     i++;
    }
    
  } while (str != "SEQEND");
  is >> n;
  assert(n == 8);
  is >> n;
  assert(n == 0);
  std::cout<<i<<std::endl;

}


void
read_entities(std::istream& is,std::ostream& os)
{Polygons poly;
  int n;
  //double x, y;
  std::string str;
  is >> n;
  assert(n == 0);
  is >> str;
  assert(str == "SECTION");
  is >> n;
  is >> str;
  assert(str == "ENTITIES");
  do {
    is >> n;
    assert(n == 0);
    is >> str;
    if(str == "POLYLINE"){
      Polygon p;
      poly.push_back(p); 
	std::cout<< "it's polyline" <<std::endl;
      read_polygon(is, poly.back());
    } else if(str == "CIRCLE"){
      read_circle(is,os);      
    } else if(str == "ENDSEC"){
      
    } else {
      std::cerr << "unknown entity" << std::endl;
      exit(0);
    }
  } while(str != "ENDSEC");
	std::cout<< "it's endsec" <<std::endl;
  is >> n;
  assert(n == 0);
  is >> str;
  assert(str == "EOF");
  std::cout << "read_entities - finished"<<std::endl;
  write_polygon(poly,os);
}

public:

void operator()(std::istream& is,std::ostream& os)
{
  skip_header(is);
  read_entities(is,os);
}

};
