#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/range_search_delaunay_2.h>
#include <list>


typedef CGAL::Exact_predicates_inexact_constructions_kernel  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt>                 Delaunay;
typedef CGAL::Delaunay_triangulation_2<Gt>::Edge_iterator  Edge_iterator;
typedef CGAL::Delaunay_triangulation_2<Gt>::Vertex_handle  Vertex_handle;

typedef Gt::Point_2   Point;
typedef Gt::Circle_2  Circle;
typedef Gt::Triangle_2  Triangle;
typedef Gt::Iso_rectangle_2  Rectangle_2;

Delaunay PS;

class check_empty {
public:
  bool    result;
  Circle  c;

  check_empty(Circle cact) : result(false), c(cact) { }

  bool get_result() const  { return result; }
  void set_result(bool nr) { result=nr; }

  bool operator()(const Point& p)
  {
    return ! c.has_on_unbounded_side(p);
  }
};

class check_empty_rectangle {
public:
  bool    result;
  Rectangle_2 r;

  check_empty_rectangle(Rectangle_2 ract) : result(false), r(ract) { }

  bool get_result() const  { return result; }
  void set_result(bool nr) { result=nr; }

  bool operator()(const Point& p)
  {
    return ! r.has_on_unbounded_side(p);
  }
};

class check_empty_triangle {
public:
  bool    result;
  Triangle t;

  check_empty_triangle(Triangle tact) : result(false), t(tact) { }

  bool get_result() const  { return result; }
  void set_result(bool nr) { result=nr; }

  bool operator()(const Point& p)
  {
    return ! t.has_on_unbounded_side(p);
  }
};


int main()
{
  Point p1(-253.357, -123.36);
  Point p2(-190.03, 216.606);
  Point p3(-343.349, 286.6);
  Point p4(141.604, 279.934);
  Point p5(276.591, -46.7012);
  Point p6(251.593, -263.347);
  Point p7(-3.38184, -343.339);
  Point p8(-380.012, -173.355);
  Point p9(-98.3726, 39.957);
  Point p10(133.271, 124.949);
  Point p11(289.923, 301.598);
  Point p12(421.577, 23.292);
  Point p13(79.9434, -93.3633);
  Point p14(-40.0449, 366.592);
  Point p15(311.587, 374.924);
  Point p16(431.576, 214.94);
  Point p17(426.576, -131.693);
  Point p18(-265.023, -285.011);
  Point p19(369.915, 89.9521);
  Point p20(368.249, -15.0376);
  Point p21(484.904, 18.2925);
  Point p22(-411.675, 283.267);
  Point p23(-250.024, 124.949);
  Point p24(-80.041, -78.3647);
  Point p25(-360.014, 31.6245);
  Point p26(-305.019, 356.593);

  // built Delaunay triangulation
  PS.insert(p1); PS.insert(p2); PS.insert(p3); PS.insert(p4);
  PS.insert(p5); PS.insert(p6); PS.insert(p7); PS.insert(p8);
  PS.insert(p9); PS.insert(p10); PS.insert(p11); PS.insert(p12);
  PS.insert(p13); PS.insert(p14); PS.insert(p15); PS.insert(p16);
  PS.insert(p17); PS.insert(p18); PS.insert(p19); PS.insert(p20);
  PS.insert(p21); PS.insert(p22); PS.insert(p23); PS.insert(p24);
  PS.insert(p25); PS.insert(p26);

  std::list<Vertex_handle> LV;

  bool correct = true;

  // circle emptiness check
  Circle cs1(Point(-23.3799, 108.284), 1124.78);
  check_empty checker(cs1);

  CGAL::range_search(PS,cs1,std::back_inserter(LV),checker,true);

  if (checker.get_result()) {
    std::cout << "circle not empty !\n";
    std::cout <<  "this is an error !\n"; correct=false;
  }
  else std::cout << "circle was empty !\n";

  Circle cs2(Point(-255.024, -100.029), 23551);
  check_empty checker2(cs2);

  CGAL::range_search(PS,cs2,std::back_inserter(LV),checker2,true);

  if (checker2.get_result()) std::cout << "circle not empty !\n";
  else {
    std::cout << "circle was empty !\n";
    std::cout <<  "this is an error !\n"; correct=false;
  }

  // triangle check
  Triangle t1(Point(-21.7134, -123.36), Point(84.9429, 74.9536), Point(209.931, -161.69));
  Triangle t2(Point(-61.7095, 164.945), Point(-88.3735, 101.618), Point(49.9463, 101.618));

  check_empty_triangle tchecker1(t1);
  CGAL::range_search(PS,t1.vertex(0),t1.vertex(1),t1.vertex(2),std::back_inserter(LV),tchecker1,true);

  if (tchecker1.get_result()) std::cout << "triangle not empty !\n";
  else {
    std::cout << "triangle was empty !\n";
    std::cout <<  "this is an error !\n"; correct=false;
  }

  check_empty_triangle tchecker2(t2);
  CGAL::range_search(PS,t2.vertex(0),t2.vertex(1),t2.vertex(2),std::back_inserter(LV),tchecker2,true);

  if (tchecker2.get_result()) {
     std::cout << "triangle not empty !\n";
     std::cout <<  "this is an error !\n"; correct=false;
  }
  else std::cout << "triangle was empty !\n";

  // rectangle check
  Rectangle_2 r1(-290.021, -175.022, -125.037, -35.0356);
  Rectangle_2 r2(-48.3774, 136.614, -23.3799, 251.603);

  check_empty_rectangle rchecker1(r1);
  CGAL::range_search(PS,r1.vertex(0),r1.vertex(1),r1.vertex(2),r1.vertex(3),std::back_inserter(LV),rchecker1,true);

  if (rchecker1.get_result()) std::cout << "rectangle not empty !\n";
  else {
    std::cout << "rectangle was empty !\n";
    std::cout <<  "this is an error !\n"; correct=false;
  }

  check_empty_rectangle rchecker2(r2);
  CGAL::range_search(PS,r2.vertex(0),r2.vertex(1),r2.vertex(2),r2.vertex(3),std::back_inserter(LV),rchecker2,true);

  if (rchecker2.get_result()) {
    std::cout << "rectangle not empty !\n";
    std::cout <<  "this is an error !\n"; correct=false;
  }
  else std::cout << "rectangle was empty !\n";

  if (correct) return 0;

  return 1;
}

