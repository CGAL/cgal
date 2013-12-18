#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/box_intersection_d.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Segment_2 Segment_2;

typedef std::pair<Point_2,std::string> Point;
typedef std::list<Point>::iterator iterator;

typedef std::pair<iterator,iterator> Segment;


void
fix_loops(std::list<std::pair<Point_2,std::string> >& points)
{
  iterator pit = points.begin();
  while(true){
  
    Point_2& p = pit->first;
    iterator qit = pit;
    ++qit;
    if(qit == points.end()){
      break;
    }

    Point_2& q = qit->first;
    
    iterator rit = qit;
    ++rit;
    if(rit == points.end()){
      break;
    }
    Point_2& r = rit->first;

    iterator sit = rit;
    ++sit;
    if(sit == points.end()){
      break;
    }
    Point_2& s = sit->first;
    Segment_2 segA(p,q), segB(r,s);
    if(CGAL::do_intersect(segA,segB)){
      //std::cerr << "intersection\n2\n" << pit->second << std::endl << qit->second << std::endl;  
      //std::cerr << "2\n" << rit->second << std::endl << sit->second << std::endl;

      // without C++11
      CGAL::cpp11::result_of<Kernel::Intersect_2(Segment_2, Segment_2)>::type
        result = intersection(segA, segB);
     
      if(result) {
        if(const Segment_2* rs = boost::get<Segment_2>(&*result)) {
          std::cerr << "intersection is a segment" << std::endl;
        } else {
          const Point_2& rp = *boost::get<Point_2>(&*result);
          //std::cerr << "point " << rp << std::endl;
          //std::cerr << "perimeter " << std::sqrt(CGAL::squared_distance(rp,q))+ std::sqrt(CGAL::squared_distance(rp,r)) + std::sqrt(CGAL::squared_distance(r,q)) << std::endl;
          std::swap(*qit,*rit);
        }
      }
    }
    ++pit;
  }
}

void contour(int count, int n)
{
  //std::cerr << "contour  " << count << " with " << n << "vertices" << std::endl;
 
  std::cin.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');
  std::list<std::pair<Point_2,std::string> > points;

  for(int i = 0; i <n; i++){
    std::string line;
    std::getline(std::cin, line);
    std::istringstream iss(line);
    double x,y,z;
    iss >> x >> y >> z;
    Point_2 p(x,z);
    if(! points.empty()){
      if(p != points.back().first){
      points.push_back(std::make_pair(p, line));
      } else {
        std::cerr << "identical consecutive points " << p << count << std::endl;
      }
    } else {
     points.push_back(std::make_pair(p, line)); 
    }
  }
  if(points.size() == 3){
    std::cerr << "ignore segment" <<  std::endl;
    std::cerr << points.front().first << " -- " << (++points.begin())->first << std::endl;
  }

  // in the data the first point is duplicated
  // to make the loop easier also duplicate the second point

  points.push_back(*(++points.begin()));
  iterator pit = points.begin();
  while(true){
    if(points.size() == 4){
      break;
    }
    iterator qit = pit;
    ++qit;
    if(qit == points.end()){
      break;
    }
    Point_2& p = pit->first;
    Point_2& q = qit->first;
    
    iterator rit = qit;
    ++rit;
    if(rit == points.end()){
      break;
    }
    Point_2& r = rit->first;

    if(CGAL::collinear(p,q,r)){
      if(! CGAL::collinear_are_strictly_ordered_along_line(p,q,r)){
        std::cerr << "not strictly collinear" << std::endl;
        std::cerr << p << std::endl;
        std::cerr << q << std::endl;
        std::cerr << r << std::endl;
        std::cerr << " in  polygon " << count << std::endl;
        points.erase(qit);
        if(pit->first == rit->first){
          std::cerr << "identical consecutive points " << pit->first << std::endl;
          points.erase(rit);
        }
        if(pit != points.begin()){
            --pit;
          }
        continue;  // and do not advance pit, because its two successors are different
      }
    }
    ++pit; 
  }
 

  if(points.size() == 4){
    std::cerr << "ignore segment" <<  std::endl;
    std::cerr << points.front().first << " -- " << (++points.begin())->first << std::endl;
  } else {
    fix_loops(points);

    // remove the added second point;
    points.pop_back(); 

    std::cout << points.size() <<std::endl;
    for(iterator it = points.begin(); it!= points.end(); ++it){
      std::cout << it->second << std::endl;
    }

  }
  
}

int main()
{
 
  int n;
  int count = 0;
  while(std::cin >> n){
    contour(count, n);
    count++;
  }
  std::cerr << "done" << std::endl;
  return 0;
}
