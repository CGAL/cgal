#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Segment_2 Segment_2;



void contour(int count, int n)
{
  std::cerr << "contour  " << count << " with " << n << "vertices" << std::endl;
 
  std::cin.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');
  std::list<std::pair<Point_2,std::string> > points;
  typedef std::list<std::pair<Point_2,std::string> >::iterator iterator;
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
  // remove the added second point;
  points.pop_back(); 


  if(points.size() == 3){
    std::cerr << "ignore segment" <<  std::endl;
    std::cerr << points.front().first << " -- " << (++points.begin())->first << std::endl;
  } else {

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
