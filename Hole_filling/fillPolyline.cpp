#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include <vector>

#include "Fill_hole.h"

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;


int main()
{
  int n;
  Point_3 p;
  std::list<Point_3> P, Q;
  std::cin >> n;
  while(n--){
    std::cin >> p;
    P.push_back(p);
  }

 std::cin >> n;
  while(n--){
    std::cin >> p;
    Q.push_back(p);
  }
 
  std::vector<CGAL::Triple<int,int,int> > facets;
  CGAL::fill_hole(P.begin(), P.end(),Q.begin(), Q.end(), std::back_inserter(facets));


  std::cout.precision(20);
  std::cout << "OFF\n" << P.size() << " " << facets.size() << " 0" << std::endl;

  for (std::list<Point_3>::iterator it = P.begin(); it != P.end(); ++it){
    std::cout << *it << std::endl;
  }
  for(int i = 0; i < facets.size(); i++){
    std::cout << "3 " << facets[i].first << " " << facets[i].second << " " << facets[i].third << std::endl;
  }


  std::cerr << "done" << std::endl;
  return 0;
}
