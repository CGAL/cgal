#include <CGAL/basic.h>
#include <CGAL/Filtred_container.h>
#include <CGAL/Random.h>

#include <list>
#include <iostream>

class Is_odd {
public:
  bool operator()(const int i) const
    {
      return (i%2)==1;
    }
};

typedef CGAL::Filtred_container<std::list<int>, Is_odd> List;

int main(int, char**)
{
  List l;

  int
    real_number_of_odds=0,
    detected_number_of_odds=0;

  for(int n=0; n<500; n++)
    {
      int i=CGAL::default_random.get_int(0,1000);
      if((i%2)==1) ++real_number_of_odds;
      l.push_back(i);
    }

  while(!l.empty())
    {
      ++detected_number_of_odds;
      l.pop_front();
    }

  std::cout << "detected = " << detected_number_of_odds << std::endl
	    << "reallly  = " << real_number_of_odds << std::endl;
  return detected_number_of_odds!=real_number_of_odds;
}
