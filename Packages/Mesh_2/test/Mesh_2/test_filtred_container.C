#include <CGAL/basic.h>
#include <CGAL/Mesh_2/Filtered_queue_container.h>
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

typedef CGAL::Mesh_2::Filtered_queue_container<int, Is_odd> Queue;

int main(int, char**)
{
  Queue q;

  int
    real_number_of_odds=0,
    detected_number_of_odds=0;

  for(int n=0; n<500; n++)
    {
      int i=CGAL::default_random.get_int(0,1000);
      if((i%2)==1) ++real_number_of_odds;
      q.add_element(i);
    }

  while(!q.empty())
    {
      ++detected_number_of_odds;
      q.remove_next_element();
    }

  std::cout << "detected = " << detected_number_of_odds << std::endl
	    << "really  = " << real_number_of_odds << std::endl;
  return detected_number_of_odds!=real_number_of_odds;
}
