#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/point_generators_2.h>

typedef CGAL::Exact_circular_kernel_2             Circular_k;

typedef CGAL::Point_2<Circular_k>                 Point_2;
typedef CGAL::Circular_arc_2<Circular_k>          Circular_arc_2;

int main()
{
	int n = 0;
	Circular_arc_2 c = Circular_arc_2(Point_2(10,0), Point_2(5,5), Point_2(0, 0));
	for(int i = 1; i <= 5; i++) {
		for(int j = 1; j <= 10; j++) {
			Point_2 p = Point_2(i, j);
			if(Circular_k().has_on_2_object()(c,p)) n++;
		}
	}
	std::cout << "There are " << n << " points in the 10x5 grid on the circular" 
	<< std::endl
	<< " arc with endpoints (0,0), (10,0) defined on the superior half plane" 
	<< std::endl;
  return 0;
};
