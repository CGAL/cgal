// example using nearest_neighbour_iterator for weighted Minkowski Distance

#include <CGAL/compiler_config.h>


#include <vector>
#include <numeric>
#include <cassert>

#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_d.h>
#include <CGAL/Binary_search_tree.h>
#include <CGAL/Weighted_Minkowski_distance.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Nearest_neighbour_general_distance.h>
#include <CGAL/Search_nearest_neighbour.h>
#include <CGAL/Random.h>

#pragma hdrstop

#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>

//example illustrating next nearest neigbbour searching
//using std::copy and illustrating distance browsing

void example2_using_nearest_neighbour_general_distance() {

typedef CGAL::Cartesian<double> R;
typedef CGAL::Point_d<R> Point;
typedef Point::FT NT;
typedef CGAL::Plane_separator<NT> Separator;
typedef CGAL::Kd_tree_traits_point<Separator,Point> Tree_traits;
typedef CGAL::Weighted_Minkowski_distance<Tree_traits> Distance;
typedef CGAL::Nearest_neighbour_general_distance<Tree_traits,
CGAL::Search_nearest_neighbour, Distance>::iterator NNN_Iterator;

int dim=12;
int point_number=5;
int nearest_neighbour_number=5;
int bucket_size=1;
double eps=0.1;
double the_power=2.0;
std::vector<double> my_weights(dim);
for (int wi=0; wi<dim; wi++) my_weights[wi]=1.0;

// Build list containing 5 random points in a 12-dim unit square
typedef std::list<Point> listd;
listd lpt;
CGAL::Random Rnd;
for (int i1=0; i1<point_number; i1++) {
        std::vector<double> vec(dim);
        for (int j=0; j<dim; j++) vec[j]=Rnd.get_double(-1.0,1.0);
        Point Random_point(dim, vec.begin(), vec.end());
        lpt.push_front(Random_point);
}

Tree_traits t(bucket_size, CGAL::SLIDING_FAIR, 2.0, true);

typedef CGAL::Binary_search_tree<Tree_traits> Tree;
Tree d(lpt.begin(), lpt.end(), t);
std::cout << "tree statistics" << std::endl;
d.statistics();

// query point q
std::vector<double> v(dim);
for (int i2=0; i2<dim; i2++) v[i2]=1.0;
v[0]=2.0;
Point q(dim,v.begin(),v.end());
std::cout << "query point is " << q << std::endl;

Distance tr(the_power,dim,my_weights);


// illustrate use of container with copy

// did not work for Visual C++


std::cout << "started test std::copy" << std::endl;

std::vector<Tree_traits::Item_with_distance> result1(point_number);
CGAL::Nearest_neighbour_general_distance<Tree_traits,
CGAL::Search_nearest_neighbour, Distance> NNN1(d,q,tr,eps);
std::copy(NNN1.begin(),NNN1.end(),result1.begin());

for (int i3=0; i3 < nearest_neighbour_number; i3++) {
        std::cout << "result[" << i3 << "]= " << *(result1[i3].first)
                  << std::endl;
        std::cout << "with distance " << result1[i3].second << std::endl;
}



// example of browsing


CGAL::Nearest_neighbour_general_distance<Tree_traits,
CGAL::Search_nearest_neighbour, Distance> NNN2(d,q,tr,eps);

// define predicate class
class GreaterThan0 {

        public:

                bool operator() (NNN_Iterator::value_type const result) const {

                return ( (*(result.first))[0] > 0.0);

        }

};

GreaterThan0 pred;

NNN_Iterator first= NNN2.begin();
NNN_Iterator last=NNN2.end();

while (first != last && !pred(*first)) ++first;

if (last != first)  {
        std::cout << "first positive neighbour is " << (*(*first).first)
        << std::endl;
  }
else  {
        std::cout << "no positive neighbour found" << std::endl;
};



}

int main() {
  example2_using_nearest_neighbour_general_distance();
  return 0;
}


