#include <CGAL/compiler_config.h>


#include <vector>
#include <numeric>
#include <cassert>

#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_d.h>
#include <CGAL/Binary_search_tree.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Nearest_neighbour_L2.h>
#include <CGAL/Random.h>
#include <CGAL/Search_nearest_neighbour.h>

#pragma hdrstop

#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>


// example illustrating the computation of a fixed number of nearest neighbours using copy_n

// as well as the computation of a fixed number of next nearest neighbours using an iterator
void example0_using_nearest_neighbour_L2() {

typedef CGAL::Cartesian<double> R;
typedef CGAL::Point_d<R> Point;
typedef Point::FT NT;
typedef CGAL::Plane_separator<NT> Separator;
typedef CGAL::Kd_tree_traits_point<Separator,Point> Traits;
typedef CGAL::Nearest_neighbour_L2<Traits,CGAL::Search_nearest_neighbour>::iterator NNN_Iterator;

// define constants
int dim=2;
int point_number=1000;
int nearest_neighbour_number=4;
int bucket_size=10;
double eps=0.1;


// create a list of random points in the unit square
typedef std::list<Point> listd;

listd lpt;
CGAL::Random Rnd;

for (int i1=0; i1<point_number; i1++) {
        std::vector<double> vec(dim);
        for (int j=0; j<dim; j++) vec[j]=Rnd.get_double(-1.0,1.0);
        Point Random_point(dim, vec.begin(), vec.end());
        lpt.push_front(Random_point);
}


Traits t(bucket_size, CGAL::SLIDING_FAIR, 4.0, true);



// store the points from the list in a binary search tree.
typedef CGAL::Binary_search_tree<Traits> Tree;
Tree d(lpt.begin(), lpt.end(), t);


// define query point
std::vector<double> v(dim);
for (int i2=0; i2<dim; i2++) v[i2]=1.0;
v[0]=2.0;
Point q(dim,v.begin(),v.end());
std::cout << "query point is " << q << std::endl;



// define iterator for browsing the binary search tree d.
NNN_Iterator NNN_Iterator1(d,q,eps);


// define vector to store the results
std::vector<Traits::Item_with_distance> result1(nearest_neighbour_number);

// compute the nearest neighbours
std::vector<Traits::Item_with_distance>::iterator it = result1.begin();
CGAL::copy_n(NNN_Iterator1, nearest_neighbour_number, it);

// compute exact distances
double distance_to_q;
std::vector<double> distance_array(point_number);

int i3=0;

for (listd::iterator pli=lpt.begin(); pli != lpt.end(); pli++) {

        distance_to_q=0.0;

        for (int i=0; i<dim; i++)

        if (sqrt(fabs(q[i]-(*pli)[i])) > distance_to_q )

        distance_to_q = sqrt(fabs(q[i]-(*pli)[i]));

        distance_array[i3]=distance_to_q;

        i3++;

};


// compare exact and approximate distances

std::partial_sort(&distance_array[0],

                    &distance_array[2*nearest_neighbour_number],

                    &distance_array[point_number]);


std::cout <<
"comparison of approximate nearest neighbour distances and real distances"
<< std::endl;


std::cout << std::endl;

for (int i4=0; i4<nearest_neighbour_number; i4++) {
        std::cout << "dist[" << i4 << "]= " << sqrt(result1[i4].second) << std::endl;
        std::cout << "distance_array[" << i4 << "]= " <<
        distance_array[i4] << std::endl;
        std::cout << "result[" << i4 << "]= "

        << *(result1[i4].first) << std::endl;

        std::cout << "approximation factor= "

        << sqrt(result1[i4].second)/distance_array[i4] - 1.0;

        std::cout << std::endl;
};

// compute the next k nearest neighbours using iterators
std::cout << "testing iterator started computing next "
<< nearest_neighbour_number << " nearest neighbours" << std::endl;
std::vector<Traits::Item_with_distance> result2(nearest_neighbour_number);
std::cout << std::endl;
for (int i5=0; i5 < nearest_neighbour_number; i5++) {
        result2[i5] = *NNN_Iterator1; ++NNN_Iterator1;
        std::cout << "dist[" << i5 << "]= " << sqrt(result2[i5].second)
        << std::endl;
        std::cout << "distance_array[" << nearest_neighbour_number+i5 << "]= "
        << distance_array[nearest_neighbour_number+i5] << std::endl;
        std::cout << "result[" << i5 << "]= "

        << *(result2[i5].first) << std::endl;

        std::cout << "approximation factor = "

        << sqrt(result1[i5].second)/distance_array[i5] - 1.0;

        std::cout << std::endl;
};


};


int main() {
  example0_using_nearest_neighbour_L2();
  return 0;
};



