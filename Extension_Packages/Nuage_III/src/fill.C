#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/algorithm.h>
#include <vector>
#include <set>
#include <algorithm>


typedef double NT;

struct K : public CGAL::Filtered_kernel<CGAL::Simple_cartesian<NT> > {};
typedef K::Point_3  Point;
typedef K::Vector_3 Vector;
typedef K::Segment_3  Segment;
typedef K::Triangle_3  Triangle;

NT
weight(const Point& p, const Point& q, const Point& r){
  NT area = std::sqrt(Triangle(p,q,r).squared_area());
  NT l1 = std::sqrt((p-q) * (p-q));
  NT l2 = std::sqrt((p-r) * (p-r));
  NT l3 = std::sqrt((q-r) * (q-r));
  if(l1>l2) std::swap(l1,l2);
  if(l2>l3) std::swap(l2,l3);
  if(l1>l2) std::swap(l1,l2);
  if(l2>l3) std::swap(l2,l3);
  
  // Taken from Piecewise-Linear Interpolation between Polygonal Slices
  // from Gill Barequet and Micha Sharir
  return 0.85 * area + 0.05 * (l1 +l2 + l3) + 0.1 * l3/l1; 
}


void
insert(std::set<CGAL::Triple<int,int,int> >& triangles, int i, int j, int k){
  std::cout << i << ", " << j << ", " << k << std::endl;
  if(i>j) std::swap(i,j);
  if(j>k) std::swap(j,k);
  if(i>j) std::swap(i,j);
  if(j>k) std::swap(j,k);
  std::cout << i << ", " << j << ", " << k << std::endl;
  triangles.insert(CGAL::make_triple(i,j,k));
}
  

void
collect(int i, int k, int n, const std::vector<NT>& O, std::set<CGAL::Triple<int,int,int> >& triangles){

  std::cout << "collect(" << i << ", " << k << ")" << std::endl;
  if((i+2) == k){
    insert(triangles, i, i+1, k);
  }else {
    int o = O[i*n+k];

    if(o != (i+1)){
      collect(i, o, n, O, triangles);
    }
    insert(triangles, i, o, k);
    if(o != (k-1)){
      collect(o, k, n, O, triangles);
    }
  }
}


 
int
main(){

  int n;
  std::cin >> n;

  std::vector<Point> points(n);
  CGAL::copy_n(std::istream_iterator<Point>(std::cin), n, points.begin());


  std::set<CGAL::Triple<int,int,int> > triangles;

  std::vector<NT> W(n*n);
  std::vector<NT> O(n*n);
  
  for(int i = 0; i <= n-2; i++){
    W[i*n + i + 1] = 0;
  }
  for(int i = 0; i <= n-3; i++){
    W[i*n + i + 2] = weight(points[i], points[i+1], points[i+3]);
  }
   
  for(int j = 3; j <= n-1; j++){
    for(int i=0; i <= n - j - 1; i++){
      int k = i + j;
      double lmin = -1;
      int lmin_index;
      for(int m = i+1; m < k; m++){
	double d = W[i*n + m] + W[m*n + k] + weight(points[i], points[m], points[k]);
	if( (lmin == -1) || (d < lmin )){
	  lmin = d;
	  lmin_index = m;
	}
      }
      W[i*n + k] = lmin;
      O[i*n + k] = lmin_index;

    }
  }

  collect(0, n-1, n, O, triangles);
  
  std::cout << "Shape {\n"
    "appearance Appearance {\n"
    "material Material { diffuseColor .9 .5 .1}}\n"
    "geometry\n"
    "IndexedFaceSet {\n"
    "coord DEF def_coords Coordinate {\n"
    "point [  \n" ;
  
  for (int i = 0; i < n; i++){
    std::cout << points[i].x() << " " << points[i].y() << " " << points[i].z() << ",\n "; 
  }
  std::cout << "]\n"
    "}\n"
    "solid FALSE\n"
    "coordIndex [ \n";
  
  for(std::set<CGAL::Triple<int,int,int> >::iterator it = triangles.begin();
      it != triangles.end();
      it++){
    std::cout << it->first << ", " << it->second << ", " << it->third  << ", -1,\n ";
  }
  std::cout << "]\n"
    "}# IndexedFaceSet\n"
    "}# Shape\n";
  
  return 1;
}


