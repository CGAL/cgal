/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 * example3.C -
 *    Simple example the CGAL KD-tree module.
 *    Example with user defined point_d.    
 *
 * Written by Sariel Har-Peled 
 *            Iddo Hanniel
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/

#include <CGAL/Cartesian.h>

#include <iostream>
#include <iterator>
#include <ctime>
#include <cassert>
#include <list>

#include <CGAL/kdtree_d.h>

template <int  DIM>
class Point_float_d 
{
private:
  double   vec[ DIM ];
  
public:
  Point_float_d()
  {
    for  ( int ind = 0; ind < DIM; ind++ )
      vec[ ind ] = 0;
  }

  int dimension() const
  {
    return  DIM;
  }
  
//not essential by specification but needed for initializing a general d-point
  void set_coord(int k, double x)
  {
    assert( 0 <= k  &&  k < DIM );
    vec[ k ] = x;
  }
  
  double  & operator[](int k)  
  {
    assert( 0 <= k  &&  k < DIM );
    return  vec[ k ];
  }

  double  operator[](int k) const
  {
    assert( 0 <= k  &&  k < DIM );
    return  vec[ k ];
  }
};

// not essential by specification but nice to have
template <int DIM>
std::ostream &operator<<(std::ostream &os, const Point_float_d<DIM> &p)
{
  std::cout << "(";
  for(int i = 0; i < DIM; i++)
    {
      std::cout << p[i] ;
      if (i < p.dimension() - 1) std::cout << ", ";
    }
  std::cout << ")";
  return os;
}

typedef Point_float_d<4>  point;
typedef CGAL::Kdtree_interface<point>  kd_interface;
typedef CGAL::Kdtree_d<kd_interface>  kd_tree;
typedef kd_tree::Box  box;
typedef std::list<point>  points_list; 

//RANDOM FUNCTIONS
// dblRand - a random number between 0..1 
#ifndef  RAND_MAX
#define  RAND_MAX    0x7fffffff
#endif

inline double dblRand( void )
{
    return  (double)rand() / (double)RAND_MAX;
}

void random_points( int  num, points_list &l, int DIM)
{
  double  x;
  
  for  (int j = 0;  j < num; j++)
    {
      point p;
      for (int i=0; i<DIM; i++)
        {
          x = dblRand()*10 ;
          p.set_coord(i,x);
        }
      l.push_front(p);
    }
}

int main()
{
  CGAL::Kdtree_d<kd_interface>  tree(3);
  
  srand( (unsigned)time(NULL) );
  
  std::cout << "Choosing randomly 30 points in the cube (0,0,0)-(10,10,10)\n" ;
  
  points_list  l , res;
  random_points( 30, l , 4);

  std::cout << "Listing of random points:\n" ;
  std::copy (l.begin(),l.end(),std::ostream_iterator<point>(std::cout,"\n") );
  std::cout << std::endl;
  
  // Building the tree for the random points
  tree.build( l );
  
  // Checking validity
  if  ( ! tree.is_valid() )
    tree.dump();
  assert( tree.is_valid() );
  
  // Searching the box r
  point p,q;
  for (int k=0;k<4;k++)
    {
      p.set_coord(k,2);
      q.set_coord(k,8);
    }
  
  box r(p, q, 4);
  tree.search( std::back_inserter( res ), r );
  
  std::cout << "Listing of the points in the box (2,2,2,2)-(8,8,8,8) : \n" ;
  std::copy (res.begin(),res.end(),
	     std::ostream_iterator<point>(std::cout,"\n") );
  std::cout << std::endl;
  
  tree.delete_all();

  return 0;
}
