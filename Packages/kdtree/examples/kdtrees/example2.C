/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 * example2.C -
 *    Simple example the CGAL KD-tree module.
 *
 * Written by Sariel Har-Peled 
 *            Iddo Hanniel
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/

//might be needed when compiling with g++-2.7.2 
//#ifdef  __GNUG__ 
//#include  <typeinfo>
//#endif  /* __GNUG__ */ 
#include <CGAL/config.h>
#include <iostream>
#include <ctime>
#include <cassert>

#include  <CGAL/Cartesian.h>
#include  <CGAL/Point_3.h>

#include <list>

#include  <CGAL/kdtree_d.h>

typedef CGAL::Point_3<CGAL::Cartesian<double> >  point;
typedef CGAL::Kdtree_interface_3d<point>  kd_interface;
typedef CGAL::Kdtree_d<kd_interface>  kd_tree;
typedef kd_tree::Box  box;
typedef std::list<point>  points_list; 





//RANDOM FUNCTIONS
// dblRand - a random number between 0..1 
#ifndef  RAND_MAX
#define  RAND_MAX    0x7fffffff
#endif  // RAND_MAX

static inline  double    dblRand( void )
{
    return  (double)rand() / (double)RAND_MAX;
}


void random_points( int  num, points_list &l )
{
  double  x,y,z;
  
  for  (int j = 0;  j < num; j++)
    {
      x = dblRand()*10 ;
      y = dblRand()*10 ;
      z = dblRand()*10 ;
      point p(x,y,z);
      l.push_front(p);
    }
  
}


int   main()
{

  CGAL::Kdtree_d<kd_interface>  tree(3);
  
  srand( (unsigned)time(NULL) );
  
  std::cout << "Choosing randomly 30 points in the cube (0,0,0)-(10,10,10)\n" ;
  
  points_list  l , res;
  random_points( 30, l);
  
  
  std::cout << "Listing of random points:\n" ;
  std::copy (l.begin(),l.end(),std::ostream_iterator<point>(std::cout,"\n") );
  std::cout << std::endl;
  
  // building the tree for the random points
  tree.build( l );
    
  //checking validity
  if  ( ! tree.is_valid() )
    tree.dump();
  assert( tree.is_valid() );
  
  
  //searching the box r
  box r(point(2,2,2), point(7,7,7) ,3);
  tree.search( std::back_inserter( res ), r );
  
  std::cout << "Listing of the points in the box (2,2,2)-(7,7,7) : \n" ;
  std::copy (res.begin(),res.end(),std::ostream_iterator<point>(std::cout,"\n") );
  std::cout << std::endl;
  
  tree.delete_all();
  
  
  
  return  0;
}








