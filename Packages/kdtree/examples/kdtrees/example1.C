/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 * example1.C -
 *    Simple example the CGAL KD-tree module.
 *
 * Written by Sariel Har-Peled 
 *            Iddo Hanniel
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/
 
//might be needed when compiling with g++-2.7.2
//#ifdef  __GNUG__ 
//#include  <typeinfo>
//#endif  /* __GNUG__ */ 
#include  <CGAL/config.h>

#include <iostream>
#include <ctime>
#include <cassert>

#include  <CGAL/Cartesian.h>
#include  <CGAL/Point_2.h>

#include <list>

#include  <CGAL/kdtree_d.h>

typedef CGAL::Point_2<CGAL::Cartesian<double> >  point;
typedef CGAL::Kdtree_interface_2d<point>  kd_interface;
typedef CGAL::Kdtree_d<kd_interface>  kd_tree;
typedef kd_tree::Box  box;
typedef std::list<point>  points_list; 

int   main()
{

    CGAL::Kdtree_d<kd_interface>  tree(2);
    points_list  l , res;

    srand( (unsigned)time(NULL) );

    std::cout << "Insering evenly 81 points  in the square (0,0)-(10,10) ...\n\n" ;
    for (int i=1; i<10; i++)
      for (int j=1; j<10; j++)
        {
          point p(i,j);
          l.push_front(p);
        }

    // building the tree for the random points
    tree.build( l );
       
    //checking validity
    if  ( ! tree.is_valid() )
        tree.dump();
    assert( tree.is_valid() );

    
    //defining and searching the box r
    double lx,ly,rx,ry;
    std::cout << "Define your query square.\nEnter left x coordinate: " ;
    std::cin >> lx ;
    std::cout << "Enter left y coordinate: ";
    std::cin >> ly;
    std::cout << "Enter right x coordinate: " ;
    std::cin >> rx ;
    std::cout << "Enter right y coordinate: ";
    std::cin >> ry;
    std::cout << std::endl; 

    box r(point(lx,ly), point(rx,ry) ,2);

    tree.search( std::back_inserter( res ), r );
    
    std::cout << "Listing of the points in the square: \n" ;
    std::copy (res.begin(),res.end(),std::ostream_iterator<point>(std::cout," \n") );
    std::cout << std::endl;

    tree.delete_all();

    

    return  0;
}





