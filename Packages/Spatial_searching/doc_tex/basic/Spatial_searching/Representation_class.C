#include <CGAL/Kd_tree_d.h>

template <int  DIM>
class Point_float_d 
{
private:
  double   vec[ DIM ];
  
public:

  // new kd tree requires definition of representation class
  class R
  {
  public:
        typedef double FT;
        typedef double RT;
        typedef Point_float_d Point_d;
  }; 

  Point_float_d()
  {
    for  ( int ind = 0; ind < DIM; ind++ )
      vec[ ind ] = 0;
  }

  int dimension() const
  {
    return  DIM;
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

typedef Point_float_d<4>  point;
typedef CGAL::Kdtree_interface<point>  kd_interface;
typedef CGAL::Kdtree_d<kd_interface>  kd_tree;


