#ifndef MY_POINTC2_H
#define MY_POINTC2_H


#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>


class MyPointC2 {

private:
  double vec[2];
  int col;

public:

  MyPointC2()
    : col(0)
  {
    *vec = 0;
    *(vec+1) = 0;
  }


  MyPointC2(const double x, const double y, int c = 0)
    : col(c)
  {
    *vec = x;
    *(vec+1) = y;
  }

  const double& x() const  { return *vec; }

  const double& y() const { return *(vec+1); }

  double & x() { return *vec; }

  double& y() { return *(vec+1); }

  int color() const { return col; }

  int& color() { return col; }


  bool operator==(const MyPointC2 &p) const
  {
    return ( *vec == *(p.vec) )  && ( *(vec+1) == *(p.vec + 1) && ( col == p.col) );
  }

  bool operator!=(const MyPointC2 &p) const
  {
      return !(*this == p);
  }

};

#endif // MY_POINTC2_H
