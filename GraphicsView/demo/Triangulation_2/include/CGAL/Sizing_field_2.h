#ifndef CGAL_SIZING_FIELD_2_H
#define CGAL_SIZING_FIELD_2_H

#include <list>

#include <CGAL/basic.h>

namespace CGAL
{

template <typename Kernel>
class Sizing_field_2 // pure virtual class
{    
public:
  typedef Point_2<Kernel> Point;
	
public:
  Sizing_field_2()
  {
  }
  virtual ~Sizing_field_2()
  {
  }

  virtual double query(const Point& p) const = 0;
};

}//namespace CGAL

#endif
