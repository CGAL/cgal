#ifndef CGAL_PM_POINT_UTILITIES_2_H
#define CGAL_PM_POINT_UTILITIES_2_H

CGAL_BEGIN_NAMESPACE

template <class Point_2_>
struct Construct_point_to_right_2
{
  typedef Point_2_                      Point_2;
  typedef typename Point_2::R::FT       FT;

  Point_2 operator()(const Point_2 & pnt) const
  { 
    return Point_2(pnt.x() + FT(1), pnt.y());
  }
};

template <class Point_2_>
struct Construct_point_to_left_2
{
  typedef Point_2_                      Point_2;
  typedef typename Point_2::R::FT       FT;

  Point_2 operator()(const Point_2 & pnt) const
  { 
    return Point_2(pnt.x() - FT(1), pnt.y());
  }
};

template<class Point_2>
bool is_left(const Point_2 & p1, const Point_2 & p2)
{
  typedef typename Point_2::R Kernel;
  return Kernel().less_x_2_object()(p1, p2); 
}


template<class Point_2>
bool is_right(const Point_2 & p1, const Point_2 & p2)
{ 
  typedef typename Point_2::R Kernel;
  return Kernel().less_x_2_object()(p2, p1);
}

template<class Point_2>
bool is_same_x(const Point_2 & p1, const Point_2 & p2)
{ 
  typedef typename Point_2::R Kernel;
  return Kernel().equal_x_object()(p1, p2);
}

template<class Point_2>
bool is_lower(const Point_2 & p1, const Point_2 & p2)
{
  typedef typename Point_2::R Kernel;
  return Kernel().less_y_2_object()(p1, p2);
}

template<class Point_2>
bool is_higher(const Point_2 & p1, const Point_2 & p2)
{ 
  typedef typename Point_2::R Kernel;
  return Kernel().less_y_2_object()(p2, p1);
}

template<class Point_2>
bool is_same_y(const Point_2 & p1, const Point_2 & p2)
{ 
  typedef typename Point_2::R Kernel;
  return Kernel().equal_y_object()(p1, p2);
}

template<class Point_2>
bool is_same(const Point_2 & p1, const Point_2 &p2)
{
  typedef typename Point_2::R Kernel;
  return (compare_x(p1, p2) == EQUAL) && (compare_y(p1, p2) == EQUAL);
}

template<class Point_2>
const Point_2 & leftmost(const Point_2 & p1, const Point_2 & p2)
{ 
  typedef typename Point_2::R Kernel;
  return (is_left(p1, p2) ? p1 : p2); 
}

template<class Point_2>
const Point_2 & rightmost(const Point_2 & p1, const Point_2 & p2)
{ 
  typedef typename Point_2::R Kernel;
  return (is_right(p1, p2) ? p1 : p2); 
}
  
template<class Point_2>
const Point_2 & lowest(const Point_2 &p1, const Point_2 & p2)
{ 
  typedef typename Point_2::R Kernel;
  return (is_lower(p1, p2) ? p1 : p2);
}
  
template<class Point_2>
const Point_2 & highest(const Point_2 &p1, const Point_2 & p2)
{
  typedef typename Point_2::R Kernel;
  return (is_higher(p1, p2) ? p1 : p2);
}

CGAL_END_NAMESPACE

#endif // CGAL_PM_POINT_UTILITIES_2_H
// EOF
