#ifndef CGAL_POINT_CONTAINER_H
#define CGAL_POINT_CONTAINER_H

#include <CGAL/basic.h>
#include <list>
#include <vector>

CGAL_BEGIN_NAMESPACE

template<class P, class C = std::list<P> >
class Point_container
{
public:
  typedef C  Container;
  typedef P  Point_2;
  typedef typename Container::iterator   Point_handle;
  typedef typename Container::size_type  size_type;

private:
  typedef Point_container<Point_2,Container> Self;

public:
  Point_container() {}

  Point_handle insert(const Point_2& p)
  {
    c.push_back(p);
    return --c.end();
  }

  void remove(Point_handle handle)
  {
    c.erase(handle);
  }

  void swap(const Self& other)
  {
    c.swap(other.c);
  }

  void clear() {
    c.clear();
  }

  size_type size() const { return c.size(); }

private:
  Container c;
};

#if 1
template<class P, unsigned long S>
class Array_point_container
{
public:
  //  typedef C  Container;
  typedef P  Point_2;
  enum { Size = S };
  typedef Point_2*      Point_handle;
  typedef unsigned long size_type;

  Array_point_container() {
    last = 0;
    c = new Point_2[Size];
  }

  Point_handle insert(const Point_2& p)
  {
    CGAL_precondition( last < Size );

    c[last] = p;
    Point_handle h = &c[last];
    last++;
    return h;
  }

  std::size_t size() const { return c.size(); }

  void clear() {
    last = 0;
  }

private:
  unsigned long last;
  Point_2* c;
};
#endif

CGAL_END_NAMESPACE

#endif // CGAL_POINT_CONTAINER_H
