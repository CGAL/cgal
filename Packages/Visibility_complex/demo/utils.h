#ifndef VISIBILITY_COMPLEX_ARR_UTILS_H
#define VISIBILITY_COMPLEX_ARR_UTILS_H

#include "types.h"

const Line_2 dual(const Point_2& p)
{
  // duality p(a,b) -> l: y=ax-b -> l: -ax+y+b=0
  return Line_2(-p.x(), 1, p.y()); // WARNING: kernel dependant!
}

const Point_2 dual(const Line_2& l)
{
  // duality l: ax+by+c=0 -> l: y = (-a/b)x - (c/b) -> p(-a/b,c/b)
  return Point_2( -l.a()/l.b(), l.c()/l.b() );
}

bool dual(const Vertex_handle vh, Point_2& p)
{
  const Line_2& line =
    vh->supporting_line();
  if(line.b() == 0) return false;
  p = dual(line);
  return true;
}

void display(CGAL::Qt_widget* widget, Vertex_handle va, Vertex_handle vb)
{
  Point_2 a, b;
  if( dual(va, a) && dual(vb, b) )
    *widget << Segment_2(a, b);
}

void print(Vertex_handle va, Vertex_handle vb)
{
  Point_2 a, b;
  if( dual(va, a) && dual(vb, b) )
    std::cerr << Segment_2(a, b) << std::endl;
}

#endif
