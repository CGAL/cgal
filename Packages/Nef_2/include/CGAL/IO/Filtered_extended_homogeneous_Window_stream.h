#ifndef FILTERED_EXTENDED_HOMOGENEOUS_WINDOW_STREAM_H
#define FILTERED_EXTENDED_HOMOGENEOUS_WINDOW_STREAM_H

#ifdef CGAL_USE_LEDA
#include <CGAL/Filtered_extended_homogeneous.h>
#include <CGAL/IO/Window_stream.h>

CGAL_BEGIN_NAMESPACE

template <class RT> 
CGAL::Window_stream& operator<<(CGAL::Window_stream& w, 
                                const Extended_point<RT>& p)
{ w.draw_filled_node(CGAL::to_double(p.x()),CGAL::to_double(p.y())); 
  return w;
}

template <class RT> 
CGAL::Window_stream& operator<<(CGAL::Window_stream& w, 
                                const Extended_segment<RT>& s)
{ w.draw_segment(CGAL::to_double(s.source().x()),
                 CGAL::to_double(s.source().y()),
                 CGAL::to_double(s.target().x()),
                 CGAL::to_double(s.target().y())); 
  return w;
}

template <class RT> 
leda_point pnt(const Extended_point<RT>& p)
{ return leda_point(CGAL::to_double(p.x()),CGAL::to_double(p.y())); }

CGAL_END_NAMESPACE
#endif // CGAL_USE_LEDA
#endif // FILTERED_EXTENDED_HOMOGENEOUS_WINDOW_STREAM_H

