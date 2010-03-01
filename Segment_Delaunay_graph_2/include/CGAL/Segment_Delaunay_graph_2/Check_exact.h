#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_CHECK_EXACT
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_CHECK_EXACT 1

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

CGAL_SEGMENT_DELAUNAY_GRAPH_2_BEGIN_NAMESPACE

#ifdef CGAL_PROFILE

template<typename T>
struct Check_exact {
  bool operator()() const { return true; }
};

template<>
struct Check_exact<double> {
  bool operator()() const { return false; }
};

template<>
struct Check_exact<float> {
  bool operator()() const { return false; }
};

template<bool b>
struct Check_exact< Interval_nt<b> > {
  bool operator()() const { return false; }
};

#endif


CGAL_SEGMENT_DELAUNAY_GRAPH_2_END_NAMESPACE

CGAL_END_NAMESPACE


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_CHECK_EXACT
