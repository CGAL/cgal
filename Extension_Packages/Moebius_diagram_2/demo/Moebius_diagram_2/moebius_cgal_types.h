#ifndef MOEBIUS_CGAL_TYPES_H
#define MOEBIUS_CGAL_TYPES_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>


#ifdef CHR_FILTERED_EXACT

#include <CGAL/Filtered_exact.h>
#include <CGAL/MP_Float.h>
typedef CGAL::Filtered_exact<double, CGAL::MP_Float> NT;

#else
#  ifdef CHR_LEDA

#include <CGAL/leda_real.h>
typedef leda_real NT;

#else

// double
typedef double NT;

#endif
#endif

#include <CGAL/Moebius_diagram_2.h>
#include <CGAL/Moebius_diagram_euclidean_traits_2.h>

typedef CGAL::Cartesian<NT> K1;
struct K : public K1 {};
typedef K::Point_2 Point;
typedef K::Conic_2 Conic;

typedef double W;
typedef CGAL::Moebius_diagram_euclidean_traits_2<K,W> Traits;

typedef CGAL::Moebius_diagram_2<Traits> MD1;
struct MD: public MD1 {};
typedef MD::Point WPoint;
typedef MD::Vertex_iterator Vertex_iterator;
typedef MD::Finite_edges_iterator Finite_edges_iterator;


#endif
