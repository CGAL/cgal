#ifndef BOPS_LINEAR_H
#define BOPS_LINEAR_H

#include <list>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                            Point;
typedef CGAL::Polygon_2<Kernel>                    Polygon;
typedef CGAL::Polygon_with_holes_2<Kernel>         Polygon_with_holes;
typedef std::list<Polygon_with_holes>              Pgn_with_holes_container;

#endif
