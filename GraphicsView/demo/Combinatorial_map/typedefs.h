#ifndef TYPEDEFS_H
#define TYPEDEFS_H


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Embedded_combinatorial_map_3.h>
#include <CGAL/Embedded_combinatorial_map_constructors.h>
#include <CGAL/Embedded_combinatorial_map_operations.h>

#include <CGAL/Timer.h>

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>



typedef double coord_type;


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef CGAL::Embedded_combinatorial_map_3<Kernel> Map;
typedef Map::Dart_handle Dart_handle;
typedef Map::Vertex      Vertex;

typedef Map::Point    Point_3;
typedef Map::Vector   Vector_3;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;

typedef CGAL::Timer Timer;

struct Scene {

  Map map;
};



#endif
