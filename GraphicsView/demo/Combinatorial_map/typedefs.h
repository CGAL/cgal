#ifndef TYPEDEFS_H
#define TYPEDEFS_H


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Combinatorial_map_3_with_embedding.h>
#include <CGAL/Combinatorial_map_with_embedding_constructors.h>
#include <CGAL/Combinatorial_map_with_embedding_operations.h>

#include <CGAL/Timer.h>

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>



typedef double coord_type;


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef CGAL::Combinatorial_map_3_with_embedding<Kernel> Map;
typedef Map::Dart_handle Dart_handle;
typedef Map::Vertex      Vertex;

typedef Kernel::Point_3    Point_3;
typedef Kernel::Vector_3   Vector_3;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;

typedef CGAL::Timer Timer;

struct Scene {

  Map map;
};



#endif
