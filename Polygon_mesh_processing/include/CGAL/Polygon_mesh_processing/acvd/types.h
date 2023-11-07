
// CGAL
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Aff_transformation_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Point_set_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;


typedef Kernel::FT                  FT;
typedef Kernel::Point_3             Point;
typedef Kernel::Vector_3            Vector;
typedef Kernel::Triangle_3          Triangle;
typedef Kernel::Plane_3             Plane;

typedef Kernel::Point_2             Point_2;
typedef Kernel::Segment_2           Segment_2;

typedef CGAL::First_of_pair_property_map<std::pair<Point, Vector>>                     Point_map;
typedef CGAL::Second_of_pair_property_map<std::pair<Point, Vector>>                    Normal_map;
typedef CGAL::Point_set_3< Point, Vector > Pointset;

// triangle fit
typedef CGAL::Polyhedron_3<Kernel>      Polyhedron;
typedef std::vector<Point>                  PointList;
typedef CGAL::Bbox_3                    Bbox;


/// @brief

typedef std::vector<float>                  DataList;
typedef std::vector<int>                    IntList;
typedef std::vector<double>                 DoubleList;
typedef std::vector<bool>                   BoolList;
typedef std::vector<Vector>                 VecList;

typedef std::set<int>                       IntSet;

typedef std::pair<int, int>                 IntPair;


#ifndef TYPE_H
#define TYPE_H

namespace qem
{
    enum class INIT_QEM_GENERATOR {RANDOM,KMEANS_FARTHEST,KMEANS_PLUSPLUS};
    enum class VERBOSE_LEVEL {LOW,MEDIUM,HIGH};
}

#endif /* TYPE_H */