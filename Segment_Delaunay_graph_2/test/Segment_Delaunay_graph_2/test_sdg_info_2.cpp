// file: test_sdg_info_2.cpp
#include <CGAL/basic.h>

// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// choose the kernel
#include <CGAL/Simple_cartesian.h>

struct Rep : public CGAL::Simple_cartesian<double> {};

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h>

typedef
CGAL::Segment_Delaunay_graph_filtered_traits_2<Rep>
Traits_x;

struct Gt : public Traits_x {};

#include "Multi_info.h"

#include "test_info.h"

typedef CGAL::Segment_Delaunay_graph_storage_traits_2<Gt>   ST_base;
typedef CGAL::Segment_Delaunay_graph_storage_traits_with_info_2
<ST_base,
 ::Multi_info<int>,
 ::Multi_info_convert_info<int>,
 ::Multi_info_merge_info<int> >
ST;

typedef CGAL::Segment_Delaunay_graph_2<Gt,ST>  SDG2;

#include "sdg_info_main.h"

// Comment for cgal_create_makefile : main()
