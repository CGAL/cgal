
#ifndef CGAL_SDG_VERBOSE
#define CGAL_SDG_DEBUG(a)
#else
#define CGAL_SDG_DEBUG(a) { a }
#endif

// standard includes
#include <iostream>
#include <fstream>
#include <cassert>
// choose the kernel
#include <CGAL/Simple_cartesian.h>

struct Rep : public CGAL::Simple_cartesian<double> {};

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h>

typedef
CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<Rep>
Traits_x;

typedef Traits_x Gt;

#include "Multi_info.h"

#include "test_info.h"

typedef Multi_info<int>               Info;
typedef Multi_info_convert_info<int>  Convert_info;
typedef Multi_info_merge_info<int>    Merge_info;

typedef CGAL::Segment_Delaunay_graph_storage_traits_with_info_2<Gt,
                                                                Info,
                                                                Convert_info,
                                                                Merge_info>
ST;

typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt,ST>  SDG2;

#include "sdg_info_main.h"

// Comment for cgal_create_makefile : main()
