// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// choose the kernel
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> K;

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h>

typedef CGAL::Segment_Delaunay_graph_filtered_traits_2<K> Gt;


#include "Multi_info.h"

#include "test_info_hierarchy.h"

typedef Multi_info<int>               Info;
typedef Multi_info_convert_info<int>  Convert_info;
typedef Multi_info_merge_info<int>    Merge_info;

typedef CGAL::Segment_Delaunay_graph_storage_traits_with_info_2<Gt,
                                                                Info,
                                                                Convert_info,
                                                                Merge_info>
ST;

typedef CGAL::Tag_false STag;

typedef CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,ST,STag>  SDG2;

#include "sdg_hierarchy_info_main.h"

// Comment for cgal_create_makefile : main()
