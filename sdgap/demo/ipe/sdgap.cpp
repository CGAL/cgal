#define CGAL_SDG_VERBOSE
#undef CGAL_SDG_VERBOSE

#ifndef CGAL_SDG_VERBOSE
#define CGAL_SDG_DEBUG(a)
#else
#define CGAL_SDG_DEBUG(a) { a }
#endif

//~ #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/CGAL_Ipelet_base.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_traits_2.h>
//#include <CGAL/Segment_Delaunay_graph_Linf_2/Bisector_Linf.h>
//#include <CGAL/Segment_Delaunay_graph_Linf_site_2.h>


//#include <CGAL/Polychain_2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>

namespace CGAL_sdgap {

  typedef CGAL::Cartesian<double>                             Kernel;
  typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<Kernel>  Gt;
  typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt>             SDG2;

  const unsigned int num_entries = 2;

  const std::string sublabel[] = {
    "Segment VD Linf axis-parallel",
    "Help"
  };

  const std::string helpmsg[] = {
    "Draw the L_inf Voronoi diagram of axis-parallel segments in Linf",
  };

  class sdgapIpelet
    : public CGAL::Ipelet_base<Kernel,num_entries> {
      public:
        sdgapIpelet()
          :CGAL::Ipelet_base<Kernel,num_entries>
             ("SVDLinfAP",sublabel,helpmsg){}
        void protected_run(int);

    };
  // --------------------------------------------------------------------

  void sdgapIpelet::protected_run(int fn)
  {
    SDG2 svd;     //Voronoi for segments

    if (fn == (num_entries-1)) {
      show_help();
      return;
    }

    std::list<Point_2> pt_list;
    std::list<Segment_2> sg_list;

    Iso_rectangle_2 bbox;

    // grab input

    bbox =
      read_active_objects(
          CGAL::dispatch_or_drop_output
          <Point_2,Polygon_2,Segment_2>(
          std::back_inserter(pt_list),
          segment_grabber(std::back_inserter(sg_list)),
          std::back_inserter(sg_list)
          )
          );

    // check input

    if (pt_list.empty() and sg_list.empty()) {
      print_error_message(("Nothing selected"));
      return;
    }

    // check that segments are all axis-parallel
    for (std::list<Segment_2>::iterator
         sit  = sg_list.begin();
         sit != sg_list.end();
         ++sit)
    {
      if (not (sit->is_horizontal() or sit->is_vertical())) {
        print_error_message(("Non axis-parallel segment"));
        return;
      }
    }


    Kernel::FT incr_len = 75;
    // slightly increase the size of the bbox
    bbox = Iso_rectangle_2(
      bbox.min()+Kernel::Vector_2(-incr_len,-incr_len),
      bbox.max()+Kernel::Vector_2(incr_len,incr_len));

    // insert input into svd

    for (std::list<Segment_2>::iterator
         sit  = sg_list.begin();
         sit != sg_list.end();
         ++sit)
    {
      svd.insert(sit->point(0),sit->point(1));
    }

    svd.insert(pt_list.begin(),pt_list.end());

    draw_dual_in_ipe(svd, bbox);

  } // end of void sdgapIpelet::protected_run(int fn)

}

CGAL_IPELET(CGAL_sdgap::sdgapIpelet)
