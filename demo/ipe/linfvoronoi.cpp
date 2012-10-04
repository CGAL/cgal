#define CGAL_SDG_VERBOSE
#undef CGAL_SDG_VERBOSE

#ifndef CGAL_SDG_VERBOSE
#define CGAL_SDG_DEBUG(a)
#else
#define CGAL_SDG_DEBUG(a) { a }
#endif

//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/CGAL_Ipelet_base.h> 
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_traits_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Object.h>
#include <CGAL/intersections.h>


namespace CGAL_linfvoronoi{

// choice of kernel
//typedef CGAL::Exact_predicates_exact_constructions_kernel   Kernel;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Cartesian<double>                             Kernel;

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<Kernel> Gt;
typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt>	           SDG2;
typedef Kernel::Segment_2                                  Segment_2;
typedef CGAL::Polychainray_2<Gt>                           pcr;
typedef CGAL::Polychainsegment_2<Gt>                       pcs;
typedef CGAL::Polychainline_2<Gt>                          pcl;
// --------------------------------------------------------------------

const std::string sublabel[] = {
  "LinfVoronoi", "Help"
};

const std::string helpmsg[] = {
  "Draw the Voronoi diagram of a set of points and segments in the L-infinity metric"
};

class linfvoronoiIpelet 
  : public CGAL::Ipelet_base<Kernel,2> {
public:
  linfvoronoiIpelet() 
    :CGAL::Ipelet_base<Kernel,2>("LinfVoronoi",sublabel,helpmsg){}
  void protected_run(int);
  struct Voronoi_from_tri{  //Class using stream to get the voronoi diagram
    std::list<pcl> pcl_list;
    std::list<pcr> pcr_list;
    std::list<pcs> pcs_list;
      
      
    void operator<<(const pcr& p){pcr_list.push_back(p);}
    void operator<<(const pcl& p){pcl_list.push_back(p);}
    void operator<<(const pcs& p){pcs_list.push_back(p);}
  };
    
  template <class T,class output_iterator>
  bool 
  cast_pcs_into_seg(const T& obj,const Iso_rectangle_2& bbox,output_iterator out_it) const{
    std::list<Segment_2> seglist;
    typedef pcs::Vertex_const_iterator VI;
    
    if (obj.size() > 1) { 
      VI source = obj.vertices_begin();
      VI target = source+1;
      for( ; target!=obj.vertices_end(); ++source, ++target)
      {
        seglist.push_back( Segment_2(*source, *target));
      }
    }
    bool ret;
    for(typename std::list<Segment_2>::iterator iteS = seglist.begin();iteS!=seglist.end();){
      typename std::list<Segment_2>::iterator itc=iteS++;
      CGAL::Object obj_cgal = CGAL::intersection(*itc,bbox);
      Segment_2 s;
      ret=CGAL::assign(s, obj_cgal);
      if (ret) *out_it++=s;
    }
    return ret;
  }
    
  //Convert infinite objects into drawable segments
  template<class iterator,class output_iterator>
  void 
  cast_pcs_into_seg(const iterator first,const iterator end,
                   const Iso_rectangle_2& bbox, output_iterator out_it) const
  {
    for (iterator it=first;it!=end;++it)
      cast_pcs_into_seg(*it,bbox,out_it);
  } 
    
  template <class T,class output_iterator>
  bool 
  cast_pcr_into_seg(const T& obj,const Iso_rectangle_2& bbox,output_iterator out_it) const{
    std::list<Segment_2> seglist;
    Ray_2 ray;
    
    typedef pcr::Vertex_const_iterator VI;
   
    VI source = obj.vertices_begin();

    if (obj.size() > 1) { 
      VI target = source+1;
      for( ; target!=obj.vertices_end(); ++source, ++target)
      {
        seglist.push_back( Segment_2(*source, *target));
      }
    }
    
    ray = Ray_2(*source, obj.get_outgoing());
    
    bool ret;
    for(typename std::list<Segment_2>::iterator iteS = seglist.begin();iteS!=seglist.end();){
      typename std::list<Segment_2>::iterator itc=iteS++;
      CGAL::Object obj_cgal = CGAL::intersection(*itc,bbox);
      Segment_2 s;
      ret=CGAL::assign(s, obj_cgal);
      if (ret) *out_it++=s;
    }
    
    CGAL::Object obj_cgal = CGAL::intersection(ray,bbox);
    Segment_2 s;
    ret=CGAL::assign(s, obj_cgal);
    if (ret) *out_it++=s;
    
    return ret;
  }
    
  //Convert infinite objects into drawable segments
  template<class iterator,class output_iterator>
  void 
  cast_pcr_into_seg(const iterator first,const iterator end,
                    const Iso_rectangle_2& bbox, output_iterator out_it) const
  {
    for (iterator it=first;it!=end;++it)
      cast_pcr_into_seg(*it,bbox,out_it);
  }    
    
    template <class T,class output_iterator>
    bool 
    cast_pcl_into_seg(const T& obj,const Iso_rectangle_2& bbox,output_iterator out_it) const{
      std::list<Segment_2> seglist;
      Ray_2 ray;
      
      typedef pcl::Vertex_const_iterator VI;
      
      VI source = obj.vertices_begin();
      //for incoming ray
      ray = Ray_2(*source, obj.get_incoming());
      CGAL::Object obj_cgal = CGAL::intersection(ray,bbox);
      Segment_2 s;
      bool ret;
      ret=CGAL::assign(s, obj_cgal);
      if (ret) *out_it++=s;
      //for segments
      if (obj.size() > 1) { 
        VI target = source+1;
        for( ; target!=obj.vertices_end(); ++source, ++target)
        {
          seglist.push_back( Segment_2(*source, *target));
        }
      }
      
      ray = Ray_2(*source, obj.get_outgoing());
      
      for(typename std::list<Segment_2>::iterator iteS = seglist.begin();iteS!=seglist.end();){
        typename std::list<Segment_2>::iterator itc=iteS++;
        CGAL::Object obj_cgal = CGAL::intersection(*itc,bbox);
        Segment_2 s;
        ret=CGAL::assign(s, obj_cgal);
        if (ret) *out_it++=s;
      }
      //for out going ray
      obj_cgal = CGAL::intersection(ray,bbox);
      ret=CGAL::assign(s, obj_cgal);
      if (ret) *out_it++=s;
      
      return ret;
    }
    
    //Convert infinite objects into drawable segments
    template<class iterator,class output_iterator>
    void 
    cast_pcl_into_seg(const iterator first,const iterator end,
                      const Iso_rectangle_2& bbox, output_iterator out_it) const
    {
      for (iterator it=first;it!=end;++it)
        cast_pcl_into_seg(*it,bbox,out_it);
    }    
    
    
  void 
  draw_dual(Voronoi_from_tri& v_recup,const Iso_rectangle_2& bbox,bool makegrp) const
  {
    std::list<Segment_2> seg_list;
    //std::cout << "inside draw_dual" << std::endl;
    cast_pcr_into_seg(v_recup.pcr_list.begin(),v_recup.pcr_list.end(),bbox,std::back_inserter(seg_list));//cast pcr into segments in bbox
    cast_pcl_into_seg(v_recup.pcl_list.begin(),v_recup.pcl_list.end(),bbox,std::back_inserter(seg_list));//cast pcl into segments in bbox
    cast_pcs_into_seg(v_recup.pcs_list.begin(),v_recup.pcs_list.end(),bbox,std::back_inserter(seg_list));//cast pcs into segments in bbox
    
    //filter degenerate segments
    for(std::list<Segment_2>::iterator iteS = seg_list.begin();iteS!=seg_list.end();){
      std::list<Segment_2>::iterator itc=iteS++;
      if (itc->is_degenerate())
        seg_list.erase(itc);
    }
    //std::cout << "before draw_in_ipe" << std::endl;
    draw_in_ipe(seg_list.begin(), seg_list.end(), makegrp);
  }
    
  template<class Triangulation>
  void 
  draw_dual_in_ipe(Triangulation& T,const Iso_rectangle_2& bbox,bool makegrp=true,bool deselect_all=false) const
  {
    //~ template<class GT,class TDS>
    //~ void draw_dual_in_ipe(const CGAL::Triangulation_2<GT,TDS>& T,const Iso_rectangle_2& bbox) const{
    //std::cout << "draw_dual_in_ipe" << std::endl;
    Voronoi_from_tri v_recup;
    T.draw_dual(v_recup);
    //std::cout << "T.draw_dual passed" << std::endl;
    draw_dual(v_recup,bbox,makegrp);
    if (deselect_all) get_IpePage()->deselectAll();
  }
    
};
// --------------------------------------------------------------------

void linfvoronoiIpelet::protected_run(int fn)
{
  SDG2 svd;     //Voronoi for segments
 
  if (fn==1) {
    show_help();
    return;
  } 
  else { 
    std::list<Point_2> pt_list;
    std::list<Segment_2> sg_list;
    std::list<Circle_2> cir_list;
  
    Iso_rectangle_2 bbox=
    read_active_objects(
      CGAL::dispatch_or_drop_output<Point_2,Polygon_2,Circle_2,Segment_2>(
        std::back_inserter(pt_list),
        segment_grabber(std::back_inserter(sg_list)),
        std::back_inserter(cir_list),
        std::back_inserter(sg_list)
      )
    );
    
    if (pt_list.empty() && sg_list.empty()){
      print_error_message(("No mark, no segment and no polygon selected"));
      return;
    }
    
    for (std::list<Segment_2>::iterator it=sg_list.begin();it!=sg_list.end();++it)
      svd.insert(it->point(0),it->point(1));
std::cout << "insert segments " << std::endl;
    if (!pt_list.empty()) {
      svd.insert(pt_list.begin(),pt_list.end());
      std::cout << "insert pnts " << std::endl;
    }
      
    Kernel::FT incr_len=50;
    //slightly increase the size of the Bbox
    bbox=Iso_rectangle_2(bbox.min()+Kernel::Vector_2(-incr_len,-incr_len),
                         bbox.max()+Kernel::Vector_2(incr_len,incr_len));
  std::cout << "before draw_dual_in_ipe " << std::endl;
    draw_dual_in_ipe(svd,bbox);   
    
   }//end of else
}//end of protected_run
}//end of namespace
CGAL_IPELET(CGAL_linfvoronoi::linfvoronoiIpelet)
