// Copyright (c) 2005-2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot, Sylvain Pion

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<boost/shared_ptr.hpp>
#include <CGAL/CGAL_Ipelet_base.h>
#include<CGAL/create_offset_polygons_2.h>
#include <boost/format.hpp>


namespace CGAL_skeleton{

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

const std::string Slab[] = {
  "Interior skeleton", "Exterior skeleton","Interior offset","Exterior offset","Interior offsets","Exterior offsets", "Help"
};

const std::string Hmsg[] = {
  "Draw the interior skeleton of one polygon",
  "Draw the exterior skeleton of one polygon",
  "Draw an interior offset of one polygon",
  "Draw an exterior offset of one polygon",
  "Draw several interior offsets of one polygon",
  "Draw several exterior offsets of one polygon"
};


class SkeletonIpelet
  : public CGAL::Ipelet_base<Kernel,7>{

  typedef boost::shared_ptr<Polygon_2>        PolygonPtr ;
  typedef std::vector<PolygonPtr>             PolygonPtrVector ;
  typedef CGAL::Straight_skeleton_2<Kernel>   Skeleton ;
  typedef boost::shared_ptr<Skeleton>         SkeletonPtr ;

  void draw_straight_skeleton(const Skeleton& skeleton,double);

public:
  SkeletonIpelet()
    :CGAL::Ipelet_base<Kernel,7>("Skeleton and offset",Slab,Hmsg){}
  void protected_run(int);
};

void SkeletonIpelet::draw_straight_skeleton(const Skeleton& skeleton,double /*max_edge*/)
{
  typedef Skeleton::Vertex_const_handle     Vertex_const_handle ;
  typedef Skeleton::Halfedge_const_handle   Halfedge_const_handle ;
  typedef Skeleton::Halfedge_const_iterator Halfedge_const_iterator ;

  Halfedge_const_handle null_halfedge ;
  Vertex_const_handle   null_vertex ;

  std::list<Segment_2> seglist;
  std::back_insert_iterator< std::list<Segment_2> > out=std::back_inserter(seglist);

  for ( Halfedge_const_iterator i = skeleton.halfedges_begin();
                                i != skeleton.halfedges_end();
                                ++i )
    if ( i->is_bisector() && ((i->id()%2)==0) ){
        out++=Segment_2(i->opposite()->vertex()->point(),i->vertex()->point());
    }
  draw_in_ipe(seglist.begin(),seglist.end());
}

void SkeletonIpelet::protected_run(int fn)
{

  if (fn==6) {
    show_help();
    return;
  }

  std::list<Polygon_2> pol_list;
  Iso_rectangle_2 bbox=
    read_active_objects( CGAL::dispatch_or_drop_output<Polygon_2>( std::back_inserter(pol_list) ) );



  if (pol_list.size()!=1){
    print_error_message("Exactly one polygon must be selected");
    return;
  }

  Polygon_2 polygon=*pol_list.begin();


  if (!polygon.is_simple()){
    print_error_message("Polygon must be simple");
    return;
  }

  if (polygon.orientation()!=CGAL::COUNTERCLOCKWISE)
    polygon.reverse_orientation();

  std::list<double> offsets;
    //~ "Interior skeleton", "Exterior skeleton","Interior offset","Exterior offset","Interior offsets","Exterior offsets", "Help"
  SkeletonPtr ss;
  double max_edge=(std::max)((bbox.xmax()-bbox.xmin()),(bbox.ymax()-bbox.ymin()));
  double dist=0.;
  int ret_val=-1;
  switch(fn){
    case 3://Exterior offset
    case 5://Exterior offsets
    case 1://Exterior skeleton
      ss = CGAL::create_exterior_straight_skeleton_2(max_edge,polygon);
      break;
    case 2://Interior offset
    case 4://Interior offsets
    case 0://Interior skeleton
      ss = CGAL::create_interior_straight_skeleton_2(polygon);
      break;
  }


  if (fn==0 || fn==1)
    draw_straight_skeleton(*ss,max_edge);
  else{
    boost::tie(ret_val,dist)=
      request_value_from_user<double>(
        (boost::format("Offset value (BBox %1%x%2%)") % (bbox.xmax()-bbox.xmin()) % (bbox.ymax()-bbox.ymin())).str()
      );
    if (ret_val == -1){
      print_error_message("Bad value provided");
      return;
    }

    if (fn==2 || fn==3)
      offsets.push_back(dist);
    else{
      for (int i=1;i<static_cast<int>(ceil(max_edge/dist/2.))+1;++i)
        offsets.push_back(i*dist);
    }

    for (std::list<double>::iterator it=offsets.begin();it!=offsets.end();++it){
      PolygonPtrVector offset_polygons = CGAL::create_offset_polygons_2<Polygon_2>(*it,*ss);
      for( PolygonPtrVector::const_iterator pi = offset_polygons.begin() ; pi != offset_polygons.end() ; ++ pi )
        draw_in_ipe(**pi);
    }

    if (offsets.size()>1)
      group_selected_objects_();
  }
}

}







CGAL_IPELET(CGAL_skeleton::SkeletonIpelet)

