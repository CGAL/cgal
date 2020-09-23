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
#include <CGAL/CGAL_Ipelet_base.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Apollonius_site_2.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>



namespace CGAL_Hull{

typedef CGAL::Exact_predicates_inexact_constructions_kernel               Kernel;

//Apollonius
typedef CGAL::Apollonius_graph_traits_2<Kernel>                           AT;
typedef CGAL::Apollonius_graph_2<AT>                                      Apollonius;
typedef Apollonius::Site_2                                                ASite;
//Delaunay
typedef CGAL::Triangulation_vertex_base_with_info_2<bool,Kernel>          Vb;
typedef CGAL::Triangulation_face_base_2<Kernel>                           Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                       Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel,Tds>                        Delaunay;
typedef Delaunay::Vertex_handle                                           Vertex_handle;
typedef Delaunay::Finite_faces_iterator                                   Face_iterator;


const std::string sublabel[] = {
  "Convex minimal", "Crust", "Help"
};


const std::string helpmsg[] = {
  "Draw the convex hull of a set of segments, circles and points",
  "Draw the result of the crust algorithm for a set of points"
};

// --------------------------------------------------------------------

class enveloppeIpelet
  : public CGAL::Ipelet_base<Kernel,3>{
public:
  enveloppeIpelet()
    : CGAL::Ipelet_base<Kernel,3>("Hulls",sublabel,helpmsg){}
  void protected_run(int);
  //compute tangency points of the convex hull of 2 circles using pythagore and Thales
  IpeVector tangency_point(double c_radius,double p_radius,const Point_2& current_pt,const Point_2& previous_pt,int f=1)
  {
    int i=(c_radius>=p_radius)?(1):(-1);
    double angle = atan2(i*(previous_pt.y()-current_pt.y()),i*(previous_pt.x()-current_pt.x()));
    if(c_radius!=p_radius){
      double dist = pow(previous_pt.x()-current_pt.x(),2)+pow(previous_pt.y()-current_pt.y(),2);
      dist = dist*pow(c_radius/(c_radius-p_radius),2);
      angle = angle - i * f * atan(sqrt(dist-pow(c_radius,2))/c_radius);
    }
    else
      angle = angle - f * CGAL_PI/2.;
    return IpeVector((current_pt.x()+c_radius*cos(angle)),(current_pt.y()+c_radius*sin(angle)));
  }

};

void enveloppeIpelet::protected_run(int fn)
{
  if (fn==2) {
    show_help();
    return;
  }

  switch(fn){
    case 0:{//Convex hull using apollonius
      std::list<Point_2> pt_list;
      std::list<Circle_2> cir_list;

      read_active_objects(
        CGAL::dispatch_or_drop_output<Point_2,Circle_2,Segment_2,Polygon_2>(
          std::back_inserter(pt_list),
          std::back_inserter(cir_list),
          point_grabber(std::back_inserter(pt_list)),
          point_grabber(std::back_inserter(pt_list))
        )
     );

      if (pt_list.empty() && cir_list.empty()) {
        print_error_message("No mark nor circle nor segment selected");
        return;
      }

      Apollonius apo;     //apollonius

      for (std::list<Point_2>::iterator it=pt_list.begin();it!=pt_list.end();++it)
        apo.insert(ASite(*it,0));
      for (std::list<Circle_2>::iterator it=cir_list.begin();it!=cir_list.end();++it)
        apo.insert(ASite(it->center(),sqrt(it->squared_radius())));

      if (apo.number_of_vertices()==1) {
        print_error_message("Need more than one mark or one circle");
        return;
      }

      Apollonius::Vertex_circulator Cvert = apo.incident_vertices(apo.infinite_vertex()); //take points incident to infinte vertex
      Apollonius::Vertex_circulator Cvert0 = Cvert;
      std::vector<ASite> Vsite0;
      do{
        Vsite0.push_back(Cvert->site());
        ++Cvert;
      }while(Cvert0!=Cvert);

      IpeVector pt_ipe1=IpeVector(0,0);
      IpeSegmentSubPath* SSPseg_ipe = new IpeSegmentSubPath;
      Vsite0.insert(Vsite0.begin(),Vsite0.back());
      Vsite0.insert(Vsite0.end(),*(Vsite0.begin()+1));
      std::vector<ASite>::iterator Vsiteite0 = Vsite0.begin();
      std::vector<ASite>::iterator Vsiteite00 = Vsite0.end()-2;
      #ifdef CGAL_USE_IPE_7
      for(std::vector<ASite>::iterator it=Vsiteite00 ; it!=Vsiteite0 ; --it){//draw precise convex hull computing tangency point to circles
        double c_rad = it->weight();
        if(c_rad!=0){
          Point_2 p_pt = (it-1)->point();  //previous neighbor
          Point_2 c_pt = it->point();
          Point_2 n_pt = (it+1)->point(); //next neighbor
          double p_rad = (it-1)->weight();
          double n_rad = (it+1)->weight();
          IpeVector pt_ipe=tangency_point(c_rad,p_rad,c_pt,p_pt);
          IpeVector pt_ipe0=tangency_point(c_rad,n_rad,c_pt,n_pt,-1);

          if(it!=Vsiteite00)
            SSPseg_ipe->appendSegment(pt_ipe1,pt_ipe0);
          SSPseg_ipe->appendArc(IpeMatrix(c_rad,0,0,c_rad,c_pt.x(),c_pt.y()),pt_ipe0,pt_ipe);
          pt_ipe1=pt_ipe;
        }
        else{
          Point_2 c_pt = it->point();
          IpeVector pt_ipe=IpeVector(c_pt.x(),c_pt.y());
          if(it!=Vsiteite00)
            SSPseg_ipe->appendSegment(pt_ipe1,pt_ipe);
          pt_ipe1=pt_ipe;
        }
      }
      SSPseg_ipe->setClosed(true);
      ipe::Shape shape;
      shape.appendSubPath(SSPseg_ipe);
      get_IpePage()->append(ipe::EPrimarySelected,CURRENTLAYER,new ipe::Path(CURRENTATTRIBUTES,shape));
      #else
      for(std::vector<ASite>::iterator it=Vsiteite00 ; it!=Vsiteite0 ; --it){//draw precise convex hull computing tangency point to circles
        double c_rad = it->weight();
        if(c_rad!=0){
          Point_2 p_pt = (it-1)->point();  //previous neighbor
          Point_2 c_pt = it->point();
          Point_2 n_pt = (it+1)->point(); //next neighbor
          double p_rad = (it-1)->weight();
          double n_rad = (it+1)->weight();
          IpeVector pt_ipe=tangency_point(c_rad,p_rad,c_pt,p_pt);
          IpeVector pt_ipe0=tangency_point(c_rad,n_rad,c_pt,n_pt,-1);

          if(it!=Vsiteite00)
            SSPseg_ipe->AppendSegment(pt_ipe1,pt_ipe0);
          SSPseg_ipe->AppendArc(IpeMatrix(c_rad,0,0,c_rad,c_pt.x(),c_pt.y()),pt_ipe0,pt_ipe);
          pt_ipe1=pt_ipe;
        }
        else{
          Point_2 c_pt = it->point();
          IpeVector pt_ipe=IpeVector(c_pt.x(),c_pt.y());
          if(it!=Vsiteite00)
            SSPseg_ipe->AppendSegment(pt_ipe1,pt_ipe);
          pt_ipe1=pt_ipe;
        }
      }
      SSPseg_ipe->SetClosed(true);
      IpePath* obj_ipe1 = new IpePath(get_IpeletHelper()->Attributes());
      obj_ipe1 -> AddSubPath(SSPseg_ipe);
      get_IpePage()->push_back(IpePgObject(IpePgObject::ESecondary,get_IpeletHelper()->CurrentLayer(),obj_ipe1));
      #endif
    }
    break;

    case 1:{          //CRUST
      std::list<Point_2> pt_list;
      read_active_objects( CGAL::dispatch_or_drop_output<Point_2>(std::back_inserter(pt_list)) );


      if (pt_list.empty()) {
        print_error_message("No mark selected");
        return;
      }

      Delaunay triang;

      for(std::list<Point_2>::iterator it = pt_list.begin();it!=pt_list.end();++it){
        Vertex_handle vI_cgal = triang.insert(*it);
        vI_cgal -> info() = true;
      }

      //add voronoi center to the delauney triangulation
      pt_list.clear();
      for (Face_iterator itt=triang.finite_faces_begin();itt!=triang.finite_faces_end();++itt)
        pt_list.push_back(triang.dual(itt));

      for (std::list<Point_2>::iterator it=pt_list.begin();it!=pt_list.end();++it){
        Vertex_handle vI_cgal = triang.insert(*it);
        vI_cgal -> info() = false;
      }

      //keep delaunay segments which endpoints are original points
      for(Delaunay::Finite_edges_iterator it= triang.finite_edges_begin();it!=triang.finite_edges_end();++it){
        if (it->first->vertex(Delaunay::cw(it->second))->info()==true &&
              it->first->vertex(Delaunay::ccw(it->second))->info()==true)
          draw_in_ipe(Segment_2( it->first->vertex(Delaunay::cw(it->second))->point(),
                                      it->first->vertex(Delaunay::ccw(it->second))->point()
                                    )
                     );
      }
      group_selected_objects_();
    }
  }
}

}
// --------------------------------------------------------------------

CGAL_IPELET(CGAL_Hull::enveloppeIpelet)


// --------------------------------------------------------------------
