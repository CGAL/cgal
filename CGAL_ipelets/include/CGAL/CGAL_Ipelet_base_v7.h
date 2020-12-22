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


#ifndef CGAL_IPELET_BASE_H
#define CGAL_IPELET_BASE_H

// Ipe headers use uint which is not standard.
#ifdef __APPLE__
typedef unsigned int uint;
#endif

#include <ipelib.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/iterator.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/grabbers.h>
#include <CGAL/iterator.h>
#include <CGAL/tuple.h>
#include<CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/assertions.h>
#include <CGAL/use.h>

#include <boost/utility.hpp>

#define CURRENTLAYER get_IpeletData()->iLayer
#define CURRENTATTRIBUTES get_IpeletData()->iAttributes

namespace CGAL{

  template <class Kernel,int nbf>
  class Ipelet_base : public ipe::Ipelet {
  private:
    const std::string* SubLab;
    const std::string* HMsg;
    std::string Name;
    ipe::IpeletData* data_;
    ipe::IpeletHelper* helper_;

  public:

    typedef ipe::Vector IpeVector; //ipe6 compatibility
    typedef ipe::Curve  IpeSegmentSubPath;//ipe6 compatibility
    typedef ipe::Matrix IpeMatrix;//ipe6 compatibility
    typedef ipe::Path   IpePath;//ipe6 compatibility
    //indicates if the selection should be primary or secondary. Exactly one primary selection should exist
    ipe::TSelect get_selection_type() const { return get_IpePage()->primarySelection()==-1 ? ipe::EPrimarySelected : ipe::ESecondarySelected;}
    //ipe6 compatibility
    void transform_selected_objects_(const IpeMatrix& tfm) const {
      for (int i=0;i<static_cast<int>(get_IpePage()->count());++i)
        if (get_IpePage()->select(i)!=ipe::ENotSelected)
          get_IpePage()->transform(i,tfm);
    }

    void delete_selected_objects_() const {
      for (unsigned i=static_cast<unsigned>(get_IpePage()->count());i>0;--i)
        if (get_IpePage()->select(i-1)!=ipe::ENotSelected)
          get_IpePage()->remove(i-1);
    }

    void group_selected_objects_() const {
      ipe::Group* grp=new ipe::Group();
      for (unsigned i=static_cast<unsigned>(get_IpePage()->count());i>0;--i)
        if (get_IpePage()->select(i-1)!=ipe::ENotSelected){
          grp->push_back( get_IpePage()->object(i-1)->clone() );
          //~ grp->push_back( get_IpePage()->object(i-1) );
          get_IpePage()->remove(i-1);
        }
      get_IpePage()->append(get_selection_type(),CURRENTLAYER,grp);
    }



    //typedefs
    typedef typename Kernel::FT                                               FT;
    typedef typename CGAL::Point_2<Kernel>                                    Point_2;
    typedef typename CGAL::Weighted_point_2<Kernel>                           Weighted_point_2;
    typedef typename Kernel::Segment_2                                        Segment_2;
    typedef typename Kernel::Ray_2                                            Ray_2;
    typedef typename Kernel::Line_2                                           Line_2;
    typedef typename Kernel::Iso_rectangle_2                                  Iso_rectangle_2;
    typedef typename Kernel::Triangle_2                                       Triangle_2;
    //~ typedef typename CGAL::Polygon_2<Kernel,std::list<Point_2> >              Polygon_2;
    typedef CGAL::Polygon_2<Kernel>                                           Polygon_2;

    typedef typename Kernel::Circle_2                                         Circle_2;
    typedef std::tuple<Circle_2,Point_2,Point_2,CGAL::Sign>           Circular_arc_2;


    Ipelet_base(const std::string NameS,const std::string SubLabS[],const std::string HMsgS[])
      :SubLab(&SubLabS[0]),HMsg(&HMsgS[0]),Name(NameS),data_(nullptr),helper_(nullptr){};


    ipe::Page* get_IpePage() const {return data_->iPage;}
    ipe::IpeletData* get_IpeletData() const {return data_;}
    ipe::IpeletHelper* get_IpeletHelper() const {return helper_;}
    int ipelibVersion() const { return ipe::IPELIB_VERSION; }
    int NumFunctions() const { return nbf; }
    virtual const char *Label() const{ return &Name[0]; }
    const char *About() const {return "https://www.cgal.org";};
    virtual const char *SubLabel(int function) const {return &SubLab[function][0];};
    virtual const char *HelpMsg(int function) const{return &HMsg[function][0];};
    bool run (int i, ipe::IpeletData* data, ipe::IpeletHelper* helper) {
      data_=data;
      helper_=helper;
      try{
        protected_run(i);
        return true;
      }
      catch(...){
        helper->messageBox("Error : Save your page in a file and submit it to \n https://www.cgal.org/bug_report.html",nullptr,ipe::IpeletHelper::EOkCancelButtons);
        return false;
      }
    };

    virtual void protected_run(int)=0;

    void show_help(bool gen=true) const{
      std::string hmsg;
      hmsg="<qt><h1>"+Name+"</h1><ul>";
      if (gen)
        for(int i=0;i<nbf-1;++i)
          hmsg=hmsg+"<li><i>"+SubLab[i]+"</i>: "+HMsg[i]+"</li>";
      else
        hmsg=hmsg+"<li>"+HMsg[0]+"</li>";
      get_IpeletHelper()->messageBox(&hmsg[0],nullptr,ipe::IpeletHelper::EOkCancelButtons);
      return;
    }


    //grabbers

    template <class output_iterator>
    struct Point_grabber:public internal::Point_grabber<Kernel,output_iterator>{
      Point_grabber(output_iterator it):internal::Point_grabber<Kernel,output_iterator>(it){}
    };

    template<class output_iterator>
    boost::function_output_iterator<Point_grabber<output_iterator> >
    point_grabber(output_iterator it){
      return boost::make_function_output_iterator(Point_grabber<output_iterator>(it));
    }


    template <class output_iterator>
    struct Segment_grabber:public internal::Segment_grabber<Kernel,output_iterator>{
      Segment_grabber(output_iterator it):internal::Segment_grabber<Kernel,output_iterator>(it){}
    };

    template<class output_iterator>
    boost::function_output_iterator<Segment_grabber<output_iterator> >
    segment_grabber(output_iterator it){
      return boost::make_function_output_iterator(Segment_grabber<output_iterator>(it));
    }

    template <class output_iterator>
    struct Wpoint_grabber:public internal::Wpoint_grabber<Kernel,output_iterator>{
      Wpoint_grabber(output_iterator it):internal::Wpoint_grabber<Kernel,output_iterator>(it){}
    };

    template<class output_iterator>
    boost::function_output_iterator<Wpoint_grabber<output_iterator> >
    wpoint_grabber(output_iterator it){
      return boost::make_function_output_iterator(Wpoint_grabber<output_iterator>(it));
    }

    //Interaction functions
    //------------------------------

    void
    print_error_message(const char* s) const
    {
      get_IpeletHelper()->message(s);
    }

    template <class T>
    std::pair<int,T>
    request_value_from_user(std::string msg) const
    {
      ipe::String str;
      std::pair<int,T> ret_val=std::make_pair(-1,T());
      if (get_IpeletHelper()-> getString(msg.c_str(),str)){
        if (!str.empty()){
          ipe::Lex lex(str);
          lex >> ret_val.second;
          ret_val.first=1;
        }
        else
          ret_val.first=0;
      }
      return ret_val;
    }

    //Conversion functions
    //------------------------------
    Point_2
    segment_endpoint(const ipe::CurveSegment& segment,ipe::Path* obj_ipe,int i) const
    {
      CGAL_precondition(i<2);
      ipe::Vector pt_ipe = obj_ipe -> matrix() * segment.cp(i);
      return Point_2((double)(pt_ipe.x),(double)(pt_ipe.y));//conversion into CGAL point
    }

    Point_2
    to_point_2(ipe::Object*  object) const
    {
      ipe::Vector pt_ipe = object-> matrix() * object-> asReference() -> position();
      return Point_2((double)(pt_ipe.x),(double)(pt_ipe.y));//conversion into CGAL point
    }

    Circle_2
    to_circle_2(ipe::Path* obj_ipe,int subpath=0) const
    {
      const ipe::Ellipse* ell_ipe = obj_ipe -> shape().subPath(subpath) -> asEllipse();
      ipe::Matrix mat_ipe = obj_ipe -> matrix() * ell_ipe -> matrix();
      FT radius = (mat_ipe*ipe::Vector(1,0)-mat_ipe.translation()).len();
      ipe::Vector pt_ipe = mat_ipe.translation();
      return Circle_2(Point_2(pt_ipe.x,pt_ipe.y),radius*radius);
    }


    //Picking functions
    //------------------------------

    bool
    is_only_rotated_or_scaled(const ipe::Matrix& m) const
    {
      return (m.a[0]==m.a[3] && m.a[1]==-m.a[2]);
    }

    bool
    is_IPE_circle(ipe::Object* object,int subpath=0) const
    {
      return ( object -> asPath() && object -> asPath() -> shape().subPath(subpath) -> asEllipse()
        && is_only_rotated_or_scaled(object ->asPath()->matrix()));
    }


public:
    //declaration
    template< class multi_output_iterator >
    bool read_one_active_object( ipe::Object* object,
      multi_output_iterator it_out) const;

public:

    template< class V,class O>
    Iso_rectangle_2
    read_active_objects (
      CGAL::Dispatch_or_drop_output_iterator<V,O> it_out,
      bool deselect_all=true,
      bool delete_selected_objects=false) const
    {
      ipe::Rect bbox_ipe;

      if (!get_IpePage()->hasSelection()) {
        return Iso_rectangle_2();
      }

      for (int i=0;i<static_cast<int>(get_IpePage()->count());++i){
        if (get_IpePage()->select(i)==ipe::ENotSelected)
          continue;

        bbox_ipe.addRect(get_IpePage()->bbox(i));

        //Test one function for segments, circles, circle arcs and polygons
        bool desel_it=read_one_active_object(get_IpePage()->object(i),it_out);
        if ( delete_selected_objects && desel_it  )
          get_IpePage()->setSelect(i,ipe::ENotSelected);
      }

      if (delete_selected_objects)
        delete_selected_objects_();

      if (deselect_all)
        get_IpePage()->deselectAll();

      Iso_rectangle_2 bbox_cgal(
        static_cast<double>(bbox_ipe.bottomLeft().x),static_cast<double>(bbox_ipe.bottomLeft().y),
        static_cast<double>(bbox_ipe.topRight().x),static_cast<double>(bbox_ipe.topRight().y)
      );

        return bbox_cgal;
    }

    //drawing functions
    //------------------------------
    void
    create_polygon_with_holes(bool delete_underlying_polygons=false) const
    {
      std::list<ipe::SubPath*> SSPqu;
      for (int i=0;i<static_cast<int>(get_IpePage()->count());++i){
        if (get_IpePage()->select(i)!=ipe::ENotSelected && get_IpePage()->object(i)->asPath()->shape().subPath(0)->closed() ){
          ipe::SubPath* ssp=new ipe::Curve(*get_IpePage()->object(i)->asPath()->shape().subPath(0)->asCurve());
          SSPqu.push_back(ssp);
        }
      }
      if (!delete_underlying_polygons)
        get_IpePage() -> deselectAll();
      ipe::Shape shape;// create new objects with current attributes
      for (std::list<ipe::SubPath*>::iterator it=SSPqu.begin();it!=SSPqu.end();++it)
        shape.appendSubPath(*it);
      if (delete_underlying_polygons)
        delete_selected_objects_();
      get_IpePage()->append(get_selection_type(),CURRENTLAYER,new ipe::Path(CURRENTATTRIBUTES,shape));
    }

    void
    center_selection_in_page() const
    {
      ipe::Vector paper_size=get_paper_size();
      ipe::Matrix tfm (1,0,0,1,paper_size.x/2.,paper_size.y/2.);
      for (int i=0;i<static_cast<int>(get_IpePage()->count());++i)
        if (get_IpePage()->select(i)!=ipe::ENotSelected )
          get_IpePage()->transform(i,tfm);
    }

    template<class iterator>
    ipe::Curve*
    create_polyline(const iterator first, const iterator last,bool setclose=false) const
    {
      if (boost::next(first)!=last){
        ipe::Curve* SSP_ipe = new ipe::Curve();
        ipe::Vector Prev_pt=ipe::Vector(CGAL::to_double(first->x()),CGAL::to_double(first->y())) ;
        for (iterator it = boost::next(first);it!=last;++it){
          ipe::Vector Cur_pt=ipe::Vector(CGAL::to_double(it->x()),CGAL::to_double(it->y()));
          SSP_ipe -> appendSegment(Prev_pt,Cur_pt);
          Prev_pt=Cur_pt;
        }
        if (setclose)
          SSP_ipe->setClosed(true);
        return SSP_ipe;
      }
      return nullptr;
    }


    template<class iterator>
    ipe::Path*
    draw_polyline_in_ipe(const iterator first, const iterator last,
                         bool setclose=false,bool deselect_all=false,
                         bool blackfill=false,
                         typename boost::enable_if<
                                    boost::is_same<
                                      typename std::iterator_traits<iterator>::value_type,
                                      Point_2
                                    >
                                  >::type* =nullptr) const
    {
      ipe::Curve* SSP_ipe=create_polyline(first,last,setclose);
      if (SSP_ipe!=nullptr){
        ipe::Shape shape;
        shape.appendSubPath(SSP_ipe);
        ipe::Path* obj_ipe=new ipe::Path(CURRENTATTRIBUTES,shape);
        if (blackfill){
          obj_ipe->setPathMode(ipe::EStrokedAndFilled);
          obj_ipe->setFill(ipe::Attribute::BLACK());
        }
        get_IpePage()->append( (deselect_all?ipe::ENotSelected:get_selection_type()),CURRENTLAYER,obj_ipe);
        return obj_ipe;
      }
      return nullptr;
    }

    void draw_in_ipe(const Circle_2& C,bool deselect_all=false) const {
      ipe::Ellipse *ellipse = new ipe::Ellipse(ipe::Matrix(sqrt(CGAL::to_double(C.squared_radius())),0,
                                                     0,sqrt(CGAL::to_double(C.squared_radius())),
                                                     to_double(C.center().x()),to_double(C.center().y())
                                           )
                                );
      ipe::Shape shape;
      shape.appendSubPath(ellipse);
      get_IpePage()->append( (deselect_all?ipe::ENotSelected:get_selection_type()),CURRENTLAYER,new ipe::Path(CURRENTATTRIBUTES,shape));
    }

    void
    draw_in_ipe(const Point_2& P,bool deselect_all=false) const
    {
      ipe::Reference *mark = new ipe::Reference(CURRENTATTRIBUTES,CURRENTATTRIBUTES.iMarkShape, ipe::Vector(CGAL::to_double(P.x()),CGAL::to_double(P.y())));
      get_IpePage()->append( (deselect_all?ipe::ENotSelected:get_selection_type()),CURRENTLAYER,mark);
    }

    void
    draw_in_ipe(const Segment_2& S,bool deselect_all=false) const
    {
      ipe::Segment seg_ipe;
      seg_ipe.iP = ipe::Vector(CGAL::to_double(S.point(0).x()),CGAL::to_double(S.point(0).y()));
      seg_ipe.iQ = ipe::Vector(CGAL::to_double(S.point(1).x()),CGAL::to_double(S.point(1).y()));
      get_IpePage()->append( (deselect_all?ipe::ENotSelected:get_selection_type()),CURRENTLAYER,new ipe::Path(CURRENTATTRIBUTES,ipe::Shape(seg_ipe)));
    }

    template<class Container>
    void
    draw_in_ipe(const CGAL::Polygon_2<Kernel,Container>& poly,bool deselect_all=false) const
    {
      std::list<Point_2> LP;
      for (typename CGAL::Polygon_2<Kernel,Container>::iterator it=poly.vertices_begin();it!= poly.vertices_end();++it)
        LP.push_back(*it);
      draw_polyline_in_ipe(LP.begin(),LP.end(),true,deselect_all,false);
    }

    void
    draw_in_ipe(const Circular_arc_2& arc,bool deselect_all=false) const
    {
      ipe::Curve* SSP_ipe = new ipe::Curve;
      ipe::Vector ipeS=ipe::Vector( CGAL::to_double(std::get<1>(arc).x()),
                                CGAL::to_double(std::get<1>(arc).y()));//convert ot ipe format
      ipe::Vector ipeT=ipe::Vector( CGAL::to_double(std::get<2>(arc).x()),
                                CGAL::to_double(std::get<2>(arc).y()));//convert ot ipe format
      SSP_ipe->appendArc(ipe::Matrix(sqrt(CGAL::to_double(std::get<0>(arc).squared_radius())),0,
                                   0,(std::get<3>(arc)==CGAL::COUNTERCLOCKWISE?1:-1)*
                                     sqrt(CGAL::to_double(std::get<0>(arc).squared_radius())),
                                   CGAL::to_double(std::get<0>(arc).center().x()),
                                   CGAL::to_double(std::get<0>(arc).center().y())),
                                   ipeS,ipeT);
      ipe::Shape shape;
      shape.appendSubPath(SSP_ipe);
      get_IpePage()->append( (deselect_all?ipe::ENotSelected:get_selection_type()),CURRENTLAYER,new ipe::Path(CURRENTATTRIBUTES,shape));
    }


    void
    draw_in_ipe(const Triangle_2& t,bool deselect_all=false) const
    {
      ipe::Curve* SSP_ipe = new ipe::Curve();
      ipe::Vector P0=ipe::Vector(t[0].x(),t[0].y());
      ipe::Vector P1=ipe::Vector(t[1].x(),t[1].y());
      ipe::Vector P2=ipe::Vector(t[2].x(),t[2].y());
      SSP_ipe->appendSegment(P0,P1);
      SSP_ipe->appendSegment(P1,P2);
      SSP_ipe->appendSegment(P2,P0);
      SSP_ipe->setClosed(true);
      ipe::Shape shape;
      shape.appendSubPath(SSP_ipe);
      get_IpePage()->append( (deselect_all?ipe::ENotSelected:get_selection_type()),CURRENTLAYER,new ipe::Path(CURRENTATTRIBUTES,shape));
    }

    void
    draw_in_ipe(const Iso_rectangle_2& r,bool deselect_all=false)
    {
      ipe::Curve* SSP_ipe = new ipe::Curve();
      ipe::Vector P0=ipe::Vector(r[0].x(),r[0].y());
      ipe::Vector P1=ipe::Vector(r[1].x(),r[1].y());
      ipe::Vector P2=ipe::Vector(r[2].x(),r[2].y());
      ipe::Vector P3=ipe::Vector(r[3].x(),r[3].y());
      SSP_ipe->appendSegment(P0,P1);
      SSP_ipe->appendSegment(P1,P2);
      SSP_ipe->appendSegment(P2,P3);
      SSP_ipe->appendSegment(P3,P0);
      SSP_ipe->setClosed(true);
      ipe::Shape shape;
      shape.appendSubPath(SSP_ipe);
      get_IpePage()->append( (deselect_all?ipe::ENotSelected:get_selection_type()),CURRENTLAYER,new ipe::Path(CURRENTATTRIBUTES,shape));
    }


    //Drawing function with bbox : global version
    template <class T>
    void
    draw_in_ipe(const T& object,const Iso_rectangle_2& bbox,bool deselect_all=false) const
    {
      Segment_2 s;
      bool success=cast_into_seg(object,bbox,&s);
      if (success)
        draw_in_ipe(s,deselect_all);
    }
  private:
    enum Type_circ_arc{SRC,TRG,OSRC,OTRG};
  public:
    void
    draw_in_ipe(const Circular_arc_2& object,const Iso_rectangle_2& bbox,bool deselect_all=false) const
    {
      std::vector<Circular_arc_2> arc_list;
      const Circle_2& circle=std::get<0>(object);
      restrict_circle_to_bbox(circle,bbox,std::back_inserter(arc_list));
      if (arc_list.empty() && bbox.has_on_bounded_side(circle.center()) ){
        draw_in_ipe(object,deselect_all);
        return;
      }

      const Point_2* source=(std::get<3>(object)==CGAL::COUNTERCLOCKWISE)?
                            &std::get<1>(object):&std::get<2>(object);
      const Point_2* target=(std::get<3>(object)==CGAL::COUNTERCLOCKWISE)?
                            &std::get<2>(object):&std::get<1>(object);
      std::multimap<double,std::pair<Type_circ_arc,const Point_2*> > map_theta;
      typedef  typename std::multimap<double,std::pair<Type_circ_arc,const Point_2*> >::iterator Map_theta_iterator;
      Map_theta_iterator s_it=map_theta.insert(
                                std::make_pair(get_theta(*source,circle),std::make_pair(OSRC,source)));
      /* Map_theta_iterator t_it=*/
                              map_theta.insert(
                                std::make_pair(get_theta(*target,circle),std::make_pair(OTRG,target)));

      for (typename std::vector<Circular_arc_2>::iterator it_arc=arc_list.begin();it_arc!=arc_list.end();++it_arc){
        const Point_2* arc_s=(std::get<3>(*it_arc)==CGAL::COUNTERCLOCKWISE)?
                             &std::get<1>(*it_arc):&std::get<2>(*it_arc);
        const Point_2* arc_t=(std::get<3>(*it_arc)==CGAL::COUNTERCLOCKWISE)?
                             &std::get<2>(*it_arc):&std::get<1>(*it_arc);
        map_theta.insert( std::make_pair(get_theta(*arc_s,circle),std::make_pair(SRC,arc_s) ) );
        map_theta.insert( std::make_pair(get_theta(*arc_t,circle),std::make_pair(TRG,arc_t) ) );
      }

      Map_theta_iterator next_s=s_it;
      ++next_s; if (next_s==map_theta.end()) next_s=map_theta.begin();
      switch (next_s->second.first){
        case TRG:
          draw_in_ipe(Circular_arc_2(circle,*source,*(next_s->second.second),CGAL::COUNTERCLOCKWISE));
          break;
        case OSRC:
          CGAL_error();
        case SRC:{
          Map_theta_iterator current=next_s;
          ++next_s; if (next_s==map_theta.end()) next_s=map_theta.begin();
          draw_in_ipe(Circular_arc_2(circle,*(current->second.second),*(next_s->second.second),CGAL::COUNTERCLOCKWISE));
          if(next_s->second.first==OTRG) return;
          break;
        }
        case OTRG:
          ++next_s; if (next_s==map_theta.end()) next_s=map_theta.end();
          if (next_s->second.first==TRG){
            draw_in_ipe(object);
            return;
          }
          else
            return;
      }

      ++next_s; if (next_s==map_theta.end()) next_s=map_theta.begin();
      Map_theta_iterator current=next_s;
      ++next_s; if (next_s==map_theta.end()) next_s=map_theta.begin();
      do{
        if (current->second.first==OTRG) return;
        CGAL_assertion(current->second.first==SRC);
        draw_in_ipe(Circular_arc_2(circle,*(current->second.second),*(next_s->second.second),CGAL::COUNTERCLOCKWISE));
        if (next_s->second.first==OTRG) return;
        CGAL_assertion(next_s->second.first==TRG);
        ++next_s; if (next_s==map_theta.end()) next_s=map_theta.begin();
        current=next_s;
        ++next_s; if (next_s==map_theta.end()) next_s=map_theta.begin();
      }while(true);

    }


    void
    draw_in_ipe(const Circle_2& object,const Iso_rectangle_2& bbox,bool deselect_all=false) const
    {
      std::vector<Circular_arc_2> arc_list;
      restrict_circle_to_bbox(object,bbox,std::back_inserter(arc_list));
      if (arc_list.empty() && bbox.has_on_bounded_side(object.center()) )
        draw_in_ipe(object,deselect_all);
      else
        draw_in_ipe(arc_list.begin(),arc_list.end(),false,deselect_all);
    }


    void
    draw_in_ipe(const Triangle_2& object,const Iso_rectangle_2& bbox,bool deselect_all=false) const
    {
       for (unsigned int i=0;i!=3;++i)
        draw_in_ipe(Segment_2(object.vertex(i),object.vertex(i+1)),bbox,deselect_all);
    }

    void
    draw_in_ipe(const Iso_rectangle_2& object,const Iso_rectangle_2& bbox,bool deselect_all=false) const
    {
      for (unsigned int i=0;i!=4;++i)
        draw_in_ipe(Segment_2(object.vertex(i),object.vertex(i+1)),bbox,deselect_all);
    }

    void
    draw_in_ipe(const Polygon_2& object,const Iso_rectangle_2& bbox,bool deselect_all=false) const
    {
      for (typename Polygon_2::Edge_const_iterator it=object.edges_begin();it!=object.edges_end();++it)
        draw_in_ipe(*it,bbox,deselect_all);
    }



    template<class GT,class TDS>
    void
    draw_in_ipe(const CGAL::Triangulation_2<GT,TDS>& tri,const Iso_rectangle_2& bbox,bool deselect_all=false) const
    {
      typedef CGAL::Triangulation_2<GT,TDS> Triangulation;
      typename Triangulation::Finite_edges_iterator first=tri.finite_edges_begin();
      typename Triangulation::Finite_edges_iterator last=tri.finite_edges_end();
      for(typename Triangulation::Finite_edges_iterator it=first;it!=last;++it)
        draw_in_ipe(tri.segment(*it),bbox);
      if (deselect_all)
        get_IpePage()->deselectAll();
    }

    void
    draw_in_ipe(const Line_2& line,bool deselect_all=false) const
    {
      ipe::Vector paper_size=get_paper_size();
      Iso_rectangle_2 bbox(0,0,paper_size.x,paper_size.y);
      draw_in_ipe(line,bbox,deselect_all);
    }

    void
    draw_in_ipe(const Ray_2& ray,bool deselect_all=false)
    {
      ipe::Vector paper_size=get_paper_size();
      Iso_rectangle_2 bbox(0,0,paper_size.x,paper_size.y);
      draw_in_ipe(ray,bbox,deselect_all);
    }

    template<class GT,class TDS>
    void
    draw_in_ipe(const CGAL::Triangulation_2<GT,TDS>& tri,bool deselect_all=false,bool make_grp=true) const
    {
      typedef CGAL::Triangulation_2<GT,TDS> Triangulation;
      typename Triangulation::Finite_edges_iterator first=tri.finite_edges_begin();
      typename Triangulation::Finite_edges_iterator last=tri.finite_edges_end();
      for(typename Triangulation::Finite_edges_iterator it=first;it!=last;++it)
        draw_in_ipe(tri.segment(*it));
      if (make_grp)
        group_selected_objects_();
      if (deselect_all)
        get_IpePage()->deselectAll();
    }

    template<class iterator>
    void
    draw_in_ipe(const iterator begin,const iterator end,bool make_grp=true,bool deselect_all=false) const
    {
      for (iterator it=begin;it!=end;++it)
        draw_in_ipe(*it);
      if (make_grp and ++iterator(begin)!=end)
        group_selected_objects_();
      if (deselect_all)
        get_IpePage()->deselectAll();
    }

    template<class iterator>
    void
    draw_in_ipe(const iterator begin,const iterator end,const Iso_rectangle_2& bbox,bool make_grp=true,bool deselect_all=false,
     typename boost::enable_if<  boost::mpl::or_< boost::is_same<typename std::iterator_traits<iterator>::value_type,Point_2> ,
                                 boost::mpl::or_< boost::is_same<typename std::iterator_traits<iterator>::value_type,Segment_2> ,
                                 boost::mpl::or_< boost::is_same<typename std::iterator_traits<iterator>::value_type,Circle_2> ,
                                 boost::mpl::or_< boost::is_same<typename std::iterator_traits<iterator>::value_type,Circular_arc_2> ,
                                                  boost::is_same<typename std::iterator_traits<iterator>::value_type,Polygon_2>
                                                > > > >
                    >::type* = nullptr) const
    {
      for (iterator it=begin;it!=end;++it)
        draw_in_ipe(*it,bbox);
      if (make_grp and ++iterator(begin)!=end)
        group_selected_objects_();
      if (deselect_all)
        get_IpePage()->deselectAll();
    }

    private:

    ipe::Vector
    get_paper_size() const {
        return data_->iDoc->cascade()->findLayout()->iPaperSize;
    }

    struct Voronoi_from_tri{  //Class using stream to get the voronoi diagram
      std::list<Ray_2> ray_list;
      std::list<Line_2> line_list;
      std::list<Segment_2> seg_list;

      void operator<<(const Ray_2& p){ray_list.push_back(p);}
      void operator<<(const Line_2& p){line_list.push_back(p);}
      void operator<<(const Segment_2& p){seg_list.push_back(p);}

    };

    template <class T,class output_iterator>
    bool
    cast_into_seg(const T& obj,const Iso_rectangle_2& bbox,output_iterator out_it) const{
      CGAL::Object obj_cgal = CGAL::intersection(obj,bbox);
      Segment_2 s;
      bool ret=CGAL::assign(s, obj_cgal);
      if (ret) *out_it++=s;
      return ret;
    }

    //Convert infinite objects into drawable segments
    template<class iterator,class output_iterator>
    void
    cast_into_seg(const iterator first,const iterator end,
                  const Iso_rectangle_2& bbox, output_iterator out_it) const
    {
      for (iterator it=first;it!=end;++it)
        cast_into_seg(*it,bbox,out_it);
    }



    void
    draw_dual_(Voronoi_from_tri& v_recup,const Iso_rectangle_2& bbox,bool makegrp) const
    {
      std::vector<Segment_2> seg_cont;
      //filter degenerate segments
      for(typename std::list<Segment_2>::iterator iteS = v_recup.seg_list.begin();iteS!=v_recup.seg_list.end();){
        typename std::list<Segment_2>::iterator itc=iteS++;
        if (itc->is_degenerate())
          v_recup.seg_list.erase(itc);
      }

      cast_into_seg(v_recup.ray_list.begin(),v_recup.ray_list.end(),bbox,std::back_inserter(seg_cont));//cast rays into segments in bbox
      cast_into_seg(v_recup.line_list.begin(),v_recup.line_list.end(),bbox,std::back_inserter(seg_cont));//cast lines into segments in bbox
      cast_into_seg(v_recup.seg_list.begin(),v_recup.seg_list.end(),bbox,std::back_inserter(seg_cont));//cast lines into segments in bbox
      draw_in_ipe(seg_cont.begin(), seg_cont.end(),makegrp);
    }

    public:
    template<class Triangulation>
    void
    draw_dual_in_ipe(Triangulation& T,const Iso_rectangle_2& bbox,bool makegrp=true,bool deselect_all=false) const
    {
    //~ template<class GT,class TDS>
    //~ void draw_dual_in_ipe(const CGAL::Triangulation_2<GT,TDS>& T,const Iso_rectangle_2& bbox) const{
      Voronoi_from_tri v_recup;
      T.draw_dual(v_recup);
      draw_dual_(v_recup,bbox,makegrp);
      if (deselect_all) get_IpePage()->deselectAll();
    }

    template<class Triangulation>
    void
    draw_skeleton_in_ipe(Triangulation& T,const Iso_rectangle_2& bbox,
                         bool makegrp=true,bool deselect_all=false) const
    {
    //~ template<class GT,class TDS>
    //~ void draw_skeleton_in_ipe(const CGAL::Triangulation_2<GT,TDS>& T,const Iso_rectangle_2& bbox) const{
      Voronoi_from_tri v_recup;
      T.draw_skeleton(v_recup);
      draw_dual_(v_recup,bbox,makegrp);
      if (deselect_all) get_IpePage()->deselectAll();
    }

    //Circle restriction
  private:
    inline double get_theta(const Point_2& point, const Circle_2& circle) const {
      return atan2(CGAL::to_double(point.y()-circle.center().y()),CGAL::to_double(point.x()-circle.center().x()));
    }

    //SK objects
    typedef CGAL::Exact_circular_kernel_2 SK;
    typedef SK::Circle_2                  Exact_circle_2;
    typedef SK::Point_2                   Exact_point_2;
    typedef SK::Circular_arc_point_2      Circular_arc_point_2;



    //s and t are given such that if center of exact_circle is inside bbox then turn COUNTERCLOCKWISE
    Circular_arc_2
    build_arc(const Exact_circle_2& exact_circle,const SK::Circular_arc_point_2& s,
              const SK::Circular_arc_point_2& t,bool sign_known=false) const
    {
      Point_2 sd=Point_2(CGAL::to_double(s.x()),CGAL::to_double(s.y()));
      Point_2 td=Point_2(CGAL::to_double(t.x()),CGAL::to_double(t.y()));
      Point_2 center(CGAL::to_double(exact_circle.center().x()),CGAL::to_double(exact_circle.center().y()));
      CGAL::Cartesian_converter<SK,Kernel> conv;
      Circle_2 approx_circle=conv(exact_circle);
      if (!sign_known){
        CGAL::Sign sign = (CGAL::orientation(sd,td,center)==CGAL::LEFT_TURN)?CGAL::POSITIVE:CGAL::NEGATIVE;
        return std::make_tuple(approx_circle,sd,td,sign);
      }
      return std::make_tuple(approx_circle,sd,td,CGAL::POSITIVE);
    }

    void
    get_pair_indices(int* array,int* pair) const {
      for (int i=0;i<8;++i)
        if (array[i]!=-1) *pair++=i;
    }

    void
    set_done(int* array,int index) const {
      for (int i =0;i<8;++i)
        if (array[i]==index) array[i]=-1;
    }


    bool
    indices_are_on_opposite_side(int* array) const {
        if (array[0]!=-1 || array[5]!=-1)
          return array[2]!=-1 || array[7]!=-1;
        if (array[1]!=-1 || array[6]!=-1)
          return array[3]!=-1 || array[4]!=-1;
        return false;
    }

    int
    count_points(int* array) const {
      int ret=0;
      for (int i =0;i<8;++i)
        if (array[i]!=-1) ++ret;
      return ret;
    }

  public:
    //
    // .-----7---------2-------.
    // |                       |
    // 3                       6
    // |     indices           |
    // 4                       1
    // |                       |
    // .-----0---------5-------.
    template <class Output_iterator>
    void
    restrict_circle_to_bbox(const Circle_2& approx_circle,
                            const Iso_rectangle_2& bbox,Output_iterator out) const
    {
      CGAL::Cartesian_converter<Kernel,SK> conv;
      Exact_circle_2 exact_circle=conv(approx_circle);

      SK::Intersect_2 inter=SK().intersect_2_object();
      std::vector< std::pair<Circular_arc_point_2,unsigned> > points;
      points.reserve(8);

      std::vector<CGAL::Object> ints;
      ints.reserve(2);
      std::pair<Circular_arc_point_2,unsigned> tmp_pt;

      int indices[8]={-1,-1,-1,-1,-1,-1,-1,-1};

      for (unsigned i=0;i!=4;++i){
        ints.clear();
        SK::Segment_2 S(conv(bbox[i]),conv(bbox[(i+1)%4]));
        inter(exact_circle,SK::Line_arc_2(S),std::back_inserter(ints));
        unsigned index=0;
        bool ok=true;
        switch (ints.size()){
          case 1:
            ok=CGAL::assign(tmp_pt,ints[0]);
            CGAL_assertion(ok); CGAL_USE(ok);
            points.push_back(tmp_pt);
            index=points.size()-1;
            indices[i]=index;
            indices[(i+1)%4+4]=index;
            break;
          case 2:
            int right_ind=i<2?0:1;
            ok=CGAL::assign(tmp_pt,ints[right_ind]);
            CGAL_assertion(ok); CGAL_USE(ok);
            points.push_back(tmp_pt);
            index=points.size()-1;
            indices[i]=index;
            ok=CGAL::assign(tmp_pt,ints[(right_ind+1)%2]);
            CGAL_assertion(ok); CGAL_USE(ok);
            points.push_back(tmp_pt);
            index=points.size()-1;
            indices[(i+1)%4+4]=index;
            break;
        }
      }

      //corner case
      for (unsigned i=0;i!=4;++i){
        if (indices[i]!=-1 && indices[i+4]!=-1){
          *out++=build_arc(exact_circle,points[ indices[i+4] ].first,points[ indices[i] ].first);
          if (points[ indices[i] ].second==1) set_done(indices,indices[i]);
          else indices[i]=-1;
          if (points[ indices[i+4] ].second==1) set_done(indices,indices[i+4]);
          else indices[i+4]=-1;
        }
      }
      int rem_pt=count_points(indices);
      if (rem_pt==4){
        //double opposite
        if (indices[0]!=-1){
          *out++=build_arc(exact_circle,points[ indices[7] ].first,points[ indices[0] ].first);
          if (indices[7]!=indices[2] && indices[0]!=indices[5])
            *out++=build_arc(exact_circle,points[ indices[5] ].first,points[ indices[2] ].first);
        }
        else{
          *out++=build_arc(exact_circle,points[ indices[6] ].first,points[ indices[3] ].first);
          if (indices[6]!=indices[1] && indices[3]!=indices[4])
            *out++=build_arc(exact_circle,points[ indices[4] ].first,points[ indices[1] ].first);
        }
        return;
      }

      if (rem_pt==2){
        int pair[2];
        get_pair_indices(indices,pair);
        if (!indices_are_on_opposite_side(indices))
          *out++=build_arc(exact_circle,points[ indices[pair[1]] ].first,points[ indices[pair[0]] ].first,true);
        else
          *out++=build_arc(exact_circle,points[ indices[pair[1]] ].first,points[ indices[pair[0]] ].first);
        return;
      }
      CGAL_assertion (rem_pt==0);
    }
  };


  //definition
  template <class Kernel,int nbf>
  template< class multi_output_iterator >
  bool
  Ipelet_base<Kernel,nbf>::read_one_active_object(ipe::Object* object,
                                                  multi_output_iterator it_out) const
  {
    if (object->asGroup()){
      bool deselect_grp=false;
      for (ipe::Group::const_iterator it=object->asGroup()->begin();
                                    it!=object->asGroup()->end();++it)
      {
        ipe::Object *obj = (*it)->clone();
        obj->setMatrix(obj->matrix()*object->matrix());
        bool cur=read_one_active_object(obj,it_out);
        deselect_grp=deselect_grp || cur;
      }
      return deselect_grp;
    }

    //detect Points
    if( object->asReference() ){
      if ( !CGAL::Is_in_tuple<Point_2,typename multi_output_iterator::Value_type_tuple>::value )
        return true;
      it_out++=to_point_2(object);
      return false;
    }
    bool to_deselect=true;
    //Test one function for segments, circles, circle arcs and polygons
    if (object->asPath()){
      //iterate on each subpath
      to_deselect=false;
      for (int i=0;i<object->asPath()->shape().countSubPaths();++i){
        if(object->asPath()-> shape().subPath(i)->asCurve()){
          std::list<Segment_2> seg_list;
          bool is_polygon=object-> asPath() -> shape().subPath(i)->closed();

          const ipe::Curve* SSP_ipe = object -> asPath() -> shape().subPath(i) -> asCurve();
          for(int j=0; j< SSP_ipe->countSegments();++j){
            if (SSP_ipe -> segment(j).type()==ipe::CurveSegment::ESegment){
              //TODO depending on if current_it  is  a polygon or not do not do the same thing
              seg_list.push_back(Segment_2(segment_endpoint(SSP_ipe -> segment(j),object -> asPath(),0),
                                  segment_endpoint(SSP_ipe -> segment(j),object -> asPath(),1)));
            }
            else{
              //retrieve circle arcs
              if(SSP_ipe -> segment(j).type()==ipe::CurveSegment::EArc &&
                 is_only_rotated_or_scaled(object->asPath()->matrix()))
              {//retreve circle arcs
                if ( !CGAL::Is_in_tuple<Circular_arc_2,typename multi_output_iterator::Value_type_tuple>::value ){
                  to_deselect=true;
                  continue;
                }
                is_polygon=false;
                ipe::Path* obj_ipe=object->asPath();
                ipe::Matrix mat=obj_ipe->matrix() * SSP_ipe->segment(j).matrix();
                ipe::Vector ipe_center=ipe::Vector(mat.a[4],mat.a[5]);
                ipe::Vector ipe_first=obj_ipe->matrix() * SSP_ipe->segment(j).cp(0);
                ipe::Vector ipe_last=obj_ipe->matrix() * SSP_ipe->segment(j).last();
                //TODO Check how object is defined
                Circular_arc_2 arc(
                            Circle_2(Point_2(ipe_center.x,ipe_center.y),(ipe_first-ipe_center).len()*(ipe_first-ipe_center).len()),
                            Point_2(ipe_first.x,ipe_first.y),
                            Point_2(ipe_last.x,ipe_last.y),
                            mat.a[0]*mat.a[3]-mat.a[1]*mat.a[2]<0?CGAL::CLOCKWISE:CGAL::COUNTERCLOCKWISE
                          );
                it_out++=arc;
              }
              else
                to_deselect=true;
            }
          }
          if (object->asPath() -> shape().subPath(i)->closed() &&
              (SSP_ipe -> segment(0).cp(0)-SSP_ipe -> segment(SSP_ipe->countSegments()-1).cp(1)).len()!=0 ){ //for close polygon, seg
            seg_list.push_back( Segment_2(
                              segment_endpoint(SSP_ipe -> segment(SSP_ipe->countSegments()-1),object-> asPath(),1),
                              segment_endpoint(SSP_ipe -> segment(0),object-> asPath(),0)
                              ));
          }
          //~ if (seg_list.empty())
            //~ to_deselect=true;

          if (is_polygon){
            if (  !CGAL::Is_in_tuple<Polygon_2,typename multi_output_iterator::Value_type_tuple>::value  )
              to_deselect=true;
            else{
              Polygon_2 polygon;
              typename std::list<Segment_2>::iterator its=seg_list.begin();
              for (;its!=seg_list.end();++its)
                polygon.push_back(its->source());
              it_out++=polygon;
            }
          }
          else{
            if (  !CGAL::Is_in_tuple<Segment_2,typename multi_output_iterator::Value_type_tuple>::value  )
              to_deselect=true;
            else{
              for (typename std::list<Segment_2>::iterator its=seg_list.begin();its!=seg_list.end();++its)
                it_out++=*its;
            }
          }
        }
        else
          if (is_IPE_circle(object,i)){
            if (  !CGAL::Is_in_tuple<Circle_2,typename multi_output_iterator::Value_type_tuple>::value )
              to_deselect=true;
            else{
              Circle_2 C=to_circle_2(object -> asPath(),i);
              it_out++=C;
            }
          }
          else
            to_deselect=true; // avoid deleting no handled objects
      }
    }
    return to_deselect;
  }

} //CGAL

#define CGAL_IPELET(T) IPELET_DECLARE ipe::Ipelet *newIpelet(){ return new T; }

#endif //CGAL_IPELET_BASE_H
