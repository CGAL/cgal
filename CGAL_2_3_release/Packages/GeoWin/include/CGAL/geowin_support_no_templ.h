// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 2001, January 30
//
// file          : include/CGAL/geowin_support_no_templ.h
// package       : GeoWin (1.2.2)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.2.2
// revision_date : 30 January 2001 
// author(s)     : Matthias Baesken, Ulrike Bartuschka, Stefan Naeher
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ======================================================================


#ifndef CGAL_GEOWIN_SUPPORT_NO_TEMPL_H
#define CGAL_GEOWIN_SUPPORT_NO_TEMPL_H


leda_point convert_to_leda(const CGAL::Point_2<CGAL::Cartesian<double> >& obj)
{
  double x = obj.x();
  double y = obj.y();
  leda_point p(x,y);
  return p;
}

leda_segment convert_to_leda(const CGAL::Segment_2<CGAL::Cartesian<double> >& obj)
{
  CGAL::Point_2<CGAL::Cartesian<double> > p1=obj.source();
  CGAL::Point_2<CGAL::Cartesian<double> > p2=obj.target();
  leda_segment seg(p1.x(), p1.y(), p2.x(), p2.y() ); 
  return seg;
}

leda_circle convert_to_leda(const CGAL::Circle_2<CGAL::Cartesian<double> >& c)
{
 CGAL::Point_2<CGAL::Cartesian<double> > p1=c.center();
 double   radius= ::sqrt(CGAL::to_double(c.squared_radius()));

 leda_point lp(p1.x(), p1.y());
 leda_circle lc(lp,radius);
 return lc;
}

leda_line convert_to_leda(const CGAL::Line_2<CGAL::Cartesian<double> >& l)
{
 CGAL::Point_2<CGAL::Cartesian<double> > p1=l.point(1);
 CGAL::Point_2<CGAL::Cartesian<double> > p2=l.point(2);
 leda_point lp1(p1.x(), p1.y());
 leda_point lp2(p2.x(), p2.y());
 leda_line lc(lp1,lp2);
 return lc;
}

leda_ray convert_to_leda(const CGAL::Ray_2<CGAL::Cartesian<double> >& r)
{
 CGAL::Point_2<CGAL::Cartesian<double> > p1=r.source();
 CGAL::Point_2<CGAL::Cartesian<double> > p2=r.point(1);

 leda_point lp1(p1.x(), p1.y());
 leda_point lp2(p2.x(), p2.y());
 leda_ray rc(lp1,lp2);
 return rc;
}

leda_polygon convert_to_leda(const CGALPolygon& p)
{
  CGALPolygon::Vertex_const_iterator it=p.vertices_begin();
  CGALPolygon::Vertex_const_iterator st=p.vertices_end();

 leda_list<leda_point> lp;

 while (it != st) {  lp.append(leda_point((*it).x(),(*it).y())); it++; }
 
 leda_polygon ph(lp);

 return ph;
}

leda_polygon convert_to_leda(const CGAL::Triangle_2<CGAL::Cartesian<double> >& t)
{
 CGAL::Point_2<CGAL::Cartesian<double> > p1= t[1];
 CGAL::Point_2<CGAL::Cartesian<double> > p2= t[2];
 CGAL::Point_2<CGAL::Cartesian<double> > p3= t[3];
 leda_point lp1(p1.x(), p1.y());
 leda_point lp2(p2.x(), p2.y());
 leda_point lp3(p3.x(), p3.y());

 leda_list<leda_point> L;
 if (t.orientation()== CGAL::LEFTTURN)
 {
  L.append(lp1); L.append(lp2); L.append(lp3);
 }
 else
 {
  L.push(lp1); L.push(lp2); L.push(lp3);  
 }
  
 return leda_polygon(L); 
}

leda_rectangle convert_to_leda(const CGAL::Iso_rectangle_2<CGAL::Cartesian<double> >& t)
{
 CGAL::Point_2<CGAL::Cartesian<double> > p1= t.min();
 CGAL::Point_2<CGAL::Cartesian<double> > p2= t.max();
 leda_point lp1(p1.x(), p1.y());
 leda_point lp2(p2.x(), p2.y());
 return leda_rectangle(lp1,lp2); 
}

ps_file& operator<<(ps_file& F,const leda_d3_point& obj);

leda_d3_point convert_to_leda(const CGAL::Point_3<CGAL::Cartesian<double> >& obj)
{
  leda_d3_point p(obj.x(), obj.y(), obj.z() );
  return p;
}


ps_file& operator<<(ps_file& F,const CGAL::Point_2<CGAL::Cartesian<double> >& o) { F << convert_to_leda(o); return F; }

ps_file& operator<<(ps_file& F,const CGAL::Segment_2<CGAL::Cartesian<double> >& o) { F << convert_to_leda(o); return F; }

ps_file& operator<<(ps_file& F,const CGAL::Circle_2<CGAL::Cartesian<double> >& o) { F << convert_to_leda(o); return F; }

ps_file& operator<<(ps_file& F,const CGAL::Line_2<CGAL::Cartesian<double> >& o) { F << convert_to_leda(o); return F; }

ps_file& operator<<(ps_file& F,const CGAL::Ray_2<CGAL::Cartesian<double> >& o) { F << convert_to_leda(o); return F; }

ps_file& operator<<(ps_file& F,const CGALPolygon& o) 
{ F << convert_to_leda(o); return F; }

ps_file& operator<<(ps_file& F,const CGAL::Triangle_2<CGAL::Cartesian<double> >& o) { F << convert_to_leda(o); return F; }

ps_file& operator<<(ps_file& F,const CGAL::Iso_rectangle_2<CGAL::Cartesian<double> >& o) 
{ F << convert_to_leda(o); return F; }

ps_file& operator<<(ps_file& F,const CGAL::Point_3<CGAL::Cartesian<double> >& o) { F << convert_to_leda(o); return F; }

#include <math.h>
#include <ctype.h>

GEOWIN_BEGIN_NAMESPACE

bool geowin_IntersectsBox(const CGAL::Point_2<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f);

bool geowin_IntersectsBox(const CGAL::Segment_2<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f);

bool geowin_IntersectsBox(const CGAL::Circle_2<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f);

bool geowin_IntersectsBox(const CGAL::Line_2<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f);

bool geowin_IntersectsBox(const CGAL::Ray_2<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f);

bool geowin_IntersectsBox(const CGALPolygon& obj, double x1,double y1,double x2, double y2,bool f);

bool geowin_IntersectsBox(const CGAL::Triangle_2<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f);

bool geowin_IntersectsBox(const CGAL::Iso_rectangle_2<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f);

bool geowin_IntersectsBox(const CGAL::Point_3<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f);


// ----------------------------------------------
// templates for the functions for the containers
// ----------------------------------------------

template<class T>
void geowin_redraw_fcn(const std::list<T>& L, 
		    leda_window& w, leda_color c, leda_color, 
		    double x1, double y1, double x2, double y2)
{
  std::list<T>::const_iterator it   = L.begin(), stop = L.end();
  while( it != stop )
    { if ( geowin_IntersectsBox(*it, x1, y1, x2, y2,true) ) {  w.set_color(c);   w << *it; } it++; }
}

GEOWIN_END_NAMESPACE


#include <LEDA/geowin.h>
#include <LEDA/geowin_init.h>

GEOWIN_BEGIN_NAMESPACE
// new geowin_init_default_type ...

template<class CONT>
class geo_scene_traits {
 leda_string description;
 
 public:
 
 typedef CONT                       CONTAINER;
 typedef typename CONT::value_type  MYTYPE;
 
 leda_string get_name() { return description; }

 leda_string (*geowin_info_fcn)(const CONTAINER& L);
 bool (*geowin_IntersectsBox)(const MYTYPE& obj, double x1,double y1,double x2, double y2,bool f);
 void (*geowin_BoundingBox)(const MYTYPE& obj, double& x1, double& x2,double& y1, double& y2);
 void (*geowin_Translate)(MYTYPE& obj, double dx, double dy);
 void (*geowin_Rotate)(MYTYPE& obj, double dx, double dy,double a);
 void (*geowin_generate_objects)(GeoWin& gw, CONTAINER& L);
 
 geo_scene_traits(leda_string (*f1)(const CONTAINER& ), \
   bool (*f2)(const MYTYPE&, double, double, double, double,bool ), \
   void (*f3)(const MYTYPE&, double&, double&, double&, double& ), \
   void (*f4)(MYTYPE&, double, double), void (*f5)(MYTYPE&, double, double, double), \
   void (*f6)(GeoWin&, CONTAINER&), leda_string scene_type_name)
 {
   description = scene_type_name;
   geowin_info_fcn = f1;
   geowin_IntersectsBox = f2;
   geowin_BoundingBox = f3;
   geowin_Translate = f4;
   geowin_Rotate = f5;
   geowin_generate_objects = f6;
 }
 
};

#if (__LEDA__ < 410)

template<class CONT> 
void geowin_init_default_type(geo_scene_traits<CONT> tr)
{ 
  CONT* t;
  GeoEditScene<CONT>* sc = make_edit_prototype(t, tr.get_name());

  sc->set_redraw_fcn(0);
  sc->set_info_fcn(tr.geowin_info_fcn);
  sc->set_box_intersection_fcn(tr.geowin_IntersectsBox);
  sc->set_get_bounding_box_fcn(tr.geowin_BoundingBox);
  sc->set_move_fcn(tr.geowin_Translate);
  sc->set_rotate_fcn(tr.geowin_Rotate);
  sc->set_generate_fcn(tr.geowin_generate_objects);
}

#endif

// 2d Points ...

const char* leda_tname(CGAL::Point_2<CGAL::Cartesian<double> >* p) {  return "CGALPoint"; }

bool geowin_IntersectsBox(const CGAL::Point_2<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f)
{
  double xw= obj.x();
  double yw= obj.y();  
  
  if (x1<=xw && x2>=xw && y1<=yw && y2>=yw) return true;
  return false;
}

void geowin_BoundingBox(const CGAL::Point_2<CGAL::Cartesian<double> >& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_2 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

void geowin_Translate(CGAL::Point_2<CGAL::Cartesian<double> >& obj, double dx, double dy)
{
  CGAL::Vector_2<CGAL::Cartesian<double> > vec(dx,dy);
  obj = obj + vec;
}

void geowin_Rotate(CGAL::Point_2<CGAL::Cartesian<double> >& obj, double x, double y, double a)
{
  leda_point p2(obj.x(), obj.y());
  p2 = p2.rotate(leda_point(x,y), a);
  obj = CGAL::Point_2<CGAL::Cartesian<double> >(p2.xcoord(), p2.ycoord()); 
}

// functions for the container

leda_string geowin_info_fcn(const std::list<CGAL::Point_2<CGAL::Cartesian<double> > >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-point");  return str;
}

void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Point_2<CGAL::Cartesian<double> > >& L)
{ 
  leda_list<leda_point> H;
  geowin_generate_objects(gw,H);

  //convert the contents
  CGAL::Point_2<CGAL::Cartesian<double> > p;
  leda_point mp;

  forall(mp,H){
   p= CGAL::Point_2<CGAL::Cartesian<double> >(mp.xcoord(), mp.ycoord());
   L.push_front(p);
  }
}


// Segments

const char* leda_tname(CGAL::Segment_2<CGAL::Cartesian<double> >* p) {  return "CGALSegment"; }

bool geowin_IntersectsBox(const CGAL::Segment_2<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f)
{
  CGAL::Point_2<CGAL::Cartesian<double> > p1,p2;
  p1=obj.source();
  p2=obj.target();
  leda_segment seg(p1.x(),p1.y(),p2.x(),p2.y() );
  
  return geowin_IntersectsBox(seg,x1,y1,x2,y2,f);
}

void geowin_BoundingBox(const CGAL::Segment_2<CGAL::Cartesian<double> >& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_2 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

void geowin_Translate(CGAL::Segment_2<CGAL::Cartesian<double> >& obj, double dx, double dy)
{
  CGAL::Vector_2<CGAL::Cartesian<double> > vec;
  vec= CGAL::Vector_2<CGAL::Cartesian<double> >(dx, dy); 
  CGAL::Aff_transformation_2<CGAL::Cartesian<double> > translate(CGAL::TRANSLATION, vec);
  obj = obj.transform(translate);
}

void geowin_Rotate(CGAL::Segment_2<CGAL::Cartesian<double> >& obj, double x, double y, double a)
{
  CGAL::Point_2<CGAL::Cartesian<double> > p1,p2;
  p1=obj.source();
  p2=obj.target();
  
  leda_segment hlp(p1.x(), p1.y(), p2.x(), p2.y() );

  hlp = hlp.rotate(leda_point(x,y), a);
  p1= CGAL::Point_2<CGAL::Cartesian<double> >( (hlp.source()).xcoord(), (hlp.source()).ycoord() );
  p2= CGAL::Point_2<CGAL::Cartesian<double> >( (hlp.target()).xcoord(), (hlp.target()).ycoord() );
  obj = CGAL::Segment_2<CGAL::Cartesian<double> >(p1,p2);
}

// functions for the container

leda_string geowin_info_fcn(const std::list<CGAL::Segment_2<CGAL::Cartesian<double> > >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-segment");  return str;
}

void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Segment_2<CGAL::Cartesian<double> > >& L)
{
  leda_list<leda_segment> H;
  geowin_generate_objects(gw,H);

  //convert the contents
  CGAL::Segment_2<CGAL::Cartesian<double> > p;
  leda_segment mp;

  forall(mp,H){
   CGAL::Point_2<CGAL::Cartesian<double> > pa(mp.source().xcoord(), mp.source().ycoord());
   CGAL::Point_2<CGAL::Cartesian<double> > pb(mp.target().xcoord(), mp.target().ycoord());   
   
   p= CGAL::Segment_2<CGAL::Cartesian<double> >(pa,pb);
   L.push_front(p);
  }
}

// Circles

void convert_from_leda(const leda_circle& c, CGAL::Circle_2<CGAL::Cartesian<double> >& cir)
{
 leda_point lp=c.center();
 double sr= c.radius()*c.radius();

 CGAL::Point_2<CGAL::Cartesian<double> > p1(lp.xcoord(), lp.ycoord());
 cir = CGAL::Circle_2<CGAL::Cartesian<double> >(p1, sr);
}

const char* leda_tname(CGAL::Circle_2<CGAL::Cartesian<double> >* p) {  return "CGALCircle"; }

bool geowin_IntersectsBox(const CGAL::Circle_2<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f)
{  
  return geowin_IntersectsBox(convert_to_leda(obj),x1,y1,x2,y2,f); 
}

void geowin_BoundingBox(const CGAL::Circle_2<CGAL::Cartesian<double> >& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_2 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

void geowin_Translate(CGAL::Circle_2<CGAL::Cartesian<double> >& obj, double dx, double dy)
{
  CGAL::Vector_2<CGAL::Cartesian<double> > vec;
  vec= CGAL::Vector_2<CGAL::Cartesian<double> >(dx, dy); 
  CGAL::Aff_transformation_2<CGAL::Cartesian<double> > translate(CGAL::TRANSLATION, vec);
  
  obj = obj.orthogonal_transform(translate);
}

void geowin_Rotate(CGAL::Circle_2<CGAL::Cartesian<double> >& obj, double x, double y, double a)
{  
  leda_circle hlp=convert_to_leda(obj);
  hlp = hlp.rotate(leda_point(x,y), a);
  convert_from_leda(hlp,obj);
}


leda_string geowin_info_fcn(const std::list<CGAL::Circle_2<CGAL::Cartesian<double> > >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-circle");  return str;
}

void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Circle_2<CGAL::Cartesian<double> > >& L)
{
  leda_list<leda_circle> H;
  geowin_generate_objects(gw,H);

  //convert the contents
  CGAL::Circle_2<CGAL::Cartesian<double> > c;
  leda_circle mp;

  forall(mp,H){
   convert_from_leda(mp,c);
   L.push_front(c);
  }
}


// Lines

void convert_from_leda(const leda_line& l, CGAL::Line_2<CGAL::Cartesian<double> >& lc)
{
 leda_point lp1=l.point1();
 leda_point lp2=l.point2();

 CGAL::Point_2<CGAL::Cartesian<double> > p1(lp1.xcoord(), lp1.ycoord());
 CGAL::Point_2<CGAL::Cartesian<double> > p2(lp2.xcoord(), lp2.ycoord());
 lc = CGAL::Line_2<CGAL::Cartesian<double> >(p1,p2);
}

const char* leda_tname(CGAL::Line_2<CGAL::Cartesian<double> >* p) {  return "CGALLine"; }

bool geowin_IntersectsBox(const CGAL::Line_2<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f)
{  
  return geowin_IntersectsBox(convert_to_leda(obj),x1,y1,x2,y2,f); 
}

void geowin_BoundingBox(const CGAL::Line_2<CGAL::Cartesian<double> >& obj, double& x1, double& x2,
	 double& y1, double& y2)
{
  geowin_BoundingBox(convert_to_leda(obj),x1,x2,y1,y2);
}

void geowin_Translate(CGAL::Line_2<CGAL::Cartesian<double> >& obj, double dx, double dy)
{
  CGAL::Vector_2<CGAL::Cartesian<double> > vec;
  vec= CGAL::Vector_2<CGAL::Cartesian<double> >(dx, dy); 
  CGAL::Aff_transformation_2<CGAL::Cartesian<double> > translate(CGAL::TRANSLATION, vec);
  
  obj = obj.transform(translate);
}

void geowin_Rotate(CGAL::Line_2<CGAL::Cartesian<double> >& obj, double x, double y, double a)
{  
  leda_line hlp=convert_to_leda(obj);
  hlp = hlp.rotate(leda_point(x,y), a);
  convert_from_leda(hlp,obj);
}

// functions for the container

leda_string geowin_info_fcn(const std::list<CGAL::Line_2<CGAL::Cartesian<double> > >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-Line");  return str;
}

void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Line_2<CGAL::Cartesian<double> > >& L)
{
  leda_list<leda_line> H;
  geowin_generate_objects(gw,H);

  //convert the contents
  CGAL::Line_2<CGAL::Cartesian<double> > obj;
  leda_line mp;

  forall(mp,H){
   convert_from_leda(mp,obj);
   L.push_front(obj);
  }
}


// Rays

void convert_from_leda(const leda_ray& r, CGAL::Ray_2<CGAL::Cartesian<double> >& rc)
{
 leda_point lp1=r.point1();
 leda_point lp2=r.point2();

 CGAL::Point_2<CGAL::Cartesian<double> > p1(lp1.xcoord(), lp1.ycoord());
 CGAL::Point_2<CGAL::Cartesian<double> > p2(lp2.xcoord(), lp2.ycoord());
 rc = CGAL::Ray_2<CGAL::Cartesian<double> >(p1,p2);
}

const char* leda_tname(CGAL::Ray_2<CGAL::Cartesian<double> >* p) {  return "CGALRay"; }

bool geowin_IntersectsBox(const CGAL::Ray_2<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f)
{  
  return geowin_IntersectsBox(convert_to_leda(obj),x1,y1,x2,y2,f); 
}

void geowin_BoundingBox(const CGAL::Ray_2<CGAL::Cartesian<double> >& obj, double& x1, double& x2,
	 double& y1, double& y2)
{
  geowin_BoundingBox(convert_to_leda(obj),x1,x2,y1,y2);
}

void geowin_Translate(CGAL::Ray_2<CGAL::Cartesian<double> >& obj, double dx, double dy)
{
  CGAL::Vector_2<CGAL::Cartesian<double> > vec;
  vec= CGAL::Vector_2<CGAL::Cartesian<double> >(dx, dy); 
  CGAL::Aff_transformation_2<CGAL::Cartesian<double> > translate(CGAL::TRANSLATION, vec);
  
  obj = obj.transform(translate); 
}

void geowin_Rotate(CGAL::Ray_2<CGAL::Cartesian<double> >& obj, double x, double y, double a)
{  
  leda_ray hlp=convert_to_leda(obj);
  hlp = hlp.rotate(leda_point(x,y), a);
  convert_from_leda(hlp,obj);
}

// functions for the container

leda_string geowin_info_fcn(const std::list<CGAL::Ray_2<CGAL::Cartesian<double> > >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-Ray");  return str;
}

void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Ray_2<CGAL::Cartesian<double> > >& L)
{
  leda_list<leda_ray> H;
  geowin_generate_objects(gw,H);

  //convert the contents
  CGAL::Ray_2<CGAL::Cartesian<double> > obj;
  leda_ray mp;

  forall(mp,H){
   convert_from_leda(mp,obj);
   L.push_front(obj);
  }
}


// Triangles

void convert_from_leda(const leda_polygon& p, CGAL::Triangle_2<CGAL::Cartesian<double> >& tr)
{
 leda_list<leda_point> L=p.vertices();
 leda_point lp1,lp2,lp3;
 lp1= L.pop(); lp2= L.pop(); lp3= L.pop();
 CGAL::Point_2<CGAL::Cartesian<double> > p1(lp1.xcoord(), lp1.ycoord());
 CGAL::Point_2<CGAL::Cartesian<double> > p2(lp2.xcoord(), lp2.ycoord());
 CGAL::Point_2<CGAL::Cartesian<double> > p3(lp3.xcoord(), lp3.ycoord());

 tr= CGAL::Triangle_2<CGAL::Cartesian<double> >(p1,p2,p3); 
}

const char* leda_tname(CGAL::Triangle_2<CGAL::Cartesian<double> >* t) {  return "CGALTriangle"; }

bool geowin_IntersectsBox(const CGAL::Triangle_2<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f)
{  
  return geowin_IntersectsBox(convert_to_leda(obj),x1,y1,x2,y2,f); 
}

void geowin_BoundingBox(const CGAL::Triangle_2<CGAL::Cartesian<double> >& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_2 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

void geowin_Translate(CGAL::Triangle_2<CGAL::Cartesian<double> >& obj, double dx, double dy)
{
  CGAL::Vector_2<CGAL::Cartesian<double> > vec;
  vec= CGAL::Vector_2<CGAL::Cartesian<double> >(dx, dy); 
  CGAL::Aff_transformation_2<CGAL::Cartesian<double> > translate(CGAL::TRANSLATION, vec);
  
  obj = obj.transform(translate);
}

void geowin_Rotate(CGAL::Triangle_2<CGAL::Cartesian<double> >& obj, double x, double y, double a)
{  
  leda_polygon hlp=convert_to_leda(obj);
  hlp = hlp.rotate(leda_point(x,y), a);
  convert_from_leda(hlp,obj);
}

leda_string geowin_info_fcn(const std::list<CGAL::Triangle_2<CGAL::Cartesian<double> > >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-triangle");  return str;
}

void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Triangle_2<CGAL::Cartesian<double> > >& L)
{
}


// Rectangles

void convert_from_leda(const leda_rectangle& r, CGAL::Iso_rectangle_2<CGAL::Cartesian<double> >& rc)
{
 leda_point lp1= r.lower_left();
 leda_point lp2= r.upper_right();
 CGAL::Point_2<CGAL::Cartesian<double> > p1(lp1.xcoord(),lp1.ycoord());
 CGAL::Point_2<CGAL::Cartesian<double> > p2(lp2.xcoord(),lp2.ycoord());

 rc = CGAL::Iso_rectangle_2<CGAL::Cartesian<double> >(p1,p2); 
}

const char* leda_tname(CGAL::Iso_rectangle_2<CGAL::Cartesian<double> >* t) {  return "CGALRectangle"; }

bool geowin_IntersectsBox(const CGAL::Iso_rectangle_2<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f)
{  
  CGAL::Point_2<CGAL::Cartesian<double> > pa;
  pa = CGAL::Point_2<CGAL::Cartesian<double> >(x1, y1);
  CGAL::Point_2<CGAL::Cartesian<double> > pb;
  pb = CGAL::Point_2<CGAL::Cartesian<double> >(x2, y2);
  
  CGAL::Iso_rectangle_2<CGAL::Cartesian<double> > r2(pa,pb);
  return CGAL::do_intersect(obj,r2);  
}

void geowin_BoundingBox(const CGAL::Iso_rectangle_2<CGAL::Cartesian<double> >& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_2 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

void geowin_Translate(CGAL::Iso_rectangle_2<CGAL::Cartesian<double> >& obj, double dx, double dy)
{
  CGAL::Vector_2<CGAL::Cartesian<double> > vec;
  vec= CGAL::Vector_2<CGAL::Cartesian<double> >(dx,dy); 
  CGAL::Aff_transformation_2<CGAL::Cartesian<double> > translate(CGAL::TRANSLATION, vec);
  
  obj = obj.transform(translate);
}

void geowin_Rotate(CGAL::Iso_rectangle_2<CGAL::Cartesian<double> >& obj, double x, double y, double a)
{  
}

leda_string geowin_info_fcn(const std::list<CGAL::Iso_rectangle_2<CGAL::Cartesian<double> > >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-rectangle");  return str;
}

void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Iso_rectangle_2<CGAL::Cartesian<double> > >& L)
{
}

// Polygons
void convert_from_leda(const leda_polygon& p, CGALPolygon& rc)
{
 leda_list<leda_point> pl= p.vertices();
 std::list<CGALPoint> sl;

 leda_point akt;
 forall(akt,pl) sl.push_back(CGALPoint(akt.xcoord(), akt.ycoord()));
 
 rc = CGALPolygon(sl.begin(),sl.end());
}

const char* leda_tname(CGALPolygon* p) 
{  return "CGALPolygon"; }

leda_window& operator >> (leda_window& w, CGALPolygon& obj)
{
  leda_polygon p;
  w >> p;
  convert_from_leda(p, obj);
  return w;
}

bool geowin_IntersectsBox(const CGALPolygon& obj, double x1,double y1,double x2, double y2,bool f)
{  
  return geowin_IntersectsBox(convert_to_leda(obj),x1,y1,x2,y2,f); 
}

void geowin_BoundingBox(const CGALPolygon& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_2 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

void geowin_Translate(CGALPolygon& obj, double dx, double dy)
{
  leda_polygon hlp=convert_to_leda(obj);
  hlp = hlp.translate(dx,dy);
  convert_from_leda(hlp,obj);   
}

void geowin_Rotate(CGALPolygon& obj, double x, double y, double a)
{  
  leda_polygon hlp=convert_to_leda(obj);
  hlp = hlp.rotate(leda_point(x,y), a);
  convert_from_leda(hlp,obj);
}

leda_string geowin_info_fcn(const std::list<CGALPolygon>& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-Polygons");  return str;
}

void geowin_generate_objects(GeoWin& gw, std::list<CGALPolygon>& L)
{
  leda_list<leda_polygon> H;
  geowin_generate_objects(gw,H);

  //convert the contents
  CGALPolygon obj;
  leda_polygon mp;

  forall(mp,H){
   convert_from_leda(mp,obj);
   L.push_front(obj);
  }
}

// 3d Points ...

void convert_from_leda(const leda_d3_point& obj, CGAL::Point_3<CGAL::Cartesian<double> >& p)
{
  p = CGAL::Point_3<CGAL::Cartesian<double> >(obj.xcoord(), obj.ycoord(), obj.zcoord() );
}

const char* leda_tname(CGAL::Point_3<CGAL::Cartesian<double> >* p) {  return "CGALPoint_3"; }

leda_window& operator << (leda_window& w, const CGALPoint_3& obj)
{
  leda_point p(obj.x(), obj.y());
  w << p;
  return w;
}

leda_window& operator >> (leda_window& w, CGALPoint_3& p)
{
  leda_point p1;
  if( w >> p1 ) p = CGALPoint_3(p1.xcoord(), p1.ycoord(), 0);
  return w;
}

bool geowin_IntersectsBox(const CGAL::Point_3<CGAL::Cartesian<double> >& obj, double x1,double y1,double x2, double y2,bool f)
{
  double xw= obj.x();
  double yw= obj.y();  
  
  if (x1<=xw && x2>=xw && y1<=yw && y2>=yw) return true;
  return false;
}

void geowin_BoundingBox(const CGAL::Point_3<CGAL::Cartesian<double> >& obj, double& x1, double& x2,
	 double& y1, double& y2)
{
  x1=obj.x(); x2=obj.x(); 
  y1=obj.y(); y2=obj.y(); 
}

void geowin_Translate(CGAL::Point_3<CGAL::Cartesian<double> >& obj, double dx, double dy)
{
  CGAL::Vector_2<CGAL::Cartesian<double> > vec;
  vec= CGAL::Vector_2<CGAL::Cartesian<double> >(dx, dy); 
  CGAL::Aff_transformation_2<CGAL::Cartesian<double> > translate(CGAL::TRANSLATION, vec);
  
  CGAL::Point_2<CGAL::Cartesian<double> > pt(obj.x(), obj.y());
  CGAL::Point_2<CGAL::Cartesian<double> > pt2 = pt.transform(translate); 
  
  obj = CGAL::Point_3<CGAL::Cartesian<double> >(pt2.x(),pt2.y(),obj.z());
}


void geowin_Rotate(CGAL::Point_3<CGAL::Cartesian<double> >& obj, double x, double y, double a)
{
}

// functions for the container

leda_string geowin_info_fcn(const std::list<CGAL::Point_3<CGAL::Cartesian<double> > >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), "3d-point");  return str;
}

void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Point_3<CGAL::Cartesian<double> > >& L)
{
  leda_list<leda_d3_point> H;
  geowin_generate_objects(gw,H);

  //convert 
  CGAL::Point_3<CGAL::Cartesian<double> > obj;
  leda_d3_point mp;

  forall(mp,H){
   convert_from_leda(mp,obj);
   L.push_front(obj);
  }
}

GEOWIN_END_NAMESPACE

#endif 




