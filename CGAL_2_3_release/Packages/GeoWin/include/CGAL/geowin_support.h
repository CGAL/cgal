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
// file          : include/CGAL/geowin_support.h
// package       : GeoWin (1.2.2)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.2.2
// revision_date : 30 January 2001 
// author(s)     : Matthias Baesken, Ulrike Bartuschka, Stefan Naeher
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ======================================================================


#ifndef CGAL_GEOWIN_SUPPORT_H
#define CGAL_GEOWIN_SUPPORT_H


#include <CGAL/Cartesian.h>
#include <CGAL/basic.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Aff_transformation_2.h>
#include <LEDA/ps_file.h>
#include <LEDA/d3_segment.h>
#include <LEDA/d3_line.h>
#include <LEDA/rat_circle.h>
#include <LEDA/float_geo_alg.h>

#include<list>
#include<vector>

#include<LEDA/d3_point.h>

#if (__LEDA__ < 410)
// fix problem with missing prefixing for d3 rays ...
#if !defined(LEDA_ROOT_INCL_ID)
#define LEDA_ROOT_INCL_ID 420047
#include <LEDA/REDEFINE_NAMES.h>
#endif

#include <LEDA/d3_ray.h>

#if LEDA_ROOT_INCL_ID == 420047
#undef LEDA_ROOT_INCL_ID
#include <LEDA/UNDEFINE_NAMES.h>
#endif

typedef d3_ray  leda_d3_ray;
#else
#include <LEDA/d3_ray.h>
#endif


typedef CGAL::Cartesian<double>  DEFREP;

typedef CGAL::Point_2< DEFREP >     CGALPoint;
typedef CGAL::Vector_2< DEFREP >    CGALVector;
typedef std::list<CGALPoint>     CGALPointlist;

typedef CGAL::Segment_2< DEFREP >   CGALSegment;
typedef std::list<CGALSegment>   CGALSegmentlist;

typedef CGAL::Circle_2< DEFREP >    CGALCircle;
typedef std::list<CGALCircle>    CGALCirclelist;

typedef CGAL::Line_2< DEFREP >      CGALLine;
typedef std::list<CGALLine>      CGALLinelist;

typedef CGAL::Ray_2< DEFREP >       CGALRay;
typedef std::list<CGALRay>       CGALRaylist;

typedef CGAL::Triangle_2< DEFREP >  CGALTriangle; 
typedef std::list<CGALTriangle>  CGALTrianglelist;

typedef CGAL::Iso_rectangle_2< DEFREP>  CGALRectangle; 
typedef std::list<CGALRectangle>     CGALRectanglelist;

typedef CGAL::Polygon_traits_2< DEFREP >                 mypolyTraits;
typedef CGAL::Polygon_2<mypolyTraits, CGALPointlist > CGALPolygon;
typedef std::list<CGALPolygon>                        CGALPolygonlist;

//3d - points
typedef CGAL::Point_3< DEFREP >      CGALPoint_3;
typedef std::list<CGALPoint_3>    CGALPoint_3_list;

//3d - segments
typedef CGAL::Segment_3< DEFREP >      CGALSegment_3;
typedef std::list<CGALSegment_3>    CGALSegment_3_list;

//3d - lines
typedef CGAL::Line_3< DEFREP >      CGALLine_3;
typedef std::list<CGALLine_3>    CGALLine_3_list;

//3d - rays
typedef CGAL::Ray_3< DEFREP >      CGALRay_3;
typedef std::list<CGALRay_3>    CGALRay_3_list;

//3d - triangles
typedef CGAL::Triangle_3< DEFREP >      CGALTriangle_3;
typedef std::list<CGALTriangle_3>    CGALTriangle_3_list;

//3d - tetrahedra
typedef CGAL::Tetrahedron_3< DEFREP >      CGALTetrahedron_3;
typedef std::list<CGALTetrahedron_3>    CGALTetrahedron_3_list;


#if defined GEOWIN_USE_NAMESPACE

#if !defined GEOWIN_NAMESPACE_NAME
#define GEOWIN_NAMESPACE_NAME CGAL
#endif

#define GEOWIN_BEGIN_NAMESPACE namespace GEOWIN_NAMESPACE_NAME {
#define GEOWIN_END_NAMESPACE }
#else
#  define GEOWIN_BEGIN_NAMESPACE
#  define GEOWIN_END_NAMESPACE
#endif


// no templates on VC++
#if defined (_MSC_VER)
#define GEOWIN_SUPPORT_NO_TEMPLATES
#endif

#if defined (GEOWIN_SUPPORT_NO_TEMPLATES)
#include <CGAL/geowin_support_no_templ.h>
#else

template<class REP>
leda_point convert_to_leda(const CGAL::Point_2<REP>& obj)
{
  double x = CGAL::to_double(obj.x());
  double y = CGAL::to_double(obj.y());
  leda_point p(x,y);
  return p;
}

template<class REP>
leda_segment convert_to_leda(const CGAL::Segment_2<REP>& obj)
{
  CGAL::Point_2<REP> p1=obj.source();
  CGAL::Point_2<REP> p2=obj.target();
  leda_segment seg(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()), CGAL::to_double(p2.x()), CGAL::to_double(p2.y()) ); 
  return seg;
}

template<class REP>
leda_circle convert_to_leda(const CGAL::Circle_2<REP>& c)
{
 CGAL::Point_2<REP> p1=c.center();
 double   radius= ::sqrt(CGAL::to_double(c.squared_radius()));

 leda_point lp(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()));
 leda_circle lc(lp,radius);

 return lc;
}

template<class REP>
leda_line convert_to_leda(const CGAL::Line_2<REP>& l)
{
 CGAL::Point_2<REP> p1=l.point(1);
 CGAL::Point_2<REP> p2=l.point(2);

 leda_point lp1(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()));
 leda_point lp2(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()));
 
 leda_line lc(lp1,lp2);

 return lc;
}

template<class REP>
leda_ray convert_to_leda(const CGAL::Ray_2<REP>& r)
{
 CGAL::Point_2<REP> p1=r.source();
 CGAL::Point_2<REP> p2=r.point(1);

 leda_point lp1(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()));
 leda_point lp2(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()));
 
 leda_ray rc(lp1,lp2);

 return rc;
}

template<class TRAITS,class CONTAINER>
leda_polygon convert_to_leda(const CGAL::Polygon_2<TRAITS,CONTAINER>& p)
{
  CGAL::Polygon_2<TRAITS,CONTAINER>::Vertex_const_iterator it=p.vertices_begin();
  CGAL::Polygon_2<TRAITS,CONTAINER>::Vertex_const_iterator st=p.vertices_end();

 leda_list<leda_point> lp;

 while (it != st) {  lp.append(leda_point(CGAL::to_double((*it).x()),CGAL::to_double((*it).y()))); it++; }
 
 leda_polygon ph(lp);

 return ph;
}

template<class REP>
leda_polygon convert_to_leda(const CGAL::Triangle_2<REP>& t)
{
 CGAL::Point_2<REP> p1= t[1];
 CGAL::Point_2<REP> p2= t[2];
 CGAL::Point_2<REP> p3= t[3];
 leda_point lp1(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()));
 leda_point lp2(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()));
 leda_point lp3(CGAL::to_double(p3.x()), CGAL::to_double(p3.y()));

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

template<class REP>
leda_rectangle convert_to_leda(const CGAL::Iso_rectangle_2<REP>& t)
{
 CGAL::Point_2<REP> p1= t.min();
 CGAL::Point_2<REP> p2= t.max();
 leda_point lp1(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()));
 leda_point lp2(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()));

 return leda_rectangle(lp1,lp2); 
}

ps_file& operator<<(ps_file& F,const leda_d3_point& obj);

template<class REP>
leda_d3_point convert_to_leda(const CGAL::Point_3<REP>& obj)
{
  leda_d3_point p(CGAL::to_double(obj.x()), CGAL::to_double(obj.y()), CGAL::to_double(obj.z()) );
  return p;
}

template<class REP>
leda_d3_segment convert_to_leda(const CGAL::Segment_3<REP>& obj)
{
  leda_d3_segment s(convert_to_leda(obj.source()), convert_to_leda(obj.target()));
  return s;
}

template<class REP>
leda_d3_line convert_to_leda(const CGAL::Line_3<REP>& obj)
{
  leda_d3_line l(convert_to_leda(obj.point(0)), convert_to_leda(obj.point(1)));
  return l;  
}

template<class REP>
leda_d3_ray convert_to_leda(const CGAL::Ray_3<REP>& obj)
{
  leda_d3_ray r(convert_to_leda(obj.source()), convert_to_leda(obj.point(1)));
  return r;  
}

// d2 projection into xy - plane ...

template<class REP>
leda_polygon convert_to_leda(const CGAL::Triangle_3<REP>& obj)
{
  CGAL::Point_3<REP> p1 = obj[0];
  CGAL::Point_3<REP> p2 = obj[1];
  CGAL::Point_3<REP> p3 = obj[2];
  leda_point pl1(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()));  
  leda_point pl2(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()));  
  leda_point pl3(CGAL::to_double(p3.x()), CGAL::to_double(p3.y()));    
  leda_list<leda_point> Lh;
  
  if (orientation(pl1,pl2,pl3)== CGAL::LEFTTURN)
  { Lh.append(pl1); Lh.append(pl2); Lh.append(pl3); }
  else
  { Lh.push(pl1); Lh.push(pl2); Lh.push(pl3); }

  leda_polygon pol(Lh);
  return pol; 
}

// d2 projection into xy - plane

template<class REP>
leda_polygon convert_to_leda(const CGAL::Tetrahedron_3<REP>& obj)
{
  CGAL::Point_3<REP> p1 = obj[0];
  CGAL::Point_3<REP> p2 = obj[1];
  CGAL::Point_3<REP> p3 = obj[2];
  CGAL::Point_3<REP> p4 = obj[3]; 
  leda_point pl1(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()));  
  leda_point pl2(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()));  
  leda_point pl3(CGAL::to_double(p3.x()), CGAL::to_double(p3.y()));    
  leda_point pl4(CGAL::to_double(p4.x()), CGAL::to_double(p4.y()));  
    
  leda_list<leda_point> Lh;
  Lh.append(pl1); Lh.append(pl2); Lh.append(pl3); Lh.append(pl4);
  
  leda_list<leda_point> Lres = CONVEX_HULL(Lh);
  
  leda_polygon pol(Lres);
  return pol;  
}

template<class REP>
ps_file& operator<<(ps_file& F,const CGAL::Point_2<REP>& o) { F << convert_to_leda(o); return F; }

template<class REP>
ps_file& operator<<(ps_file& F,const CGAL::Segment_2<REP>& o) { F << convert_to_leda(o); return F; }

template<class REP>
ps_file& operator<<(ps_file& F,const CGAL::Circle_2<REP>& o) { F << convert_to_leda(o); return F; }

template<class REP>
ps_file& operator<<(ps_file& F,const CGAL::Line_2<REP>& o) { F << convert_to_leda(o); return F; }

template<class REP>
ps_file& operator<<(ps_file& F,const CGAL::Ray_2<REP>& o) { F << convert_to_leda(o); return F; }

template<class TRAITS,class CONTAINER>
ps_file& operator<<(ps_file& F,const CGAL::Polygon_2<TRAITS,CONTAINER>& o) 
{ F << convert_to_leda(o); return F; }

template<class REP>
ps_file& operator<<(ps_file& F,const CGAL::Triangle_2<REP>& o) { F << convert_to_leda(o); return F; }

template<class REP>
ps_file& operator<<(ps_file& F,const CGAL::Iso_rectangle_2<REP>& o) 
{ F << convert_to_leda(o); return F; }

template<class REP>
ps_file& operator<<(ps_file& F,const CGAL::Point_3<REP>& o) { F << convert_to_leda(o); return F; }

template<class REP>
ps_file& operator<<(ps_file& F,const CGAL::Segment_3<REP>& o) 
{ leda_d3_segment seg = convert_to_leda(o); 
  F << seg.project_xy();
  return F; 
}

template<class REP>
ps_file& operator<<(ps_file& F,const CGAL::Line_3<REP>& o) 
{ leda_d3_line l = convert_to_leda(o);
  leda_line m;
  
  if (l.project_xy(m)){ // projection is a line ...
    F << m;
  }
  else { // ... a point ...
    F << leda_point(l.point1().xcoord(), l.point1().ycoord());
  }
  
  return F; 
}

template<class REP>
ps_file& operator<<(ps_file& F,const CGAL::Ray_3<REP>& o) 
{ leda_d3_ray r = convert_to_leda(o);
  leda_ray m;
  
  if (r.project_xy(m)){ // projection is a line ...
    F << m;
  }
  else { // ... a point ...
    F << leda_point(r.point1().xcoord(), r.point1().ycoord());
  }
  
  return F; 
}

template<class REP>
ps_file& operator<<(ps_file& F,const CGAL::Triangle_3<REP>& o) 
{ F << convert_to_leda(o); return F; }

template<class REP>
ps_file& operator<<(ps_file& F,const CGAL::Tetrahedron_3<REP>& obj) 
{
  CGAL::Point_3<REP> p1 = obj[0];
  CGAL::Point_3<REP> p2 = obj[1];
  CGAL::Point_3<REP> p3 = obj[2];
  CGAL::Point_3<REP> p4 = obj[3]; 
  leda_point pl1(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()));  
  leda_point pl2(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()));  
  leda_point pl3(CGAL::to_double(p3.x()), CGAL::to_double(p3.y()));    
  leda_point pl4(CGAL::to_double(p4.x()), CGAL::to_double(p4.y()));  
  F << leda_segment(pl1,pl2); F << leda_segment(pl1,pl3); F << leda_segment(pl1,pl4);
  F << leda_segment(pl2,pl3); F << leda_segment(pl2,pl4); F << leda_segment(pl3,pl4);
  return F;
}


static void geowin_generate_circle_segments(leda_list<leda_segment>& LS, leda_circle C, int n)
{
  leda_list<leda_rat_point> L;
  leda_point p = C.point1(), q = C.point2(), r = C.point3();
  leda_rat_point rp(p), rq(q), rr(r);
  leda_rat_circle R(rp,rq,rr);
  
  double d = (2*LEDA_PI)/n;
  double eps = 0.001;
  double a = 0;
  
  for(int i=0; i < n; i++)
    { 
      leda_rat_point pp = R.point_on_circle(a,eps);
      L.append(pp);
      a += d;
    }
    
  // now generate the segments desribing the circle ...
  list_item lit = L.first();
  
  while(lit && L.succ(lit))
  {
    LS.append(leda_segment(L[lit].to_point(), L[L.succ(lit)].to_point()));
    lit = L.succ(lit);
  }
  LS.append(leda_segment(L.tail().to_point(), L.head().to_point()));  
}


#include <math.h>
#include <ctype.h>

GEOWIN_BEGIN_NAMESPACE

template<class REP>
bool geowin_IntersectsBox(const CGAL::Point_2<REP>& obj, double x1,double y1,double x2, double y2,bool f);

template<class REP>
bool geowin_IntersectsBox(const CGAL::Segment_2<REP>& obj, double x1,double y1,double x2, double y2,bool f);

template<class REP>
bool geowin_IntersectsBox(const CGAL::Circle_2<REP>& obj, double x1,double y1,double x2, double y2,bool f);

template<class REP>
bool geowin_IntersectsBox(const CGAL::Line_2<REP>& obj, double x1,double y1,double x2, double y2,bool f);

template<class REP>
bool geowin_IntersectsBox(const CGAL::Ray_2<REP>& obj, double x1,double y1,double x2, double y2,bool f);

template<class TRAITS, class CONTAINER>
bool geowin_IntersectsBox(const CGAL::Polygon_2<TRAITS,CONTAINER>& obj, double x1,double y1,double x2, double y2,bool f);

template<class REP>
bool geowin_IntersectsBox(const CGAL::Triangle_2<REP>& obj, double x1,double y1,double x2, double y2,bool f);

template<class REP>
bool geowin_IntersectsBox(const CGAL::Iso_rectangle_2<REP>& obj, double x1,double y1,double x2, double y2,bool f);

template<class REP>
bool geowin_IntersectsBox(const CGAL::Point_3<REP>& obj, double x1,double y1,double x2, double y2,bool f);

template<class REP>
bool geowin_IntersectsBox(const CGAL::Segment_3<REP>& obj, double x1,double y1,double x2, double y2,bool f);

template<class REP>
bool geowin_IntersectsBox(const CGAL::Line_3<REP>& obj, double x1,double y1,double x2, double y2,bool f);

template<class REP>
bool geowin_IntersectsBox(const CGAL::Ray_3<REP>& obj, double x1,double y1,double x2, double y2,bool f);

template<class REP>
bool geowin_IntersectsBox(const CGAL::Triangle_3<REP>& obj, double x1,double y1,double x2, double y2,bool f);

template<class REP>
bool geowin_IntersectsBox(const CGAL::Tetrahedron_3<REP>& obj, double x1,double y1,double x2, double y2,bool f);


// ----------------------------------------------
// templates for the functions for the containers
// ----------------------------------------------

template<class T>
void geowin_redraw_fcn(const std::list<T>& L, 
		    leda_window& w, leda_color c, leda_color, 
		    double x1, double y1, double x2, double y2)
{
#if defined (__GNUC__)
  typename std::list<T>::const_iterator it   = L.begin(), stop = L.end();
#else
  std::list<T>::const_iterator it   = L.begin(), stop = L.end();
#endif
  
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

template<class REP>
const char* leda_tname(CGAL::Point_2<REP>* p) {  return "CGALPoint"; }

template<class REP>
bool geowin_IntersectsBox(const CGAL::Point_2<REP>& obj, double x1,double y1,double x2, double y2,bool f)
{
  double xw= CGAL::to_double(obj.x());
  double yw= CGAL::to_double(obj.y());  
  
  if (x1<=xw && x2>=xw && y1<=yw && y2>=yw) return true;
  return false;
}

template<class REP>
void geowin_BoundingBox(const CGAL::Point_2<REP>& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_2 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

template<class REP>
void geowin_Translate(CGAL::Point_2<REP>& obj, double dx, double dy)
{
  typedef typename REP::RT RT;
  
  CGAL::Vector_2<REP> vec;
  vec= CGAL::Vector_2<REP>(RT(dx), RT(dy)); 
  obj = obj + vec;
}

template<class REP>
void geowin_Rotate(CGAL::Point_2<REP>& obj, double x, double y, double a)
{
  typedef typename REP::RT RT;

  leda_point p2(CGAL::to_double(obj.x()), CGAL::to_double(obj.y()));
  p2 = p2.rotate(leda_point(x,y), a);
  obj = CGAL::Point_2<REP>(RT(p2.xcoord()), RT(p2.ycoord())); 
}

// functions for the container

template<class REP>
leda_string geowin_info_fcn(const std::list<CGAL::Point_2<REP> >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-point");  return str;
}


template<class REP>
void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Point_2<REP> >& L)
{
  typedef typename REP::RT RT;
  
  leda_list<leda_point> H;
  geowin_generate_objects(gw,H);

  //convert the contents
  CGAL::Point_2<REP> p;
  leda_point mp;

  forall(mp,H){
   p= CGAL::Point_2<REP>(RT(mp.xcoord()), RT(mp.ycoord()));
   L.push_front(p);
  }
}

// 3d output ...
template<class T>
void cgal_Point_2_d3(const T& L, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GRAPH<leda_d3_point,int> G;
 typename T::const_iterator iter = L.begin();
 
 for(;iter != L.end(); iter++) {
   leda_point p = convert_to_leda(*iter);
   G.new_node(leda_d3_point(p.xcoord(), p.ycoord(),0));
 }
 H.join(G);
}


// Segments

template<class REP>
const char* leda_tname(CGAL::Segment_2<REP>* p) {  return "CGALSegment"; }

template<class REP>
bool geowin_IntersectsBox(const CGAL::Segment_2<REP>& obj, double x1,double y1,double x2, double y2,bool f)
{
  CGAL::Point_2<REP> p1,p2;
  p1=obj.source();
  p2=obj.target();
  leda_segment seg( CGAL::to_double(p1.x()),CGAL::to_double(p1.y()),CGAL::to_double(p2.x()),CGAL::to_double(p2.y()) );
  
  return geowin_IntersectsBox(seg,x1,y1,x2,y2,f);
}

template<class REP>
void geowin_BoundingBox(const CGAL::Segment_2<REP>& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_2 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

template<class REP>
void geowin_Translate(CGAL::Segment_2<REP>& obj, double dx, double dy)
{
  typedef typename REP::RT RT;

  CGAL::Vector_2<REP> vec;
  vec= CGAL::Vector_2<REP>(RT(dx), RT(dy)); 
  CGAL::Aff_transformation_2<REP> translate(CGAL::TRANSLATION, vec);
  
  obj = obj.transform(translate);
}

template<class REP>
void geowin_Rotate(CGAL::Segment_2<REP>& obj, double x, double y, double a)
{
  typedef typename REP::RT RT;
  
  CGAL::Point_2<REP> p1,p2;
  p1=obj.source();
  p2=obj.target();
  
  leda_segment hlp( CGAL::to_double(p1.x()), CGAL::to_double(p1.y()), CGAL::to_double(p2.x()), CGAL::to_double(p2.y()) );

  hlp = hlp.rotate(leda_point(x,y), a);
  p1= CGAL::Point_2<REP>( RT((hlp.source()).xcoord()), RT((hlp.source()).ycoord()) );
  p2= CGAL::Point_2<REP>( RT((hlp.target()).xcoord()), RT((hlp.target()).ycoord()));
  obj = CGAL::Segment_2<REP>(p1,p2);
}

// functions for the container

template<class REP>
leda_string geowin_info_fcn(const std::list<CGAL::Segment_2<REP> >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-segment");  return str;
}

template<class REP>
void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Segment_2<REP> >& L)
{
  typedef typename REP::RT RT;

  leda_list<leda_segment> H;
  geowin_generate_objects(gw,H);

  //convert the contents
  CGAL::Segment_2<REP> p;
  leda_segment mp;

  forall(mp,H){
   CGAL::Point_2<REP> pa(RT(mp.source().xcoord()), RT(mp.source().ycoord()));
   CGAL::Point_2<REP> pb(RT(mp.target().xcoord()), RT(mp.target().ycoord()));   
   
   p= CGAL::Segment_2<REP>(pa,pb);
   L.push_front(p);
  }
}


// 3d output ...
template<class T>
void cgal_Segment_2_d3(const T& L, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GRAPH<leda_d3_point,int> G;
 typename T::const_iterator iter = L.begin();
 
 for(;iter != L.end(); iter++) {
   leda_segment s = convert_to_leda(*iter);
   leda_node v1 = G.new_node(leda_d3_point(s.source().xcoord(),s.source().ycoord(),0));
   leda_node v2 = G.new_node(leda_d3_point(s.target().xcoord(),s.target().ycoord(),0));   
   leda_edge e1 = G.new_edge(v1,v2);
   leda_edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);
 }
 H.join(G);
}


// Circles

template<class REP>
void convert_from_leda(const leda_circle& c, CGAL::Circle_2<REP>& cir)
{
 typedef typename REP::RT RT;
 
 leda_point lp=c.center();
 double sr= c.radius()*c.radius();

 CGAL::Point_2<REP> p1(RT(lp.xcoord()), RT(lp.ycoord()));
 cir = CGAL::Circle_2<REP>(p1, RT(sr));
}

template<class REP>
const char* leda_tname(CGAL::Circle_2<REP>* p) {  return "CGALCircle"; }

template<class REP>
bool geowin_IntersectsBox(const CGAL::Circle_2<REP>& obj, double x1,double y1,double x2, double y2,bool f)
{  
  return geowin_IntersectsBox(convert_to_leda(obj),x1,y1,x2,y2,f); 
}

template<class REP>
void geowin_BoundingBox(const CGAL::Circle_2<REP>& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_2 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

template<class REP>
void geowin_Translate(CGAL::Circle_2<REP>& obj, double dx, double dy)
{
  typedef typename REP::RT RT;

  CGAL::Vector_2<REP> vec;
  vec= CGAL::Vector_2<REP>(RT(dx), RT(dy)); 
  CGAL::Aff_transformation_2<REP> translate(CGAL::TRANSLATION, vec);
  
  obj = obj.orthogonal_transform(translate);
}

template<class REP>
void geowin_Rotate(CGAL::Circle_2<REP>& obj, double x, double y, double a)
{  
  leda_circle hlp=convert_to_leda(obj);
  hlp = hlp.rotate(leda_point(x,y), a);
  convert_from_leda(hlp,obj);
}


template<class REP>
leda_string geowin_info_fcn(const std::list<CGAL::Circle_2<REP> >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-circle");  return str;
}

template<class REP>
void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Circle_2<REP> >& L)
{
  leda_list<leda_circle> H;
  geowin_generate_objects(gw,H);

  //convert the contents
  CGAL::Circle_2<REP> c;
  leda_circle mp;

  forall(mp,H){
   convert_from_leda(mp,c);
   L.push_front(c);
  }
}

// 3d output ...
template<class T>
void cgal_Circle_2_d3(const T& L, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GRAPH<leda_d3_point,int> G;
 typename T::const_iterator iter = L.begin();
 
 for(;iter != L.end(); iter++) {
   leda_circle c = convert_to_leda(*iter);
   leda_list<leda_segment> LS;
   geowin_generate_circle_segments(LS,c,30);
   leda_segment siter;
   forall(siter,LS){
     leda_node v1 = G.new_node(leda_d3_point(siter.source().xcoord(),siter.source().ycoord(),0));
     leda_node v2 = G.new_node(leda_d3_point(siter.target().xcoord(),siter.target().ycoord(),0));   
     leda_edge e1 = G.new_edge(v1,v2);
     leda_edge e2 = G.new_edge(v2,v1);
     G.set_reversal(e1,e2);     
   }   
 }
 H.join(G);
}


// Lines

template<class REP>
void convert_from_leda(const leda_line& l, CGAL::Line_2<REP>& lc)
{
 typedef typename REP::RT RT;

 leda_point lp1=l.point1();
 leda_point lp2=l.point2();

 CGAL::Point_2<REP> p1(RT(lp1.xcoord()), RT(lp1.ycoord()));
 CGAL::Point_2<REP> p2(RT(lp2.xcoord()), RT(lp2.ycoord()));
 lc = CGAL::Line_2<REP>(p1,p2);
}

template<class REP>
const char* leda_tname(CGAL::Line_2<REP>* p) {  return "CGALLine"; }

template<class REP>
bool geowin_IntersectsBox(const CGAL::Line_2<REP>& obj, double x1,double y1,double x2, double y2,bool f)
{  
  return geowin_IntersectsBox(convert_to_leda(obj),x1,y1,x2,y2,f); 
}

template<class REP>
void geowin_BoundingBox(const CGAL::Line_2<REP>& obj, double& x1, double& x2,
	 double& y1, double& y2)
{
  geowin_BoundingBox(convert_to_leda(obj),x1,x2,y1,y2);
}

template<class REP>
void geowin_Translate(CGAL::Line_2<REP>& obj, double dx, double dy)
{
  typedef typename REP::RT RT;

  CGAL::Vector_2<REP> vec;
  vec= CGAL::Vector_2<REP>(RT(dx), RT(dy)); 
  CGAL::Aff_transformation_2<REP> translate(CGAL::TRANSLATION, vec);
  
  obj = obj.transform(translate);
}

template<class REP>
void geowin_Rotate(CGAL::Line_2<REP>& obj, double x, double y, double a)
{  
  leda_line hlp=convert_to_leda(obj);
  hlp = hlp.rotate(leda_point(x,y), a);
  convert_from_leda(hlp,obj);
}

// functions for the container

template<class REP>
leda_string geowin_info_fcn(const std::list<CGAL::Line_2<REP> >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-Line");  return str;
}

template<class REP>
void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Line_2<REP> >& L)
{
  leda_list<leda_line> H;
  geowin_generate_objects(gw,H);

  //convert the contents
  CGAL::Line_2<REP> obj;
  leda_line mp;

  forall(mp,H){
   convert_from_leda(mp,obj);
   L.push_front(obj);
  }
}

// 3d output ...
template<class T>
void cgal_Line_2_d3(const T& L, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GRAPH<leda_d3_point,int> G;
 typename T::const_iterator iter = L.begin();
 
 for(;iter != L.end(); iter++) {
   leda_line li = convert_to_leda(*iter);
   leda_point pm((li.point1().xcoord()+li.point2().xcoord())/2,(li.point1().ycoord()+li.point2().ycoord())/2);
   leda_vector v = li.point1() - li.point2();
   v= v * 50;
   leda_point p1=pm+v, p2=pm-v;
   leda_node v1 = G.new_node(leda_d3_point(p1.xcoord(),p1.ycoord(),0));
   leda_node v2 = G.new_node(leda_d3_point(p2.xcoord(),p2.ycoord(),0));   
   leda_edge e1 = G.new_edge(v1,v2);
   leda_edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);   
 }
 H.join(G);
}

// Rays

template<class REP>
void convert_from_leda(const leda_ray& r, CGAL::Ray_2<REP>& rc)
{
 typedef typename REP::RT RT;

 leda_point lp1=r.point1();
 leda_point lp2=r.point2();

 CGAL::Point_2<REP> p1(RT(lp1.xcoord()), RT(lp1.ycoord()));
 CGAL::Point_2<REP> p2(RT(lp2.xcoord()), RT(lp2.ycoord()));
 rc = CGAL::Ray_2<REP>(p1,p2);
}

template<class REP>
const char* leda_tname(CGAL::Ray_2<REP>* p) {  return "CGALRay"; }

template<class REP>
bool geowin_IntersectsBox(const CGAL::Ray_2<REP>& obj, double x1,double y1,double x2, double y2,bool f)
{  
  return geowin_IntersectsBox(convert_to_leda(obj),x1,y1,x2,y2,f); 
}

template<class REP>
void geowin_BoundingBox(const CGAL::Ray_2<REP>& obj, double& x1, double& x2,
	 double& y1, double& y2)
{
  geowin_BoundingBox(convert_to_leda(obj),x1,x2,y1,y2);
}

template<class REP>
void geowin_Translate(CGAL::Ray_2<REP>& obj, double dx, double dy)
{
  typedef typename REP::RT RT;

  CGAL::Vector_2<REP> vec;
  vec= CGAL::Vector_2<REP>(RT(dx), RT(dy)); 
  CGAL::Aff_transformation_2<REP> translate(CGAL::TRANSLATION, vec);
  
  obj = obj.transform(translate); 
}

template<class REP>
void geowin_Rotate(CGAL::Ray_2<REP>& obj, double x, double y, double a)
{  
  leda_ray hlp=convert_to_leda(obj);
  hlp = hlp.rotate(leda_point(x,y), a);
  convert_from_leda(hlp,obj);
}

// functions for the container

template<class REP>
leda_string geowin_info_fcn(const std::list<CGAL::Ray_2<REP> >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-Ray");  return str;
}

template<class REP>
void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Ray_2<REP> >& L)
{
  leda_list<leda_ray> H;
  geowin_generate_objects(gw,H);

  //convert the contents
  CGAL::Ray_2<REP> obj;
  leda_ray mp;

  forall(mp,H){
   convert_from_leda(mp,obj);
   L.push_front(obj);
  }
}

// 3d output ...
template<class T>
void cgal_Ray_2_d3(const T& L, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GRAPH<leda_d3_point,int> G;
 typename T::const_iterator iter = L.begin();
 
 for(;iter != L.end(); iter++) {
   leda_ray r = convert_to_leda(*iter);
   leda_vector v = r.point2() - r.point1();
   v= v * 50;
   leda_point p1=r.source(), p2=p1 + v;
   leda_node v1 = G.new_node(leda_d3_point(p1.xcoord(),p1.ycoord(),0));
   leda_node v2 = G.new_node(leda_d3_point(p2.xcoord(),p2.ycoord(),0));   
   leda_edge e1 = G.new_edge(v1,v2);
   leda_edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);   
 }
 H.join(G);
}

// Triangles

template<class REP>
void convert_from_leda(const leda_polygon& p, CGAL::Triangle_2<REP>& tr)
{
 typedef typename REP::RT RT;
 
 leda_list<leda_point> L=p.vertices();
 leda_point lp1,lp2,lp3;
 lp1= L.pop(); lp2= L.pop(); lp3= L.pop();
 CGAL::Point_2<REP> p1(RT(lp1.xcoord()), RT(lp1.ycoord()));
 CGAL::Point_2<REP> p2(RT(lp2.xcoord()), RT(lp2.ycoord()));
 CGAL::Point_2<REP> p3(RT(lp3.xcoord()), RT(lp3.ycoord()));

 tr= CGAL::Triangle_2<REP>(p1,p2,p3); 
}

template<class REP>
const char* leda_tname(CGAL::Triangle_2<REP>* t) {  return "CGALTriangle"; }

template<class REP>
bool geowin_IntersectsBox(const CGAL::Triangle_2<REP>& obj, double x1,double y1,double x2, double y2,bool f)
{  
  return geowin_IntersectsBox(convert_to_leda(obj),x1,y1,x2,y2,f); 
}

template<class REP>
void geowin_BoundingBox(const CGAL::Triangle_2<REP>& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_2 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

template<class REP>
void geowin_Translate(CGAL::Triangle_2<REP>& obj, double dx, double dy)
{
  typedef typename REP::RT RT;

  CGAL::Vector_2<REP> vec;
  vec= CGAL::Vector_2<REP>(RT(dx), RT(dy)); 
  CGAL::Aff_transformation_2<REP> translate(CGAL::TRANSLATION, vec);
  
  obj = obj.transform(translate);
}

template<class REP>
void geowin_Rotate(CGAL::Triangle_2<REP>& obj, double x, double y, double a)
{  
  leda_polygon hlp=convert_to_leda(obj);
  hlp = hlp.rotate(leda_point(x,y), a);
  convert_from_leda(hlp,obj);
}

template<class REP>
leda_string geowin_info_fcn(const std::list<CGAL::Triangle_2<REP> >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-triangle");  return str;
}

template<class REP>
void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Triangle_2<REP> >& L)
{
}

// 3d output ...
template<class T>
void cgal_Triangle_2_d3(const T& L, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GRAPH<leda_d3_point,int> G;
 typename T::const_iterator iter = L.begin();
 
 for(;iter != L.end(); iter++) {
   leda_point p0 = convert_to_leda((*iter).vertex(0));
   leda_point p1 = convert_to_leda((*iter).vertex(1)); 
   leda_point p2 = convert_to_leda((*iter).vertex(2)); 
    
   leda_node v0 = G.new_node(leda_d3_point(p0.xcoord(),p0.ycoord(),0));
   leda_node v1 = G.new_node(leda_d3_point(p1.xcoord(),p1.ycoord(),0));  
   leda_node v2 = G.new_node(leda_d3_point(p2.xcoord(),p2.ycoord(),0)); 

   leda_edge e1 = G.new_edge(v0,v1);
   leda_edge e2 = G.new_edge(v1,v0);       
   leda_edge e3 = G.new_edge(v1,v2);
   leda_edge e4 = G.new_edge(v2,v1);
   leda_edge e5 = G.new_edge(v2,v0);
   leda_edge e6 = G.new_edge(v0,v2);
   
   G.set_reversal(e1,e2);   
   G.set_reversal(e3,e4);
   G.set_reversal(e5,e6);
 }
 H.join(G);
}

// Rectangles

template<class REP>
void convert_from_leda(const leda_rectangle& r, CGAL::Iso_rectangle_2<REP>& rc)
{
 typedef typename REP::RT RT;
 
 leda_point lp1= r.lower_left();
 leda_point lp2= r.upper_right();
 CGAL::Point_2<REP> p1(RT(lp1.xcoord()),RT(lp1.ycoord()));
 CGAL::Point_2<REP> p2(RT(lp2.xcoord()),RT(lp2.ycoord()));

 rc = CGAL::Iso_rectangle_2<REP>(p1,p2); 
}

template<class REP>
const char* leda_tname(CGAL::Iso_rectangle_2<REP>* t) {  return "CGALRectangle"; }

template<class REP>
bool geowin_IntersectsBox(const CGAL::Iso_rectangle_2<REP>& obj, double x1,double y1,double x2, double y2,bool f)
{  
  typedef typename REP::RT RT;
  
  CGAL::Point_2<REP> pa;
  pa = CGAL::Point_2<REP>(RT(x1), RT(y1));
  CGAL::Point_2<REP> pb;
  pb = CGAL::Point_2<REP>(RT(x2), RT(y2));
  
  CGAL::Iso_rectangle_2<REP> r2(pa,pb);
  
  return CGAL::do_intersect(obj,r2);  
}

template<class REP>
void geowin_BoundingBox(const CGAL::Iso_rectangle_2<REP>& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_2 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

template<class REP>
void geowin_Translate(CGAL::Iso_rectangle_2<REP>& obj, double dx, double dy)
{
  typedef typename REP::RT RT;

  CGAL::Vector_2<REP> vec;
  vec= CGAL::Vector_2<REP>(RT(dx), RT(dy)); 
  CGAL::Aff_transformation_2<REP> translate(CGAL::TRANSLATION, vec);
  
  obj = obj.transform(translate);
}

template<class REP>
void geowin_Rotate(CGAL::Iso_rectangle_2<REP>& obj, double x, double y, double a)
{  
}

template<class REP>
leda_string geowin_info_fcn(const std::list<CGAL::Iso_rectangle_2<REP> >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-rectangle");  return str;
}

template<class REP>
void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Iso_rectangle_2<REP> >& L)
{
}

// 3d output ...
template<class T>
void cgal_Iso_rectangle_2_d3(const T& L, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GRAPH<leda_d3_point,int> G;
 typename T::const_iterator iter = L.begin();
 
 for(;iter != L.end(); iter++) {
   leda_point p0 = convert_to_leda((*iter).vertex(0));
   leda_point p1 = convert_to_leda((*iter).vertex(1)); 
   leda_point p2 = convert_to_leda((*iter).vertex(2)); 
   leda_point p3 = convert_to_leda((*iter).vertex(3));    
    
   leda_node v0 = G.new_node(leda_d3_point(p0.xcoord(),p0.ycoord(),0));
   leda_node v1 = G.new_node(leda_d3_point(p1.xcoord(),p1.ycoord(),0));  
   leda_node v2 = G.new_node(leda_d3_point(p2.xcoord(),p2.ycoord(),0));
   leda_node v3 = G.new_node(leda_d3_point(p3.xcoord(),p3.ycoord(),0));    

   leda_edge e1 = G.new_edge(v0,v1);
   leda_edge e2 = G.new_edge(v1,v0);       
   leda_edge e3 = G.new_edge(v1,v2);
   leda_edge e4 = G.new_edge(v2,v1);
   leda_edge e5 = G.new_edge(v2,v3);
   leda_edge e6 = G.new_edge(v3,v2);
   leda_edge e7 = G.new_edge(v3,v0);
   leda_edge e8 = G.new_edge(v0,v3);
      
   G.set_reversal(e1,e2);   
   G.set_reversal(e3,e4);
   G.set_reversal(e5,e6);
   G.set_reversal(e7,e8);   
 }
 H.join(G);
}


// Polygons
template<class TRAITS, class CONTAINER>
void convert_from_leda(const leda_polygon& p, CGAL::Polygon_2<TRAITS,CONTAINER>& rc)
{
 typedef typename TRAITS::Point_2 POINT;
 typedef typename POINT::RT  RT;

 leda_list<leda_point> pl= p.vertices();
 std::list<POINT> sl;

 leda_point akt;
 forall(akt,pl) sl.push_back(POINT(RT(akt.xcoord()), RT(akt.ycoord())));
 
 rc = CGAL::Polygon_2<TRAITS,CONTAINER>(sl.begin(),sl.end());
}

template<class TRAITS, class CONTAINER>
const char* leda_tname(CGAL::Polygon_2<TRAITS,CONTAINER>* p) 
{  return "CGALPolygon"; }

template<class TRAITS, class CONTAINER>
leda_window& operator >> (leda_window& w, CGAL::Polygon_2<TRAITS,CONTAINER>& obj)
{
  leda_polygon p;
  w >> p;
  convert_from_leda(p, obj);
  return w;
}

template<class TRAITS,class CONTAINER>
bool geowin_IntersectsBox(const CGAL::Polygon_2<TRAITS,CONTAINER>& obj, double x1,double y1,double x2, double y2,bool f)
{  
  return geowin_IntersectsBox(convert_to_leda(obj),x1,y1,x2,y2,f); 
}

template<class TRAITS,class CONTAINER>
void geowin_BoundingBox(const CGAL::Polygon_2<TRAITS,CONTAINER>& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_2 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

template<class TRAITS, class CONTAINER>
void geowin_Translate(CGAL::Polygon_2<TRAITS,CONTAINER>& obj, double dx, double dy)
{
 typedef typename TRAITS::Point_2 POINT;
 typedef typename POINT::RT  RT;
 typedef typename POINT::R   REP;
 
  CGAL::Vector_2<REP> vec;
  vec= CGAL::Vector_2<REP>(RT(dx), RT(dy)); 
  CGAL::Aff_transformation_2<REP> translate(CGAL::TRANSLATION, vec);
  
  obj = CGAL::transform(translate,obj);
}

template<class TRAITS, class CONTAINER>
void geowin_Rotate(CGAL::Polygon_2<TRAITS,CONTAINER>& obj, double x, double y, double a)
{  
  leda_polygon hlp=convert_to_leda(obj);
  hlp = hlp.rotate(leda_point(x,y), a);
  convert_from_leda(hlp,obj);
}

template<class TRAITS,class CONTAINER>
leda_string geowin_info_fcn(const std::list<CGAL::Polygon_2<TRAITS,CONTAINER> >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-Polygons");  return str;
}

template<class TRAITS,class CONTAINER>
void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Polygon_2<TRAITS,CONTAINER> >& L)
{
  leda_list<leda_polygon> H;
  geowin_generate_objects(gw,H);

  //convert the contents
  CGAL::Polygon_2<TRAITS,CONTAINER> obj;
  leda_polygon mp;

  forall(mp,H){
   convert_from_leda(mp,obj);
   L.push_front(obj);
  }
}


// 3d output ...
template<class T>
void cgal_Polygon_2_d3(const T& L, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GRAPH<leda_d3_point,int> G;
 typename T::const_iterator iter = L.begin();
 leda_list<leda_segment> LS;
 
 for(;iter != L.end(); iter++) {
   leda_polygon p = convert_to_leda(*iter);
   LS = p.segments();
   leda_segment siter;
   forall(siter,LS) {
    leda_node v1 = G.new_node(leda_d3_point(siter.source().xcoord(),siter.source().ycoord(),0));
    leda_node v2 = G.new_node(leda_d3_point(siter.target().xcoord(),siter.target().ycoord(),0));   
    leda_edge e1 = G.new_edge(v1,v2);
    leda_edge e2 = G.new_edge(v2,v1);
    G.set_reversal(e1,e2);
   }       
 }
 H.join(G);
}


// 3d Points ...

template<class REP>
void convert_from_leda(const leda_d3_point& obj, CGAL::Point_3<REP>& p)
{
  typedef typename REP::RT RT;

  p = CGAL::Point_3<REP>(RT(obj.xcoord()), RT(obj.ycoord()), RT(obj.zcoord()) );
}

template<class REP>
const char* leda_tname(CGAL::Point_3<REP>* p) {  return "CGALPoint_3"; }

template<class REP>
leda_window& operator << (leda_window& w, const CGAL::Point_3<REP>& obj)
{
  leda_point p(CGAL::to_double(obj.x()), CGAL::to_double(obj.y()));
  w << p;
  return w;
}

template<class REP>
leda_window& operator >> (leda_window& w, CGAL::Point_3<REP>& p)
{
  typedef typename REP::RT RT;

  leda_point p1;
  if( w >> p1 ) p = CGAL::Point_3<REP>(RT(p1.xcoord()), RT(p1.ycoord()), RT(0));
  return w;
}

template<class REP>
bool geowin_IntersectsBox(const CGAL::Point_3<REP>& obj, double x1,double y1,double x2, double y2,bool f)
{
  double xw= CGAL::to_double(obj.x());
  double yw= CGAL::to_double(obj.y());  
  
  if (x1<=xw && x2>=xw && y1<=yw && y2>=yw) return true;
  return false;
}

template<class REP>
void geowin_BoundingBox(const CGAL::Point_3<REP>& obj, double& x1, double& x2,
	 double& y1, double& y2)
{
  x1=CGAL::to_double(obj.x()); x2=CGAL::to_double(obj.x()); 
  y1=CGAL::to_double(obj.y()); y2=CGAL::to_double(obj.y()); 
}

template<class REP>
void geowin_Translate(CGAL::Point_3<REP>& obj, double dx, double dy)
{
  typedef typename REP::RT  RT;
 
  CGAL::Vector_2<REP> vec;
  vec= CGAL::Vector_2<REP>(RT(dx), RT(dy)); 
  CGAL::Aff_transformation_2<REP> translate(CGAL::TRANSLATION, vec);
  
  CGAL::Point_2<REP> pt(obj.x(), obj.y());
  CGAL::Point_2<REP> pt2 = pt.transform(translate); 
  
  obj = CGAL::Point_3<REP>(pt2.x(),pt2.y(),obj.z());
}

template<class REP>
void geowin_Rotate(CGAL::Point_3<REP>& obj, double x, double y, double a)
{
}

// functions for the container

template<class REP>
leda_string geowin_info_fcn(const std::list<CGAL::Point_3<REP> >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), "3d-point");  return str;
}

template<class REP>
void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Point_3<REP> >& L)
{
  leda_list<leda_d3_point> H;
  geowin_generate_objects(gw,H);

  //convert 
  CGAL::Point_3<REP> obj;
  leda_d3_point mp;

  forall(mp,H){
   convert_from_leda(mp,obj);
   L.push_back(obj);
  }
}

// 3d output ...
template<class T>
void cgal_Point_3_d3(const T& L, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GRAPH<leda_d3_point,int> G;
 typename T::const_iterator iter = L.begin();
 for(;iter != L.end(); iter++) {
   leda_d3_point p = convert_to_leda(*iter);
   G.new_node(p);
 }
 H.join(G);
}


// 3d Segments

template<class REP>
const char* leda_tname(CGAL::Segment_3<REP>* p) {  return "CGALSegment_3"; }

template<class REP>
leda_window& operator << (leda_window& w, const CGAL::Segment_3<REP>& obj)
{
  leda_d3_segment seg3 = convert_to_leda(obj);
  leda_segment seg(seg3.source().xcoord(),seg3.source().ycoord(), seg3.target().xcoord(),seg3.target().ycoord());
  w << seg;
  return w;
}

template<class REP>
leda_window& operator >> (leda_window& w, CGAL::Segment_3<REP>& p)
{
  typedef typename REP::RT RT;

  leda_segment seg;
  if( w >> seg ) {
    leda_point ps = seg.source(), pt = seg.target();
    CGAL::Point_3<REP> cps(RT(ps.xcoord()), RT(ps.ycoord()), RT(0));
    CGAL::Point_3<REP> cpt(RT(pt.xcoord()), RT(pt.ycoord()), RT(0));    
    
    p = CGAL::Segment_3<REP>(cps, cpt);
  }
  return w;
}

template<class REP>
bool geowin_IntersectsBox(const CGAL::Segment_3<REP>& obj, double x1,double y1,double x2, double y2,bool f)
{
  leda_d3_segment seg3 = convert_to_leda(obj);
  leda_segment seg(seg3.source().xcoord(),seg3.source().ycoord(), seg3.target().xcoord(),seg3.target().ycoord() );
  return geowin_IntersectsBox(seg,x1,y1,x2,y2,f);
}

template<class REP>
void geowin_BoundingBox(const CGAL::Segment_3<REP>& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_3 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

template<class REP>
void geowin_Translate(CGAL::Segment_3<REP>& obj, double dx, double dy)
{
  typedef typename REP::RT RT;

  CGAL::Vector_3<REP> vec;
  vec= CGAL::Vector_3<REP>(RT(dx), RT(dy), RT(0)); 
  
  CGAL::Aff_transformation_3<REP> translate(CGAL::TRANSLATION, vec);
  
  obj = obj.transform(translate);
}

template<class REP>
void geowin_Rotate(CGAL::Segment_3<REP>& obj, double x, double y, double a)
{
}

// functions for the container

template<class REP>
leda_string geowin_info_fcn(const std::list<CGAL::Segment_3<REP> >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " d3-segment");  return str;
}

template<class REP>
void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Segment_3<REP> >& L)
{
}

// 3d output ...
template<class T>
void cgal_Segment_3_d3(const T& L, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GRAPH<leda_d3_point,int> G;
 typename T::const_iterator iter = L.begin();
 
 for(;iter != L.end(); iter++) {
   leda_d3_segment s = convert_to_leda(*iter);
   leda_node v1 = G.new_node(s.source());
   leda_node v2 = G.new_node(s.target());   
   leda_edge e1 = G.new_edge(v1,v2);
   leda_edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);
 }
 H.join(G);
}

// 3d Lines

template<class REP>
const char* leda_tname(CGAL::Line_3<REP>* p) {  return "CGALLine_3"; }

template<class REP>
leda_window& operator << (leda_window& w, const CGAL::Line_3<REP>& obj)
{
  leda_d3_line l3 = convert_to_leda(obj);
  leda_d3_point p1 = l3.point1();
  leda_d3_point p2 = l3.point2();
  
  if ((p1.xcoord()==p2.xcoord()) && (p1.ycoord()==p2.ycoord())) 
    w << leda_point(p1.xcoord(), p1.ycoord());
  else {
   leda_line l(leda_point(l3.point1().xcoord(),l3.point1().ycoord()), leda_point(l3.point2().xcoord(), l3.point2().ycoord()));
   w << l;
  } 
   
  return w;
}

template<class REP>
leda_window& operator >> (leda_window& w, CGAL::Line_3<REP>& l)
{
  typedef typename REP::RT RT;

  leda_line lh;
  if( w >> lh ) {
    leda_point ps = lh.point1(), pt = lh.point2();
    CGAL::Point_3<REP> cps(RT(ps.xcoord()), RT(ps.ycoord()), RT(0));
    CGAL::Point_3<REP> cpt(RT(pt.xcoord()), RT(pt.ycoord()), RT(0)); 
    
    if (cps == cpt) cpt = CGAL::Point_3<REP>(RT(pt.xcoord()), RT(pt.ycoord()), RT(10));   
    
    l = CGAL::Line_3<REP>(cps, cpt);
  }
  return w;
}

template<class REP>
bool geowin_IntersectsBox(const CGAL::Line_3<REP>& obj, double x1,double y1,double x2, double y2,bool f)
{
  leda_d3_line l3 = convert_to_leda(obj);
  leda_segment seg(l3.point1().xcoord(),l3.point1().ycoord(), l3.point2().xcoord(),l3.point2().ycoord() );
  return geowin_IntersectsBox(seg,x1,y1,x2,y2,f);
}

template<class REP>
void geowin_BoundingBox(const CGAL::Line_3<REP>& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Point_3<REP> cps = obj.point(0), cpt = obj.point(1);
  CGAL::Segment_3<REP> hobj(cps, cpt);
  CGAL::Bbox_3 bb= hobj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

template<class REP>
void geowin_Translate(CGAL::Line_3<REP>& obj, double dx, double dy)
{
  typedef typename REP::RT RT;

  CGAL::Vector_3<REP> vec;
  vec= CGAL::Vector_3<REP>(RT(dx), RT(dy), RT(0)); 
  
  CGAL::Aff_transformation_3<REP> translate(CGAL::TRANSLATION, vec);
  
  obj = obj.transform(translate);
}

template<class REP>
void geowin_Rotate(CGAL::Line_3<REP>& obj, double x, double y, double a)
{
}

// functions for the container

template<class REP>
leda_string geowin_info_fcn(const std::list<CGAL::Line_3<REP> >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " d3-line");  return str;
}

template<class REP>
void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Line_3<REP> >& L)
{
}


// 3d output ...
template<class T>
void cgal_Line_3_d3(const T& L, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GRAPH<leda_d3_point,int> G;
 typename T::const_iterator iter = L.begin();
 
 for(;iter != L.end(); iter++) {
   leda_d3_line l = convert_to_leda(*iter);

   leda_d3_point pm((l.point1().xcoord()+l.point2().xcoord())/2, \
               (l.point1().ycoord()+l.point2().ycoord())/2, \
	       (l.point1().zcoord()+l.point2().zcoord())/2);
	       
   leda_vector v = l.point1() - l.point2();
   v= v * 50;
   leda_d3_point p1=pm+v, p2=pm-v;
   leda_node v1 = G.new_node(p1);
   leda_node v2 = G.new_node(p2);   
   leda_edge e1 = G.new_edge(v1,v2);
   leda_edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);   
 }
 H.join(G);
}


// 3d Ray

template<class REP>
const char* leda_tname(CGAL::Ray_3<REP>* p) {  return "CGALRay_3"; }

template<class REP>
leda_window& operator << (leda_window& w, const CGAL::Ray_3<REP>& obj)
{
  leda_d3_ray r3 = convert_to_leda(obj);
  leda_d3_point p1 = r3.point1();
  leda_d3_point p2 = r3.point2();
  
  if ((p1.xcoord()==p2.xcoord()) && (p1.ycoord()==p2.ycoord())) 
    w << leda_point(p1.xcoord(), p1.ycoord());
  else {
   leda_ray r(leda_point(p1.xcoord(),p1.ycoord()), leda_point(p2.xcoord(), p2.ycoord()));
   w << r;
  } 
   
  return w;
}

template<class REP>
leda_window& operator >> (leda_window& w, CGAL::Ray_3<REP>& r)
{
  typedef typename REP::RT RT;

  leda_ray rh;
  if( w >> rh ) {
    leda_point ps = rh.point1(), pt = rh.point2();
    CGAL::Point_3<REP> cps(RT(ps.xcoord()), RT(ps.ycoord()), RT(0));
    CGAL::Point_3<REP> cpt(RT(pt.xcoord()), RT(pt.ycoord()), RT(0));    
  
    if (cps == cpt) cpt = CGAL::Point_3<REP>(RT(pt.xcoord()), RT(pt.ycoord()), RT(10));   
    
    r = CGAL::Ray_3<REP>(cps, cpt);
  }
  return w;
}

template<class REP>
bool geowin_IntersectsBox(const CGAL::Ray_3<REP>& obj, double x1,double y1,double x2, double y2,bool f)
{
  leda_d3_ray r3 = convert_to_leda(obj);
  leda_segment seg(r3.point1().xcoord(),r3.point1().ycoord(), r3.point2().xcoord(),r3.point2().ycoord() );
  return geowin_IntersectsBox(seg,x1,y1,x2,y2,f);
}

template<class REP>
void geowin_BoundingBox(const CGAL::Ray_3<REP>& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Point_3<REP> cps = obj.point(0), cpt = obj.point(1);
  CGAL::Segment_3<REP> hobj(cps, cpt);
  CGAL::Bbox_3 bb= hobj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

template<class REP>
void geowin_Translate(CGAL::Ray_3<REP>& obj, double dx, double dy)
{
  typedef typename REP::RT RT;

  CGAL::Vector_3<REP> vec;
  vec= CGAL::Vector_3<REP>(RT(dx), RT(dy), RT(0)); 
  
  CGAL::Aff_transformation_3<REP> translate(CGAL::TRANSLATION, vec);
  
  obj = obj.transform(translate);
}

template<class REP>
void geowin_Rotate(CGAL::Ray_3<REP>& obj, double x, double y, double a)
{
}

// functions for the container

template<class REP>
leda_string geowin_info_fcn(const std::list<CGAL::Ray_3<REP> >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " d3-ray");  return str;
}

template<class REP>
void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Ray_3<REP> >& L)
{
}

// 3d output ...
template<class T>
void cgal_Ray_3_d3(const T& L, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GRAPH<leda_d3_point,int> G;
 typename T::const_iterator iter = L.begin();
 
 for(;iter != L.end(); iter++) {
   leda_d3_ray r = convert_to_leda(*iter);
   leda_vector v = r.point2() - r.point1();
   v= v * 50;
   leda_d3_point p1=r.source(), p2=p1 + v;
   leda_node v1 = G.new_node(p1);
   leda_node v2 = G.new_node(p2);   
   leda_edge e1 = G.new_edge(v1,v2);
   leda_edge e2 = G.new_edge(v2,v1);
   G.set_reversal(e1,e2);   
 }
 H.join(G);
}


// 3d Triangles

template<class REP>
const char* leda_tname(CGAL::Triangle_3<REP>* t) {  return "CGALTriangle_3"; }

template<class REP>
leda_window& operator << (leda_window& w, const CGAL::Triangle_3<REP>& obj)
{  
  w << convert_to_leda(obj);
  return w;
}

template<class REP>
leda_window& operator >> (leda_window& w, CGAL::Triangle_3<REP>& obj)
{
  typedef typename REP::RT RT;
  
  CGAL::Triangle_2<REP> th;
  if (w >> th){
    CGAL::Point_2<REP> pa1(th.vertex(0)), pa2(th.vertex(1)), pa3(th.vertex(2));
  
    CGAL::Point_3<REP> p1(pa1.hx(), pa1.hy(), RT(0), pa1.hw());
    CGAL::Point_3<REP> p2(pa2.hx(), pa2.hy(), RT(0), pa2.hw());     
    CGAL::Point_3<REP> p3(pa3.hx(), pa3.hy(), RT(0), pa3.hw());     
    
    obj = CGAL::Triangle_3<REP>(p1, p2, p3);
  }
  return w;
}

template<class REP>
bool geowin_IntersectsBox(const CGAL::Triangle_3<REP>& obj, double x1,double y1,double x2, double y2,bool f)
{  
  return geowin_IntersectsBox(convert_to_leda(obj),x1,y1,x2,y2,f); 
}

template<class REP>
void geowin_BoundingBox(const CGAL::Triangle_3<REP>& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_3 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

template<class REP>
void geowin_Translate(CGAL::Triangle_3<REP>& obj, double dx, double dy)
{
  typedef typename REP::RT RT;

  CGAL::Vector_3<REP> vec;
  vec= CGAL::Vector_3<REP>(RT(dx), RT(dy), RT(0)); 
  
  CGAL::Aff_transformation_3<REP> translate(CGAL::TRANSLATION, vec);
  
  obj = obj.transform(translate);
}

template<class REP>
void geowin_Rotate(CGAL::Triangle_3<REP>& obj, double x, double y, double a)
{  
}

template<class REP>
leda_string geowin_info_fcn(const std::list<CGAL::Triangle_3<REP> >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " 3d-triangle");  return str;
}

template<class REP>
void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Triangle_3<REP> >& L)
{
}


// 3d output ...
template<class T>
void cgal_Triangle_3_d3(const T& L, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GRAPH<leda_d3_point,int> G;
 typename T::const_iterator iter = L.begin();
 
 for(;iter != L.end(); iter++) {
   leda_d3_point p0 = convert_to_leda((*iter).vertex(0));
   leda_d3_point p1 = convert_to_leda((*iter).vertex(1)); 
   leda_d3_point p2 = convert_to_leda((*iter).vertex(2)); 
    
   leda_node v0 = G.new_node(p0);
   leda_node v1 = G.new_node(p1);  
   leda_node v2 = G.new_node(p2); 

   leda_edge e1 = G.new_edge(v0,v1);
   leda_edge e2 = G.new_edge(v1,v0);       
   leda_edge e3 = G.new_edge(v1,v2);
   leda_edge e4 = G.new_edge(v2,v1);
   leda_edge e5 = G.new_edge(v2,v0);
   leda_edge e6 = G.new_edge(v0,v2);
   
   G.set_reversal(e1,e2);   
   G.set_reversal(e3,e4);
   G.set_reversal(e5,e6);
 }
 H.join(G);
}

// 3d tetrahedron ...

template<class REP>
const char* leda_tname(CGAL::Tetrahedron_3<REP>* t) {  return "CGALTetrahedron_3"; }

template<class REP>
leda_window& operator << (leda_window& w, const CGAL::Tetrahedron_3<REP>& t)
{  
  w << CGAL::Segment_3<REP>( t[0], t[1]);
  w << CGAL::Segment_3<REP>( t[1], t[2]);
  w << CGAL::Segment_3<REP>( t[2], t[0]);
  w << CGAL::Segment_3<REP>( t[0], t[3]);
  w << CGAL::Segment_3<REP>( t[1], t[3]);
  w << CGAL::Segment_3<REP>( t[2], t[3]);
  
  return w;
}

template<class REP>
leda_window& operator >> (leda_window& w, CGAL::Tetrahedron_3<REP>& obj)
{
  typedef typename REP::RT   RT;
  
  CGAL::Triangle_3<REP> q;
  w >> q;
  double x0 = CGAL::to_double( q[0].x());
  double y0 = CGAL::to_double( q[0].y());
  double x1, y1;
  w.read_mouse_seg(x0, y0, x1, y1);
  CGAL::Point_3<REP> p = CGAL::Point_3<REP>( RT(x1), RT(y1), RT(0));
  w << CGAL::Segment_3<REP>( q[0], p);
  w << CGAL::Segment_3<REP>( q[1], p);
  w << CGAL::Segment_3<REP>( q[2], p);
  obj =  CGAL::Tetrahedron_3<REP>( q[0], q[1], q[2], p);
  return w;
}

template<class REP>
bool geowin_IntersectsBox(const CGAL::Tetrahedron_3<REP>& obj, double x1,double y1,double x2, double y2,bool f)
{  
  return geowin_IntersectsBox(convert_to_leda(obj),x1,y1,x2,y2,f); 
}

template<class REP>
void geowin_BoundingBox(const CGAL::Tetrahedron_3<REP>& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_3 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

template<class REP>
void geowin_Translate(CGAL::Tetrahedron_3<REP>& obj, double dx, double dy)
{
  typedef typename REP::RT RT;

  CGAL::Vector_3<REP> vec;
  vec= CGAL::Vector_3<REP>(RT(dx), RT(dy), RT(0)); 
  
  CGAL::Aff_transformation_3<REP> translate(CGAL::TRANSLATION, vec);
  
  obj = obj.transform(translate);
}

template<class REP>
void geowin_Rotate(CGAL::Tetrahedron_3<REP>& obj, double x, double y, double a)
{  
}

template<class REP>
leda_string geowin_info_fcn(const std::list<CGAL::Tetrahedron_3<REP> >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %s", L.size(), " 3d-tetrahedra");  return str;
}

template<class REP>
void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Tetrahedron_3<REP> >& L)
{
}

// 3d output ...
template<class T>
void cgal_Tetrahedron_3_d3(const T& L, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GRAPH<leda_d3_point,int> G;
 typename T::const_iterator iter = L.begin();
 
 for(;iter != L.end(); iter++) {
   leda_d3_point p0 = convert_to_leda((*iter).vertex(0));
   leda_d3_point p1 = convert_to_leda((*iter).vertex(1)); 
   leda_d3_point p2 = convert_to_leda((*iter).vertex(2)); 
   leda_d3_point p3 = convert_to_leda((*iter).vertex(3));  
    
   leda_node v0 = G.new_node(p0);
   leda_node v1 = G.new_node(p1);  
   leda_node v2 = G.new_node(p2); 
   leda_node v3 = G.new_node(p3);    

   leda_edge e1 = G.new_edge(v0,v1);
   leda_edge e2 = G.new_edge(v1,v0);       
   leda_edge e3 = G.new_edge(v1,v2);
   leda_edge e4 = G.new_edge(v2,v1);
   leda_edge e5 = G.new_edge(v2,v0);
   leda_edge e6 = G.new_edge(v0,v2);
   
   leda_edge e7 = G.new_edge(v0,v3);
   leda_edge e8 = G.new_edge(v3,v0);       
   leda_edge e9 = G.new_edge(v1,v3);
   leda_edge e10 = G.new_edge(v3,v1);
   leda_edge e11 = G.new_edge(v2,v3);
   leda_edge e12 = G.new_edge(v3,v2);
      
   G.set_reversal(e1,e2);   
   G.set_reversal(e3,e4);
   G.set_reversal(e5,e6);
   G.set_reversal(e7,e8);   
   G.set_reversal(e9,e10);
   G.set_reversal(e11,e12);   
 }
 H.join(G);
}


GEOWIN_END_NAMESPACE

#endif

#endif // CGAL_GEOWIN_SUPPORT




