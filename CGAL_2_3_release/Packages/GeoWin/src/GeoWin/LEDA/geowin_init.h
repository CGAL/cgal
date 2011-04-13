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
// file          : src/GeoWin/LEDA/geowin_init.h
// package       : GeoWin (1.2.2)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.2.2
// revision_date : 30 January 2000 
// author(s)     : Matthias Baesken, Ulrike Bartuschka, Stefan Naeher
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ======================================================================



#ifndef LEDA_GEOWIN_ADDITIONAL_H
#define LEDA_GEOWIN_ADDITIONAL_H

#if !defined(LEDA_ROOT_INCL_ID)
#define LEDA_ROOT_INCL_ID 400951
#include <LEDA/REDEFINE_NAMES.h>
#endif

#include <LEDA/geowin.h>
#include <LEDA/d3_rat_point.h>


GEOWIN_BEGIN_NAMESPACE

bool geowin_IntersectsBox(const point&, double, double, double, double, bool); 
bool geowin_IntersectsBox(const segment&, double, double, double, double, bool );
bool geowin_IntersectsBox(const ray&, double, double, double, double, bool);
bool geowin_IntersectsBox(const line&, double, double, double, double, bool);
bool geowin_IntersectsBox(const circle&, double, double, double, double, bool);
bool geowin_IntersectsBox(const polygon&,  double, double, double, double, bool);
bool geowin_IntersectsBox(const gen_polygon&,  double, double, double, double, bool);
bool geowin_IntersectsBox(const d3_point&, double, double, double, double, bool);
bool geowin_IntersectsBox(const rectangle&, double, double, double, double, bool);
#if (__LEDA__ >= 420)
bool geowin_IntersectsBox(const triangle&, double, double, double, double, bool);
#endif


void geowin_BoundingBox(const point&, double&, double&, double&, double&);
void geowin_BoundingBox(const segment&, double&, double&, double&, double&);
void geowin_BoundingBox(const ray&, double&, double&, double&, double&);
void geowin_BoundingBox(const line&, double&, double&, double&, double&);
void geowin_BoundingBox(const circle&, double&, double&, double&, double&);
void geowin_BoundingBox(const polygon&, double&, double&, double&, double&);
void geowin_BoundingBox(const gen_polygon&, double&, double&, double&, double&);
void geowin_BoundingBox(const d3_point&, double&, double&, double&, double&);
void geowin_BoundingBox(const rectangle&, double&, double&, double&, double&);
#if (__LEDA__ >= 420)
void geowin_BoundingBox(const triangle&, double&, double&, double&, double&);
#endif


void geowin_Translate(point& obj, double dx, double dy);
void geowin_Translate(segment& obj, double dx, double dy);
void geowin_Translate(ray& obj, double dx, double dy);
void geowin_Translate(line& obj, double dx, double dy);
void geowin_Translate(circle& obj, double dx, double dy);
void geowin_Translate(polygon& obj, double dx, double dy);
void geowin_Translate(gen_polygon& obj, double dx, double dy);
void geowin_Translate(d3_point& p, double dx, double dy);
void geowin_Translate(rectangle& obj, double dx, double dy);
#if (__LEDA__ >= 420)
void geowin_Translate(triangle& obj, double dx, double dy);
#endif

// move point functions ...
void geowin_Translatepoint_seg(segment& obj, double dx, double dy, int pnr);
void geowin_Translatepoint_ray(ray& obj, double dx, double dy, int pnr);
void geowin_Translatepoint_line(line& obj, double dx, double dy, int pnr);
void geowin_Translatepoint_circle(circle& obj, double dx, double dy, int pnr);
void geowin_Translatepoint_poly(polygon& obj, double dx, double dy, int pnr);
void geowin_Translatepoint_gpoly(gen_polygon& obj, double dx, double dy, int pnr);


void geowin_Rotate(point& obj, double x, double y, double a);
void geowin_Rotate(segment& obj, double x, double y, double a);
void geowin_Rotate(ray& obj, double x, double y, double a);
void geowin_Rotate(line& obj, double x, double y, double a);
void geowin_Rotate(circle& obj, double x, double y, double a);
void geowin_Rotate(polygon& obj, double x, double y, double a);
void geowin_Rotate(gen_polygon& obj, double x, double y, double a);
void geowin_Rotate(d3_point& p, double x, double y, double a);
void geowin_Rotate(rectangle& obj, double x, double y, double a);
#if (__LEDA__ >= 420)
void geowin_Rotate(triangle& obj, double x, double y, double a);
#endif


void geowin_generate_objects(GeoWin& gw, list<point>& L);
void geowin_generate_objects(GeoWin& gw, list<segment>& S);
void geowin_generate_objects(GeoWin& gw, list<ray>& L);
void geowin_generate_objects(GeoWin& gw, list<line>& L);
void geowin_generate_objects(GeoWin& gw, list<circle>& C);
void geowin_generate_objects(GeoWin& gw, list<polygon>& Pl);
void geowin_generate_objects(GeoWin& gw, list<gen_polygon>& Pl);
void geowin_generate_objects(GeoWin& gw, list<d3_point>& Pl);
void geowin_generate_objects(GeoWin& gw, list<rectangle>& Pl);
#if (__LEDA__ >= 420)
void geowin_generate_objects(GeoWin& gw, list<triangle>& Pl);
#endif


#if !defined(NO_RAT_ALGORITHMS)

void geowin_Translate(rat_point& obj, double dx, double dy);
void geowin_Translate(rat_segment& obj, double dx, double dy);
void geowin_Translate(rat_ray& obj, double dx, double dy);
void geowin_Translate(rat_line& obj, double dx, double dy);
void geowin_Translate(rat_circle& obj, double dx, double dy);
void geowin_Translate(rat_polygon& obj, double dx, double dy);
void geowin_Translate(rat_gen_polygon& obj, double dx, double dy);
void geowin_Translate(d3_rat_point& p, double dx, double dy);
void geowin_Translate(rat_rectangle& obj, double dx, double dy);
#if (__LEDA__ >= 420)
void geowin_Translate(rat_triangle& obj, double dx, double dy);
#endif

// move point functions ...
void geowin_Translatepoint_rat_seg(rat_segment& obj, double dx, double dy, int pnr);
void geowin_Translatepoint_rat_ray(rat_ray& obj, double dx, double dy, int pnr);
void geowin_Translatepoint_rat_line(rat_line& obj, double dx, double dy, int pnr);
void geowin_Translatepoint_rat_circle(rat_circle& obj, double dx, double dy, int pnr);
void geowin_Translatepoint_rat_poly(rat_polygon& obj, double dx, double dy, int pnr);
void geowin_Translatepoint_rat_gpoly(rat_gen_polygon& obj, double dx, double dy, int pnr);


void geowin_Rotate(rat_point& obj, double x, double y, double a);
void geowin_Rotate(rat_segment& obj, double x, double y, double a);
void geowin_Rotate(rat_ray& obj, double x, double y, double a);
void geowin_Rotate(rat_line& obj, double x, double y, double a);
void geowin_Rotate(rat_circle& obj, double x, double y, double a);
void geowin_Rotate(rat_polygon& obj, double x, double y, double a);
void geowin_Rotate(rat_gen_polygon& obj, double x, double y,double a);
void geowin_Rotate(d3_rat_point& p, double x, double y, double a);
void geowin_Rotate(rat_rectangle& obj, double x, double y, double a);
#if (__LEDA__ >= 420)
void geowin_Rotate(rat_triangle& obj, double x, double y, double a);
#endif


bool geowin_IntersectsBox(const rat_point& p, double, double, double, double, bool);
bool geowin_IntersectsBox(const rat_segment& s, double, double, double, double, bool);
bool geowin_IntersectsBox(const rat_ray& r, double, double, double, double, bool); 
bool geowin_IntersectsBox(const rat_line& l, double, double, double, double, bool); 
bool geowin_IntersectsBox(const rat_circle& c, double, double, double, double, bool); 
bool geowin_IntersectsBox(const rat_polygon& p, double, double, double, double, bool);
bool geowin_IntersectsBox(const rat_gen_polygon& p, double, double, double, double, bool);
bool geowin_IntersectsBox(const d3_rat_point&, double, double, double, double, bool);
bool geowin_IntersectsBox(const rat_rectangle&, double, double, double, double, bool);
#if (__LEDA__ >= 420)
bool geowin_IntersectsBox(const rat_triangle&, double, double, double, double, bool);
#endif


void geowin_BoundingBox(const rat_point&, double&, double&, double&, double&);
void geowin_BoundingBox(const rat_segment&, double&, double&, double&, double&);
void geowin_BoundingBox(const rat_ray&, double&, double&, double&, double&);
void geowin_BoundingBox(const rat_line&, double&, double&, double&, double&);
void geowin_BoundingBox(const rat_circle&, double&, double&, double&, double&);
void geowin_BoundingBox(const rat_polygon&, double&, double&, double&, double&);
void geowin_BoundingBox(const rat_gen_polygon&, double&, double&, double&, double&);
void geowin_BoundingBox(const d3_rat_point&, double&, double&, double&, double&);
void geowin_BoundingBox(const rat_rectangle&, double&, double&, double&, double&);
#if (__LEDA__ >= 420)
void geowin_BoundingBox(const rat_triangle&, double&, double&, double&, double&);
#endif


void geowin_generate_objects(GeoWin& gw, list<rat_point>& L);
void geowin_generate_objects(GeoWin& gw, list<rat_segment>& L);
void geowin_generate_objects(GeoWin& gw, list<rat_ray>& L);
void geowin_generate_objects(GeoWin& gw, list<rat_line>& L);
void geowin_generate_objects(GeoWin& gw, list<rat_circle>& L);
void geowin_generate_objects(GeoWin& gw, list<rat_polygon>& L);
void geowin_generate_objects(GeoWin& gw, list<rat_gen_polygon>& L);
void geowin_generate_objects(GeoWin& gw, list<d3_rat_point>& L); 
void geowin_generate_objects(GeoWin& gw, list<rat_rectangle>& L); 
#if (__LEDA__ >= 420)
void geowin_generate_objects(GeoWin& gw, list<rat_triangle>& L);
#endif

#endif




template<class T>
void geowin_fredraw_fcn(const list<T>& L, window& w, color c1, color c2, 
		     bool filled,  double x1, double y1, double x2, double y2)
{
  T t;  
  forall(t, L) 
    if ( geowin_IntersectsBox(t, x1, y1, x2, y2,filled) )
    { if (filled)
         w.set_fill_color(c2);
      else
         w.set_fill_color(invisible);
      w.set_color(c1);
      w << t;
     }
}


template<class T>
string geowin_info_fcn(const list<T>& L)
{ int n = L.length();
  if (n == 1)
    return string("\\black\\tt %d %s", n, leda_tname((T*)0));
  else
    return string("\\black\\tt %d %ss",n, leda_tname((T*)0));
}

/*{\Manpage {GeoWin} {} {Geometry Windows - Initialization}}*/

 /*{\Mtext
    \medskip
    The following non-member functions can be used to register additional containers storing
    geometric objects for usage in $GeoWin$. Note that they have to be called before
    creating a $GeoWin$.
    
    |template<class T>| 
    |GeoEditScene<T>* geowin_init_default_type( T* t, string str);|

    |template<class T>| 
    |GeoEditScene<T>* geowin_init_basic_type(T* t,string str);|

    |t| is a pointer to the container and |str| is the name that will be assoziated with
    this container in |GeoWin|.
    Note that the container of |T| type has to provide
    \begin{itemize}
    \item STL-style iteration
    \item |value_type| - type of the values the container stores
    \item |iterator|
    \item operations |begin| and |end| returning iterators, that can the used for traversal
          of the container
    \item operation |void push_back(const T::value_type&)| for inserting elements at the end of
          the container
    \item operation |iterator insert(iterator, const T::value_type&)| for inserting elements
    \item operation |void erase(iterator it)| for deleting elements in the container at position |it|
    \item operation |bool empty()| returning |true| if the container is empty, |false| otherwise
    \end{itemize} 
    
    The following functions have to be provided as well for the geometric objects (type |G|)
    and the container of type |CONT|:
    \begin{itemize}
    \item operators for stream and LEDA window I/O and for writing to |ps_file| \\
     |ostream& operator<<(ostream&, const G&);| \\
     |istream& operator>>(istream&, G&);| \\
     |window& operator<<(window&, const G&);| \\
     |window& operator>>(window&, G&);| \\
     |ps_file& operator<<(ps_file&, const G&);| \\
    \item translate and rotate operations for translating/rotating objects\\
     |void geowin_Translate(G& o, double dx, double dy);| \\
     |void geowin_Rotate(G& o, double x, double y, double phi);| \\
    \item simple geometric queries for getting the bounding box of an object
          and queriing, if an object intersects a box \\
     |bool geowin_BoundingBox(const G& o, double& x0, double& y0, double& x1, double& y1);| \\
     |bool geowin_IntersectsBox(const G& o, double x0, double y0, double x1, double y1);|
     \item other functions: object generation and info function \\
      |template<class CONT>|
      |void geowin_generate_objects(GeoWin& gw, CONT& C);|

      |template<class CONT>|
      |string geowin_info_fcn(const CONT& C);|
   \end{itemize}
   
   Note that for the |geowin_init_basic_type|-function the functions
   |geowin_Translate|, |geowin_Rotate|, |geowin_generate_objects| and |geowin_info_fcn| are not needed.
 }*/
 
#if !defined(__SUNPRO_CC)
template<class TRAITS> 
GeoEditScene<__typename TRAITS::CONTAINER>* geowin_init_default_type(TRAITS tr)
{
  typedef typename TRAITS::CONTAINER  CONT;
  
  CONT* t;
  string str = tr.get_name();
  GeoBaseScene<CONT>* bsc= make_base_prototype(t, str + string("Base"));
  GeoEditScene<CONT>* sc = make_edit_prototype(t, str);

  sc->set_info_fcn(tr.geowin_info_fcn);
 
  bsc->set_box_intersection_fcn(tr.geowin_IntersectsBox);
  bsc->set_get_bounding_box_fcn(tr.geowin_BoundingBox);  
  
  sc->set_box_intersection_fcn(tr.geowin_IntersectsBox);
  sc->set_get_bounding_box_fcn(tr.geowin_BoundingBox);
  sc->set_move_fcn(tr.geowin_Translate);
  sc->set_rotate_fcn(tr.geowin_Rotate);
  sc->generate_fcn = tr.geowin_generate_objects;
  
  return sc;
}
#endif

template<class T> 
GeoEditScene<T>* geowin_init_default_type( T* t, string str)
{ 
  //links a string (the scene-typename) with the real type
  GeoBaseScene< T >* bsc= make_base_prototype(t, str+string("Base"));
  GeoEditScene< T >* sc = make_edit_prototype(t, str);

#if !defined GEOWIN_USE_NAMESPACE
  string (*f)(const T&) = geowin_info_fcn;
  sc->set_info_fcn(f);
  
  bsc->set_box_intersection_fcn(geowin_IntersectsBox);
  bsc->set_get_bounding_box_fcn(geowin_BoundingBox);
    
  sc->set_box_intersection_fcn(geowin_IntersectsBox);
  sc->set_get_bounding_box_fcn(geowin_BoundingBox);
  sc->set_move_fcn(geowin_Translate);
  sc->set_rotate_fcn(geowin_Rotate);
  sc->generate_fcn = geowin_generate_objects;
#else
  string (*f)(const T&) = GEOWIN_NAMESPACE_NAME::geowin_info_fcn;
  sc->set_info_fcn(f);

  bsc->set_box_intersection_fcn(GEOWIN_NAMESPACE_NAME::geowin_IntersectsBox);
  bsc->set_get_bounding_box_fcn(GEOWIN_NAMESPACE_NAME::geowin_BoundingBox);
  
  sc->set_box_intersection_fcn(GEOWIN_NAMESPACE_NAME::geowin_IntersectsBox);
  sc->set_get_bounding_box_fcn(GEOWIN_NAMESPACE_NAME::geowin_BoundingBox);
  sc->set_move_fcn(GEOWIN_NAMESPACE_NAME::geowin_Translate);
  sc->set_rotate_fcn(GEOWIN_NAMESPACE_NAME::geowin_Rotate);
  sc->generate_fcn = GEOWIN_NAMESPACE_NAME::geowin_generate_objects;
#endif

  return sc;
} 

template<class T,class F> 
GeoEditScene<T>* geowin_init_default_type( T* t, string str, F d3_f)
{ 
  //links a string (the scene-typename) with the real type
  GeoBaseScene< T >* bsc = make_base_prototype(t, str+string("Base"));
  bsc->set_objects_init_d3_window(d3_f);
 
  GeoEditScene< T >* sc = make_edit_prototype(t, str);
  sc->set_objects_init_d3_window(d3_f);

#if !defined GEOWIN_USE_NAMESPACE
  string (*f)(const T&) = geowin_info_fcn;
  sc->set_info_fcn(f);
  
  bsc->set_box_intersection_fcn(geowin_IntersectsBox);
  bsc->set_get_bounding_box_fcn(geowin_BoundingBox);  
  
  sc->set_box_intersection_fcn(geowin_IntersectsBox);
  sc->set_get_bounding_box_fcn(geowin_BoundingBox);
  sc->set_move_fcn(geowin_Translate);
  sc->set_rotate_fcn(geowin_Rotate);
  sc->generate_fcn = geowin_generate_objects;
#else
  string (*f)(const T&) = GEOWIN_NAMESPACE_NAME::geowin_info_fcn;
  sc->set_info_fcn(f);

  bsc->set_box_intersection_fcn(GEOWIN_NAMESPACE_NAME::geowin_IntersectsBox);
  bsc->set_get_bounding_box_fcn(GEOWIN_NAMESPACE_NAME::geowin_BoundingBox);  
  
  sc->set_box_intersection_fcn(GEOWIN_NAMESPACE_NAME::geowin_IntersectsBox);
  sc->set_get_bounding_box_fcn(GEOWIN_NAMESPACE_NAME::geowin_BoundingBox);
  sc->set_move_fcn(GEOWIN_NAMESPACE_NAME::geowin_Translate);
  sc->set_rotate_fcn(GEOWIN_NAMESPACE_NAME::geowin_Rotate);
  sc->generate_fcn = GEOWIN_NAMESPACE_NAME::geowin_generate_objects;
#endif

  return sc;
} 

template<class T> 
GeoEditScene<T>* geowin_init_basic_type(T* t,string str)
{
  GeoBaseScene< T >* bsc= make_base_prototype(t, str+string("Base"));
  GeoEditScene< T >* sc = make_edit_prototype(t, str);

#if !defined GEOWIN_USE_NAMESPACE  
  bsc->set_box_intersection_fcn(geowin_IntersectsBox);
  bsc->set_get_bounding_box_fcn(geowin_BoundingBox);

  sc->set_box_intersection_fcn(geowin_IntersectsBox);
  sc->set_get_bounding_box_fcn(geowin_BoundingBox);
#else
  bsc->set_box_intersection_fcn(GEOWIN_NAMESPACE_NAME::geowin_IntersectsBox);
  bsc->set_get_bounding_box_fcn(GEOWIN_NAMESPACE_NAME::geowin_BoundingBox);

  sc->set_box_intersection_fcn(GEOWIN_NAMESPACE_NAME::geowin_IntersectsBox);
  sc->set_get_bounding_box_fcn(GEOWIN_NAMESPACE_NAME::geowin_BoundingBox);
#endif

  return sc;
}

template<class T,class F> 
GeoEditScene<T>* geowin_init_basic_type(T* t,string str, F d3_f)
{
  GeoBaseScene< T >* bsc = make_base_prototype(t, str+string("Base"));
  bsc->set_objects_init_d3_window(d3_f);
 
  GeoEditScene< T >* sc = make_edit_prototype(t, str);
  sc->set_objects_init_d3_window(d3_f);
 
#if !defined GEOWIN_USE_NAMESPACE  
  bsc->set_box_intersection_fcn(geowin_IntersectsBox);
  bsc->set_get_bounding_box_fcn(geowin_BoundingBox);

  sc->set_box_intersection_fcn(geowin_IntersectsBox);
  sc->set_get_bounding_box_fcn(geowin_BoundingBox);
#else
  bsc->set_box_intersection_fcn(geowin_IntersectsBox);
  bsc->set_get_bounding_box_fcn(geowin_BoundingBox);

  sc->set_box_intersection_fcn(GEOWIN_NAMESPACE_NAME::geowin_IntersectsBox);
  sc->set_get_bounding_box_fcn(GEOWIN_NAMESPACE_NAME::geowin_BoundingBox);
#endif

  return sc;
}



GEOWIN_END_NAMESPACE

#if LEDA_ROOT_INCL_ID == 400951
#undef LEDA_ROOT_INCL_ID
#include <LEDA/UNDEFINE_NAMES.h>
#endif


#endif // _GEOWIN_ADDITIONAL_

