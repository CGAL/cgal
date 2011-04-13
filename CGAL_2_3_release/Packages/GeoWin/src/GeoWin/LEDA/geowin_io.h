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
// file          : src/GeoWin/LEDA/geowin_io.h
// package       : GeoWin (1.2.2)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.2.2
// revision_date : 30 January 2000 
// author(s)     : Matthias Baesken, Ulrike Bartuschka, Stefan Naeher
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ======================================================================



#ifndef LEDA_GEOWIN_IO_H
#define LEDA_GEOWIN_IO_H

#if !defined(LEDA_ROOT_INCL_ID)
#define LEDA_ROOT_INCL_ID 450096
#include <LEDA/REDEFINE_NAMES.h>
#endif

#include <LEDA/geowin.h>
#include <LEDA/geowin_init_d3.h>

// ----------------------------------------------------------------------
//   importers and exporters for LEDA types
// ----------------------------------------------------------------------

template<class T>
void geowin_load_leda_objects(list<T>& L, ifstream& is)
{
  for(;;)
  { char c;
    while (is.get(c) && isspace(c));    if (!is) break;   is.putback(c);
    T obj;  is >> obj;  L.append(obj);
  }     
}


bool geowin_identify_file(ifstream& is, string& type_name, string& sc_name)
{
     // read header ...
     string identifier;
     double version;
  
     char* t_name = new char[200];
     char* s_name = new char[200];  
     is >> identifier >> version;   
     
     if( !(identifier == string("GEOWIN")) ) return false;

     // read rest of line before type name ...
     is.getline(t_name,200);  
     // read type name
     is.getline(t_name,200);    
     type_name = string(t_name);
     // read scene name
     is.getline(s_name,200);  
     sc_name = string(s_name);

     delete t_name; delete s_name;
     return true;   
}

void geowin_get_type_names(point test_obj,string& rat_type_name,string& float_type_name)
{ rat_type_name="RationalPoints"; float_type_name="Points"; }

void geowin_get_type_names(segment test_obj,string& rat_type_name,string& float_type_name)
{ rat_type_name="RationalSegments"; float_type_name="Segments"; }

void geowin_get_type_names(line test_obj,string& rat_type_name,string& float_type_name)
{ rat_type_name="RationalLines"; float_type_name="Lines"; }

void geowin_get_type_names(ray test_obj,string& rat_type_name,string& float_type_name)
{ rat_type_name="RationalRays"; float_type_name="Rays"; }

void geowin_get_type_names(circle test_obj,string& rat_type_name,string& float_type_name)
{ rat_type_name="RationalCircles"; float_type_name="Circles"; }

void geowin_get_type_names(triangle test_obj,string& rat_type_name,string& float_type_name)
{ rat_type_name="RationalTriangles"; float_type_name="Triangles"; }

void geowin_get_type_names(rectangle test_obj,string& rat_type_name,string& float_type_name)
{ rat_type_name="RationalRectangles"; float_type_name="Rectangles"; }

void geowin_get_type_names(polygon test_obj,string& rat_type_name,string& float_type_name)
{ rat_type_name="RationalSimplePolygons"; float_type_name="SimplePolygons"; }

void geowin_get_type_names(gen_polygon test_obj,string& rat_type_name,string& float_type_name)
{ rat_type_name="RationalGeneralizedPolygons"; float_type_name="GeneralizedPolygons"; }

void geowin_get_type_names(d3_point test_obj,string& rat_type_name,string& float_type_name)
{ rat_type_name="D3-RationalPoints"; float_type_name="D3-Points"; }

// TF ... float kernel type, TR ... rat kernel type

template<class TF, class TR>
class geowin_import_leda_objects : public geowin_import {
  virtual void operator()(geo_scene sc, string filename) { 
     ifstream is(filename);    
     if (! is) { error_handler(0,"No such file or directory."); return; }
     
     // what kind of scene is it ??
     if (sc==NULL) return;
     
     TF test_obj;
     
     string rat_type_name, float_type_name;
     geowin_get_type_names(test_obj, rat_type_name, float_type_name);
     
     string tn = sc->get_type_name();
     
     if ((tn != rat_type_name) && (tn != float_type_name)) return;
  
     list<TF>   f_objs;
     list<TR>   r_objs;
     
     string type_name, sc_name;
     
     bool b = geowin_identify_file(is, type_name, sc_name);
     if (! b) return;
     
     if ((type_name != rat_type_name) && (type_name != float_type_name)) return;
     
     if (type_name == rat_type_name) geowin_load_leda_objects(r_objs, is);
     else geowin_load_leda_objects(f_objs, is);
     
     GeoWin* gw = get_geowin(sc);
     
     if (tn == rat_type_name){ 
       if (type_name == rat_type_name) geowin_set_objects(*gw, sc, r_objs);
       else {  // convert ...
         TF iter; 
	 forall(iter, f_objs) r_objs.append(TR(iter));
	 geowin_set_objects(*gw, sc, r_objs);
       }
     }
     else { 
       if (type_name == float_type_name) geowin_set_objects(*gw, sc, f_objs);
       else {  // convert ...
         TR riter;
	 forall(riter, r_objs) f_objs.append(riter.to_float());
         geowin_set_objects(*gw, sc, f_objs);
       }     
     }
  }
};



// d3 objects

template<class TF, class TR>
class geowin_import_leda_d3_objects : public geowin_import {
  virtual void operator()(geo_scene sc, string filename) { 
     ifstream is(filename);    
     if (! is) { error_handler(0,"No such file or directory."); return; }  
  
     // what kind of scene is it ??
     if (sc==NULL) return;
     
     TF test_obj;
     
     string rat_type_name, float_type_name;
     geowin_get_type_names(test_obj, rat_type_name, float_type_name);
     
     string tn = sc->get_type_name();
     
     if ((tn != rat_type_name) && (tn != float_type_name)) return;
  
     list<TF>   f_objs;
     list<TR>   r_objs;
     
     string type_name, sc_name;
     
     bool b = geowin_identify_file(is, type_name, sc_name);
     if (! b) return;
     
     if ((type_name != rat_type_name) && (type_name != float_type_name)) return;
     
     if (type_name == rat_type_name) geowin_load_leda_objects(r_objs, is);
     else geowin_load_leda_objects(f_objs, is);
     
     GeoWin* gw = get_geowin(sc);
     
     if (tn == rat_type_name){ 
       if (type_name == rat_type_name) geowin_set_objects(*gw, sc, r_objs);
       else {  // convert ...
         TF iter; 
	 forall(iter, f_objs) {
	   r_objs.append(leda_convert_to(iter));
	 }  
	 geowin_set_objects(*gw, sc, r_objs);
       }
     }
     else { 
       if (type_name == float_type_name) geowin_set_objects(*gw, sc, f_objs);
       else {  // convert ...
         TR riter;
	 forall(riter, r_objs) f_objs.append(riter.to_float());
         geowin_set_objects(*gw, sc, f_objs);
       }     
     }
  }
};


// ---------------------------------------------------------------------
//   import CGAL objects
// ---------------------------------------------------------------------

string geowin_get_cgal_type_name(point test_obj){ return string("CGALPointList"); }
string geowin_get_cgal_type_name(segment test_obj){ return string("CGALSegmentList"); }
string geowin_get_cgal_type_name(line test_obj){ return string("CGALLineList"); }
string geowin_get_cgal_type_name(ray test_obj){ return string("CGALRayList"); }
string geowin_get_cgal_type_name(circle test_obj){ return string("CGALCircleList"); }

string geowin_get_cgal_type_name(triangle test_obj){ return string("CGALTriangleList"); }
string geowin_get_cgal_type_name(rectangle test_obj){ return string("CGALRectangleList"); }

string geowin_get_cgal_type_name(polygon test_obj){ return string("CGALPolygonList"); }
string geowin_get_cgal_type_name(gen_polygon test_obj){ return string("CGALPolygonList"); }
string geowin_get_cgal_type_name(d3_point test_obj){ return string("CGALPoint_3_List"); }


string geowin_cgal_eq(string tn)
{
  if ((tn=="Points") || (tn=="RationalPoints"))     return string("CGALPointList");
  if ((tn=="Segments") || (tn=="RationalSegments")) return string("CGALSegmentList");
  if ((tn=="Lines") || (tn=="RationalLines"))       return string("CGALLineList");
  if ((tn=="Rays") || (tn=="RationalRays"))         return string("CGALRayList");
  if ((tn=="Circles") || (tn=="RationalCircles"))   return string("CGALCircleList");
  
  if ((tn=="Triangles") || (tn=="RationalTriangles"))                     return string("CGALTriangleList");
  if ((tn=="Rectangles") || (tn=="RationalRectangles"))                   return string("CGALRectangleList");  
  
  if ((tn=="SimplePolygons") || (tn=="RationalSimplePolygons"))           return string("CGALPolygonList");
  if ((tn=="GeneralizedPolygons") || (tn=="RationalGeneralizedPolygons")) return string("CGALPolygonList");
  if ((tn=="D3-Points") || (tn=="D3-RationalPoints"))                     return string("CGALPoint_3_List");
  
  return string("");
}

// --------------------------------------------------------------
// import routines for CGAL objects
// --------------------------------------------------------------


void geowin_read_in(ifstream& is, point& obj, int repres)
{
  double x, y, w;
  if (repres==0) { // Cartesian
    is >> x; is >> y; 
    obj = point(x,y);
  }
  else { // Homogeneous
    is >> x; is >> y; is >> w;
    obj = point(x,y,w);  
  }
}

void geowin_read_in(ifstream& is, segment& obj, int repres)
{
  double x1, y1, w1, x2, y2, w2;
  if (repres==0) { // Cartesian
    is >> x1; is >> y1; is >> x2; is >> y2;
    obj = segment(x1,y1,x2,y2);
  }
  else { // Homogeneous
    is >> x1; is >> y1; is >> w1; is >> x2; is >> y2; is >> w2;    
    obj = segment(x1/w1,y1/w1,x2/w2,y2/w2);  
  }
}

void geowin_read_in(ifstream& is, circle& obj, int repres)
{
  double xm, ym, sr;
  int ori;
  if (repres==0) { // Cartesian
    is >> xm; is >> ym; is >> sr; is >> ori;
    obj = circle(xm,ym, sqrt(sr));
    if (ori == -1) obj = obj.reverse();
  }
  else { // Homogeneous
    //is >> x1; is >> y1; is >> w1; is >> x2; is >> y2; is >> w2;    
    //obj = segment(x1/w1,y1/w1,x2/w2,y2/w2);  
  }
}


void geowin_read_in(ifstream& is, line& obj, int repres)
{
  double a, b, c;

 
    is >> a; is >> b; is >> c; 
    double x1, y1, x2, y2;
    
    if (b==0) { x1 = -c/a; x2=x1; y1 = 0; y2 = 100; }
    else {
      x1 = 10; x2 = 100;
      y1 = (-c - a*x1)/b;
      y2 = (-c - a*x2)/b;
    }
    obj = line(point(x1,y1), point(x2,y2));
}


void geowin_read_in(ifstream& is, triangle& obj, int repres)
{  
  point p1, p2, p3;
  //cout << "read triangle !\n";
  geowin_read_in(is, p1, repres);
  geowin_read_in(is, p2, repres);
  geowin_read_in(is, p3, repres);
  obj = triangle(p1, p2, p3);
}

void geowin_read_in(ifstream& is, rectangle& obj, int repres)
{  
  point p1, p2;
 
  geowin_read_in(is, p1, repres);
  geowin_read_in(is, p2, repres); 
  obj = rectangle(p1, p2);
}


void geowin_read_in(ifstream& is, polygon& obj, int repres)
{
  int anz,i;
 
  is >> anz;
  
  list<point> poly_points;
  point pact;
  
  for(i=0; i<anz; i++) {
    geowin_read_in(is, pact, repres);
    poly_points.append(pact);
  }
  
  obj = polygon(poly_points);
}

void geowin_read_in(ifstream& is, gen_polygon& obj, int repres)
{
  int anz,i;
 
  is >> anz;
  
  list<point> poly_points;
  point pact;
  
  for(i=0; i<anz; i++) {
    geowin_read_in(is, pact, repres);
    poly_points.append(pact);
  }
  
  obj = gen_polygon(poly_points);
}


template<class TF, class TR>
void geowin_load_cgal_objects(int repres, ifstream& is, list<TF>& f_objs, list<TR>& r_objs)
{
  for(;;)
  { char c;
    while (is.get(c) && isspace(c));    if (!is) break;   is.putback(c);
    
    TF obj;
    geowin_read_in(is, obj, repres);
  
    f_objs.append(obj);  
    r_objs.append(TR(obj)); // this is still a problem with the d3 kernel types ...
  }       
}


template<class TF, class TR>
class geowin_import_cgal_objects : public geowin_import {
  virtual void operator()(geo_scene sc, string filename) { 
     ifstream is(filename);    
     if (! is) { error_handler(0,"No such file or directory."); return; }  
  
     // what kind of scene is it ??
     if (sc==NULL) return;
     
     TF test_obj;
     
     string rat_type_name, float_type_name;
     geowin_get_type_names(test_obj, rat_type_name, float_type_name);
     
     string cgal_type_name;
     cgal_type_name = geowin_get_cgal_type_name(test_obj);
     
     string tn = sc->get_type_name();
     
     if (geowin_cgal_eq(tn) != cgal_type_name) return;
  
     list<TF>   f_objs;
     list<TR>   r_objs;
     
     string type_name, sc_name;
     
     //cout << "identify file!\n"; cout.flush();
     bool b = geowin_identify_file(is, type_name, sc_name);
     if (! b) return;
     
     if (type_name != cgal_type_name) return;
     
     // Homogeneous or Cartesian representation ??
     // ask user
     GeoWin* gw = get_geowin(sc);
     window& win = gw->get_window();
     
     string rep[2];
     rep[0] = string("Cartesian");
     rep[1] = string("Homogeneous");
     int repres = win.read_panel("Cartesian or Homogeneous Representation?", 2, rep);
     
     geowin_load_cgal_objects(repres, is, f_objs, r_objs);
     
     if (tn == rat_type_name) geowin_set_objects(*gw, sc, r_objs);
     else  geowin_set_objects(*gw, sc, f_objs);
  }
};

// -------------------------------------------------------------------------------------------
// import routines for ESRI shapefile objects ...
// -------------------------------------------------------------------------------------------


void geowin_esri_read_int(ifstream& I, unsigned int& number);
void geowin_esri_read_int_big(ifstream& I, unsigned int& number);

bool geowin_esri_read_header(ifstream& I, string& shapetype, unsigned int& length)
{
#if defined(__SUNPRO_CC) || defined(__KCC) || defined(_MSC_VER)
 char ch;
#else 
 unsigned char ch;
#endif
 
 int i;
 unsigned int code;
 geowin_esri_read_int_big(I, code);
 
 if (code != 9994) return false; 
 
 // 20 unused bytes
 for (i=0;i<20;i++) I.get(ch);  
 
 // file length
 geowin_esri_read_int_big(I, length); length=length*2;
 
 // now little endian format ...
 unsigned int version; 
 geowin_esri_read_int(I, version);
 
 //shape type ...
 unsigned int shape_type; 
 geowin_esri_read_int(I, shape_type);
 
 bool ret = false;
 
 switch(shape_type){ 
   case 0: {  cout << "Null shape!\n"; shapetype = "Null shape"; ret = true; break; }
   case 1: {  cout << "Point!\n"; shapetype = "Point"; ret = true; break; }
   case 3: {  cout << "PolyLine!\n"; shapetype = "PolyLine"; ret = true; break; }
   case 5: {  cout << "Polygon!\n"; shapetype = "Polygon"; ret = true; break; }
   case 8: {  cout << "MultiPoint!\n"; shapetype = "MultiPoint"; ret = true; break; }
   case 11:{  cout << "PointZ!\n"; shapetype = "PointZ"; ret = true; break; }
   case 13:{  cout << "PolyLineZ!\n"; shapetype = "PolyLineZ"; ret = true; break; }
   case 15:{  cout << "PolygonZ!\n"; shapetype = "PolygonZ"; ret = true; break; }
   case 18:{  cout << "MultiPointZ!\n"; shapetype = "MultiPointZ"; ret = true; break; }
   case 21:{  cout << "PointM!\n"; shapetype = "PointM"; ret = true; break; }
   case 23:{  cout << "PolyLineM!\n"; shapetype = "PolyLineM"; ret = true; break; }
   case 25:{  cout << "PolygonM!\n"; shapetype = "PolygonM"; ret = true; break; }
   case 28:{  cout << "MultiPointM!\n"; shapetype = "MultiPointM"; ret = true; break;}
   case 31:{  cout << "MultiPatch!\n"; shapetype = "MultiPatch"; ret = true; break; }
   default: { cout << "Unknown shape type!\n"; ret = false; break; }
 }
 
 // read in bounding box parameters
 int cnt;
 for (cnt=0;cnt < 8; cnt++){
  for(i=0;i<8;i++) I.get(ch); 
 }
  
 return ret;
}

bool geowin_read_esri_record(ifstream& I, unsigned int& clength, unsigned int& prev, bool run)
{  
  unsigned int nr; 
  if (I.eof()) return false;
  
  // get record number ...
  geowin_esri_read_int_big(I, nr);
    
  if (run) { if (nr<prev) return false; }
  
  // get content length ...
  geowin_esri_read_int_big(I, clength); clength=clength*2;
  prev = nr;  
  return true;
}


// help functions for esri read

void geowin_esri_read_bbox(ifstream& I, double& xmin, double& ymin, double& xmax, double& ymax)
{
 double *ptr;
 
#if defined(__SUNPRO_CC) || defined(__KCC) || defined(_MSC_VER)
 char HF[8], ch;
#else
 unsigned char HF[8], ch;
#endif 
 
 int i;
 
#if defined(LITTLE_ENDIAN_MACHINE)
 //cout << "little endian!\n";
 for(i=0;i<8;i++) { I.get(ch); HF[i] = ch; }
 ptr = (double*) HF;  
 xmin = *ptr; 
 for(i=0;i<8;i++) { I.get(ch); HF[i] = ch; } 
 ymin = *ptr; 
 for(i=0;i<8;i++) { I.get(ch); HF[i] = ch; } 
 xmax = *ptr;   
 for(i=0;i<8;i++) { I.get(ch); HF[i] = ch; }
 ymax = *ptr; 
#else 
 //cout << "big endian!\n";
 for(i=0;i<8;i++) { I.get(ch); HF[7-i] = ch; }
 ptr = (double*) HF;  
 xmin = *ptr; 
 for(i=0;i<8;i++) { I.get(ch); HF[7-i] = ch; } 
 ymin = *ptr; 
 for(i=0;i<8;i++) { I.get(ch); HF[7-i] = ch; } 
 xmax = *ptr;   
 for(i=0;i<8;i++) { I.get(ch); HF[7-i] = ch; } 
 ymax = *ptr; 
#endif  
}

// read integer (little endian format ...)
void geowin_esri_read_int(ifstream& I, unsigned int& number)
{
 int ih1 = I.get(); 
 int ih2 = I.get(); 
 int ih3 = I.get();
 int ih4 = I.get();
 number = ih1 + 256*ih2 + 65536*ih3 + 16777216*ih4; 
}

// read integer (big endian format ...)
void geowin_esri_read_int_big(ifstream& I, unsigned int& number)
{
 int ih1 = I.get(); 
 int ih2 = I.get(); 
 int ih3 = I.get();
 int ih4 = I.get();
 number = ih4 + 256*ih3 + 65536*ih2 + 16777216*ih1; 
}


void geowin_esri_read_double(ifstream& I, double& number)
{
 double *ptr;
 
#if defined(__SUNPRO_CC) || defined(__KCC) || defined(_MSC_VER)
 char HF[8], ch;
#else 
 unsigned char HF[8], ch;
#endif 
 
 int i;

#if defined(LITTLE_ENDIAN_MACHINE) 
 for(i=0;i<8;i++) { I.get(ch); HF[i] = ch; }
 ptr = (double*) HF;  number = *ptr; 
#else
 for(i=0;i<8;i++) { I.get(ch); HF[7-i] = ch; }
 ptr = (double*) HF;  number = *ptr; 
#endif 
}

void geowin_esri_read_point(ifstream& I, point& p)
{
#if defined(__SUNPRO_CC) || defined(__KCC) || defined(_MSC_VER)
 char HF[8], ch;
#else 
 unsigned char HF[8], ch;
#endif 
 
 int i;

#if defined(LITTLE_ENDIAN_MACHINE)
 int nb;
 for(i=0;i<8;i++) { nb = I.get();  HF[i] = nb; }
 double* ptr = (double*) HF;  
 double x = *ptr;
 for(i=0;i<8;i++) { nb = I.get();  HF[i] = nb; }
 ptr = (double*) HF;   
 double y = *ptr;
#else 
 for(i=0;i<8;i++) { I.get(ch); HF[7-i] = ch; }
 double* ptr = (double*) HF;  
 double x = *ptr;
 for(i=0;i<8;i++) { I.get(ch); HF[7-i] = ch; }
 ptr = (double*) HF;   
 double y = *ptr;
#endif 

 p = point(x,y);
}

bool geowin_read_esri_point(unsigned int length, ifstream& I, point& pact)
{
 if (I.eof()) return false;
 geowin_esri_read_point(I,pact); 
 return true;
}


bool geowin_read_esri_polygon(unsigned int length, ifstream& I, list<polygon>& poly_act)
{
 double xmin, ymin, xmax, ymax;
 // read bounding box
 geowin_esri_read_bbox(I, xmin, ymin, xmax, ymax);
 //cout << xmin << "  " << ymin << "  " << xmax << "  " << ymax << "\n";
 
 if (I.eof()) return false;
 
 unsigned int num_parts, num_points, act;
 geowin_esri_read_int(I, num_parts); 
 geowin_esri_read_int(I, num_points);
 
 //cout << "#parts:" << num_parts << "  #points:" << num_points << "\n";
 
 array<point> poly_points(num_points);
 array<unsigned int> parts(num_parts);
 
 point p;
 unsigned int z, i;
 
 // read parts information ...
 for(z=0;z<num_parts;z++){
   geowin_esri_read_int(I,act); parts[z]=act;
 }
 
 // read polygon information
 for(z=0;z<num_points;z++){
   geowin_esri_read_point(I,p); poly_points[z]=p;
 }
 
 // build polygon(s) ...
 list<polygon> rings;
 list<point>   act_poly;
 point pact;
 //unsigned int anz =0;
 unsigned int w1,w2;
 
 for(z=0;z<num_parts;z++){
   w1 = parts[z];
   if (z==num_parts-1) w2 = num_points; else w2 = parts[z+1];
   act_poly.clear();
   for(i=w1;i<w2;i++) act_poly.push(poly_points[i]);
   if (! act_poly.empty()) act_poly.pop();
   rings.append(polygon(act_poly));
 }

 poly_act.conc(rings);
 return true;
}

bool geowin_read_esri_poly_line(unsigned int length, ifstream& I, list<segment>& sact)
{
 double xmin, ymin, xmax, ymax;
 geowin_esri_read_bbox(I, xmin, ymin, xmax, ymax);
 
 if (I.eof()) return false;
 
 unsigned int num_parts, num_points, act;
 geowin_esri_read_int(I, num_parts); geowin_esri_read_int(I, num_points);
 
 array<point>        points(num_points);
 array<unsigned int> parts(num_parts); 
 unsigned int z,i;
 point p;

 // read parts information ...
 for(z=0;z<num_parts;z++){ geowin_esri_read_int(I,act); parts[z]=act; }
 
 // read poly line information
 for(z=0;z<num_points;z++){ geowin_esri_read_point(I,p); points[z]=p; } 
 
  // build poly line ...
 unsigned int w1,w2;
 
 for(z=0;z<num_parts;z++){
   w1 = parts[z];
   if (z==num_parts-1) w2 = num_points; else w2 = parts[z+1];
   for(i=w1;i<w2;i++) if (i != w1) sact.append(segment(points[i-1],points[i]));
 }
 
 return true;
}


template<class TF, class TR>
class geowin_import_esri_objects : public geowin_import {
  virtual void operator()(geo_scene sc, string filename) { 
 
#if defined(__win32__)
     ifstream is(filename, ios::binary);
#else
     ifstream is(filename);
#endif
    
     if (! is) { error_handler(0,"No such file or directory."); return; }  
  
     if (sc==NULL) return;
     GeoWin* gw = get_geowin(sc);
     
     TF test_obj;
     
     string rat_type_name, float_type_name;
     geowin_get_type_names(test_obj, rat_type_name, float_type_name);    
     
     string tn = sc->get_type_name();
     
     string shapetype;
     unsigned int length;
     //bool b = 
     geowin_esri_read_header(is, shapetype, length);
     bool run = false;
     unsigned int number; // record number     
     
     //polygons
     if (shapetype == "Polygon"){
      if ((float_type_name == "SimplePolygons") || (float_type_name == "GeneralizedPolygons")) {
      
       list<polygon> LP;
       list<gen_polygon> GP;
       
       polygon pol;
       gen_polygon gpol;       
       
       while (true) {
         unsigned int clength;
         bool act = geowin_read_esri_record(is, clength, number, run);
	     run = true;
	 
	     if (! act) break;
	     else { // read in a polygon	  
	      //get shape type again ...
	      unsigned int sht; 
              geowin_esri_read_int(is,sht);
	  
	      // add here a check for sht
	      list<polygon> Lact;
	      bool w = geowin_read_esri_polygon(clength, is, Lact);
	      if (! w) break;
	      GP.append(gen_polygon(Lact));
	      LP.conc(Lact);
	     }
       }
       
       // set scene contents ...
       if (float_type_name == "SimplePolygons") {
        if (tn == rat_type_name) { 
         list<rat_polygon> LRP;
	     forall(pol,LP) LRP.append(rat_polygon(pol));
             geowin_set_objects(*gw, sc, LRP);
        }
        else  geowin_set_objects(*gw, sc, LP);
       }
       else { // GeneralizedPolygons
        if (tn == rat_type_name) { 
         list<rat_gen_polygon> LRP;
	 forall(gpol,GP) LRP.append(rat_gen_polygon(gpol));
         geowin_set_objects(*gw, sc, LRP);
        }
        else  geowin_set_objects(*gw, sc, GP);       
       }   
      } 	 
     }
     
     //points
     if (shapetype == "Point"){
      if (float_type_name == "Points") {
        
        list<point> LP;
	
        while (true) {
         unsigned int clength;
         bool act = geowin_read_esri_record(is, clength, number, run);
	     run = true;
	     if (! act) break;
	     else { // read in a point	  
	      unsigned int sht; 
              geowin_esri_read_int(is,sht);
	  
	      point pact;
	      bool w = geowin_read_esri_point(clength, is, pact);
	      if (! w) break;
	      LP.append(pact);
	     }
        }  	

        if (tn == rat_type_name) { 
         list<rat_point> LRP;
	     point pt;
	     forall(pt,LP) LRP.append(rat_point(pt));
             geowin_set_objects(*gw, sc, LRP);
        }
        else  geowin_set_objects(*gw, sc, LP);        
      }
     }
     
     //poly lines
     if (shapetype == "PolyLine"){
      if (float_type_name == "Segments") {      

        list<segment> SG;
	
        while (true) {
         unsigned int clength;
         bool act = geowin_read_esri_record(is, clength, number, run);
	     run = true;
	     if (! act) break;
	     else { // read in segments  
	      unsigned int sht; 
              geowin_esri_read_int(is,sht);
	  
	      list<segment> sacts;
	      bool w = geowin_read_esri_poly_line(clength, is, sacts);
	      if (! w) break;
	      SG.conc(sacts);
	     }
        }  	

        if (tn == rat_type_name) { 
         list<rat_segment> SGR;
	     segment st;
	     forall(st,SG) SGR.append(rat_segment(st));
             geowin_set_objects(*gw, sc, SGR);
        }
        else  geowin_set_objects(*gw, sc, SG);        
           
      }
     }
  }
};

// -------------------------------------------------------------------------------------------
// end ESRI - shapefile import
// -------------------------------------------------------------------------------------------


// -----------------------------   some additional input objects ----------------------------

// construct circles by providing 3 points ...

template<class T>
class geowin_input_circle_by_points : public GeoInputObject<T> {
public:
 geowin_input_circle_by_points() {}

 void operator()(GeoWin& gw, list<T>& LC)
 {
  gw.message("input 3 points constructing a circle");
  window& w = gw.get_window();
  __typename T::point_type p1, p2, p3;
  
  w >> p1; w << p1;
  w >> p2; w << p2;
  w >> p3; w << p3;
  
  // now lets construct the circle    
  T crc(p1,p2,p3);
  LC.append(crc);
  gw.message("");
 }

};

// construct a circle providing the diameter ...

template<class T>
class geowin_input_circle_by_diameter : public GeoInputObject<T> {
public:
 geowin_input_circle_by_diameter() {}
 
 void operator()(GeoWin& gw, list<T>& LC)
 {
  gw.message("input the diameter of a circle");
  window& w = gw.get_window();
  
  segment seg;
  w >> seg;
  point mp = midpoint(seg.start(), seg.end());
  
  __typename T::point_type p1(mp), p2(seg.end());
  
  // now lets construct the circle    
  T crc(p1,p2);
  LC.append(crc);  
  gw.message("");
 }
  
};


// -------------------------------------------------------------------------------------------


#if LEDA_ROOT_INCL_ID == 450096
#undef LEDA_ROOT_INCL_ID
#include <LEDA/UNDEFINE_NAMES.h>
#endif


#endif











