// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//

// release       :
// release_date  :
//
// file          : include/CGAL/IO/Geomview_stream.h
// source        : web/geomview.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_GEOMVIEW_STREAM_H
#define CGAL_GEOMVIEW_STREAM_H

#include <CGAL/Cartesian.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/IO/Color.h>

#include <stdio.h>
#include <unistd.h>
#include <signal.h>
#include <string.h>
#include <strstream.h>
#include <iomanip.h>

#include <sys/types.h>
#include <sys/uio.h>


class CGAL_Geomview_stream {
public:
    CGAL_Geomview_stream(const CGAL_Bbox_3 &bbox = CGAL_Bbox_3(0,0,0, 1,1,1),
                         const char *machine = (char*)NULL,
                         const char *login = (char*)NULL);

    // kept for backward compatibility
    CGAL_Geomview_stream(const char *machine,
                         const char *login,
                         const CGAL_Bbox_3 &bbox = CGAL_Bbox_3(0,0,0, 1,1,1));

    ~CGAL_Geomview_stream();

    void clear();
    void look_recenter() const;

    void set_bg_color(const CGAL_Color &c);

    CGAL_Geomview_stream &operator<<(const CGAL_Color &c);
    CGAL_Color get_vertex_color() const;
    CGAL_Color get_edge_color() const;
    CGAL_Color get_face_color() const;

    CGAL_Color set_vertex_color(const CGAL_Color&);
    CGAL_Color set_edge_color(const CGAL_Color&);
    CGAL_Color set_face_color(const CGAL_Color&);

    double vcr() const;
    double vcg() const;
    double vcb() const;

    double ecr() const;
    double ecg() const;
    double ecb() const;

    double fcr() const;
    double fcg() const;
    double fcb() const;

    double get_vertex_radius() const;
    double set_vertex_radius(double r);

    int get_line_width() const;
    int set_line_width(int w);

    CGAL_Geomview_stream &operator<<(const char *cptr);

    CGAL_Geomview_stream &operator<<(int i);
    CGAL_Geomview_stream &operator<<(double d);

    bool get_trace() const;
    bool set_trace(bool b);

    void trace(const char *cptr) const;
    void trace(double d) const;
    void trace(int i) const;

    void set_binary_mode();
    void set_ascii_mode();
    bool in_binary_mode() const;
    bool in_ascii_mode() const;

    CGAL_Geomview_stream &operator<< (
                          CGAL_Geomview_stream& (*fct)(CGAL_Geomview_stream&));

    CGAL_Geomview_stream &operator>>(char *expr);

    int bbox_count;
    int triangle_count;
    int segment_count;
    int point_count;
    int tetrahedron_count;

    char sexpr[1024];
private:
    void setup_geomview(const char *machine, const char *login);
    void frame(const CGAL_Bbox_3 &bbox);
    void pickplane(const CGAL_Bbox_3 &bbox);

    CGAL_Color col, vertex_color, edge_color, face_color;
    bool _trace;
    int in;       // file descriptor for input pipe
    int out;      // file descriptor for output pipe
    int pid;      // the geomview process identification
    int bflag ;   // bool that makes operator<< write binary format
    double _radius;
    int _line_width;
};


inline
CGAL_Geomview_stream&
binary(CGAL_Geomview_stream &os)
{
    os.set_binary_mode();
    return os;
}

inline
CGAL_Geomview_stream&
ascii(CGAL_Geomview_stream &os)
{
    os.set_ascii_mode();
    return os;
}


#ifdef CGAL_POINT_2_H
#ifndef CGAL_GV_OUT_POINT_2_H
#define CGAL_GV_OUT_POINT_2_H

template < class R >
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream &gv,
           const CGAL_Point_2<R> &p)
{
    ostrstream os;
    os << "p" << gv.point_count++ << ends ;
    char *id = os.str();

    gv << ascii
       << "(geometry " << id
       << " {appearance {linewidth 5 material {edgecolor "
       << gv.vcr() << gv.vcg() << gv.vcb()
       << "}}{SKEL 1 1 " ;

    gv << CGAL_to_double(p.x())
       << CGAL_to_double(p.y())
       << 0.0
       << "1 0\n"
       << "}})" ;

    return gv;
}
#endif // CGAL_GV_OUT_POINT_2_H
#endif // CGAL_POINT_2_H
#ifdef CGAL_POINT_3_H
#ifndef CGAL_GV_OUT_POINT_3_H
#define CGAL_GV_OUT_POINT_3_H

template < class R >
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream &gv,
           const CGAL_Point_3<R> &p)
{
    ostrstream os;
    os << "p" << gv.point_count++ << ends ;
    char *id = os.str();

    gv << ascii
       << "(geometry " << id
       << " {appearance {linewidth 5 material {edgecolor "
       << gv.vcr() << gv.vcg() << gv.vcb()
       << "}}{SKEL 1 1 "

       << CGAL_to_double(p.x())
       << CGAL_to_double(p.y())
       << CGAL_to_double(p.z())
       << "1 0\n"
       << "}})" ;

    return gv;
}
#endif // CGAL_GV_OUT_POINT_3_H
#endif // CGAL_POINT_3_H
#ifdef CGAL_SEGMENT_2_H
#ifndef CGAL_GV_OUT_SEGMENT_2_H
#define CGAL_GV_OUT_SEGMENT_2_H
template < class R >
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream &gv,
           const CGAL_Segment_2<R> &segment)
{
    ostrstream os;
    os << "seg" << gv.segment_count++ << ends ;
    char *id = os.str();

    gv << ascii
       << "(geometry " << id << " {appearance {linewidth "
       << gv.get_line_width() << "}{VECT "
       << 1 <<  2 << 1    // 1 polyline, two vertices, 1 color
       << 2               // the first polyline contains 2 vertices
       << 1               // and it has 1 color

    // here are start and end points
       << CGAL_to_double(segment.source().x())
       << CGAL_to_double(segment.source().y())
       << 0.0
       << CGAL_to_double(segment.target().x())
       << CGAL_to_double(segment.target().y())
       << 0.0

    // and the color of the segment and its opaqueness
        << gv.ecr()  << gv.ecg() << gv.ecb()  << 1.0

    // close the text bracket
       << "}})" ;

    return gv;
}
#endif // CGAL_GV_OUT_SEGMENT_2_H
#endif //CGAL_SEGMENT_2_H
#ifdef CGAL_SEGMENT_3_H
#ifndef CGAL_GV_OUT_SEGMENT_3_H
#define CGAL_GV_OUT_SEGMENT_3_H
template < class R >
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream &gv,
           const CGAL_Segment_3<R> &segment)
{
    ostrstream os;
    os << "seg" << gv.segment_count++ << ends ;
    char *id = os.str();

    gv << ascii
       << "(geometry " << id
       << " {appearance {linewidth "
       << gv.get_line_width()
       << " }{VECT "
       << 1 << 2 << 1   // 1 polyline, two vertices, 1 color
       << 2             // the first polyline contains 2 vertices
       << 1             // and it has 1 color

    // here are start and end points
       << CGAL_to_double(segment.source().x())
       << CGAL_to_double(segment.source().y())
       << CGAL_to_double(segment.source().z())
       << CGAL_to_double(segment.target().x())
       << CGAL_to_double(segment.target().y())
       << CGAL_to_double(segment.target().z())

    // and the color of the segment and its opaqueness
        << gv.ecr()  << gv.ecg()  << gv.ecb()  << 1.0

    // close the text bracket
       << "}})";

    return gv;
}
#endif // CGAL_GV_OUT_SEGMENT_3_H
#endif // CGAL_SEGMENT_3_H
#ifdef CGAL_TRIANGLE_2_H
#ifndef CGAL_GV_OUT_TRIANGLE_2_H
#define CGAL_GV_OUT_TRIANGLE_2_H

template < class R >
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream &gv,
           const CGAL_Triangle_2<R> &t)
{
    ostrstream os;
    os << "tr" << gv.triangle_count++ << ends;
    char *id = os.str();

    gv << ascii
       << "(geometry " << id << " {appearance {+edge material {edgecolor "
       << gv.ecr()  << gv.ecg()  << gv.ecb() <<  " } shading constant}{ "
       << binary
    // it's a planar polygon
       << "OFF BINARY\n"

    // it has 3 vertices, 1 face and 3 edges
       << 3 << 1 << 3;

    for(int i=0; i<3; i++){
        gv << CGAL_to_double(t[i].x())
           << CGAL_to_double(t[i].y())
           << 0.0;
    }

    // the face
    gv << 3 << 0 << 1 << 2 << 4 << gv.fcr() << gv.fcg() << gv.fcb() << 1.0
       << "}})"
       << ascii;

    return gv;
}
#endif // CGAL_GV_OUT_TRIANGLE_2_H
#endif // CGAL_TRIANGLE_2_H
#ifdef CGAL_TRIANGLE_3_H
#ifndef CGAL_GV_OUT_TRIANGLE_3_H
#define CGAL_GV_OUT_TRIANGLE_3_H

template < class R >
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream &gv,
           const CGAL_Triangle_3<R> &t)
{
    ostrstream os;
    os << "tr" << gv.triangle_count++ << ends;
    char *id = os.str();

    gv << ascii
       << "(geometry " << id << " {appearance {+edge material {edgecolor "
       << gv.ecr()  << gv.ecg()  << gv.ecb() <<  "} shading constant}{ "
       << binary
    // it's a planar polygon
       << "OFF BINARY\n"

    // it has 3 vertices, 1 face and 3 edges
       << 3 << 1 << 3;

    for(int i=0; i<3; i++){
        gv << CGAL_to_double(t[i].x())
           << CGAL_to_double(t[i].y())
           << CGAL_to_double(t[i].z());
    }

    // the face
    gv << 3 << 0 << 1 << 2 << 4 << gv.fcr() << gv.fcg() << gv.fcb() << 1.0
       << "}})"
       << ascii;

    return gv;
}
#endif // CGAL_GV_OUT_TRIANGLE_3_H
#endif // CGAL_TRIANGLE_3_H
#ifdef CGAL_TETRAHEDRON_3_H
#ifndef CGAL_GV_OUT_TETRAHEDRON_3_H
#define CGAL_GV_OUT_TETRAHEDRON_3_H

template < class R >
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream &gv,
           const CGAL_Tetrahedron_3<R> &t)
{
    ostrstream os;
    os << "tetra" << gv.tetrahedron_count++ << ends ;
    char *id = os.str();

    gv << ascii
       << "(geometry " << id  << "  {appearance {}{ "
       << binary
       << "OFF BINARY\n"

    // it has 4 vertices, 4 face and 6 edges
       << 4 << 4 << 6 ;

    // the vertices
    for(int i=0; i<4; i++){
        gv << CGAL_to_double(t[i].x())
           << CGAL_to_double(t[i].y())
           << CGAL_to_double(t[i].z());
    }

    // the faces
    double r = gv.fcr(),
           g = gv.fcg(),
           b = gv.fcb();
    gv << 3 << 0 << 1 << 2 << 4 << r << g << b << 1.0
       << 3 << 3 << 0 << 1 << 4 << r << g << b << 1.0
       << 3 << 3 << 1 << 2 << 4 << r << g << b << 1.0
       << 3 << 3 << 0 << 2 << 4 << r << g << b << 1.0
       << "}})" << ascii;
    return gv;
}
#endif // CGAL_GV_OUT_TETRAHEDRON_3_H
#endif  // CGAL_TETRAHEDRON_3_H
#ifdef CGAL_BBOX_2_H
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream &gv,
           const CGAL_Bbox_2 &bbox);
#endif // CGAL_BBOX_2_H
#ifdef CGAL_BBOX_3_H
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream &gv,
           const CGAL_Bbox_3 &bbox);
#endif // CGAL_BBOX_3_H
#ifdef CGAL_TETRAHEDRALIZATION_3_H
#ifndef CGAL_GV_OUT_CGAL_TETRAHEDRALIZATION_3_H
#define CGAL_GV_OUT_CGAL_TETRAHEDRALIZATION_3_H

template < class Tr >
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream& os, CGAL_Tetrahedralization_3<Tr>& T)
{
    CGAL_Tetrahedralization_3<Tr>::Simplex_iterator
        it = T.simplices_begin(),
        end = T.simplices_end();

    while(it != end){
        if(! T.is_infinite(*it)){
            os << *it;
        }
        ++it;
  }
  return os;
}
#endif // CGAL_GV_OUT_CGAL_TETRAHEDRALIZATION_3_H
#endif // CGAL_TETRAHEDRALIZATION_3_H
#ifdef CGAL_TETRAHEDRALIZATION_SIMPLEX_H
#ifndef CGAL_GV_OUT_CGAL_TETRAHEDRALIZATION_SIMPLEX_H
#define CGAL_GV_OUT_CGAL_TETRAHEDRALIZATION_SIMPLEX_H

template < class V >
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream &gv,
           const CGAL_Tetrahedralization_simplex<V>* s)
{
    ostrstream os;
    os << "Simplex_" << (unsigned long int)s  << ends ;
    char *id = os.str();

    gv << ascii
       << "(geometry " << id  << "  {appearance {+edge material {edgecolor "
       << gv.ecr()  << gv.ecg()  << gv.ecb() <<  "} shading constant}{ "
       << binary
       << "OFF BINARY\n";
    // it has 4 vertices, 4 face and 6 edges
    gv << 4 << 4 << 6 ;

    // the vertices
    for(int i=0; i<4; i++){
        gv << CGAL_to_double(s->vertex(i)->point().x())
           << CGAL_to_double(s->vertex(i)->point().y())
           << CGAL_to_double(s->vertex(i)->point().z()) ;
    }

    // the faces
    double r = gv.fcr(),
           g = gv.fcg(),
           b = gv.fcb();

    gv << 3 << 3 << 1 << 2 << 4 << r << g << b << 1.0
       << 3 << 3 << 0 << 2 << 4 << r << g << b << 1.0
       << 3 << 3 << 0 << 1 << 4 << r << g << b << 1.0
       << 3 << 0 << 1 << 2 << 4 << r << g << b << 1.0
       << "}})"
       << ascii;

    return gv;
}
#endif // CGAL_GV_OUT_CGAL_TETRAHEDRALIZATION_SIMPLEX_H
#endif  // CGAL_TETRAHEDRALIZATION_SIMPLEX_H
#ifdef CGAL_TETRAHEDRALIZATION_VERTEX_H
#ifndef CGAL_GV_OUT_CGAL_TETRAHEDRALIZATION_VERTEX_H
#define CGAL_GV_OUT_CGAL_TETRAHEDRALIZATION_VERTEX_H

template < class P >
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream &gv,
           const CGAL_Tetrahedralization_vertex<P>* v)
{
    ostrstream os;
    os << "Vertex_" << (unsigned long int)v  << ends ;
    char *id = os.str();
    double x = CGAL_to_double(v->point().x());
    double y = CGAL_to_double(v->point().y());
    double z = CGAL_to_double(v->point().z());
    double radius = gv.get_vertex_radius();
    gv << ascii
       << "(geometry " << id  << "  {appearance {}{ "
       << binary
       << "OFF BINARY\n"

    // it has 4 vertices, 6 face and 12 edges
       << 8 << 6 << 12;

    // the vertices
    gv << x-radius << y-radius << z-radius
       << x+radius << y-radius << z-radius
       << x+radius << y+radius << z-radius
       << x-radius << y+radius << z-radius
       << x-radius << y+radius << z+radius
       << x+radius << y+radius << z+radius
       << x+radius << y-radius << z+radius
       << x-radius << y-radius << z+radius;

    // the faces
    double r = gv.vcr(),
           g = gv.vcg(),
           b = gv.vcb();
    gv << 4 << 0 << 1 << 6 << 7 << 4 << r << g << b << 1.0
       << 4 << 1 << 2 << 5 << 6 << 4 << r << g << b << 1.0
       << 4 << 2 << 3 << 4 << 5 << 4 << r << g << b << 1.0
       << 4 << 3 << 0 << 7 << 4 << 4 << r << g << b << 1.0
       << 4 << 0 << 1 << 2 << 3 << 4 << r << g << b << 1.0
       << 4 << 4 << 5 << 6 << 7 << 4 << r << g << b << 1.0
       << "}})" << ascii;

    return gv;
}
#endif // CGAL_GV_OUT_CGAL_TETRAHEDRALIZATION_VERTEX_H
#endif  // CGAL_TETRAHEDRALIZATION_VERTEX_H

#ifdef CGAL_DELAUNAY_TETRAHEDRALIZATION_3_H
#ifndef CGAL_GV_OUT_CGAL_DELAUNAY_TETRAHEDRALIZATION_3_H
#define CGAL_GV_OUT_CGAL_DELAUNAY_TETRAHEDRALIZATION_3_H
template < class I >
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream& gv,
           const CGAL_Delaunay_tetrahedralization_3<I> &DT)
{
  // return gv << (const CGAL_Tetrahedralization_3<I>&)DT;
  CGAL_Tetrahedralization_3<I>::Simplex_iterator
      it = DT.simplices_begin(),
      end = DT.simplices_end();

  while(it != end){
    if(! DT.is_infinite(*it)){
      gv << *it;
    }
    ++it;
  }
  return gv;
}
#endif // CGAL_GV_OUT_CGAL_DELAUNAY_TETRAHEDRALIZATION_3_H
#endif // CGAL_DELAUNAY_TETRAHEDRALIZATION_3_H


char*
CGAL_nth(char* s, int count);

bool
CGAL_is_prefix(const char* p, const char* w);


#ifdef CGAL_POINT_3_H
template < class R >
void
parse_point(char* pickpoint,
            CGAL_Point_3<R> &point)
{
    strstream ss;
    ss << pickpoint << ends ;

    double x, y, z, w;
    char parenthesis;
    ss >> parenthesis >> x >> y >> z >> w;
    point  = CGAL_Point_3<R>(x, y, z, w);
}
#endif // CGAL_POINT_3_H


#ifdef CGAL_POINT_3_H
#ifndef CGAL_GV_IN_POINT_3_H
#define CGAL_GV_IN_POINT_3_H

template < class R >
CGAL_Geomview_stream&
operator>>(CGAL_Geomview_stream &gv,
           CGAL_Point_3<R> &point)
{
    char gclpick[100];
    strcpy(gclpick, "(pick world pickplane * nil nil nil nil nil nil nil)");

    gv << ascii
       << "(pickable pickplane yes) (ui-target pickplane yes)"
       << "(interest " << gclpick <<")";

    gv >> gv.sexpr;  // this reads a gcl expression

    char* pickpoint = CGAL_nth(gv.sexpr, 3);
    // this gives something as: (0.0607123 0.0607125 4.76837e-07 0.529628)
    parse_point(pickpoint,point);

    // we echo the input
    gv << point;

    // we are done and tell geomview to stop sending pick events
    gv << "(uninterest " << gclpick  << ") (pickable pickplane no)" ;

    return gv ;
}
#endif // CGAL_GV_IN_POINT_3_H
#endif // CGAL_POINT_3_H
#ifdef CGAL_TETRAHEDRALIZATION_SIMPLEX_H
#ifndef CGAL_GV_IN_CGAL_TETRAHEDRALIZATION_SIMPLEX_H
#define CGAL_GV_IN_CGAL_TETRAHEDRALIZATION_SIMPLEX_H

#include <stdlib.h>

template < class V >
CGAL_Geomview_stream&
operator>>(CGAL_Geomview_stream &gv,
           CGAL_Tetrahedralization_simplex<V>*& s)
{
    char* id;
    char gclpick[100];
    strcpy(gclpick, "(pick world * nil nil nil nil nil nil nil nil)");

    gv << ascii ;
    gv << "(interest " << gclpick << ")" ;

    while(true) {
        gv >> gv.sexpr;  // this reads a gcl expression
        id = CGAL_nth(gv.sexpr, 2);

        if(! CGAL_is_prefix("Simplex_", id)){
            cerr << "You did not click on a simplex" << endl;
            continue;
        } else {
            break;
        }
    }
    gv << "(ui-target " << id << " yes)";
    gv << "(uninterest " << gclpick << ")";
    id+=8;                   // remove first 7 chars
    unsigned long int ui = strtoul(id, (char **)NULL, 10);
    s = (CGAL_Tetrahedralization_simplex<V>*)ui;

    return gv;
}
#endif // CGAL_GV_IN_CGAL_TETRAHEDRALIZATION_SIMPLEX_H
#endif  // CGAL_TETRAHEDRALIZATION_SIMPLEX_H
#ifdef CGAL_TETRAHEDRALIZATION_VERTEX_H
#ifndef CGAL_GV_IN_CGAL_TETRAHEDRALIZATION_VERTEX_H
#define CGAL_GV_IN_CGAL_TETRAHEDRALIZATION_VERTEX_H

#include <stdlib.h>

template < class P >
CGAL_Geomview_stream&
operator>>(CGAL_Geomview_stream &gv,
           CGAL_Tetrahedralization_vertex<P>*& v)
{
    char* id;
    char gclpick[100];
    strcpy(gclpick, "(pick world * nil nil nil nil nil nil nil nil)");

    gv << ascii
       << "(interest " << gclpick << ")" ;

    while(true) {
        gv >> gv.sexpr;  // this reads a gcl expression
        id = CGAL_nth(gv.sexpr, 2);

        if(! CGAL_is_prefix("Vertex_", id)){
            cerr << "You did not click on a vertex" << endl;
            continue;
        } else {
            break;
        }
    }
    gv << "(ui-target " << id << " yes)";
    gv << "(uninterest " << gclpick << ")";

    id+=7;                   // cut first 6 chars
    unsigned long int ui = strtoul(id, (char **)NULL, 10);
    v = (CGAL_Tetrahedralization_vertex<P>*)ui;

    return gv;
}
#endif // CGAL_GV_IN_CGAL_TETRAHEDRALIZATION_VERTEX_H
#endif  // CGAL_TETRAHEDRALIZATION_VERTEX_H


#ifdef CGAL_POLYHEDRON_3
#include <CGAL/IO/Polyhedron_geomview_ostream.h>
#endif

#endif // CGAL_GEOMVIEW_STREAM_H
