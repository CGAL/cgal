// ============================================================================
//
// Copyright (c) 1997,1998,1999,2000 The CGAL Consortium
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
// package       : Geomview 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Sylvain Pion
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_GEOMVIEW_STREAM_H
#define CGAL_GEOMVIEW_STREAM_H

#include <CGAL/Cartesian.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/IO/Color.h>

#include <string>
#include <strstream> // deprecated
// #include <sstream>

CGAL_BEGIN_NAMESPACE

class Geomview_stream {
public:
    Geomview_stream(const Bbox_3 &bbox = Bbox_3(0,0,0, 1,1,1),
		    const char *machine = NULL,
		    const char *login = NULL);

    ~Geomview_stream();

    void clear();
    void look_recenter();

    void set_bg_color(const Color &c);

    Color get_vertex_color() const;
    Color get_edge_color() const;
    Color get_face_color() const;

    Color set_vertex_color(const Color&);
    Color set_edge_color(const Color&);
    Color set_face_color(const Color&);

    double vcr() const;
    double vcg() const;
    double vcb() const;

    double ecr() const;
    double ecg() const;
    double ecb() const;

    double fcr() const;
    double fcg() const;
    double fcb() const;

    double get_vertex_radius() const
    {
	return radius_;
    }
    double set_vertex_radius(double r)
    {
	double old = radius_;
	radius_ = r;
	return old;
    }

    int get_line_width() const
    {
	return line_width_;
    }
    int set_line_width(int w)
    {
        int old = line_width_;
        line_width_ = w;
        return old;
    }

    Geomview_stream &operator<<(const Color &c);
    Geomview_stream &operator<<(std::string s);
    Geomview_stream &operator<<(int i);
    Geomview_stream &operator<<(double d);

    bool set_trace(bool b)
    {
	bool old = trace_flag;
	trace_flag = b;
	return old;
    }
    bool get_trace() const
    {
	return trace_flag;
    }

    void trace(const std::string s) const
    {
        if (get_trace())
            std::cerr << s;
    }
    void trace(double d) const
    {
        if (get_trace())
            std::cerr << d << ' ';
    }
    void trace(int i) const
    {
        if (get_trace())
            std::cerr << i << ' ';
    }

    void set_binary_mode()         { binary_flag = true; }
    void set_ascii_mode()          { binary_flag = false; }
    bool in_binary_mode() const    { return binary_flag; }
    bool in_ascii_mode() const     { return ! binary_flag; }

    Geomview_stream &operator<<( Geomview_stream& (*fct)(Geomview_stream&));

    Geomview_stream &operator>>(char *expr);

    int bbox_count;
    int triangle_count;
    int sphere_count;
    int segment_count;
    int point_count;
    int tetrahedron_count;

private:
    void setup_geomview(const char *machine, const char *login);
    void frame(const Bbox_3 &bbox);
    void pickplane(const Bbox_3 &bbox);

    Color vertex_color, edge_color, face_color;
    bool trace_flag;  // bool that makes operator<<() write a trace on cerr.
    bool binary_flag; // bool that makes operator<< write binary format
    double radius_;   // radius of vertices
    int line_width_;  // width of edges
    int in, out;      // file descriptors for input and output pipes
    int pid;          // the geomview process identification
};

inline
Geomview_stream&
binary(Geomview_stream &gv)
{
    gv.set_binary_mode();
    return gv;
}

inline
Geomview_stream&
ascii(Geomview_stream &gv)
{
    gv.set_ascii_mode();
    return gv;
}

#if defined CGAL_POINT_2_H && \
   !defined CGAL_GV_OUT_POINT_2_H
#define CGAL_GV_OUT_POINT_2_H
template < class R >
Geomview_stream&
operator<<(Geomview_stream &gv, const Point_2<R> &p)
{
    gv << ascii
       << "(geometry P" << gv.point_count++
       << " {appearance {linewidth 5 material {edgecolor "
       << gv.vcr() << gv.vcg() << gv.vcb()
       << "}}{SKEL 1 1 " ;

    gv << CGAL::to_double(p.x())
       << CGAL::to_double(p.y())
       << 0.0
       << "1 0\n"
       << "}})" ;

    return gv;
}
#endif

#if defined CGAL_POINT_3_H && \
   !defined CGAL_GV_OUT_POINT_3_H
#define CGAL_GV_OUT_POINT_3_H
template < class R >
Geomview_stream&
operator<<(Geomview_stream &gv, const Point_3<R> &p)
{
    gv << ascii
       << "(geometry P" << gv.point_count++
       << " {appearance {linewidth 5 material {edgecolor "
       << gv.vcr() << gv.vcg() << gv.vcb()
       << "}}{SKEL 1 1 "

       << CGAL::to_double(p.x())
       << CGAL::to_double(p.y())
       << CGAL::to_double(p.z())
       << "1 0\n"
       << "}})" ;

    return gv;
}
#endif

#if defined CGAL_SEGMENT_2_H && \
   !defined CGAL_GV_OUT_SEGMENT_2_H
#define CGAL_GV_OUT_SEGMENT_2_H
template < class R >
Geomview_stream&
operator<<(Geomview_stream &gv, const Segment_2<R> &segment)
{
    gv << ascii
       << "(geometry Seg" << gv.segment_count++
       << " {appearance {linewidth "
       << gv.get_line_width() << "}{VECT "
       << 1 <<  2 << 1    // 1 polyline, two vertices, 1 color
       << 2               // the first polyline contains 2 vertices
       << 1               // and it has 1 color

    // here are start and end points
       << CGAL::to_double(segment.source().x())
       << CGAL::to_double(segment.source().y())
       << 0.0
       << CGAL::to_double(segment.target().x())
       << CGAL::to_double(segment.target().y())
       << 0.0

    // and the color of the segment and its opaqueness
        << gv.ecr()  << gv.ecg() << gv.ecb()  << 1.0

    // close the text bracket
       << "}})" ;

    return gv;
}
#endif

#if defined CGAL_SEGMENT_3_H && \
   !defined CGAL_GV_OUT_SEGMENT_3_H
#define CGAL_GV_OUT_SEGMENT_3_H
template < class R >
Geomview_stream&
operator<<(Geomview_stream &gv, const Segment_3<R> &segment)
{
    gv << ascii
       << "(geometry Seg" << gv.segment_count++
       << " {appearance {linewidth "
       << gv.get_line_width()
       << " }{VECT "
       << 1 << 2 << 1   // 1 polyline, two vertices, 1 color
       << 2             // the first polyline contains 2 vertices
       << 1             // and it has 1 color

    // here are start and end points
       << CGAL::to_double(segment.source().x())
       << CGAL::to_double(segment.source().y())
       << CGAL::to_double(segment.source().z())
       << CGAL::to_double(segment.target().x())
       << CGAL::to_double(segment.target().y())
       << CGAL::to_double(segment.target().z())

    // and the color of the segment and its opaqueness
        << gv.ecr()  << gv.ecg()  << gv.ecb()  << 1.0

    // close the text bracket
       << "}})";

    return gv;
}
#endif

#if defined CGAL_TRIANGLE_2_H && \
   !defined CGAL_GV_OUT_TRIANGLE_2_H
#define CGAL_GV_OUT_TRIANGLE_2_H
template < class R >
Geomview_stream&
operator<<(Geomview_stream &gv, const Triangle_2<R> &t)
{
    gv << ascii
       << "(geometry Tr" << gv.triangle_count++
       << " {appearance {+edge material {edgecolor "
       << gv.ecr()  << gv.ecg()  << gv.ecb() <<  " } shading constant}{ "
       << binary
    // it's a planar polygon
       << "OFF BINARY\n"

    // it has 3 vertices, 1 face and 3 edges
       << 3 << 1 << 3;

    for(int i=0; i<3; i++){
        gv << CGAL::to_double(t[i].x())
           << CGAL::to_double(t[i].y())
           << 0.0;
    }

    // the face
    gv << 3 << 0 << 1 << 2 << 4 << gv.fcr() << gv.fcg() << gv.fcb() << 1.0
       << "}})"
       << ascii;

    return gv;
}
#endif

#if defined CGAL_TRIANGLE_3_H && \
   !defined CGAL_GV_OUT_TRIANGLE_3_H
#define CGAL_GV_OUT_TRIANGLE_3_H
template < class R >
Geomview_stream&
operator<<(Geomview_stream &gv, const Triangle_3<R> &t)
{
    gv << ascii
       << "(geometry Tr" << gv.triangle_count++
       << " {appearance {+edge material {edgecolor "
       << gv.ecr()  << gv.ecg()  << gv.ecb() <<  "} shading constant}{ "
       << binary
    // it's a planar polygon
       << "OFF BINARY\n"

    // it has 3 vertices, 1 face and 3 edges
       << 3 << 1 << 3;

    for(int i=0; i<3; i++){
        gv << CGAL::to_double(t[i].x())
           << CGAL::to_double(t[i].y())
           << CGAL::to_double(t[i].z());
    }

    // the face
    gv << 3 << 0 << 1 << 2 << 4 << gv.fcr() << gv.fcg() << gv.fcb() << 1.0
       << "}})"
       << ascii;

    return gv;
}
#endif

#if defined CGAL_TETRAHEDRON_3_H && \
   !defined CGAL_GV_OUT_TETRAHEDRON_3_H
#define CGAL_GV_OUT_TETRAHEDRON_3_H
template < class R >
Geomview_stream&
operator<<(Geomview_stream &gv, const Tetrahedron_3<R> &t)
{
    gv << ascii
       << "(geometry Tetra" << gv.tetrahedron_count++
       << " {appearance {}{ "
       << binary
       << "OFF BINARY\n"

    // it has 4 vertices, 4 face and 6 edges
       << 4 << 4 << 6 ;

    // the vertices
    for(int i=0; i<4; i++){
        gv << CGAL::to_double(t[i].x())
           << CGAL::to_double(t[i].y())
           << CGAL::to_double(t[i].z());
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
#endif

#if defined CGAL_SPHERE_3_H && \
   !defined CGAL_GV_OUT_SPHERE_3_H
#define CGAL_GV_OUT_SPHERE_3_H
template < class R >
Geomview_stream&
operator<<(Geomview_stream &gv, const Sphere_3<R> &S)
{
    gv << ascii
       << "(geometry Sph" << gv.sphere_count++
       << " {appearance {+edge material {edgecolor "
       << gv.ecr()  << gv.ecg()  << gv.ecb() <<  "} shading constant}{ "
       << "SPHERE\n"
       << CGAL_CLIB_STD::sqrt(CGAL::to_double(S.squared_radius())) << "\n"
       << CGAL::to_double(S.center().x()) << " "
       << CGAL::to_double(S.center().y()) << " "
       << CGAL::to_double(S.center().z()) << "}})";

    return gv;
}
#endif

Geomview_stream&
operator<<(Geomview_stream &gv, const Bbox_2 &bbox);

Geomview_stream&
operator<<(Geomview_stream &gv, const Bbox_3 &bbox);

char*
nth(char* s, int count);

#ifdef CGAL_POINT_3_H
template < class R >
void
parse_point(const char* pickpoint, Point_3<R> &point)
{
    // std::stringstream ss;
    std::strstream ss;
    ss << pickpoint << std::ends;

    double x, y, z, w;
    char parenthesis;
    ss >> parenthesis >> x >> y >> z >> w;
    point  = Point_3<R>(x, y, z, w);
}
#endif

#if defined CGAL_POINT_3_H && !defined CGAL_GV_IN_POINT_3_H
#define CGAL_GV_IN_POINT_3_H
template < class R >
Geomview_stream&
operator>>(Geomview_stream &gv, Point_3<R> &point)
{
    const char *gclpick =
	"(pick world pickplane * nil nil nil nil nil nil nil)";

    gv << ascii
       << "(pickable pickplane yes) (ui-target pickplane yes)"
       << "(interest " << gclpick <<")";

    char sexpr[1024];
    gv >> sexpr;  // this reads a gcl expression

    const char* pickpoint = nth(sexpr, 3);
    // this gives something as: (0.0607123 0.0607125 4.76837e-07 0.529628)
    parse_point(pickpoint, point);

    // we echo the input
    gv << point; // FIXME : have an echo_mode boolean ?

    // we are done and tell geomview to stop sending pick events
    gv << "(uninterest " << gclpick  << ") (pickable pickplane no)" ;

    return gv ;
}
#endif

CGAL_END_NAMESPACE

#endif // CGAL_GEOMVIEW_STREAM_H
