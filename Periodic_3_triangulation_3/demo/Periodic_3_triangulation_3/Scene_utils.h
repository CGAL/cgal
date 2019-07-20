#ifndef SCENE_UTILS_H
#define SCENE_UTILS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>

// Making available the Periodic_3_Delaunay_triangulation_3 to be
// drawn in the Scene.
typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<EPIC> K;

typedef CGAL::Periodic_3_Delaunay_triangulation_3<K> P3DT;

typedef K::RT NT;
typedef K::FT FT;
typedef K::Vector_3 Vector;

typedef P3DT::Cell_handle Cell_handle;
typedef P3DT::Vertex_handle Vertex_handle;
typedef P3DT::Vertex_iterator Vertex_iterator;
typedef P3DT::Facet Facet;
typedef P3DT::Offset Offset;

typedef P3DT::Point Point;
typedef P3DT::Segment Segment;
typedef P3DT::Triangle Triangle;
typedef P3DT::Tetrahedron Tetrahedron;

typedef P3DT::Periodic_point_iterator Point_iterator;
typedef P3DT::Periodic_segment_iterator Segment_iterator;
typedef P3DT::Periodic_triangle_iterator Triangle_iterator;
typedef P3DT::Periodic_tetrahedron_iterator Tetrahedron_iterator;

typedef P3DT::Iterator_type Iterator_type;

typedef CGAL::Random_points_in_cube_3<Point> RandPts;

// Segment comparison needed by the segment set
template <class Segment>
struct Compare_segment {
  bool operator()(Segment s1, Segment s2) const {
    if (s1.source() == s2.source()) return (s1.target() < s2.target());
    return (s1.source() < s2.source());
  }
};

// Projected triangle needed by the sorting of the triangles in z
// direction to draw them transparently. This is used by
// gl_draw_location and gl_draw_conflict
class Projected_triangle {
  double m_z;
  Triangle m_t;
public:
  Projected_triangle() {}
  Projected_triangle(double z, Triangle tr) : m_z(z), m_t(tr) {}
  Projected_triangle(const Projected_triangle &pt) : m_z(pt.z()), m_t(pt.t()) {}

  Triangle t() { return m_t; }
  const Triangle t() const { return m_t; }

  double& z() { return m_z; }
  double z() const { return m_z; }

  static bool closer(const Projected_triangle& t1, const Projected_triangle& t2) {
    return t1.z() < t2.z();
  }
};

// Segment clipper used to clip segments with the unit cube/square.
// Implements the Cohen-Sutherland line clipping algorithm
struct Segment_clipper {

  int cs_outcode(const Point& p) {
    int outcode = 0;
    outcode |= (p.x()<0.0 ?  1:0);
    outcode |= (p.x()>1.0 ?  2:0);
    outcode |= (p.y()<0.0 ?  4:0);
    outcode |= (p.y()>1.0 ?  8:0);
    outcode |= (p.z()<0.0 ? 16:0);
    outcode |= (p.z()>1.0 ? 32:0);
    return outcode;
  }

  bool operator()(Point& p1, Point& p2) {
      int p1_outcode = cs_outcode(p1);
      int p2_outcode = cs_outcode(p2);
      FT x(0),y(0),z(0), t(0);
      bool accept = false;
      bool done = false;
      do {
        if (p1_outcode == 0 && p2_outcode == 0) { accept = true; done = true; }
        else if ((p1_outcode & p2_outcode) != 0) { done = true; }
        else {
          if (p1_outcode == 0) {
            std::swap(p1_outcode, p2_outcode);
            std::swap(p1,p2);
          }
          if ((p1_outcode & 1) != 0) {
            t = (0-p1.x())/(p2.x()-p1.x());
            x = 0;
            y =p1.y()+t*(p2.y()-p1.y());
            z =p1.z()+t*(p2.z()-p1.z());
          } else {
            if ((p1_outcode & 2) != 0) {
              t = (1-p1.x())/(p2.x()-p1.x());
              x = 1;
              y =p1.y()+t*(p2.y()-p1.y());
              z =p1.z()+t*(p2.z()-p1.z());
            } else {
              if ((p1_outcode & 4) != 0) {
                t = (0-p1.y())/(p2.y()-p1.y());
                x =p1.x()+t*(p2.x()-p1.x());
                y = 0;
                z =p1.z()+t*(p2.z()-p1.z());
              } else {
                if ((p1_outcode & 8) != 0) {
                  t = (1-p1.y())/(p2.y()-p1.y());
                  x =p1.x()+t*(p2.x()-p1.x());
                  y = 1;
                  z =p1.z()+t*(p2.z()-p1.z());
                } else {
                  if ((p1_outcode & 16) != 0) {
                    t = (0-p1.z())/(p2.z()-p1.z());
                    x =p1.x()+t*(p2.x()-p1.x());
                    y =p1.y()+t*(p2.y()-p1.y());
                    z = 0;
                  } else {
                    if ((p1_outcode & 32) != 0) {
                      t = (1-p1.z())/(p2.z()-p1.z());
                      x =p1.x()+t*(p2.x()-p1.x());
                      y =p1.y()+t*(p2.y()-p1.y());
                      z = 1;
                    }
                  }
                }
              }
            }
          }
          p1 = Point(x,y,z);
          p1_outcode = cs_outcode(p1);
        }
      }
      while (!done);
    return accept;
  }
};

#endif // SCENE_UTILS_H
