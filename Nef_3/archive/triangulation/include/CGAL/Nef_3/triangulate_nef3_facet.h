#ifndef TRIAGULATE_NEF3_FACET_H
#define TRIAGULATE_NEF3_FACET_H

#include <CGAL/Nef_3/pm_from_nef3_facet.h>
#include <CGAL/partition_y_monotone_2.h>
#include <CGAL/triangulate_monotone_polygon_2.h>

#ifdef _DEBUG_WINDOW
#include <CGAL/IO/Pm_Window_stream.h>
extern CGAL::Window_stream W;
#endif

#undef _DEBUG
#define _DEBUG 11
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

typedef enum { XY, NXY, ZX, NZX, YZ, NYZ, UNINIT } Plane_label;

template <typename Point_3, typename Point_2>
class Project_point_3_to_point_2
{
 public:
  typedef Point_3 argument_type;
  typedef Point_2 result_type;
  Project_point_3_to_point_2() : label(UNINT) {}
  Project_point_3_to_point_2( const Project_point_3_to_point_2& projector) :
    label( projector.label) {}
  Project_point_3_to_point_2( Plane_label l) : label(l) {}
  Point_2 operator()( const Point_3& p3) {
    CGAL_assertion( label != UNINIT);
    Point_2 p2;
    if( label == YZ)
      p2 = Point_2( p3.hy(), p3.hz(), p3.hw());
    else if( label == NYZ) // -YZ
      p2 = Point_2( p3.hz(), p3.hy(), p3.hw());
    else if( label == ZX) // ZX
      p2 = Point_2( p3.hz(), p3.hx(), p3.hw());
    else if( label == NZX) // -ZX
      p2 = Point_2( p3.hx(), p3.hz(), p3.hw());
    else if( label == XY) // XY
      p2 = Point_2( p3.hx(), p3.hy(), p3.hw());
    else if( label == NXY) // -XY
      p2 = Point_2( p3.hy(), p3.hx(), p3.hw());
    CGAL_NEF_TRACEN("projecting "<<p3<<" on "<<p2);
    return p2;
  }
  template <typename Plane_3>
  Point_3 operator()( const Point_2& p2, const Plane_3& plane) {
    CGAL_assertion( label != UNINIT);
    Point_3 p3;
    if( label == YZ)
      p3 = Point_3( 0, p2.hx(), p2.hy(), p2.hw());
    else if( label == NYZ) // -YZ
      p3 = Point_3( 0, p2.hy(), p2.hx(), p2.hw());
    else if( label == ZX) // ZX
      p3 = Point_3( p2.hy(), 0, p2.hx(), p2.hw());
    else if( label == NZX) // -ZX
      p3 = Point_3( p2.hx(), 0, p2.hy(), p2.hw());
    else if( label == XY) // XY
      p3 = Point_3( p2.hx(), p2.hy(), 0, p2.hw());
    else if( label == NXY) // -XY
      p3 = Point_3( p2.hy(), p2.hx(), 0, p2.hw());
    CGAL_NEF_TRACEN("(inv)projecting "<<p2<<" on "<<plane.projection(p3));
    return plane.projection(p3);
  }
 private:
  Plane_label label;
};

template <typename Diagonal_iterator, typename Planar_map>
void divide_pm_by_diagonals( Diagonal_iterator begin,
                             Diagonal_iterator beyond,
                             Planar_map& pm) {
  typedef typename Planar_map::X_curve X_curve;
  typedef typename Planar_map::Vertex_handle Vertex_handle;

  for( Diagonal_iterator i = begin; i != beyond; ++i) {
    Vertex_handle v1(i->first.current_circulator()->source());
    Vertex_handle v2(i->second.current_circulator()->source());
    X_curve segment( v1->point(), v2->point());
    CGAL_NEF_TRACEN( "Diagonal { " << v1->point() << ", " << v2->point() << " }");
    pm.insert_at_vertices( segment, v1, v2);
  }
}

template <typename SNC_structure,
          typename Halffacet_handle,
          typename OutputTriangleIterator,
          typename Traits>
void triangulate_nef3_facet( Halffacet_handle facet,
                             OutputTriangleIterator triangles,
                             const Traits& traits) {

  typedef typename Traits::Point_2 Point_2;

  typedef typename Traits::Planar_map_2 Planar_map_2;
  typedef typename Traits::Partition_traits_2 Partition_traits_2;
  typedef typename Traits::Monotone_polygon_triangulation_traits_2
    Monotone_polygon_triangulation_traits_2;

  typedef typename Traits::Circulator_project Circulator_project;
  typedef std::list<Circulator_project> Cycles_list;
  typedef typename Cycles_list::const_iterator Cycle_iterator;

  typedef typename Traits::Diagonal Diagonal;
  typedef std::list<Diagonal> Diagonal_list;
  typedef typename Diagonal_list::const_iterator Diagonal_iterator;

#define USING(x) typedef typename Planar_map_2::x x;
  USING(Vertex_handle);
  USING(Face_handle);
  USING(Face_iterator);
  USING(Holes_iterator);
  USING(Ccb_halfedge_circulator);
#undef USING

  typedef typename SNC_structure::SNC_decorator SNC_decorator;
  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Vector_3 Vector_3;
  typedef typename SNC_structure::Plane_3 Plane_3;
  typedef typename SNC_structure::Triangle_3 Triangle_3;

  typedef Project_point_3_to_point_2< Point_3, Point_2> Projector;

  // create a planar map to support the face creation

  Planar_map_2 pm;
  Face_handle polygon;
  Plane_label label;
  SNC_decorator D;
  Plane_3 fplane(D.plane(facet));
  CGAL_assertion(!fplane.is_degenerate());
  Vector_3 pf(fplane.orthogonal_vector()), pxy(0,0,1), pyz(1,0,0), pzx(0,1,0);
  if( !CGAL_NTS is_zero(pf*pyz) )
     label = (CGAL_NTS is_positive(pf*pyz) ? YZ : NYZ);
  else if( !CGAL_NTS is_zero(pf*pzx))
    label = (CGAL_NTS is_positive(pf*pzx) ? ZX : NZX);
  else {
    CGAL_assertion( !CGAL_NTS is_zero(pf*pxy) );
    label = (CGAL_NTS is_positive(pf*pxy) ? XY : NXY);
  }
  Projector projector(label);
  polygon = pm_from_nef3_facet<SNC_structure>( facet, pm, projector);

  CGAL_NEF_TRACEN("triangulating facet on plane "<<fplane);
  //CGAL_NEF_TRACEN("the facet will be proyected on the plane "<<(axis<0?"-":"")<<(axis==1?"YZ":(axis==2?"XZ":"XY")));

#ifdef _DEBUG_WINDOW
  W.init(-1000,1000,-1000);
  W.clear();
  W << BLUE << pm;
  Point_2 pause;
  //W >> pause;
#endif

  // mark the faces that correspond to holes

  Unique_hash_map< Face_handle, bool> hole_mark(false);
  for( Holes_iterator hi = polygon->holes_begin();
       hi != polygon->holes_end();
       ++hi) {

    Ccb_halfedge_circulator c(*hi), cend(c);
    CGAL_For_all( c, cend)
      if( c->face() != c->twin()->face())
        break;
    if( c != cend)
      continue; // TO VERIFY: the hole does not bound an area

    Face_handle hole_face = (*hi)->twin()->face();
    hole_mark[hole_face] = true;
  }

  // get a list of the cycles bounding the input polygon on the planar map

  Cycles_list cycles;
  cycles.push_back(Circulator_project(polygon->outer_ccb()));
  for( Holes_iterator hi = polygon->holes_begin();
       hi != polygon->holes_end(); ++hi)
    cycles.push_back(Circulator_project((*hi)->ccb()));

  // divide the polygon into y-monotone pieces

  Diagonal_list diagonals;
  partition_y_monotone_2( cycles.begin(), cycles.end(),
                          std::back_inserter(diagonals),
                          Partition_traits_2());
  divide_pm_by_diagonals( diagonals.begin(), diagonals.end(), pm);

  // mark the faces corresponding to y-monoton pieces

  Unique_hash_map< Face_handle, bool> piece_mark(false);
  for( Face_iterator fi = pm.faces_begin(); fi != pm.faces_end(); ++fi) {
    if( fi->is_unbounded())
      continue;
    if( hole_mark[fi])
      continue;
    piece_mark[fi] = true;
  }

#ifdef _DEBUG_WINDOW
  W.clear();
  W << GREEN << pm;
  W >> pause;
#endif

  // triangulate each monotone part ...

  for( Face_iterator fi = pm.faces_begin(); fi != pm.faces_end(); ++fi) {
    if( !piece_mark[fi])
      continue;
    Face_handle piece = fi;
    CGAL_assertion( std::distance( piece->holes_begin(),
                                   piece->holes_end()) == 0); // no holes
    Diagonal_list diagonals;
    Circulator_project circulator(piece->outer_ccb());
    triangulate_monotone_polygon_2( circulator, std::back_inserter(diagonals),
                                    Monotone_polygon_triangulation_traits_2());
    divide_pm_by_diagonals( diagonals.begin(), diagonals.end(), pm);
  }

#ifdef _DEBUG_WINDOW
  W.clear();
  W << YELLOW << pm;
#endif

  // and extract triangles from the topological map

  for( Face_iterator fi = pm.faces_begin(); fi != pm.faces_end(); ++fi) {
    if( fi->is_unbounded())
      continue;
    if( hole_mark[fi])
      continue;
    Face_handle triangle = fi;
    CGAL_assertion( std::distance( triangle->holes_begin(),
                                   triangle->holes_end()) == 0); // no holes
    Ccb_halfedge_circulator c(triangle->outer_ccb());
    CGAL_assertion_code( Ccb_halfedge_circulator cend(c));
    Point_2 p2[3];
    int i = 0;
    p2[i++] = c->source()->point();
    c++;
    CGAL_warning( c != cend);
    p2[i++] = c->source()->point();
    c++;
    CGAL_warning( c != cend);
    p2[i++] = c->source()->point();
    CGAL_assertion( ++c == cend);
    Point_3 p3[3];
    for( int i = 0; i < 3; ++i)
      p3[i] = projector(p2[i], fplane);
    CGAL_NEF_TRACEN(Triangle_3( p3[0], p3[1], p3[2]));
    *triangles++ = Triangle_3( p3[0], p3[1], p3[2]);
  }

  return;
}

} //namespace CGAL

#endif // TRIAGULATE_NEF3_FACET_H
