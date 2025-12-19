#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>

#include <memory>

#include <CGAL/draw_polygon_with_holes_2.h>

#include <vector>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                    Point ;
typedef CGAL::Polygon_2<K>            Polygon_2 ;
typedef CGAL::Polygon_with_holes_2<K> PolygonWithHoles ;

typedef std::shared_ptr<PolygonWithHoles> PolygonWithHolesPtr ;
typedef std::shared_ptr<Polygon_2> PolygonPtr ;

typedef std::vector<PolygonWithHolesPtr> PolygonWithHolesPtrVector;
typedef std::vector<PolygonPtr> PolygonPtrVector;

PolygonWithHolesPtrVector
exterior_offset_of_disjoint_polygons_with_holes(double lOffset, const std::vector<PolygonWithHoles>& pwhs)
{
  std::vector<Point> outer_vertices;
  for (const PolygonWithHoles& pwh : pwhs)
    outer_vertices.insert(outer_vertices.end(),
                          pwh.outer_boundary().container().begin(),
                          pwh.outer_boundary().container().end());
  std::optional<double> margin = compute_outer_frame_margin(outer_vertices.begin(),
                                                              outer_vertices.end(),
                                                              lOffset);

  if ( margin )
  {
    double lm = *margin;
    CGAL::Bbox_2 bbox = bbox_2(outer_vertices.begin(), outer_vertices.end());

    double fxmin = bbox.xmin() - lm ;
    double fxmax = bbox.xmax() + lm ;
    double fymin = bbox.ymin() - lm ;
    double fymax = bbox.ymax() + lm ;

    Polygon_2 frame ;
    frame.push_back( Point(fxmin,fymin) );
    frame.push_back( Point(fxmax,fymin) );
    frame.push_back( Point(fxmax,fymax) );
    frame.push_back( Point(fxmin,fymax) );

    std::vector<Polygon_2> outer_as_holes;
    outer_as_holes.reserve(pwhs.size());
    for (const PolygonWithHoles& pwh : pwhs)
      outer_as_holes.emplace_back(pwh.outer_boundary().container().rbegin(),
                                  pwh.outer_boundary().container().rend());

    PolygonWithHoles pwh(frame, outer_as_holes.begin(), outer_as_holes.end());
    PolygonPtrVector off_polys = CGAL::create_interior_skeleton_and_offset_polygons_2(lOffset,pwh);

    // filter outer frame
    std::swap(off_polys[0], off_polys.back());
    off_polys.pop_back();

    for (PolygonPtr ptr : off_polys)
      ptr->reverse_orientation();

    // offset of holes
    for (const PolygonWithHoles& pwh : pwhs)
    {
      for (PolygonWithHoles::Hole_const_iterator hit=pwh.holes_begin();
                                                 hit!=pwh.holes_end();
                                                 ++hit)
      {
        Polygon_2 h = *hit;
        h.reverse_orientation();
        PolygonPtrVector off_hole = CGAL::create_interior_skeleton_and_offset_polygons_2(lOffset,h);
        off_polys.insert(off_polys.end(), off_hole.begin(), off_hole.end());
      }
    }

    return CGAL::arrange_offset_polygons_2<PolygonWithHoles>(off_polys);
  }

  return PolygonWithHolesPtrVector();
}

int main()
{
  std::vector<PolygonWithHoles> pwhs;

  for (int i=0; i<4; ++i)
  {
    Polygon_2 outer;
    outer.push_back( Point(i+0+i*10, 0) );
    outer.push_back( Point(i+0+(i+1)*10, 0) );
    outer.push_back( Point(i+0+(i+1)*10, 10) );
    outer.push_back( Point(i+0+i*10, 10) );
    pwhs.emplace_back(outer);

    Polygon_2 hole;
    hole.push_back( Point(i+3+i*10,3) ) ;
    hole.push_back( Point(i+6+i*10,3) ) ;
    hole.push_back( Point(i+6+i*10,6) ) ;
    hole.push_back( Point(i+3+i*10,6) ) ;
    pwhs[i].add_hole( hole ) ;
  }

  double lOffset = 1.1 ;

  PolygonWithHolesPtrVector offset_poly_with_holes = exterior_offset_of_disjoint_polygons_with_holes(lOffset,pwhs);
  std::cout << offset_poly_with_holes.size() << " offset polygons" << std::endl;

  for (PolygonWithHolesPtr ptr : offset_poly_with_holes)
    CGAL::draw(*ptr);

  return 0;
}
