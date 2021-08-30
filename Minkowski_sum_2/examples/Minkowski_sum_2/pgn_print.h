#ifndef PGN_PRINT_H
#define PGN_PRINT_H

#include <iostream>
#include <vector>

#include <CGAL/General_polygon_2.h>

//-----------------------------------------------------------------------------
// Pretty-print a CGAL polygon.
//
template <typename Polygon_2> void print_polygon(const Polygon_2 & pgn)
{
  std::cout << "[ " << pgn.size() << " vertices: (";
  typename Polygon_2::Vertex_const_iterator  vit;
  for (vit = pgn.vertices_begin(); vit != pgn.vertices_end(); ++vit)
    std::cout << "(" << *vit << ')';
  std::cout << ") ]" << std::endl;
}

template <typename Traits>
void print_polygon(const CGAL::General_polygon_2<Traits> & pgn)
{
  std::cout << "[ " << pgn.size() << " curves:" << std::endl;
  typename CGAL::General_polygon_2<Traits>::Curve_const_iterator  cit;
  for (cit = pgn.curves_begin(); cit != pgn.curves_end(); ++cit)
    std::cout << *cit;
  std::cout << " ]" << std::endl;
}

//-----------------------------------------------------------------------------
// Pretty-print a polygon with holes.
//
template <typename Polygon_with_holes>
void print_polygon_with_holes(const Polygon_with_holes & pwh)
{
  if (! pwh.is_unbounded()) {
    std::cout << "{ Outer boundary = ";
    print_polygon (pwh.outer_boundary());
  }
  else std::cout << "{ Unbounded polygon." << std::endl;

  unsigned int  k = 1;
  typename Polygon_with_holes::Hole_const_iterator hit;
  std::cout << "  " << pwh.number_of_holes() << " holes:" << std::endl;
  for (hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit, ++k) {
    std::cout << "    Hole #" << k << " = ";
    print_polygon(*hit);
  }
  std::cout << " }" << std::endl;
}

//-----------------------------------------------------------------------------
// Pretty-print a polygon set.
//
template <typename Polygon_set>
void print_polygon_set(const Polygon_set & pgn_set)
{
  typedef typename Polygon_set::Polygon_with_holes_2 Polygon_with_holes;
  typedef std::vector<Polygon_with_holes>            Pgn_with_holes_container;

  Pgn_with_holes_container res(pgn_set.number_of_polygons_with_holes());
  pgn_set.polygons_with_holes(res.begin());
  std::cout << "The result contains " << res.size() << " components:"
            << std::endl;
  typename Pgn_with_holes_container::const_iterator  it;
  for (it = res.begin(); it != res.end(); ++it) {
    std::cout << "--> ";
    print_polygon_with_holes(*it);
  }
}

#endif
