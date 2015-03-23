#ifndef CGAL_MINKOWSKI_SUM_HOLE_FILTER_2_H
#define CGAL_MINKOWSKI_SUM_HOLE_FILTER_2_H

#include <CGAL/basic.h>
#include <vector>

namespace CGAL {

/*! \class
 * This class applies filter to a polygon with holes,
 * by removing all of its holes that cannot possibly contribute
 * to the Minkowski sum boundary.
 */
template <class Kernel_, class Container_>
class Hole_filter_2
{

private:

  typedef Kernel_ Kernel;
  typedef Container_ Container;

  typedef CGAL::Polygon_with_holes_2<Kernel, Container> Polygon_with_holes_2;

public:

  void operator()(const Polygon_with_holes_2 &pgn1, const Polygon_with_holes_2 &pgn2,
                  Polygon_with_holes_2 &filtered_pgn1) const
  {
    filtered_pgn1 = pgn1;

    Bbox_2 boundary_bbox;
    typename Polygon_with_holes_2::Hole_iterator it;

    std::vector<typename Polygon_with_holes_2::Hole_iterator> to_erase;
    boundary_bbox = pgn2.outer_boundary().bbox();

    it = filtered_pgn1.holes_begin();
    while (it != filtered_pgn1.holes_end())
    {
      Bbox_2 hole_bbox = (*it).bbox();

      if (hole_bbox.ymax()-hole_bbox.ymin() < boundary_bbox.ymax()-boundary_bbox.ymin() ||
          hole_bbox.xmax()-hole_bbox.xmin() < boundary_bbox.xmax()-boundary_bbox.xmin())
      {
        to_erase.push_back(it);
      }
      ++it;
    }

    typename std::vector<typename Polygon_with_holes_2::Hole_iterator>::iterator it2 = to_erase.begin();
    while (it2 != to_erase.end())
    {
      filtered_pgn1.erase_hole(*it2);
      ++it2;
    }
  }
};

} // namespace CGAL

#endif
