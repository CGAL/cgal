#ifndef CGAL_MESH_3_TRIANGULATION_MESHER_LEVEL_TRAITS_3_H
#define CGAL_MESH_3_TRIANGULATION_MESHER_LEVEL_TRAITS_3_H

#include <list>

namespace CGAL {

  namespace Mesh_3 {
    namespace details {

      template <typename Tag, typename Tr>
      struct Type_of_points
      {
        typedef typename Tr::Point Point;
      };

      template <typename Tr>
      struct Type_of_points<Tag_true, Tr>
      {
        typedef typename Tr::Weighted_point Point;
      };

    } // end namespace Mesh_3::details
  } // end namespace Mesh_3

template <typename Tr>
struct Triangulation_mesher_level_traits_3
{
  typedef Tr Triangulation;

  typedef typename  Mesh_3::details::Type_of_points<typename Tr::Weighted_tag,
                                                    Tr>::Point Point;

  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Facet Facet;

  class Zone {
    typedef std::list<Cell_handle> Cells;
    typedef std::list<Facet> Facets;
  public:
    typedef typename Cells::iterator Cells_iterator;
    typedef typename Facets::iterator Facets_iterator;
    typedef typename Cells::const_iterator Cells_const_iterator;
    typedef typename Facets::const_iterator Facets_const_iterator;

    typedef typename Tr::Locate_type Locate_type;

    Locate_type locate_type;
    Cell_handle cell;
    int i, j;

    Cells cells;
    Facets boundary_facets;
    Facets internal_facets;
  };

  Zone get_conflicts_zone(Tr& t, const Point& p) const
  {
    Zone zone;

    zone.cell = t.locate(p, zone.locate_type, zone.i, zone.j);

    if( zone.locate_type == Tr::VERTEX ) return zone;

    t.find_conflicts(p, zone.cell,
                     std::back_inserter(zone.boundary_facets),
                     std::back_inserter(zone.cells),
                     std::back_inserter(zone.internal_facets));
    return zone;
  }

  Vertex_handle insert(Tr&t, const Point& p, Zone& zone)
  {
    if( zone.locate_type == Tr::VERTEX 
	) return zone.cell->vertex(zone.i);

    const Facet& f = *(zone.boundary_facets.begin());
    return t.insert_in_hole(p, zone.cells.begin(), zone.cells.end(),
                            f.first, f.second);
  }

}; // end Triangulation_mesher_level_traits_3

}; // end namespace CGAL

#endif // CGAL_MESH_3_TRIANGULATION_MESHER_LEVEL_TRAITS_3_H
