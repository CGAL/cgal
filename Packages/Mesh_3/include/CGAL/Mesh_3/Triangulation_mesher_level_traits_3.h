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

    Cells cells;
    Facets boundary_facets;
    Facets internal_facets;
  };

  Zone get_conflicts_zone(Tr& t, const Point& p) const
  {
    Zone zone;

    const Cell_handle ch = t.locate(p);

    t.find_conflicts(p, ch,
                     std::back_inserter(zone.boundary_facets),
                     std::back_inserter(zone.cells),
                     std::back_inserter(zone.internal_facets));
    return zone;
  }

  Vertex_handle insert(Tr&t, const Point& p, Zone& zone)
  {
    const Facet& f = *(zone.boundary_facets.begin());
    return t.insert_in_hole(p, zone.cells.begin(), zone.cells.end(),
                            f.first, f.second);
  }

}; // end Triangulation_mesher_level_traits_3

}; // end namespace CGAL

#endif // CGAL_MESH_3_TRIANGULATION_MESHER_LEVEL_TRAITS_3_H
