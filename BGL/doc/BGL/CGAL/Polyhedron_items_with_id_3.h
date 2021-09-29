
namespace CGAL {

/*!
\ingroup BGLGraphExternalIndices

The class `Polyhedron_items_with_id_3` is a model of the `PolyhedronItems_3`
concept. It provides definitions for vertices with points, halfedges,
and faces with plane equations, all of them with an additional integer
field which can be used to index the items in a \bgl algorithm.
The polyhedron traits class must provide the respective types for
the point and the plane equation.
Vertices and facets both contain a halfedge handle to an incident
halfedge.

\cgalModels `PolyhedronItems_3`

\cgalHeading{Operations}

Supported as required by the `PolyhedronItems_3` concept.

\cgalHeading{Additional Methods in All Three Items}

\code
int id() const; // Returns the index.
int& id(); // Returns a reference to the index stored in the item.
\endcode

\sa `CGAL::Polyhedron_3<Traits>`
\sa `CGAL::Polyhedron_min_items_3`
\sa `CGAL::HalfedgeDS_min_items`
\sa `CGAL::HalfedgeDS_items_2`
\sa `CGAL::HalfedgeDS_vertex_max_base_with_id<Refs>`
\sa `CGAL::HalfedgeDS_halfedge_max_base_with_id<Refs>`
\sa `CGAL::HalfedgeDS_face_max_base_with_id<Refs>`

*/

class Polyhedron_items_with_id_3 {
public:

  template < class Refs, class Traits>
  struct Vertex_wrapper {
    typedef typename Traits::Point_3 Point;
    typedef CGAL::Tag_true Supports_vertex_point;
  };

  template < class Refs, class Traits>
  struct Face_wrapper {
    typedef typename Traits::Plane_3 Plane;
    typedef CGAL::Tag_true Supports_face_plane;
  };
  /// \name Creation
  /// @{

  /*!
    %Default constructor.
  */
  Polyhedron_items_with_id_3();

  /// @}
}; /* end Polyhedron_items_with_id_3 */


/*!
\ingroup BGLGraphExternalIndices

Given a `CGAL::Polyhedron_3`,
for each simplex type (vertex, halfedge, facet) associates an index from
0 to the number of simplices minus 1 to each simplex of `hds`.
All simplex types must provide an `id()` method return a reference to a variable
that can be assigned a `std::size_t`. An item class
suited for this use is `CGAL::Polyhedron_items_with_id_3`.
*/
template<class HalfedgeDS_with_id>
void set_halfedgeds_items_id ( Polyhedron_with_id& P );

} /* end namespace CGAL */
