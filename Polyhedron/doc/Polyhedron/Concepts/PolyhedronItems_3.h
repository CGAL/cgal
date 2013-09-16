/*!
\ingroup PkgPolyhedronConcepts
\cgalConcept

The `PolyhedronItems_3` concept refines the `HalfedgeDSItems`
concept. In addition to the requirements stated there, a model for
this concept must fulfill the following requirements for the local
\ref HalfedgeDSItems::Vertex "PolyhedronItems_3::Vertex_wrapper<Refs,Traits>::Vertex" type and
\ref HalfedgeDSItems::Face "PolyhedronItems_3::Face_wrapper<Refs,Traits>::Face" type in order to
support the point for vertices and the optional plane equation for
facets. Note that the items class uses face instead of facet. Only the
polyhedral surface renames faces to facets.

\cgalRefines `HalfedgeDSItems`

\cgalHasModel `CGAL::Polyhedron_items_3`
\cgalHasModel `CGAL::Polyhedron_min_items_3`

\sa `CGAL::Polyhedron_3<Traits>`
\sa `HalfedgeDSItems` 
\sa `CGAL::HalfedgeDS_items_2` 
\sa `CGAL::HalfedgeDS_vertex_base<Refs>` 
\sa `CGAL::HalfedgeDS_halfedge_base<Refs>` 
\sa `CGAL::HalfedgeDS_face_base<Refs>` 

\cgalHeading{Example}

We define our own items class based on the available 
`CGAL::HalfedgeDS_face_base` base class for faces. We derive the 
the `Halfedge_wrapper` without further modifications from the 
`CGAL::HalfedgeDS_items_2`, replace the `Face_wrapper` 
definition with our new definition, and also replace the 
`Vertex_wrapper` with a definition that uses `Point_3` instead 
of `Point_2` as point type. The result is a model for the 
`PolyhedronItems_3` concept similar to the available 
`CGAL::Polyhedron_items_3` class. See also there for another 
illustrative example. 

\code{.cpp} 

#include <CGAL/HalfedgeDS_bases.h> 

struct My_items : public CGAL::HalfedgeDS_items_2 { 
template < class Refs, class Traits> 
struct Vertex_wrapper { 
typedef typename Traits::Point_3 Point; 
typedef CGAL::HalfedgeDS_vertex_base< Refs, CGAL::Tag_true, Point> Vertex; 
}; 
template < class Refs, class Traits> 
struct Face_wrapper { 
typedef typename Traits::Plane_3 Plane; 
typedef CGAL::HalfedgeDS_face_base< Refs, CGAL::Tag_true, Plane> Face; 
}; 
}; 

\endcode 

*/

class PolyhedronItems_3 {
public:
  /// \cgalConcept
  ///
  /// The `PolyhedronItems_3::Vertex` concept refines the
  /// `HalfedgeDSItems::Vertex` concept. This concept adds the requirement
  /// for the support of the geometry.
  ///
  /// \cgalRefines `HalfedgeDSItems::Vertex`
  class Vertex {
  public:
    /// \name Types in PolyhedronItems_3::Vertex_wrapper<Refs,Traits>::Vertex
    /// @{

    /// point type stored in vertices. A `HalfedgeDS
    /// has no dimension, so this type is named Point` and not `Point_3`
    typedef unspecified_type Point;

    /// \f$\equiv\f$ `CGAL::Tag_true`. A point is always required.
    typedef unspecified_type Supports_vertex_point;
    /// @}

    /// \name Creation
    /// @{

    /// Default constructor.
    Vertex();

    /// Constructor from a point
    Vertex(Point);
    /// @}

    /// \name Operations
    /// @{
    Point& point();
    const Point& point() const;
    /// @}
  };

  /// \cgalConcept
  ///
  /// The `PolyhedronItems_3::Face` concept refines the
  /// `HalfedgeDSItems::Face` concept. This concept adds the requirement
  /// for the support of the geometry.
  ///
  /// \cgalRefines `HalfedgeDSItems::Face`
  class Face {
  public:
    /*!
      \name Types in PolyhedronItems_3::Vertex_wrapper<Refs,Traits>::Face
      Types for (optionally) associated geometry in faces. If it is not
      supported the respective type has to be defined, although it can be
      an arbitrary dummy type, such as `void*` or `Tag_false`.
    */
    /// @{
    
    /// plane type stored in faces. A `HalfedgeDS` has no
    /// dimension, so this type is named `Plane` and not
    /// `Plane_3`.
    typedef unspecified_type Plane;


    /// either `CGAL::Tag_true` or `CGAL::Tag_false`.
    typedef unspecified_type Supports_face_plane;
    /// @}

    /// \name Only required when Supports_face_plane == Tag_true
    /// @{
    Plane& plane();
    const Plane& plane() const;
    /// @}
  };
}; /* end PolyhedronItems_3 */

