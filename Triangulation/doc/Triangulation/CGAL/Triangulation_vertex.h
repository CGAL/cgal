
namespace CGAL {

/*!
\ingroup PkgTriangulationsVertexCellClasses

The class `Triangulation_vertex` is a model of the concept `TriangulationVertex`. It is
used by default for representing vertices in the class
`Triangulation<TriangulationTraits_, TriangulationDataStructure_>`.

A `Triangulation_vertex` stores a point and an incident full cell.


\tparam TriangulationTraits_ must be a model of the concept `TriangulationTraits`. It
provides geometric types and predicates for use in the
`Triangulation<TriangulationTraits_, TriangulationDataStructure_>` class. It is of interest here for its
declaration of the `Point` type.

\tparam Data is an optional type of data to be stored in the vertex class. The
class template `Triangulation_vertex` accepts that no second parameter be specified. In
this case, `Data` defaults to `CGAL::No_vertex_data`.
`CGAL::No_vertex_data` can be explicitely specified to allow to access the
third parameter.

\tparam TriangulationDSVertex_ must be a model of the concept `TriangulationDSVertex`. The
class template `Triangulation_vertex` accepts that no third parameter be specified. It
also accepts the tag `CGAL::Default` as third parameter. In both cases,
`TriangulationDSVertex_` defaults to `CGAL::Triangulation_ds_vertex<>`.

\cgalModels `TriangulationVertex` Additionally, the class
`Triangulation_vertex` provides the following types, constructors
and methods:

\sa `Triangulation_full_cell<TriangulationTraits_, Data, TriangulationDSFullCell_>`
\sa `Triangulation_data_structure<Dimensionality, TriangulationDSVertex_, TriangulationDSFullCell_>`
\sa `Triangulation<TriangulationTraits_, TriangulationDataStructure_>`
\sa `Delaunay_triangulation<DelaunayTriangulationTraits_, TriangulationDataStructure_>`
*/
template< typename TriangulationTraits_, typename Data, typename TriangulationDSVertex_ >
class Triangulation_vertex {
public:

/// \name Types
/// @{

/*!
The point type.
*/
typedef TriangulationTraits_::Point_d Point;

/*!
The type of the additional data stored in the
vertex. If you read a `Triangulation_vertex` from a stream (a file) or write a `Triangulation_vertex` to a stream, then streaming operators `<<` and `>>` must be provided for this
type.
*/
typedef Data Data;

/// @}

/// \name Creation
/// @{

/*!
Constructs a vertex with incident full cell
`c`. The vertex is embedded at point `p` and the parameter `t` is
passed to the `Data` constructor.
*/
template< typename T>
Triangulation_vertex(Full_cell_handle c,
const Point & p, const T & t);

/*!
Same as above, but without incident full cell.
*/
template< typename T> Triangulation_vertex(const Point
& p, const T & t);

/*!
Same as above, but with default-constructed `Point` and `Data`.
*/
Triangulation_vertex();

/// @}

/// \name Data access
/// @{

/*!
Returns a const reference to the stored data.
*/
const Data & data() const;

/*!
Returns a non-const reference to the stored data.
*/
Data & data();

/*!
Inputs the non-combinatorial information given by the vertex, i.e.,
the point and other possible information. The data of type `Data` is
also read.
*/
std::istream & operator>>(std::istream & is, Triangulation_vertex & v);

/*!
Outputs the non-combinatorial information given by the vertex, i.e.,
the point and other possible information. The data of type `Data` is
also written.
*/
std::ostream & operator<<(std::ostream & os, const Triangulation_vertex & v);

/// @}

}; /* end Triangulation_vertex */
} /* end namespace CGAL */
