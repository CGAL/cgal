namespace CGAL {

/*!
\ingroup PkgPointSet2Ref

An instance `PS` of the data type `Point_set_2` is a <I>Delaunay Triangulation</I>
of its vertex set. The class `Point_set_2` is inherited from the \cgal Delaunay triangulation,
and provides additional nearest neighbor query operations and range searching operations.

The `Point_set_2` class of \cgal depends on template parameters standing for the
geometric traits classes used by the point set and by the Delaunay triangulation (Gt)
and for the triangulation data structure (Tds).

*/
template< typename Gt,
          typename Tds = Triangulation_data_structure_2<
            Triangulation_vertex_base_2<Gt> >
          >
class Point_set_2 : public Delaunay_triangulation_2<Gt,Tds> {
public:

/// \name Types
/// @{

/*!
the point type
*/
typedef Gt::Point_2 Point;

/*!
the segment type
*/
typedef Gt::Segment_2 Segment;

/*!
the circle type
*/
typedef Gt::Circle_2 Circle;

/*!
the representation field number type.
*/
typedef Gt::FT Numb_type;

/*!
the underlying triangulation type.
*/
typedef unspecified_type Triangulation;

/*!
the size type of the underlying triangulation.
*/
typedef Triangulation::size_type size_type;

/*!
the vertex type of the underlying triangulation.
*/
typedef Triangulation::Vertex Vertex;

/*!
the edge type of the underlying triangulation.
*/
typedef Triangulation::Edge Edge;

/*!
handles to vertices.
*/
typedef Triangulation::Vertex_handle Vertex_handle;

/// @}

/// \name Creation
/// @{

/*!
creates an empty `Point_set_2`.
*/
Point_set_2();

/*!
creates a `Point_set_2` `PS` of the points in the range
[`first`,`last`).
*/
template<class InputIterator>
Point_set_2(InputIterator first, InputIterator last);

/// @}

/// \name Operations
/// @{

/*!
if `PS` contains a vertex `v` with `v.point() == p`
the result is a handle to `v` otherwise the result is `NULL`.
*/
Vertex\_handle lookup(Point p);

/*!
computes a handle to a vertex `v` of `PS` that is closest to `p`.
If `PS` is empty, `NULL` is returned.

*/
Vertex\_handle nearest_neighbor(Point p);

/*!
computes a handle to a vertex `w` of `PS` that is closest to `v`.
If `v` is the only vertex in `PS` , `NULL` is returned.

*/
Vertex\_handle nearest_neighbor(Vertex_handle v);

/*!
computes the `k` nearest neighbors of  `p` in `PS`, and places the
handles to the corresponding vertices as a sequence of objects of type
Vertex_handle in a container of value type of `res`
which points to the first object in the sequence. The function
returns an output iterator pointing to the position beyond the end
of the sequence.
*/
template<class OutputIterator>
OutputIterator nearest_neighbors(Point p, size_type k, OutputIterator res);

/*!
computes the `k` nearest neighbors of `v`, and places them as a sequence of objects of type
Vertex_handle in a container of value type of `res`
which points to the first object in the sequence. The function
returns an output iterator pointing to the position beyond the end
of the sequence.
*/
template<class OutputIterator>
OutputIterator nearest_neighbors(Vertex_handle v, size_type k,OutputIterator res);

/*!
computes handles to all vertices contained in the closure of disk `C`.
The computed vertex handles will be placed as a sequence of objects in a container of value type
of `res`
which points to the first object in the sequence. The function
returns an output iterator pointing to the position beyond the end
of the sequence.

*/
template<class OutputIterator>
OutputIterator range_search(const Circle& C, OutputIterator res);

/*!
computes handles to all vertices contained in the closure of the triangle `(a,b,c)`.

\pre `a`, `b`, and `c` must not be collinear.
The computed vertex handles will be placed as a sequence of objects in a container of value type
of `res`
which points to the first object in the sequence. The function
returns an output iterator pointing to the position beyond the end
of the sequence.

*/
template<class OutputIterator>
OutputIterator range_search(const Point& a, const Point& b, const Point& c,OutputIterator res);

/*!
computes handles to all vertices contained in the closure of the iso-rectangle `(a1,b1,c1,d1)`.

\pre `a1` is the upper left point, `b1` the lower left, `c1` the lower right and `d1` the upper right point of the iso rectangle.
The computed vertex handles will be placed as a sequence of objects in a container of value type
of `res`
which points to the first object in the sequence. The function
returns an output iterator pointing to the position beyond the end
of the sequence.

*/
template<class OutputIterator>
OutputIterator range_search(const Point& a1, const Point& b1, const Point& c1,const Point&
d1,OutputIterator res);

/// @}

}; /* end Point_set_2 */
} /* end namespace CGAL */
