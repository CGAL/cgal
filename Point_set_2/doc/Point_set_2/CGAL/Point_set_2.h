namespace CGAL {

/*!
\ingroup PkgPointSet2

An instance \f$ PS\f$ of the data type `Point_set_2` is a <I>Delaunay Triangulation</I> 
of its vertex set. The class `Point_set_2` is inherited from the \cgal Delaunay triangulation, 
and provides additional nearest neighbor query operations and range searching operations. 

The `Point_set_2` class of \cgal depends on template parameters standing for the 
geometric traits classes used by the point set and by the Delaunay triangulation (Gt) 
and for the triangulation data structure (Tds). 

*/
template< typename Gt, typename Tds >
class Point_set_2 {
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
typedef Hidden_type Triangulation; 

/*! 
the size type of the underlying triangulation. 
*/ 
typedef Hidden_type Triangulation::size_type; 

/*! 
the vertex type of the underlying triangulation. 
*/ 
typedef Hidden_type Triangulation::Vertex; 

/*! 
the edge type of the underlying triangulation. 
*/ 
typedef Hidden_type Triangulation::Edge; 

/*! 
handles to vertices. 
*/ 
typedef Hidden_type Triangulation::Vertex_handle ; 

/// @} 

/// \name Creation 
/// @{

/*! 
creates an empty `Point_set_2` . 
*/ 
Point_set_2(); 

/*! 
creates a `Point_set_2` `PS` of the points in the range 
[\f$ first\f$,\f$ last\f$). 
*/ 
template<class InputIterator> 
Point_set_2(InputIterator first, InputIterator last); 

/// @} 

/// \name Operations 
/// @{

/*! 
if `PS` contains a Vertex \f$ v\f$ with \f$ |pos(v)| = p\f$ 
the result is a handle to \f$ v\f$ otherwise the result is \f$ NULL\f$. 
*/ 
Vertex\_handle lookup(Point p); 

/*! 
computes a handle to a vertex \f$ v\f$ of `PS` that is closest to \f$ p\f$. 
If `PS` is empty, \f$ NULL\f$ is returned. 

*/ 
Vertex\_handle nearest_neighbor(Point p); 

/*! 
computes a handle to a vertex \f$ w\f$ of `PS` that is closest to \f$ v\f$. 
If \f$ v\f$ is the only vertex in `PS` , \f$ NULL\f$ is returned. 

*/ 
Vertex\_handle nearest_neighbor(Vertex_handle v); 

/*! 
computes the \f$ k\f$ nearest neighbors of \f$ p\f$ in `PS`, and places the 
handles to the corresponding vertices as a sequence of objects of type 
Vertex_handle in a container of value type of `res` 
which points to the first object in the sequence. The function 
returns an output iterator pointing to the position beyond the end 
of the sequence. 
*/ 
template<class OutputIterator> 
OutputIterator nearest_neighbors(Point p, size_type k, OutputIterator res); 

/*! 
computes the \f$ k\f$ nearest neighbors of \f$ v\f$, and places them as a sequence of objects of type 
Vertex_handle in a container of value type of `res` 
which points to the first object in the sequence. The function 
returns an output iterator pointing to the position beyond the end 
of the sequence. 
*/ 
template<class OutputIterator> 
OutputIterator nearest_neighbors(Vertex_handle v, size_type k,OutputIterator res); 

/*! 
computes handles to all vertices contained in the closure of disk \f$ C\f$. 
The computed vertex handles will be placed as a sequence of objects in a container of value type 
of `res` 
which points to the first object in the sequence. The function 
returns an output iterator pointing to the position beyond the end 
of the sequence. 

*/ 
template<class OutputIterator> 
OutputIterator range_search(const Circle& C, OutputIterator res); 

/*! 
computes handles to all vertices contained in the closure of the triangle \f$ (a,b,c)\f$. 

\pre \f$ a\f$, \f$ b\f$, and \f$ c\f$ must not be collinear. 
The computed vertex handles will be placed as a sequence of objects in a container of value type 
of `res` 
which points to the first object in the sequence. The function 
returns an output iterator pointing to the position beyond the end 
of the sequence. 

*/ 
template<class OutputIterator> 
OutputIterator range_search(const Point& a, const Point& b, const Point& c,OutputIterator res); 

/*! 
computes handles to all vertices contained in the closure of the iso-rectangle \f$ (a1,b1,c1,d1)\f$. 

\pre \f$ a1\f$ is the upper left point, \f$ b1\f$ the lower left, \f$ c1\f$ the lower right and \f$ d1\f$ the upper right point of the iso rectangle. 
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
