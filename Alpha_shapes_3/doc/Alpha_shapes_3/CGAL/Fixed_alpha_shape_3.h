
namespace CGAL {

/*!
\ingroup PkgAlphaShapes3

The class `Fixed_alpha_shape_3` represents one (fixed) 
alpha shape of points in the 3D space for a real 
\f$ \alpha\f$. It maintains an underlying triangulation 
of the class `Dt` which 
represents connectivity and order among its faces. Each 
\f$ k\f$-dimensional face of the `Dt` is associated with 
a classification that specifies its status in the alpha complex, alpha being fixed. 

Note that this class can be used at the same time to build a <I>basic</I> or 
a <I>weighted</I> Alpha Shape. 

The modifying functions `insert` and `remove` will overwrite
the one inherited from the underlying triangulation class `Dt`.
At the moment, only the static version is implemented.

\cgalHeading{I/O}

The I/O operators are defined for `iostream`, and for 
the window stream provided by \cgal. The format for the iostream 
is an internal format. 

*/
template< typename Dt >
class Fixed_alpha_shape_3 : public Dt {
public:

/// \name Types 
/// @{

/*!
the alpha shape traits type. It has to derive from a triangulation
traits class. For example `Dt::Point` is a Point class.
*/ 
typedef unspecified_type Gt; 

/*!
the number type of alpha. 
*/ 
typedef Gt::FT FT; 

/*!
Enum to classify the simplices of the underlying 
triangulation with respect to a given alpha value. 

Each k-dimensional simplex of the triangulation 
can be classified as EXTERIOR, SINGULAR, REGULAR 
or INTERIOR. 
A \f$ k\f$ simplex is REGULAR if it is on the boundary 
of the alpha complex and belongs to a \f$ k+1\f$ simplex in this complex 
and it is SINGULAR if it is a boundary simplex that is not included in a \f$ k+1\f$ simplex of the complex. 

*/ 
enum Classification_type {EXTERIOR, SINGULAR, REGULAR, INTERIOR}; 

/// @} 

/// \name Creation 
/// @{

/*!
Introduces an empty fixed alpha shape and sets the alpha value to `alpha`. 
*/ 
Fixed_alpha_shape_3(FT alpha = 0); 

/*!
Builds a fixed alpha shape from the triangulation `dt`, 
and sets the alpha value to `alpha`. 
\attention This operation destroys the triangulation. 
*/ 
Fixed_alpha_shape_3(Dt& dt,FT alpha = 0); 

/*!
Builds a fixed alpha shape for the points in the range 
`[first,last)` and sets the alpha value to `alpha`. 
\tparam InputIterator must be an input iterator with value type `Point` (the type point of the underlying triangulation.) 
*/ 
template < class InputIterator > 
Fixed_alpha_shape_3( 
InputIterator first, 
InputIterator last, 
const FT& alpha = 0); 

/// @} 

/// \name Modifiers 
/// @{

/*!

Inserts point p in the underlying triangulation and returns the corresponding vertex. 
The optional argument `start` is used as a starting place for the search. 
The classification types of the new simplices are computed and that of the simplices incident 
to the new ones are updated. 

*/ 
Vertex_handle insert (Point p,Cell_handle start = Cell_handle()); 

/*!

Removes the vertex v from the underlying triangulation. 
The classification types of new simplices and their incident faces are set or reset.

*/ 
void remove (Vertex_handle v); 

/*!
Clears the structure. 
*/ 
void 
clear(); 

/// @} 

/// \name Query Functions 
/// @{

/*!
Returns the \f$ \alpha\f$-value. 
*/ 
const FT& 
get_alpha(void) const; 

/*!
Classifies the cell `c` of the underlying triangulation in the alpha complex. 
*/ 
Classification_type 
classify(Cell_handle c) const; 

/*!
Classifies the facet `f` of the underlying triangulation in the alpha complex. 
*/ 
Classification_type classify(Facet f) const; 

/*!
Classifies the facet of the cell `f` opposite to the vertex with index 
`i` of the underlying triangulation in the alpha complex. 
*/ 
Classification_type classify(Cell_handle f, int i) const; 

/*!
Classifies the edge `e` of the underlying triangulation in the alpha complex. 
*/ 
Classification_type classify(const Edge& e) const; 

/*!
Classifies the vertex `v` of the underlying triangulation in the alpha complex. 
*/ 
Classification_type classify(Vertex_handle v) const; 

/*!
Writes the cells which are of type `type` in the alpha complex 
to the sequence 
pointed to by the output iterator `it`. Returns past the end 
of the output sequence. 
*/ 
template<class OutputIterator> 
OutputIterator get_alpha_shape_cells(OutputIterator it, Classification_type type); 

/*!
Writes the facets which are of type `type` in the alpha complex 
to the sequence pointed to by the output iterator `it`. Returns past the end 
of the output sequence. 
*/ 
template<class OutputIterator> 
OutputIterator get_alpha_shape_facets(OutputIterator it, Classification_type type); 

/*!
Writes the edges which are of type `type` in the alpha complex 
to the sequence 
pointed to by the output iterator `it`. Returns past the end 
of the output sequence. 
*/ 
template<class OutputIterator> 
OutputIterator get_alpha_shape_edges(OutputIterator it, Classification_type type); 

/*!
Writes the vertices which are of type `type` in the alpha complex 
to the sequence pointed to by the output iterator `it`. Returns past the end 
of the output sequence. 
*/ 
template<class OutputIterator> 
OutputIterator get_alpha_shape_vertices(OutputIterator it, Classification_type type); 

/// @}

}; /* end Fixed_alpha_shape_3 */

/*!
Inserts the fixed alpha shape `A` into the stream `os`. 


An overlaoad of `operator<<` must be available for `GT::Point`.
\relates Fixed_alpha_shape_3 
*/ 
ostream& operator<<(ostream& os, const Fixed_alpha_shape_3<Dt>& A); 


} /* end namespace CGAL */
