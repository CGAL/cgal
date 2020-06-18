
namespace CGAL {

/*!
\ingroup PkgTDS3Classes

The class `Triangulation_utils_3` defines operations on the indices of vertices
and neighbors within a cell.

\cgalFigureAnchor{Triangulation3figutils}
<center>
  <table border=0>
  <tr>
    <td>
    \image html utils.png
    \image latex utils.png
    </td>
    <td>
    \image html utils-next-around-edge.png
    \image latex utils-next-around-edge.png
    </td>
  </tr>
  </table>
</center>
\cgalFigureCaptionBegin{Triangulation3figutils}
Operations on indices.
\cgalFigureCaptionEnd
*/

struct Triangulation_utils_3 {

/// \name Operations
/// @{

/*!
In dimension 3, index of the neighbor `n` that is next to the current cell,
when turning positively around an oriented edge whose endpoints are
indexed `i` and `j`. According to the usual numbering of
vertices and neighbors in a given cell, it is also the index of the vertex
opposite to this neighbor `n`. (see \cgalFigureRef{Triangulation3figutils}).
\pre `( i < 4 ) && ( j < 4 ) && ( i != j )`.
*/
static unsigned int next_around_edge(unsigned int i, unsigned int j);

/*!
In dimension 3, index of the `j`'th vertex in counterclockwise order on the face opposite to vertex with `i` of the cell.
\pre `( i < 4 ) && ( j < 3 )`.
*/
  static int vertex_triple_index(const int i, const int j);

/*!
Has a meaning only in dimension 2.

Computes the index of the vertex that is next to the vertex numbered
`i` in counterclockwise direction. (see
\cgalFigureRef{Triangulation3figutils}).
\pre `i<3`.
*/
static unsigned int ccw(unsigned int i);

/*!
Same for clockwise.
*/
static unsigned int cw(unsigned int i);

/// @}

}; /* end Triangulation_utils_3 */
} /* end namespace CGAL */
