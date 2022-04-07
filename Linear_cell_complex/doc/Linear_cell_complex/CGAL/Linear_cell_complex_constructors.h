namespace CGAL {

/*!
\ingroup PkgLinearCellComplexConstructions

Imports an embedded plane graph read from `ais` into `lcc`, a model of the `LinearCellComplex` concept. Objects are added in `lcc`, existing darts are not modified. Returns a dart created during the import.
\pre \link GenericMap::dimension `LCC::dimension`\endlink \f$ \geq\f$ 2 and \link LinearCellComplex::ambient_dimension `LCC::ambient_dimension`\endlink==2.

\cgalHeading{File format}

The file format must be the following. First the number of vertices and the number of edges of the planar graph. Then, for each vertex of the planar graph, the coordinates of the \f$ i^{\mbox{th}}\f$ vertex (two numbers for \f$ x\f$ and \f$ y\f$ coordinates). The first vertex index is 0. Then for each edge of the planar graph, the two indices of the two vertices (two numbers between 0 and the number of vertices minus 1).

Here a small example:
\verbatim
5 6
1.0 3.0 0.0 2.0 2.0 2.0 0.0 0.0 2.0 0.0
0 1 0 2 1 2 1 3 2 4 3 4
\endverbatim

\image html lcc_import_graph.svg "Example of import_graph reading the above file as istream, middle for combinatorial map as combinatorial data-structure, right for generalized maps."
\image latex lcc_import_graph.svg "Example of import_graph reading the above file as istream, middle for combinatorial map as combinatorial data-structure, right for generalized maps."

<B>Left</B>: A planar graph embedded in the plane with <I>P0</I>=(1.0,3.0), <I>P1</I>=(0.0,2.0), <I>P2</I>=(2.0,2.0), <I>P3</I>=(0.0,0.0), <I>P4</I>=(2.0,0.0).
<B>Middle</B>: the 2D linear cell complex reconstructed if combinatorial maps are the combinatorial data-structure.
<B>Right</B>: the 2D linear cell complex reconstructed if generalized maps are the combinatorial data-structure.

\sa `CGAL::import_from_triangulation_3<LCC,Triangulation>`
\sa `CGAL::import_from_polyhedron_3<LCC,Polyhedron>`

*/
template<class LCC>
typename LCC::Dart_handle import_from_plane_graph(LCC& lcc,
std::istream& ais);

} /* namespace CGAL */

