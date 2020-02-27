namespace CGAL {

///@{
/*!
\ingroup PkgIOTriangulation3

The triangulation streamed in `is`, of original type `Tr_src`,  is written into `tr`, of type `Tr_tgt`. As the vertex and cell
 types might be different and incompatible, the creation of new cells and vertices 
is made thanks to the functors `convert_vertex` and `convert_cell`, that convert 
vertex and cell types. For each vertex `v_src` in `in`, the corresponding 
vertex `v_tgt` in `tr` is a copy of the vertex returned by `convert_vertex(v_src)`. 
The same operations are done for cells with the functor convert_cell, except cells
in `tr` are created using the default constructor, and then filled with the data
contained in the stream.


 - A model of `ConvertVertex` must provide two `operator()`s that are responsible
 for converting the source vertex `v_src` into the target vertex:
  - `Tr_tgt::Vertex operator()(const Tr_src::Vertex& v_src) const;` This operator is 
used to create the vertex from `v_src`.
  - `void operator()(const Tr_src::Vertex& v_src, Tr_tgt::Vertex& v_tgt) const;` This
 operator is meant to be used in case heavy data should be transferred to `v_tgt`.
 - A model of ConvertCell must provide an `operator()` that is responsible for 
converting the source cell `c_src` into the target cell:
  - `void operator()(const Tr_src::Cell& c_src, Tr_tgt::Cell& c_tgt) const;` This operator
 is meant to be used in case data should be transferred to `c_tgt`.

\note The triangulation contained in `is` can be obtained with the `operator>>` of a `Triangulation_3`.
*/
template <typename Tr_src, 
          typename Tr_tgt,
          typename ConvertVertex,
          typename ConvertCell>
std::istream& file_input(std::istream& is, Tr_tgt &tr,
                         ConvertVertex convert_vertex = ConvertVertex(),
                         ConvertCell convert_cell = ConvertCell());
///@}
}
