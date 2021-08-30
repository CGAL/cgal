namespace CGAL {

/*!
\ingroup PkgMatrixSearchRef

\brief computes the maximum (as specified by `compare_strictly`)
entry for each row of `m` and writes the corresponding column
to `t`, i.e., `t[i]` is set to the index of the column
containing the maximum element in row `i`.\ The maximum \f$ m_r\f$
of a row \f$ r\f$ is the leftmost element for which
`compare_strictly`\f$ (m_r,\,x)\f$ is false for all elements \f$ x\f$ in
\f$ r\f$.

The function `monotone_matrix_search()` computes the maxima
for all rows of a totally monotone matrix.

More precisely, monotony for matrices is defined as follows.

Let \f$ K\f$ be a totally ordered set, \f$ M \in K^{(n,\, m)}\f$
a matrix over \f$ K\f$ and for \f$ 0 \le i < n\f$:
\f[
rmax_M(i) :\in \left\{ \min_{0 \le j < m} j \: \left|\:
M[i,\, j] = \max_{0 \le k < m} M[i,\, k] \right.\right\}
\f]
the (leftmost) column containing the maximum entry in row
\f$ i\f$. \f$ M\f$ is called monotone, iff
\f[
\forall\, 0 \le i_1 < i_2 < n\: :\: rmax_M(i_1) \le
rmax_M(i_2)\; .
\f]
\f$ M\f$ is totally monotone, iff all of its submatrices are
monotone (or equivalently: iff all \f$ 2 \times 2\f$ submatrices are
monotone).



\pre `t` points to a structure of size at least `m.number_of_rows()`

\tparam Matrix is a model of `MonotoneMatrixSearchTraits`.
\tparam RandomAccessIC is a model of `RandomAccessIterator` with `int` as value type.
If `compare_strictly` is defined, it is an adaptable
binary function: `Matrix::Value` \f$ \times\f$
`Matrix::Value` \f$ \rightarrow\f$ `bool` describing a strict
(non-reflexive) total ordering on `Matrix::Value`.

\sa `MonotoneMatrixSearchTraits`
\sa `all_furthest_neighbors_2()`
\sa `maximum_area_inscribed_k_gon_2()`
\sa `maximum_perimeter_inscribed_k_gon_2()`
\sa `extremal_polygon_2()`

\cgalHeading{Implementation}

The implementation uses an algorithm by Aggarwal
et al.\cgalCite{akmsw-gamsa-87}. The runtime is linear in the number
of rows and columns of the matrix.

*/
template < class Matrix, class RandomAccessIC,
class Compare_strictly > void monotone_matrix_search( const
Matrix& m, RandomAccessIC t, const Compare_strictly&
compare_strictly = less< Matrix::Value >());

} /* namespace CGAL */

