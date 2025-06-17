/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `AbstractChainComplex` describes the requirements for (topological) chain complexes associated to abstract complexes used in the concept `CGAL::HDVF`.
 
It provides methods to:
 
 - get the dimension of the complex, the number of cells in each dimension
 - get the boundary and co-boundary of cell(s)
 - get the vertices of a given cell
 - output the complex in text format
 
 Cells are indexed along each dimension and thus identified by their index together with their dimension.
 

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Abstract_simplicial_chain_complex<Ring>`}
\cgalHasModelsBare{`CGAL::Sub_chain_complex_mask<Ring, AbstractChainComplex>`}
\cgalHasModelsEnd

*/

class AbstractChainComplex
{
public:
/// \name Types
/// @{

/*!
 Type of column-major chains (returned by the boundary operator)
 */
typedef OSM::Sparse_chain<CoefficientType, OSM::COLUMN> CChain;

/*!
 Type of row-major chains (returned by the co-boundary operator)
 */
typedef OSM::Sparse_chain<CoefficientType, OSM::ROW> RChain ;

/*!
 Type of column-major sparse matrices (used to store the boundary operator)
 */
typedef OSM::Sparse_matrix<CoefficientType, OSM::COLUMN> CMatrix;

/// @}

/// \name Operators
/// @{

/*!
Affectation operator.
 
The operator creates a "fresh" copy of `complex`.
 */
    AbstractChainComplex& operator= (const AbstractChainComplex& complex);
    
/// @}
    
/// \name Access functions
/// @{

/*!
Returns the dimension of the complex, that is, the largest dimension of cells.
 */
int dim();

/*!
Returns the number of cells of dimension `q`.
If `q` is negative of larger than the dimension of the complex, returns 0.
 */
size_t nb_cells(int q);

/*!
\brief Returns all boundary matrices.
 
The function returns constant reference to a vector of column-major sparse matrices. The `q`-th element of this vector is the matrix of \f$\partial_q\f$, which gives the boundary of cells of dimension `q`(as a linear combination of `q`-1 cells).

 */
const vector<CMatrix> & get_bnd_matrices() const;

/*!
 \brief Returns the boundary matrix of dimension `q` (ie. the matrix of \f$\partial_q}\f$).
 
The function returns a column-major sparse matrices.
 */
const CMatrix & get_bnd_matrix(int q) const;

/*!
 \brief Returns the boundary of the cell of index `id_cell` in dimension `q`.

This boundary is a finite linear combination of cells of dimension `q`-1. It is encoded as a column-major chain (which maps each cell with a non-zero coefficient to this coefficient). This boundary  is thus the `id_cell`-th column of the boundary matrix in dimension `q`.
 
*/
CChain d(size_t id_cell, int q);

    
/*!
 \brief Returns the co-boundary of the cell of index `id_cell` in dimension `q`.

This boundary is a finite linear combination of cells of dimension `q`+1. It is encoded as a row-major chain (which maps each cell with a non-zero coefficient to this coefficient). This co-boundary  is thus the `id_cell`-th row of the boundary matrix in dimension `q`+1.
*/
RChain cod(size_t id_cell, int q);

/*!
 \brief Returns the vertices of a given cell (that is, the indices of its faces of dimension 0).
 
*/
std::vector<size_t> bottom_faces(size_t id_cell, int q) const;
    
/*!
 * \brief Returns the cofaces of a given chain in dimension `q`.
 *
 * The resulting chain lies in dimension `q`+1 and is null if this dimension exceeds the dimension of the complex.
*/
template <typename CoefficientT, int ChainTypeF>
CChain cofaces_chain (SparseChain<CoefficientT, ChainTypeF> chain, int q) const;

/// @}

/// \name Output functions
/// @{

/*!
Outputs the chain complex in text mode.
 
By default, outputs the complex to `std::cout`.
*/
std::ostream& print_complex(std::ostream& out = std::cout) const;
    
/// @}

};
