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
\cgalHasModelsBare{`CGAL::Abstract_simplicial_chain_complex<CoefficientType>`}
\cgalHasModelsEnd

*/

class AbstractChainComplex
{
    /// \name Types
    /// @{
    
    /*!
     Type of column-major chains (returned by the boundary operator)
     */
    typedef OSM::Chain<CoefficientType, OSM::COLUMN> CChain;
    
    /*!
     Type of row-major chains (returned by the co-boundary operator)
     */
    typedef OSM::Chain<CoefficientType, OSM::ROW> RChain ;
    
    /*!
     Type of column-major sparse matrices (used to store the boundary operator)
     */
    typedef OSM::SparseMatrix<CoefficientType, OSM::COLUMN> CMatrix;
    
    /// @}
    
    /// \name Operations
    /// @{
    
    /*!
        Returns the dimension of the complex (ie. largest dimension of cells)
     */
    int dim() ;
    
    /*!
        Returns the number of cells of dimension `dim`.
     
        If `dim` is negative of larger than the dimension of the complex, returns 0.
     */
    int nb_cells(int dim) ;
    
    /*!
        Returns all boundary matrices.
     
        The function returns a vector of column-major sparse matrices. The `q`-th element of this vector is the matrix of \f$\partial d_q\f$, which gives the boundary of cells of dimension `q`(as a linear combination of `q`-1 cells).
     */
    const vector<CMatrix> & get_bnd_matrices() const ;
    
    /*!
        Returns the boundary of the cell of index `id_cell` in dimension `dim`.
     
        This boundary is a finite linear combination of cells of dimension `dim`-1. It is encoded as a column-major chain (which maps each cell with a non-zero coefficient to this coefficient). The boundary of the this cell is thus the `id_cell`-th column of the boundary matrix in dimension `dim`.
     */
    CChain d(int id_cell, int dim) ;
    
    /// @}
    
};
