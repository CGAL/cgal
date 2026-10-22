/*!
 \ingroup PkgHDVFConcepts
 \cgalConcept

 The concept `AbstractChainComplex` describes the requirements for (topological) chain complexes associated to abstract complexes used in the concept `HDVF`.

 It provides methods to:

 - get the dimension of the complex, the number of cells in each dimension
 - get the boundary and co-boundary of cell(s)
 - get the vertices of a given cell
 - output the complex in text format

 Cells are indexed along each dimension and thus identified by their index together with their dimension.


 \cgalHasModelsBegin
 \cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Abstract_simplicial_chain_complex<IntegralDomainWithoutDivision>`}
 \cgalHasModelsEnd

 */

class AbstractChainComplex
{
public:
    /// \name Types
    /// @{

    /*!
     Type of coefficients ring used to compute homology, model of `IntegralDomainWithoutDivision`
     */
    typedef unspecified_type Coefficient_ring;

    /**
     Type of structure of sparse matrix used.

     A model of `SparseMatrixStructure` providing a model of `SparseChain` and an associated model of `SparseMatrix`.

     \warning Check if the concept `SparseMatrixStructure`and documentation are correct...
     */
    typedef unspecified_type Sparse_matrix_structure;

    /*!
     Type of column-major chains (returned by the boundary operator)
     */
    typedef Sparse_matrix_structure:: template Sparse_chain_type<Coefficient_ring, CGAL::OSM::COLUMN> Column_chain;

    /*!
     Type of row-major chains (returned by the co-boundary operator)
     */
    typedef Sparse_matrix_structure:: template Sparse_chain_type<Coefficient_ring, CGAL::OSM::ROW> Row_chain ;

    /*!
     Type of column-major sparse matrices (used to store the boundary operator)
     */
    typedef Sparse_matrix_structure:: template Sparse_matrix_type<Coefficient_ring, CGAL::OSM::COLUMN> Column_matrix;

    /// @}

    /// \name Verification operations
    /// @{

    /** \brief Checks if `q` belongs to the range of dimensions of cells in the complex. */
    bool is_valid_cell_dimension(int q);
    /** \brief Checks if the cell of index `i` and dimension `q` belongs to the range of dimensions of cells in the complex. */
    bool is_valid_cell(size_t i, int q);

    /// @}


    /// \name Access functions
    /// @{

    /*!
     \brief Returns the dimension of the complex.
     *
     * Returns the dimension of the simplicial complex (i.e. largest dimension of cells).
     *
     * \return The dimension of the complex.
     */
    int dimension();

    /** \brief Get (unique) object Id.
     * For comparison of constant references to the complex.
     */
    size_t id () const;

    /*!
     \brief Returns the number of cells in a given dimension.
     *
     * \param q Dimension along which the number of cells is returned.
     *
     * \return Number of cells in dimension `q`. If `q` is not in the interval `[0,d]` with `d` the dimension of the complex, the function returns an empty column chain.
     */
    size_t number_of_cells(int q);

    /*!
     \brief Returns a constant reference to the vector of boundary matrices (along each dimension).
     *
     * Returns a constant reference to the vector of boundary matrices along each dimension. The q-th element of this vector is a column-major sparse matrix containing the boundaries of q-cells (i.e. rows encode q-1 cells and columns q cells).
     *
     * \return A constant reference to the vector of column-major boundary matrices along each dimension.
     */
    const vector<Column_matrix> & boundary_matrices() const;

    /*!
     \brief Returns a constant reference over the boundary matrix of dimension `q` (i.e.\ column-major matrix of \f$\partial_q\f$).
     *
     * It is a column-major sparse matrix containing the boundaries of q-cells (i.e. rows encode q-1 cells and columns q cells).
     *
     * \param q Dimension of the boundary matrix (i.e. columns will contain the boundary of dimension q cells).
     *
     * \return A constant reference over the column-major sparse matrix containing the matrix of the boundary operator of dimension q.
     */
    const Column_matrix & boundary_matrix(int q) const;

    /*!
     * \brief Returns the boundary of the cell `id_cell` in dimension `q`.
     *
     * Returns a copy of the column-major chain stored in the boundary matrix of dimension dim: boundary of the cell `id_cell` in dimension `q`.
     *
     * \param id_cell %Index of the cell.
     * \param q Dimension of the cell.
     *
     * \return The column-major chain containing the boundary of the cell id_cell in dimension `q`. If `q` is  not in the interval `[0,d]` with `d` the dimension of the complex, the function returns an empty column chain.
     */
    Column_chain d(size_t id_cell, int q);


    /*!
     \brief Returns the co-boundary of the cell `id_cell` in dimension `q`.
     *
     * Returns a row-major chain containing the co-boundary of the cell `id_cell` in dimension `q` (so actually a row of the boundary matrix).
     *
     * \warning As the boundary matrix is stored column-major, this entails crossing the full matrix to extract the row coefficients (O(number of non empty columns))
     *
     * \param id_cell %Index of the cell.
     * \param q Dimension of the cell.
     *
     * \return The row-major chain containing the co-boundary of the cell id_cell in dimension `q`. If `q` is not in the interval `[0,d]` with `d` the dimension of the complex, the function returns an empty row chain.
     */
    Row_chain cod(size_t id_cell, int q);

    /*!
     \brief Returns the vertices of a given cell (that is, the indices of its faces of dimension 0).

     */
    std::vector<size_t> bottom_faces(size_t id_cell, int q) const;

    /*!
     * \brief Returns the cofaces of a given chain in dimension `q`.
     *
     * The resulting chain, whatever the storage format of the input, is column-major, lies in dimension `q`+1 and is null if this dimension exceeds the dimension of the complex.
     *
     * \tparam CoefficientR `CoefficientRing` of the chain.
     * \tparam StorageF `StorageFormat` of the chain.
     */
    template <typename CoefficientR, int StorageF>
    Column_chain cofaces_chain (SparseChain<CoefficientR, StorageF> chain, int q) const;

    /// @}

    /// \name Input, output
    /// @{
    ///
    /*!
     \brief Inserts the chain complex in text mode in the stream.

     Displays the number of cells in each dimension and the boundary matrix in each dimension.

     \param out The output stream.
     \param complex Constant reference over the complex considered.
     */
    friend std::ostream& operator<<(std::ostream& out, const AbstractChainComplex& complex) const;
    /// @}
};




