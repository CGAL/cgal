/*!
 \ingroup PkgHDVFConcepts
 \cgalConcept

 The concept `SparseMatrixStructure` describes the requirements for a class providing a model of `SparseChain` and an associated model of `SparseMatrix`.

 \cgalHasModelsBegin
 \cgalHasModelsBare{`CGAL::OSM::Sparse_matrix<SparseChain>`}
 \cgalHasModelsBare{`CGAL::OSM::Sub_sparse_matrix<SparseChain>`}
 \cgalHasModelsEnd

 \sa `IntegralDomainWithoutDivision`
 \sa `SparseMatrix`
 \sa `SparseChain`
 */

class SparseMatrixStructure {

public:
    /// \name Types
    /// @{

    /*!
     * \brief Type of sparse chains used for storage (a model of `SparseChain`).
     */
    typedef unspecified_type Sparse_chain_type;

    /**
     * \brief Associated type of sparse matrices (a model of `SparseMatrix`).
     *
     * This model uses `Sparse_chain_type` as internal storage (thus `SparseMatrix::Sparse_chain_type` are `Sparse_chain_type` similar).
     */
    typedef unspecified_type Sparse_matrix_type;

    /// @}

};
