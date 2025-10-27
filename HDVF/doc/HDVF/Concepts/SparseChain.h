/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `SparseChain` describes the requirements for sparse vectors (called *sparse chains* in homology) optimized for topological computations. More precisely, `SparseChain` provides all the operations on chains required by the `SparseMatrix` concept.

 `SparseChain` encodes non zero coefficients of a sparse chain.

 `SparseChain`  can be either row or column vectors. The following constants, called `StorageFormat`, encode the direction of sparse chains (and sparse matrices).
 - `CGAL::OSM::COLUMN` for column-major chains and matrices (which is the default),
 - `CGAL::OSM::ROW` for row-major chains and matrices.


\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::OSM::Sparse_chain<IntegralDomainWithoutDivision, StorageFormat>`}
\cgalHasModelsEnd

 \sa `IntegralDomainWithoutDivision`
 \sa `SparseMatrix`
*/

class SparseChain {

public:
    /// \name Types
    /// @{

    /*!
     * \brief Type of coefficients stored in the matrix (a model of `IntegralDomainWithoutDivision`).
     */
    typedef unspecified_type Coefficient_ring;

    /*!
     * \brief Matrix and chain storage format (either OSM::ROW or OSM::COLUMN).
     */
    typedef int Storage_format;

    /*!
     * \brief SparseChain iterator type.
     */
    typedef unspecified_type iterator;

    /*!
     * \brief SparseChain constant iterator type.
     */
    typedef unspecified_type const_iterator;

    /// @}

    /// \name Creation, filling
    /// @{

    /*!
     * \brief Creates new empty sparse chain.
     *
     * Creates a sparse chain encoding an empty linear combination of cells.
     */
    SparseChain();

    /*!
     * \brief Creates new empty sparse chain (ie. zero-chain) of given size.
     *
     * Constructor with size, initializes an empty sparse chain encoding a linear combination of cells with all coefficients null.
     */
    SparseChain(size_t chainSize);

    /*!
     * \brief Creates new sparse chain by copy.
     *
     * Copy constructor, initialize a sparse chain from an existing sparse chain of same `ChaintypeFlag`.
     */
    SparseChain(const SparseChain &other);

    /*!
     * \brief Assigns to other chain.
     *
     * Assign to other chain coefficient-wise, equivalent to copying it.
     *
     * SparseChain must have the same `CoefficientRing`.
     */
    SparseChain& operator=(const SparseChain &other);

    /// @}

    /// \name Matrix informations and iterators
    /// @{

    /*!
     * \brief Returns the dimension of the basis (that is, size of the chain).
     */
    size_t dimension() const;

    /*!
     * \brief Iterator to the beginning of the chain.
     *
     * The function returns an iterator to the first non zero index.
     */
    iterator begin() noexcept;

    /*!
     * \brief Constant iterator to the beginning of the chain.
     *
     * The function returns a constant iterator to the first non zero index.
     */
    const_iterator begin() const noexcept;

    /*!
     * \brief Constant iterator to the beginning of the chain.
     *
     * The function returns a constant iterator to the first non zero index.
     */
    const_iterator cbegin() const noexcept;

    /*!
     * \brief Past-the-endiIterator of the chain.
     *
     * The function returns an iterator past the end of the chain.
     */
    iterator end() noexcept;

    /*!
     * \brief Constant past-the-enditerator of the chain.
     *
     * The function returns a constant iterator past the end of the chain.
     */
    const_iterator end() const noexcept;

    /*!
     * \brief Constant past-the-end iterator of the chain.
     *
     * The function returns a constant iterator past the endi of the chain.
     */
    const_iterator cend() const noexcept;

    /// @}

    /// \name Output
    /// @{

    /*!
     * \brief Inserts `chain` in the output stream.
     */
    friend std::ostream& operator<<(std::ostream &stream, const SparseChain& chain);

    /// @}

    /// \name Linear algebra operators
    /// @{

    /*!
     * \brief Adds two chains together.
     *
     * Add two chains and return the result in a new matrix.
     * Chains must have the same `CoefficientRing` and the same `StorageFormat`.
     */
    friend SparseChain operator+(const SparseChain &first, const SparseChain &second);

    /*!
     * \brief Subtracts a chain from another chain.
     *
     * Subtract two chains and return the result in a new matrix.
     * Chains must have the same `CoefficientRing` and the same `StorageFormat`.
     */
    friend SparseChain operator-(const SparseChain &first, const SparseChain &second);

    /*!
     * \brief Applies multiplication on each coefficient.
     */
    friend SparseChain operator*(const CoefficientRing& lambda, const SparseChain &chain);

    /*!
     * \brief Applies multiplication on each coefficient.
     */
    friend SparseChain operator*(const SparseChain &_chain, const CoefficientRing& lambda);

    /*!
     * \brief Performs matrix multiplication between two chains (COLUMN x ROW) and return a COLUMN matrix.
     *
     * Generate a column-based matrix from the matrix multiplication and return it.
     * Chains must have the same `CoefficientRing`.
     */
    friend Sparse_matrix<CoefficientRing, COLUMN> operator*(const Sparse_chain<CoefficientRing, COLUMN> &column, const Sparse_chain<CoefficientRing, ROW> &row);

    /*!
     * \brief Performs matrix multiplication between two chains (COLUMN x ROW) and return a ROW matrix.
     *
     * Generate a row-based matrix from the matrix multiplication and return it.
     * Chains must have the same `CoefficientRing`.
     */
    friend Sparse_matrix<CoefficientRing, ROW> operator%(const Sparse_chain<CoefficientRing, COLUMN> &column, const Sparse_chain<CoefficientRing, ROW> &row);

    /*!
     * \brief Performs dot product between two chains (ROW x COLUMN).
     *
     * Chains must have the same `CoefficientRing`.
     */
    friend CoefficientRing operator*(const Sparse_chain<CoefficientRing, ROW> &row, const Sparse_chain<CoefficientRing, COLUMN> &column);

    /*!
     * \brief Adds a chain to `this`.
     *
     * Add a chain to `this`.
     * Chains must have the same `CoefficientRing` and the same `StorageFormat`.
     */
    SparseChain& operator+=(const SparseChain &_other);

    /*!
     * \brief Subtracts a chain from `this`.
     *
     * Subtract a chain from `this`.
     * Chains must have the same `CoefficientRing` and the same `StorageFormat`.
     */
    SparseChain& operator-=(const SparseChain &_other);

    /*!
     * \brief Applies multiplication on each coefficient of `this`.
     */
    SparseChain& operator*=(const CoefficientRing& lambda);

    /*!
     * \brief Transposes a SparseChain.
     *
     * The result is a chain with `StorageFormat` switched between COLUMN and ROW.
     */
    SparseChain transpose();

    /// @}

    /// \name Access and modifications
    /// @{

    /*!
     * \brief Compares two chains.
     */
    bool operator==(const SparseChain& other_chain);

    /*!
     * \brief Gets the value of a coefficient of the chain.
     */
    CoefficientRing operator[](size_t index);

    /*!
     * \brief Gets the value of a coefficient of the chain.
     */
    CoefficientRing get_coefficient(size_t index) const ;

    /**
     * \brief Sets a given coefficient of the chain.
     */
    void set_coefficient(size_t index, CoefficientRing d);

    /*!
     * \brief Checks if a coefficient of the chain is null.
     */
    bool is_null(size_t index) const;

    /*!
     * \brief Checks if the chain is null.
     */
    bool is_null() const;

    /*!
     * \brief Gets a sub-chain from the chain.
     *
     * Return a new chain where all coefficients of indices provided in the vector are removed.
     */
    friend SparseChain operator/(const SparseChain &chain, const std::vector<size_t> &indices);

    /*!
     * \brief Gets a sub-chain from the chain.
     *
     * Return a new chain where the coefficients at a given index is removed.
     */
    friend SparseChain operator/(const SparseChain &chain, size_t indices);

    /*!
     * \brief Restricts the chain to a sub-chain by removing indices.
     *
     * Removes all indices provided in the vector from the chain. Return a reference to the modified chain.
     */
    SparseChain& operator/=(const std::vector<size_t> &indexes);

    /**
     * \brief Restricts the chain to a sub-chain by removing a given index.
     *
     * Removes the index provided from the chain. Return a reference to the modified chain.
     */
    SparseChain& operator/=(size_t index);

    /**
     * \brief Removes all coefficients from the chain.
     *
     * The function comes to set all coefficients to zero.
     */
    void nullify();

    /*!
     * \brief Checks if chain is a column.
     */
    bool is_column() const;

    /*!
     * \brief Checks if chain is a row.
     */
    bool is_row() const;
    /// @}
};
