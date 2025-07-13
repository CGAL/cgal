/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `SparseChain` describes the requirements for sparse vectors (called *sparse chains* in homology) optimized for topological computations. More precisely, `SparseChains` provide all the operations on chains required by the `SparseMatrix` concept.

 `SparseChains` encode non zero coefficients of (sparse) chains.
 
 `SparseChain`  can be either row or column vectors. The following constants, called `ChainTypeFlag`, encode the direction of sparse chains (and sparse matrices).
 - `OSM::COLUMN` for column-major chains and matrices (which is the default),
 - `OSM::ROW` for row-major chains and matrices.
 
 
\cgalHasModelsBegin
\cgalHasModelsBare{`OSM::Sparse_chain<Ring, ChainTypeFlag>`}
\cgalHasModelsEnd

 \sa `Ring`
 \sa `SparseMatrix`
*/

class SparseChain {
    
public:
    /// \name Types
    /// @{
    
    /*!
     * \brief Type of coefficients stored in the matrix (a model of `Ring`).
     */
    typedef Ring CoefficientType;
    
    /*!
     * \brief Matrix and chain type (either OSM::ROW or OSM::COLUMN).
     */
    typedef int ChainTypeFlag;

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
     * \brief Create new empty SparseChain.
     */
    SparseChain();
    
    /*!
     * \brief Create new empty SparseChain of given size.
     *
     * Constructor with size, initialize an empty SparseChain.
     */
    SparseChain(size_t chainSize);
    
    /*!
     * \brief Create new SparseChain by copy.
     *
     * Copy constructor, initialize a SparseChain from an existing SparseChain of same `ChaintypeFlag`.
     */
    SparseChain(const Sparse_chain &otherToCopy);
    
    /*!
     * \brief Assign to other chain.
     *
     * Assign to other chain coefficient-wise, equivalent to copying it.
     *
     * SparseChain must have the same `CoefficientType`.
     */
    SparseChain& operator=(const SparseChain &otherToCopy);
    
    /// @}
    
    /// \name Matrix informations and iterators
    /// @{
    
    /*!
     * \brief Dimension of the basis (that is, size of the chain).
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
     * \brief Iterator to the end of the chain.
     *
     * The function returns an iterator to the ending of the chain.
     */
    iterator end() noexcept;
    
    /*!
     * \brief Constant iterator to the end of the chain.
     *
     * The function returns a constant iterator to the ending of the chain.
     */
    const_iterator end() const noexcept;
    
    /*!
     * \brief Constant iterator to the end of the chain.
     *
     * The function returns a constant iterator to the ending of the chain.
     */
    const_iterator cend() const noexcept;
    
    /// @}
    
    /// \name Output
    /// @{
    
    /*!
     * \brief Displays a SparseChain in the output stream.
     */
    friend std::ostream& operator<<(std::ostream &stream, const Sparse_chain& chain);
    
    /// @}
    
    /// \name Linear algebra operators
    /// @{
    
    /*!
     * \brief Adds two chains together.
     * 
     * Adds two chains together and return the result in a new matrix.
     * Chains must have the same `CoefficientType` and the same `ChainTypeFlag`.
     */
    friend Sparse_chain operator+(const Sparse_chain &first, const Sparse_chain &second);
    
    /*!
     * \brief Subtract two chains together.
     * 
     * Subtract two chains together and return the result in a new matrix.
     * Chains must have the same `CoefficientType` and the same `ChainTypeFlag`.
     */
    friend Sparse_chain operator-(const Sparse_chain &first, const Sparse_chain &second);
    
    /*!
     * \brief Apply factor on each coefficients.
     */
    friend Sparse_chain operator*(const CoefficientType& lambda, const Sparse_chain &chain);
    
    /*!
     * \brief Apply factor on each coefficients.
     */
    friend Sparse_chain operator*(const Sparse_chain &_chain, const CoefficientType& lambda);
    
    /*!
     * \brief Perform matrix multiplication between two chains (COLUMN x ROW) and return a COLUMN matrix.
     *
     * Generate a column-based matrix from the matrix multiplication and return it.
     * Chains must have the same `CoefficientType`.
     */
    friend Sparse_matrix<CoefficientType, COLUMN> operator*(const Sparse_chain<CoefficientType, COLUMN> &column, const Sparse_chain<CoefficientType, ROW> &row);
    
    /*!
     * \brief Perform matrix multiplication between two chains (COLUMN x ROW) and return a ROW matrix.
     *
     * Generate a row-based matrix from the matrix multiplication and return it.
     * Chains must have the same `CoefficientType`.
     */
    friend Sparse_matrix<CoefficientType, ROW> operator%(const Sparse_chain<CoefficientType, COLUMN> &column, const Sparse_chain<CoefficientType, ROW> &row);
    
    /*!
     * \brief Perform dot product between two chains (ROW x COLUMN).
     *
     * Chains must have the same `CoefficientType`.
     */
    friend CoefficientType operator*(const Sparse_chain<CoefficientType, ROW> &row, const Sparse_chain<CoefficientType, COLUMN> &column);
    
    /*!
     * \brief Add a chain to `this`.
     *
     * Add a chain to `this`.
     * Chains must have the same `CoefficientType` and the same `ChainTypeFlag`.
     */
    Sparse_chain& operator+=(const Sparse_chain &_other);
    
    /*!
     * \brief Subtract a chain from `this`.
     *
     * Subtract a chain from `this`.
     * Chains must have the same `CoefficientType` and the same `ChainTypeFlag`.
     */
    Sparse_chain& operator-=(const Sparse_chain &_other);
    
    /*!
     * \brief Apply factor on each coefficients of `this`.
     */
    Sparse_chain& operator*=(const CoefficientType& lambda);
    
    /*!
     * \brief Transpose a SparseChain.
     *
     * The result is a chain with `ChainTypeFlag` switched between COLUMN and ROW.
     */
    Sparse_chain transpose();
    
    /// @}
    
    /// \name Access and modifications
    /// @{
    
    /*!
     * \brief Compare two chains.
     */
    bool operator==(const Sparse_chain& other_chain);
    
    /*!
     * \brief Get the value of a coefficient of the chain.
     */
    CoefficientType operator[](size_t index);

    /*!
     * \brief Get the value of a coefficient of the chain.
     */
    CoefficientType get_coef(size_t index) const ;
    
    /**
     * \brief Set a given coefficient of the chain.
     */
    void set_coef(size_t index, CoefficientType d);
    
    /*!
     * \brief Checks if a coefficient of the chain is null.
     */
    bool is_null(size_t index) const;
    
    /*!
     * \brief Checks if the chain is null.
     */
    bool is_null() const;
    
    /*!
     * \brief Get a subchain from the chain.
     *
     * Return a new chain where all coefficients of indices provided in the vector are removed.
     */
    friend Sparse_chain operator/(const Sparse_chain &chain, const std::vector<size_t> &indices);
    
    /*!
     * \brief Get a subchain from the chain.
     *
     * Return a new chain where the coefficients at a given index is removed.
     */
    friend Sparse_chain operator/(const Sparse_chain &chain, size_t indices);
    
    /*!
     * \brief Restrict the chain to a sub-chain by removing indices.
     *
     * Removes all indices provided in the vector from the chain. Return a reference to the modified chain.
     */
    Sparse_chain& operator/=(const std::vector<size_t> &indexes);
    
    /**
     * \brief Restrict the chain to a sub-chain by removing a given index.
     *
     * Removes the index provided from the chain. Return a reference to the modified chain.
     */
    Sparse_chain& operator/=(size_t index);
    
    /**
     * \brief Remove all coefficients from the chain.
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
