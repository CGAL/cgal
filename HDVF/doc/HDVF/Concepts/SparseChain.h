/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `SparseChain` describes the requirements for sparse vectors (called *sparse chains* in homology) optimized for topological computations.
More precisely, `SparseChain` provides all the operations on chains required by the `SparseMatrix` concept.

 `SparseChain` encodes non zero coefficients of a sparse vector.

 `SparseChain`  can be either row or column vectors. The following integer constants, called `Storage_format`, encode the direction of sparse chains (and sparse matrices).
 - `CGAL::OSM::COLUMN` for column-major chains and matrices (which is the default),
 - `CGAL::OSM::ROW` for row-major chains and matrices.


\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::OSM::Sparse_chain<CoefficientRing, StorageFormat>`}
\cgalHasModelsEnd

 \sa `IntegralDomainWithoutDivision`
 \sa `SparseMatrix`
*/

class SparseChain {

public:
    /// \name Member data
    /// @{

    /*!
     * \brief Chain storage format (either `CGAL::OSM::ROW` or `CGAL::OSM::COLUMN`).
     */
    int Storage_format;

    /// \name Types
    /// @{

    /*!
     * \brief Type of coefficients stored in the matrix (a model of `IntegralDomainWithoutDivision`).
     */
    typedef unspecified_type Coefficient_ring;

    /** \brief Type of derived sparse matrices. */
    template <typename CR, int SF>
    using Sparse_matrix_derived = OSM::Sparse_matrix_structure<OSM::Sparse_chain>::Sparse_matrix_type<CR, SF>;

    /*!
     * \brief Sparse chain iterator type.
     */
    typedef unspecified_type iterator;

    /*!
     * \brief Sparse chain constant iterator type.
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
     * \brief Creates new empty sparse chain (ie.\ zero-chain) of given size.
     *
     * Constructor with size, initializes an empty sparse chain encoding a linear combination of cells with all coefficients null.
     */
    SparseChain(size_t chainSize);

    /**
     * \brief Assignment.
     *
     *The assignment operator converts the storage format if needed (complexity \f$\mathscr O(k)\f$ with \f$k\f$ the number of non zero coefficients of `other`).
     *
     * \pre Chains must have the same `Coefficient_ring`.
     *
     * \param other Chain to copy into `this`.
     */
    Sparse_chain& operator=(const Sparse_chain &other);

    /// @}

    /// \name Matrix informations and iterators
    /// @{

    /*!
     * \brief Size of the chain
     *
     * \return Dimension of the basis (that is, the size of the chain).
     */
    size_t dimension() const;

    /*!
     * \brief Iterator to the beginning of the chain.
     *
     * \warning The chain is stored in an unordered map for speed reason.
     *
     * \return The function returns an iterator to the first non zero index.
     */
    iterator begin() noexcept;

    /*!
     * \brief Constant iterator to the beginning of the chain.
     *
     * \warning The chain is stored unordered for speed reason.
     *
     * \return The function returns a constant iterator to the first non zero index.
     */
    const_iterator begin() const noexcept;

    /*!
     * \brief Constant iterator to the beginning of the chain.
     *
     * \warning The chain is stored unordered for speed reason.
     *
     * \return The function returns a constant iterator to the first non zero index.
     */
    const_iterator cbegin() const noexcept;

    /*!
     * \brief Iterator to the end of the chain.
     *
     * \warning The chain is stored unordered for speed reason.
     *
     * \return The function returns an iterator past-the-end of the chain.
     */
    iterator end() noexcept;

    /*!
     * \brief Constant iterator to the end of the chain.
     *
     * \warning The chain is stored unordered for speed reason.
     *
     * \return The function returns a constant iterator past-the end of the chain.
     */
    const_iterator end() const noexcept;

    /*!
     * \brief Constant iterator to the end of the chain.
     *
     * \warning The chain is stored unordered for speed reason.
     *
     * \return The function returns a constant iterator to the ending of the chain.
     */
    const_iterator cend() const noexcept;

    /// @}

    /// \name Output
    /// @{

    /*!
     * \brief Inserts `chain` in the output stream with a compact "user friendly" format.
     *
     * \param stream The output stream.
     * \param chain The chain to display.
     *
     * \return A reference to the modified stream.
     */
    friend std::ostream& operator<<(std::ostream &stream, const SparseChain& chain);

    /**
     * @brief  Write a sparse chain to an output stream  using the `.osm` file format.
     *
     * Outputs a sparse chain  (the chain can be reloaded using `read_chain`), see \ref OSM_FILE_FORMAT "OSM file format documentation".
     */
    template <typename _CT, int _SF>
    friend std::ostream& write_chain (const Sparse_chain<_CT, _SF>& chain, std::ostream& out);

    /**
     * @brief  Write a sparse chain to a file  using the `.osm` file format.
     *
     * Outputs a sparse chain  (the chain can be reloaded using `read_chain`), see \ref OSM_FILE_FORMAT "OSM file format documentation".
     */
    template <typename _CT, int _SF>
    friend void write_chain (const Sparse_chain<_CT, _SF>& chain, std::string filename);


    /**
     * \brief Reads a sparse chain from an input stream using the `.osm` file format.
     *
     * Loads a sparse chain  (the chain can be exported using `write_chain`), see \ref OSM_FILE_FORMAT "OSM file format documentation".
     */
   template <typename _CT, int _SF>
   friend std::istream& read_chain (Sparse_chain<_CT, _SF>& chain, std::istream& in);

    /**
     * \brief Reads a sparse chain from a file using the `.osm` file format.
     *
     * Loads a sparse chain  (the chain can be exported using `write_chain`), see \ref OSM_FILE_FORMAT "OSM file format documentation".
     */
    template <typename _CT, int _SF>
    friend void read_chain (Sparse_chain<_CT, _SF>& chain, std::string filename);

    /// @}

    /// \name Linear algebra operators
    /// @{

    /*!
     * \brief Adds two chains.
     *
     * Adds `other` chain to current chain and returns the result in a new matrix.
     *
     * \pre Chains must have the same `Coefficient_ring` and the same `Storage_format`.
     *
     * \warning Will raise an error if the two chains are not the same `Coefficient_ring`.
     * \warning Will raise a compilation error if the two chains don't have the same `Storage_format`.
     *
     * \param other The other chain.
     *
     * \return A new chain representing the result.
     */
    Sparse_chain operator+(const Sparse_chain &other);

    /*!
     * \brief Subtracts a chain from current chain.
     *
     * Subtracts `other` chain from current chain and returns the result in a new matrix.
     *
     * \pre Chains must have the same `Coefficient_ring` and the same `Storage_format`.
     *
     * \warning Will raise an error if the two chains are not the same `Coefficient_ring`.
     * \warning Will raise a compilation error if the two chains don't have the same `Storage_format`.
     *
     * \param other The other chain.
     *
     * \return A new chain representing the result.
     */
    Sparse_chain operator-(const Sparse_chain &other);

    /*!
     * \brief Multiplies each coefficient of `chain` by `lambda`.
     *
     * \param lambda The factor to apply.
     * \param chain The  chain.
     *
     * \return A new chain representing the result.
     */
    friend SparseChain operator*(const Coefficient_ring& lambda, const SparseChain &chain);

    /*!
     * \brief Multiplies each coefficient of `chain` by `lambda`.
     *
     * \param lambda The factor to apply.
     *
     * \return A new chain representing the result.
     */
    Sparse_chain operator*(const Coefficient_ring& lambda);

    /*!
     * \brief Performs matrix multiplication between two chains (COLUMN x ROW) and return a COLUMN matrix.
     *
     * Generate a column-based matrix from the matrix multiplication and return it.
     *
     * \pre Chains must have the same `Coefficient_ring`.
     * \pre `column.is_column()`  and `row.is_row()`  must be `true`
     *
     * \warning Will raise an error if chains do not have the same `Coefficient_ring`.
     *
     * \param column The column chain.
     * \param row The row chain.
     *
     * \return The result of the matrix multiplication, column-based.
     */
    friend Sparse_matrix_derived<Coefficient_ring, COLUMN> operator*(const SparseChain &column, const SparseChain &row);

    /*!
     * \brief Performs matrix multiplication between two chains (COLUMN x ROW) and return a ROW matrix.
     *
     * Generate a row-based matrix from the matrix multiplication and return it.
     *
     * \pre Chains must have the same `Coefficient_ring`.
     * \pre `column.is_column()`  and `row.is_row()`  must be `true`
     *
     * \warning Will raise an error if chains do not have the same `Coefficient_ring`.
     *
     * \param column The column chain.
     * \param row The row chain.
     *
     * \return The result of the matrix multiplication, row-based.
     */
    friend SparseMatrix operator%(const Sparse_chain &column, const Sparse_chain &row);

    /*!
     * \brief Performs dot product between two chains (ROW x COLUMN).
     *
     * \pre Chains must have the same `Coefficient_ring`.
     * \pre `row.is_row()` and `column.is_column()` must be `true`
     *
     * \warning Will raise an error if the chains do not have the same `Coefficient_ring`.
     *
     * \param row The row chain.
     * \param column The column chain.
     *
     * \return The result of type `Coefficient_ring`.
     */
    friend Coefficient_ring operator*(const Sparse_chain &row, const Sparse_chain &column);

    /*!
     * \brief Adds a chain to `this`.
     *
     * Add a chain to `this`.
     *
     * \pre Chains must have the same `Coefficient_ring` and the same `Storage_format`.
     *
     * \warning Will raise an error if the two chains are not the same `Coefficient_ring`.
     * \warning Will raise an error if the two chains don't have the same `Storage_format`.
     *
     * \param other The other chain.
     *
     * \return The modified chain representing the result.
     */
    SparseChain& operator+=(const SparseChain &other);

    /*!
     * \brief Sustracts a chain to `this`.
     *
     * Subtract a chain to `this`.
     *
     * \pre Chains must have the same `Coefficient_ring` and the same `Storage_format`.
     *
     * \warning Will raise an error if the two chains are not the same `Coefficient_ring`.
     * \warning Will raise an error if the two chains don't have the same `Storage_format`.
     *
     * \param other The other chain.
     *
     * \return The modified chain representing the result.
     */
    SparseChain& operator-=(const SparseChain &other);

    /*!
     * \brief Multiplies each coefficient of `this` by `lambda`.
     *
     * If `lambda` is null, this function comes to nullify the chain.
     *
     * \param lambda The factor to apply.
     *
     * \return The modified chain representing the result.
     */
    SparseChain& operator*=(const Coefficient_ring& lambda);

    /*!
     * \brief Transposes a sparse chain.
     *
     * The result is a chain with `Storage_format` switched between COLUMN and ROW.
     *
     * \return A new chain where the `Storage_format` is changed.
     */
    SparseChain transpose();

    /// @}

    /// \name Test, access and modifications
    /// @{

    /*!
     * \brief Checks if chain is a column.
     */
    bool is_column() const;

    /*!
     * \brief Checks if chain is a row.
     */
    bool is_row() const;

    /*!
     * \brief Checks if a coefficient is null (that is, not store in the sparse chain).
     *
     * \param index The index to check.
     * \return `true` if the data at given index is null, `false`.
     */
    bool is_null(size_t index) const;

    /*!
     * \brief Checks if the chain is null.
     *
     * \return `true` if the chain is null (that is, contains no non zero coefficients), `false` otherwise.
     */
    bool is_null() const;

    /*!
     * \brief Tests if two chains are equal.
     *
     * If the storage formats of chains differ, return `false`, otherwise, tests if all coefficients are equal.
     *
     * \param chain First chain to compare.
     * \param other Second chain to compare.
     */
    friend bool operator==(const Sparse_chain& chain, const Sparse_chain& other);

    /*!
     * \brief Gets the value of a coefficient index of the chain.
     *
     * \warning The chain will perform boundary check.
     *
     * \param index The coefficient index.
     *
     * \return The value of the coefficient.
     */
    Coefficient_ring get_coefficient(size_t index) const ;

    /**
     * \brief Sets coefficient index of the chain to value `d`.
     *
     * Sets the value of the coefficient in the chain at `index`.
     *
     * \warning The chain will perform boundary check.
     *
     * \param index The coefficient index.
     * \param d Value of the coefficient.
     */
    void set_coefficient(size_t index, Coefficient_ring d);

    /*!
     * \brief Gets a subchain from the chain.
     *
     * Returns a copy of `chain` where all coefficients of indices provided in the vector are removed.
     *
     * \note Will return a copy of the chain if `indices` is empty.
     *
     * \param indices The indices to remove.
     *
     * \return A new chain representing the result.
     */
    friend SparseChain operator/(const SparseChain &chain, const std::vector<size_t> &indices);

    /*!
     * \brief Gets a subchain from the chain.
     *
     * Returns a copy of `chain` where the coefficient at a given index is removed.
     *
     * \param index The index to remove.
     */
    friend SparseChain operator/(const SparseChain &chain, size_t index);

    /*!
     * \brief Restricts the chain to a sub-chain by removing indices.
     *
     * Removes all indices provided in the vector from the chain. Returns a reference to the modified chain.
     *
     * \note Will not alter the chain if `indices` is empty.
     *
     * \param indices The indices to remove.
     *
     * \return A reference to the modified chain.
     */
    SparseChain& operator/=(const std::vector<size_t> &indices);

    /**
     * \brief Restricts the chain to a sub-chain by removing a given index.
     *
     * Removes the index provided from the chain.
     *
     * \param index The index to remove.
     *
     * \return A reference to the modified chain.
     */
    SparseChain& operator/=(size_t index);

    /**
     * \brief Removes all coefficients from the chain.
     *
     * The function comes to set all coefficients to zero.
     */
    void nullify();
    /// @}
};
