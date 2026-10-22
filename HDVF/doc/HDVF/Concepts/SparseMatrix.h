/*!
 \ingroup PkgHDVFConcepts
 \cgalConcept

 The concept `SparseMatrix` describes the requirements for sparse matrices optimized for topological computations.
 Traditionally, sparse matrices data structures encode non zero coefficients of (sparse) matrices in order to optimize
 either matrices memory footprint, or linear algebra operations (which usually comes to optimize iterators over non zero coefficients
 and access to coefficients). However, topological operations require slightly different features:

 - fast access to row or columns of matrices (which are actually the images under the application encoded by the matrix)
 - fast block operations (especially along row or columns)

 The `SparseMatrix` concept describes requirements for such sparse matrix. It relies on the model of `SparseChain`
 which encodes sparse row or column vectors. Matrices are either column major or row major (hence they either store column sparse chains or row sparse chains).

 The following constants, called `Storage_format`, encode the major direction of both sparse chains and sparse matrices.
 - `CGAL::OSM::COLUMN` for column-major chains and matrices (which is the default),
 - `CGAL::OSM::ROW` for row-major chains and matrices.

 For instance, given the \f$5\times 4\f$ matrix:
 \f[
 A = \left(\begin{array}{cccc}
 1 & \cdot & \cdot & \cdot \\
 -1 & \cdot & 2 & \cdot\\
 \cdot & \cdot & 1 & \cdot \\
 \cdot & \cdot & \cdot & \cdot \\
 \cdot & \cdot & \cdot & \cdot \\
 \end{array}\right)
 \f]
 where \f$\cdot\f$ means \f$0\f$.
 - A column-major representation of \f$A\f$ is: \f$[0\mapsto c_0, 2\mapsto c_2]\f$ with the column-chains: \f$c_0 = [0\mapsto 1, 1\mapsto -1]\f$ and \f$c_2 = [1\mapsto 2, 2\mapsto 1]\f$.
 - A row-major representation of \f$A\f$ is: \f$[0\mapsto c_0, 1\mapsto c_1, 2\mapsto c_2]\f$ with the row-chains: \f$c_0 = [0\mapsto 1]\f$ and \f$c_1 = [0\mapsto -1, 2\mapsto 2]\f$ and \f$c_2 = [2\mapsto 1]\f$.

 \cgalHasModelsBegin
 \cgalHasModelsBare{`CGAL::OSM::Sparse_matrix_core<IntegralDomainWithoutDivision, Storage_format, SparseChain>`}
 \cgalHasModelsBare{`CGAL::OSM::Sub_sparse_matrix_core<IntegralDomainWithoutDivision, Storage_format, SparseChain>`}
 \cgalHasModelsEnd

 \sa `IntegralDomainWithoutDivision`
 \sa `SparseChain`
 */

class SparseMatrix {

public:
    /// \name Member data
    /// @{

    /*!
     * \brief Matrix storage format (either `CGAL::OSM::ROW` or `CGAL::OSM::COLUMN`).
     */
    int Storage_format;

    /// @}

    /// \name Types
    /// @{

    /*!
     * \brief Type of coefficients stored in the matrix (a model of `IntegralDomainWithoutDivision`).
     */
    typedef unspecified_type Coefficient_ring;

    /*!
     * \brief Type of `SparseChain` used for storage.
     */
    typedef unspecified_type Sparse_chain_type;

    /*!
     * \brief A range of indices of non empty chains.
     */
    typedef unspecified_type Non_zero_chain_indices;

    /*!
     * \brief Type of the chains stored in the matrix.
     */
    typedef unspecified_type Matrix_chain;




    /// \name Creation, Filling
    /// @{

    /*!
     * \brief Creates an empty new sparse matrix object.
     *
     * Default constructor, initialize an empty matrix of type `Storage_format` with coefficients of type `Coefficient_ring`.
     * The default matrix size is 0x0.
     */
    SparseMatrix() ;

    /*!
     * \brief Creates a new sparse matrix object with given rows/columns sizes.
     *
     * Constructor with sizes, initialize an empty matrix of type `Storage_format` with coefficients of type `Coefficient_ring` and a given size along rows/columns.
     */
    SparseMatrix(const size_t rowCount, const size_t columnCount) ;

    /**
     * \brief Cleans a sparse matrix (set all coefficients to zero).
     *
     * Empties all structures of the sparse matrix.
     *
     */
    void nullify();

    /** \brief Sets an identity matrix.
     *
     * \param nrows Number of rows of the matrix.
     * \param ncols Number of columns of the matrix.
     */
    void eye (size_t nrows, size_t  ncols);

    /** \brief Sets an identity matrix.
     *
     * \param dimensions Pair containing the number of rows and number of columns of the matrix.
     */
    void eye (std::pair<size_t, size_t> dimensions);

    /// @}

    /// \name Matrix informations and iterators
    /// @{

    /**
     * \brief Tests if a sparse matrix is null.
     *
     * The function returns `true` if the sparse matrix is null (technically, the sparse matrix is empty, and thus encodes a matrix in which all coefficients are zero) and `false` otherwise.
     */
    bool is_null();

    /**
     * \brief Tests if a `Sparse_matrix_core` is empty.
     *
     * The function return `true` is the `Sparse_matrix_core` is empty (that is, if the number of columns or rows is null) and `false` otherwise.
     */
    bool is_empty();

    /** \brief Returns the storage format of the matrix. */
    int storage_format();
    /**
     * \brief Gets the matrix size.
     *
     * The matrix size as a row/column pair.
     */
    std::pair<size_t, size_t> dimensions() const;

    /**
     * \brief Iterator to the index of the first non null chain.
     *
     * Return an iterator to the index of the first non null chain (the iterator visits indices of non null chains along the major dimension of the matrix).
     *
     * \return The iterator to the index of the first non null chain.
     */
    inline Non_zero_chain_indices::iterator begin() const noexcept;

    /**
     * \brief Iterator past-the-end of chain indices.
     *
     * \return The iterator past-the-end of chain indices.
     */
    inline Non_zero_chain_indices::iterator end() const noexcept;

    /**
     * \brief Reverse iterator to the index of the last non null chain.
     *
     * Return a reverse iiterator to the index of the last non null chain (the iterator visits indices of non null chains, in decreading order, along the major dimension of the matrix).
     *
     * \return The reverse iterator to the index of the last non null chain.
     */
    inline Non_zero_chain_indices::iterator reverse_begin() const noexcept;

    /**
     * \brief Reverse iterator past-the-end of chain indices.
     *
     * \return The reverse iterator past-the-end of chain indices.
     */
    inline Non_zero_chain_indices::iterator reverse_end() const noexcept;

    /// @}

    /// \name Input, output
    /// @{

    /**
     * \brief Writes a sparse matrix in the output stream with a compact "user friendly" format.
     *
     * \param _stream The output stream.
     * \param matrix The matrix to display.
     *
     * \return A reference to the modified stream.
     */
    friend std::ostream& operator<<(std::ostream &_stream, const SparseMatrix &matrix);

    /**
     * \brief  Write a sparse matrix to an output stream using the `.osm` file format.
     *
     * Output a sparse matrix  (the matrix can be reloaded using `read_matrix`), see \ref OSM_FILE_FORMAT "OSM file format documentation".
     */
    friend std::ostream& write_matrix (const SparseMatrix& matrix, std::ostream& out);

    /**
     * \brief  Write a sparse matrix to a file using the `.osm` file format.
     *
     * Output a sparse matrix  (the matrix can be reloaded using `read_matrix`), see \ref OSM_FILE_FORMAT "OSM file format documentation".
     */
    friend std::ostream& write_matrix (const SparseMatrix& matrix, std::string filename);

    /**
     * @brief  Read a sparse matrix from an input stream using the `.osm` file format.
     *
     * Read a sparse matrix (the input must respect the `.osm` file format), .see \ref OSM_FILE_FORMAT "OSM file format documentation".
     */

    template <typename _CT>
    friend std::istream& read_matrix (SparseMatrix& matrix, std::istream& in);

    /**
     * @brief  Read a sparse matrix from a file using the `.osm` file format.
     *
     * Read a sparse matrix (the input must respect the `.osm` file format), .see \ref OSM_FILE_FORMAT "OSM file format documentation".
     */

    template <typename _CT>
    friend std::istream& read_matrix (SparseMatrix& matrix, std::string filename);

    /// @}

    /// \name Linear algebra operators
    /// @{

    /**
     * \brief Perform addition of matrices and assign the result to `matrix`.
     *
     * Perform standard linear algebra addition. Matrices must have the same `Coefficient_ring` but can have different `Storage_format`. The addition is optimized for each combination of `Storage_format`.
     *
     * @param matrix The first matrix.
     * @param other The second matrix.
     * @return A reference to the modified matrix `matrix` containing the result.
     */
    friend SparseMatrix& operator+=(SparseMatrix &matrix, const SparseMatrix &other);

    /**
     * \brief  Perform addition of matrices and assign the result to `matrix`.
     *
     * Perform standard linear algebra addition. Matrices must have the same `Coefficient_ring` but can have different `Storage_format`. The addition is optimized for each combination of `Storage_format`.
     *
     * @param first The first matrix.
     * @param second The second matrix.
     * @return A reference to the modified matrix `matrix` containing the result.
     */
    friend SparseMatrix operator+(const SparseMatrix &first, const SparseMatrix &second);


    /**
     * \brief Subtracts two matrices into a new matrix.
     *
     * Subtracts each coefficient of the matrix `other` and returns a new matrix (of the same type as `matrix`) representing the result (when possible, prefer `-=` for efficiency).
     *
     * \pre Matrices must have the same `Coefficient_ring` but can have different `Storage_format`.
     *
     * \warning Will raise an error if the other matrix is not the same `Coefficient_ring`.
     *
     * \param matrix The first matrix.
     * \param other The second matrix.
     *
     * \return A new matrix representing the result.
     */
    SparseMatrix& operator-=(SparseMatrix &matrix, const SparseMatrix &other);

    /**
     * \brief Subtracts two matrices into a new matrix.
     *
     * Subtracts each coefficient of the matrix `other` and returns a new matrix (of the same type as `this`) representing the result (when possible, prefer `-=` for efficiency).
     *
     * \pre Matrices must have the same `Coefficient_ring` but can have different `Storage_format`.
     *
     * \warning Will raise an error if the other matrix is not the same `Coefficient_ring`.
     *
     * @param first The first matrix.
     * @param second The second matrix.
     *
     * \return A new matrix representing the result.
     */
    friend SparseMatrix operator-(const SparseMatrix &first, const SparseMatrix &second);

    /**
     * \brief Computes the negative of a matrix (unary operator).
     *
     * \return The resulting matrix.
     */
    friend SparseMatrix operator-(const SparseMatrix& matrix);

    /**
     * \brief Applies factor on each coefficients into a new matrix.
     *
     * This method creates a new matrix obtained by multiplying the matrix by a scalar factor `lambda`.
     *
     * If `lambda` is zero, the function comes to nullify the matrix (when possible, prefer `*=` for efficiency).
     *
     * \param lambda The factor to apply.
     * \param matrix The matrix.
     *
     * \return A new matrix representing the result.
     */
    friend SparseMatrix operator*(const Coefficient_ring& lambda, const SparseMatrix &matrix);

    /**
     * \brief Applies factor on each coefficients into a new matrix.
     *
     * This method creates a new matrix obtained by multiplying the matrix by a scalar factor `lambda`.
     *
     * If `lambda` is zero, the function comes to nullify the matrix (when possible, prefer `*=` for efficiency).
     *
     * \param lambda The factor to apply.
     *
     * \return A new matrix representing the result.
     */
    friend SparseMatrix operator*(const SparseMatrix &matrix, const Coefficient_ring& lambda);

    /**
     * \brief Applies factor on each coefficient and assign.
     *
     * This method multiplies the matrix by a scalar factor `lambda`.
     */
    SparseMatrix& operator*=(const Coefficient_ring& lambda);

    /**
     * @brief  Perform multiplication of matrices and assign the result to `matrix`.
     *
     * Perform standard linear algebra multiplication. Matrices must have the same `Coefficient_ring` but can have different `Storage_format`. The multiplication is optimized for each combination of `Storage_format`.
     *
     * @param matrix The first matrix.
     * @param other The second matrix.
     * @return A reference to the modified matrix `matrix` containing the result.
     */
    friend SparseMatrix& operator*=(SparseMatrix &matrix, const SparseMatrix &other);

    /**
     * @brief  Perform multiplication of matrices and return a new column-major matrix.
     *
     * Perform standard linear algebra multiplication of matrices and return a new column-major matrix (when possible, prefer `*=` for efficiency). Matrices must have the same `Coefficient_ring` but can have different `Storage_format`. The multiplication is optimized (using standard definition or block multiplications) for each combination of `Storage_format`. However, efficiency depends on `Storage_format` (when possible, prefer row-major by column-major multiplications).
     *
     * @param first The first matrix.
     * @param second The second matrix.
     * @return The result of the matrix multiplication, column-based.
     */
    friend SparseMatrix operator*(const SparseMatrix &first, const SparseMatrix &second);

    /**
     * @brief  Perform multiplication of matrices and return a new row-major matrix.
     *
     * Perform standard linear algebra multiplication of matrices and return a new row-major matrix (when possible, prefer `*=` for efficiency). Matrices must have the same `Coefficient_ring` but can have different `Storage_format`. The multiplication is optimized (using standard definition or block multiplications) for each combination of `Storage_format`. However, efficiency depends on `Storage_format`.
     *
     * @param first The first matrix.
     * @param second The second matrix.
     * @return The result of the matrix multiplication, row-based.
     */
    friend SparseMatrix operator%(const SparseMatrix &first, const SparseMatrix &second);

    /**
     * @brief  Perform multiplication of a sparse matrix (column or row major) and a column-chain. The function returns a new column-major chain.
     *
     * Perform standard linear algebra multiplication of a matrix and a column-chain (ie.\ matrix / column vector multiplication) and returns a new column-major chain. Both arguments must have the same `Coefficient_ring` but the matrix can have any `Storage_format` (and the multiplication is optimized for each of them).
     *
     * @param matrix The matrix.
     * @param column The column-major chain.
     * @return The result of the matrix multiplication, column-based.
     *
     * \pre the matrix and the chain must have the same `Coefficient_ring`.
     * \pre `column.is_column()` must be `true`
     */
    friend SparseChain operator*(const SparseMatrix& matrix, const SparseChain& column);

    /**
     * @brief  Perform multiplication of a row-chain and a sparse matrix (column or row major). The function returns a new row-major chain.
     *
     * Perform standard linear algebra multiplication of a row-chain and a matrix (ie.\ row vector / matrix multiplication) and returns a new row-major chain. Both arguments must have the same `Coefficient_ring` but the matrix can have any `Storage_format` (and the multiplication is optimized for each of them).
     *
     * @param row The row-major chain.
     * @param matrix The matrix.
     * @return The result of the matrix multiplication, row-based.
     *
     * \pre the matrix and the chain must have the same `Coefficient_ring`.
     * \pre `row.is_row()` must be `true`
     */
    friend SparseChain operator*(const SparseChain& row, const SparseMatrix& matrix);

    /**
     * \brief Transposes a matrix.
     *
     * \return A new matrix where the `Storage_format` has been swapped between COLUMN and ROW and data chains have been transposed.
     */
    SparseMatrix transpose();


    /// @}

    /// \name Access and blocks operations
    /// @{

    /**
     * \brief Gets the value of a chain from a const matrix.
     *
     * For COLUMN matrices, returns the `index`-th column, for ROW matrices, returns the `index`-th row.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param index The chain index.
     *
     * \return The chain stored at given index.
     */
    Matrix_chain operator[](size_t index) const;

    /**
     * \brief Sets a given coefficient in `matrix`.
     *
     * Assign the scalar `d` to the coefficient on row `i` and column `j`.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param matrix Reference on the matrix to modify.
     * \param i The row index.
     * \param j The column index.
     * \param d The value.
     */
    friend void set_coefficient(SparseMatrix& matrix, size_t i, size_t j, const Coefficient_ring d);

    /**
     * \brief Gets a given coefficient.
     *
     * Returns the coefficient on row `i` and column `j` of the matrix.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param matrix Constant reference on the matrix.
     * \param i The row index.
     * \param j The column index.
     *
     * \return The value of the given coefficient.
     */
    friend Coefficient_ring get_coefficient(const SparseMatrix& matrix, size_t i, size_t j);

    /**
     * @brief Gets the value of the column at a given `index` from the matrix (whatever the `Storage_format` of the matrix).
     *
     * For column-matrices, the operation comes to copy non zero coefficients of the column; for row-matrices a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     * Whenever possible, prefer `cget_column()`.
     *
     * \warning The matrix will perform boundary check.
     *
     * @param matrix The  matrix considered.
     * @param index The coefficient index.
     * @return The column at given index.
     */
    friend SparseChain get_column(const SparseMatrix &matrix, size_t index);

    /**
     * @brief Get sthe value of the row at a given `index` from the matrix (whatever the `Storage_format` of the matrix).
     *
     * For row-matrices, the operation comes to copy non zero coefficients of the row; for column-matrices a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     * Whenever possible, prefer `cget_row()`.
     *
     * \warning The matrix will perform boundary check.
     *
     * @param matrix The  matrix considered.
     * @param index The coefficient index.
     * @return The row at given index.
     */
    friend SparseChain get_row(const SparseMatrix &matrix, size_t index);

    /**
     * \brief Gets a const reference over a column from a column matrix.
     *
     * Constant time get.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param matrix The matrix considered.
     * \param index The column index.
     *
     * \return A constant reference over the column stored at given index.
     * \pre `matrix.is_column()` must be `true`.
     */
    const SparseChain& cget_column(const SparseMatrix& matrix, size_t index);

    /**
     * \brief Gets a constant reference over a row from a row matrix
     *
     * Constant time get.
     *
     * \warning The matrix will perform boundary check.
     *
     * \param matrix The matrix considered.
     * \param index The row index.
     *
     * \return A const reference over the row stored at given index.
     * \pre `matrix.is_row()` must be `true`.
     */
    const SparseChain& cget_row(const SparseMatrix& matrix, size_t index);

    /**
     * \brief Set the value of the column at a given `index` from the matrix to `column` (whatever the `Storage_format` of the matrix).
     *
     * \note For column-matrices, the operation copies non zero coefficients of the `column` chain, for row-matrices a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     *
     * \warning The matrix will perform boundary check.
     *
     * @param matrix The  matrix.
     * @param index The column index.
     * @param column The new column value.
     * \pre `column.is_column()` must be `true`
     */
    void set_column(SparseMatrix &matrix, size_t index, const SparseChain &column);

    /**
     * @brief Set the value of the row at a given `index` from the matrix to `chain` (whatever the `StorageFormat` of the matrix).
     *
     * \note For row-matrices, the operation copies non zero coefficients of the `row` chain, for column-matrices a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     *
     * \warning The matrix will perform boundary check.
     *
     * @param matrix The  matrix.
     * @param index The row index.
     * @param row The new row value.
     * \pre `row.is_row()` must be `true`
     */
    void set_row(SparseMatrix &matrix, size_t index, const SparseChain &row);

    /**
     * @brief Build and return a new matrix by copying a `matrix` and removing  the chain at a given `index`. The result is thus a block of the initial matrix. When possible, for efficiency, prefer `/=`.
     *
     * \note Will return a copy of the matrix if given vector is empty.
     *
     * \warning The matrix will perform boundary check.
     *
     * @param matrix The initial matrix.
     * @param index The index of the chain to remove.
     *
     * @return A new matrix containing the result.
     */
    friend SparseMatrix operator/(const SparseMatrix &matrix, size_t index);

    /**
     * @brief Build and return a new matrix by copying a `matrix` and removing  the chains with indices in the `indices` vector. The result is thus a block of the initial matrix. When possible, for efficiency, prefer `/=`.
     *
     * \note Will return a copy of the matrix if given vector is empty.
     *
     * \warning The matrix will perform boundary check.
     *
     * @param matrix The initial matrix.
     * @param indices The indices of the chains to remove.
     *
     * @return A new matrix containing the result.
     */
    friend SparseMatrix operator/(const SparseMatrix &matrix, const std::vector<size_t> &indices);

    /**
     * @brief Remove  the chain at a given `index`. The result is thus a block of the initial matrix.
     *
     * \warning The matrix will perform boundary check.
     *
     * @param index The index of the chain to remove.
     *
     * @return A reference on the modified matrix.
     */
    SparseMatrix& operator/=(size_t index);

    /**
     * @brief Build and return a new matrix by copying a `matrix` and removing  all chains with indices in the `indices` vector. The result is thus a block of the initial matrix. When possible, for efficiency, prefer `/=`.
     *
     * \note Will return a copy of the matrix if given vector is empty.
     *
     * \warning The matrix will perform boundary check.
     *
     * @param indices The indices of the chains to remove.
     *
     * @return A new matrix containing the result.
     */
    SparseMatrix& operator/=(const std::vector<size_t> &indices);

    /**
     * \brief Removes a column from the matrix.
     *
     * Removes column of index `index` whatever the `Storage_format` of the matrix. For column matrices, it just comes to the `/=` operator and for row matrices, it entails a traversal of the matrix.
     *
     * \param matrix Reference on the matrix to modify.
     * \param index The index to remove.
     *
     * \return A reference over the modified matrix.
     */
    friend SparseMatrix& remove_column(SparseMatrix& matrix, size_t index);

    /**
     * \brief Removes a row from the matrix.
     *
     * Removes row of index `index` whatever the `Storage_format` of the matrix. For row matrices, it just comes to the `/=` operator and for column matrices, it entails a traversal of the matrix.
     *
     * \param matrix Reference on the matrix to modify.
     * \param index The index to remove.
     *
     * \return A reference over the modified matrix.
     */
    friend SparseMatrix& remove_row(SparseMatrix& matrix, size_t index);

    /**
     * \brief Removes a coefficient from the matrix.
     *
     * Removes coefficient at row `i` and column `j`.
     *
     * \param matrix Reference on the matrix to modify.
     * \param i Index of the row
     * \param j Index of the column
     *
     * \return The modified matrix representing the result.
     */
    friend SparseMatrix& remove_coefficient(SparseMatrix& matrix, size_t i, size_t j);

    /**
     * \brief Swaps columns in a matrix.
     *
     * Swaps columns at indices `i` and `j` whatever the `StorageFormat` of the matrix. For column matrices, the operation runs in constant time, while it entails a traversal of the matrix for row matrices.
     *
     * \param matrix Reference on the matrix to modify.
     * \param i Index of the first column.
     * \param j Index of the second column.
     *
     * \return A reference over the modified matrix.
     */
        friend Sparse_matrix_core& swap_columns(Sparse_matrix_core& matrix, size_t i, size_t j);

    /**
    * \brief Swaps rows in a matrix.
    *
    * Swaps rows at indices `i` and `j` whatever the `Storage_format` of the matrix. For row matrices, the operation runs in constant time, while it entails a traversal of the matrix for column matrices.

    \param matrix Reference on the matrix to modify.
    \param i Index of the first row.
    \param j Index of the second row.

    \return A reference over the modified matrix.
    */
   friend Sparse_matrix_core& swap_rows(Sparse_matrix_core& matrix, size_t i, size_t j);



    /// @}

};
