/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `SparseMatrix` describes the requirements for sparse matrices optimized for topological computations.
Traditionally, sparse matrices data structures encode non zero coefficients of (sparse) matrices in order to optimize either matrices memory footprint, or linear algebra operations (which usually comes to optimize iterators over non zero coefficients and access to coefficients). However, topological operations require slightly different features:
 - fast access to row or columns of matrices (which are actually the images under the application encoded by the matrix)
 - fast block operations (especially along row or columns)

 The `SparseMatrix` concept describes requirements for such sparse matrix. It relies on the model of `SparseChain` which encodes sparse row or column vectors. Matrices are either column major or row major (hence they either store column sparse chains or row sparse chains).

 The following constants, called `StorageFormat`, encode the major direction of both sparse chains and sparse matrices.
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
\cgalHasModelsBare{`CGAL::OSM::Sparse_matrix<IntegralDomainWithoutDivision, StorageFormat>`}
\cgalHasModelsEnd

 \sa `IntegralDomainWithoutDivision`
 \sa `SparseChain`
*/

class SparseMatrix {

public:
    /// \name Types
    /// @{

    /*!
     * \brief Type of coefficients stored in the matrix (a model of `IntegralDomainWithoutDivision`).
     */
    typedef unspecified_type Coefficient_ring;

    /*!
     * \brief Matrix and chain storage format (either `ROW` or `COLUMN`).
     */
    typedef int Storage_format;

    /*!
     * \brief A data structure storing indices of non empty chains.
     */
    typedef unspecified_type Non_zero_chain_indices;

    /*!
     * \brief Type of the chains stored in the matrix.
     */
    typedef unspecified_type Matrix_chain;

    /// @}


    /// \name Creation, Filling
    /// @{

    /*!
     * \brief Creates an empty new sparse matrix object.
     *
     * Default constructor, initialize an empty matrix of type `StorageFormat` with coefficients of type `Coefficient_ring`.
     * The default matrix size is 0x0.
     */
    SparseMatrix() ;

    /*!
     * \brief Creates a new sparse matrix object with given rows/columns sizes.
     *
     * Constructor with sizes, initialize an empty matrix of type `StorageFormat` with coefficients of type `Coefficient_ring` and a given size along rows/columns.
     */
    SparseMatrix(const size_t rowCount, const size_t columnCount) ;

    /**
     * \brief Creates a new sparse matrix from another sparse matrix object (with possibly a different `StorageFormat`).
     *
     * Copy constructor, initialize a sparse matrix of same sizes, containing the same coefficients (but not necessarly of the same `StorageFormat`).
     * If types are different, the constructor performs conversion.
     */
    SparseMatrix(const SparseMatrix& otherToCopy);

    /*!
     * \brief Assigns to other matrix.
     *
     * Assign to other matrix coefficient-wise, equivalent to copying it.
     *
     * Matrices must have the same type.
     */
    SparseMatrix& operator=(const SparseMatrix& _otherToCopy);

    /**
     * \brief Cleans a sparse matrix (set all coefficients to zero).
     *
     * Empty all structures of the sparse matrix.
     *
     */
    void nullify();

    /// @}

    /// \name Matrix informations and iterators
    /// @{

    /**
     * \brief Tests if a sparse matrix is null.
     *
     * The function returns `true` if the sparse matrix is null (that is, empty) and `false` otherwise.
    */
    bool is_null();

    /**
     * \brief Gets the matrix size.
     *
     * The matrix size as a row/column pair.
     */
    std::pair<size_t, size_t> dimensions() const;

    /**
     * \brief Iterator to the beginning of the chain indices (visited by increasing indices).
     *
     * The function returns an iterator to the beginning of the (non zero) chain indices.
     */
    inline Non_zero_chain_indices::iterator begin() const noexcept;

    /**
     * \brief Past-the-end iterator of the chain indices (visited by increasing indices).
     *
     * The function returns an iterator past the end of the chain indices.
     */
    inline Non_zero_chain_indices::iterator end() const noexcept;

    /**
     * \brief Reverse iterator to the beginning of the chain indices (visited by decreasing indices).
     *
     * The function returns a reverse iterator to the beginning of the (non zero) chain indices.
     */
    inline Non_zero_chain_indices::iterator reverse_begin() const noexcept;

    /**
     * \brief Past-the-end reverse iterator of the chain indices (visited by decreasing indices).
     *
     * The function returns a past-the-end reverse iterator of the chain indices.
     */
    inline Non_zero_chain_indices::iterator reverse_end() const noexcept;

    /// @}

    /// \name Output
    /// @{

    /**
     * \brief Inserts `matrix` in the output stream.
     */
    friend std::ostream& operator<<(std::ostream &_stream, const SparseMatrix &matrix);

    /**
     * \brief Inserts `matrix` in an output stream.
     */
    friend std::ostream& write_matrix (const SparseMatrix& matrix, std::ostream& out);

    /**
     * \brief Extracts a sparse matrix from an input stream.
     */

    template <typename _CT>
    friend std::istream& read_matrix (SparseMatrix& matrix, std::istream& in);
    /// @}

    /// \name Linear algebra operators
    /// @{

    /**
     * \brief Adds a matrix and assign.
     *
     * Adds each coefficient of the matrix together and stores the result in `this`.
     * Matrices must have the same `Coefficient_ring` but can have different `StorageFormat`.
     */
    friend SparseMatrix& operator+=(SparseMatrix &matrix, const SparseMatrix &other);

    /**
     * \brief Adds two matrices together.
     *
     * Adds each coefficient of the matrices together and returns a new matrix (of the same type as `first`) representing the result (when possible, prefer `+=` for efficiency).
     * Matrices must have the same `Coefficient_ring` but can have different `StorageFormat`.
     */
    friend SparseMatrix operator+(const SparseMatrix &first, const SparseMatrix &second);


    /**
     * \brief Subtracts a matrix and assign.
     *
     * Subtracts each coefficient of the matrix together and stores the result in `matrix`.
     * Matrices must have the same `Coefficient_ring` but can have different `StorageFormat`.
     */
    SparseMatrix& operator-=(SparseMatrix &matrix, const SparseMatrix &other);

    /**
     * \brief Subtracts two matrices together.
     *
     * Subtracts each coefficient of the matrix together and returns a new matrix (of the same type as `first`) representing the result (when possible, prefer `-=` for efficiency).
     * Matrices must have the same `Coefficient_ring` but can have different `StorageFormat`.
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
     * This method creates a new matrix obtained by multiplying the matrix by a scalar factor `lambda`. If `lambda` is zero, the function comes to nullify the matrix (when possible, prefer `*=` for efficiency).
     */
    friend SparseMatrix operator*(const Coefficient_ring& lambda, const SparseMatrix &matrix);

    /**
     * \brief Applies factor on each coefficients into a new matrix.
     *
     * This method creates a new matrix obtained by multiplying the matrix by a scalar factor `lambda`. If `lambda` is zero, the function comes to nullify the matrix (when possible, prefer `*=` for efficiency).
     */
    friend SparseMatrix operator*(const SparseMatrix &matrix, const Coefficient_ring& lambda);

    /**
     * \brief Applies factor on each coefficient and assign.
     *
     * This method multiplies the matrix by a scalar factor `lambda`.
     */
    SparseMatrix& operator*=(const Coefficient_ring& lambda);

    /**
     * \brief Multiplies a matrix and assign.
     *
     * Multiply each coefficient of the matrix together and stores the result in `matrix`.
     * Matrices must have the same `Coefficient_ring` but can have different `StorageFormat`.
     */
    friend SparseMatrix& operator*=(SparseMatrix &matrix, const SparseMatrix &other);

    /**
     * \brief Performs multiplication between matrices and returns a new column-major matrix.
     *
     * Perform standard linear algebra product between matrices and returns a new column-major matrix (when possible, prefer `*=` for efficiency). Matrices must have the same `Coefficient_ring` but can have different `StorageFormat`. Efficiency of the product depends of `StorageFormat` (when possible, prefer row-major by column-major products).
     */
    friend SparseMatrix operator*(const SparseMatrix &first, const SparseMatrix &second);

    /**
     * \brief Performs multiplication between matrices and returns a new row-major matrix.
     *
     * Perform standard linear algebra product between matrices and returns a new row-major matrix (when possible, prefer `*=` for efficiency). Matrices must have the same `Coefficient_ring` but can have different `StorageFormat`. Efficiency of the product depends of `StorageFormat`.
     */
    friend SparseMatrix operator%(const SparseMatrix &first, const SparseMatrix &second);

    /**
     * \brief Performs multiplication between a matrix and a column-based chain.
     *
     * Generates a column-based chain from the matrix multiplication and returns it.
     *
     * \pre the matrix and the chain must have the same `Coefficient_ring`.
     * \pre `column.is_column()` must be `true`
     */
    friend SparseChain operator*(const SparseMatrix& matrix, const SparseChain& column);

    /**
     * \brief Performs multiplication between a row-based chain and a matrix.
     *
     * Generates a row-based chain from the matrix multiplication and returns it.
     *
     * \pre the matrix and the chain must have the same `Coefficient_ring`.
     * \pre `row.is_row()` must be `true`
     */
    friend SparseChain operator*(const SparseChain& row, const SparseMatrix& matrix);

    /**
     * \brief Transposes a matrix.
     *
     * The function returns a new matrix where the `StorageFormat` is changed.
     */
    SparseMatrix transpose();


    /// @}

    /// \name Access and blocks operations
    /// @{

    /**
     * \brief Gets a chain from a const matrix.
     *
     * If the matrix is column-major, returns the `i`th column, and if the matrix is row-major, returns the `i`th row.
     */
    Matrix_chain operator[](size_t index) const;

    /**
     * \brief Sets a given coefficient.
     *
     * Assign the scalar `d` to the coefficient on row `i` and column `j`.
     */
    friend void set_coefficient(SparseMatrix& matrix, size_t i, size_t j, const Coefficient_ring d);

    /**
     * \brief Gets a given coefficient.
     *
     * Returns the coefficient on row `i` and column `j` of the matrix.
     */
    friend Coefficient_ring get_coefficient(const SparseMatrix& matrix, size_t i, size_t j);

    /**
     * \brief Gets the value of the column at a given `index` from the matrix (whatever the `StorageFormat` of the matrix).
     *
     * For column-matrices, it is equivalent to `operator[]`, for row-matrices a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     *
     *\returns a column chain.
     */
    friend SparseChain get_column(const SparseMatrix &matrix, size_t index);

    /**
     * \brief Gets the value of the row at a given `index` from the matrix (whatever the `StorageFormat` of the matrix).
     *
     * For row-matrices, it is equivalent to `operator[]`, for column-matrices a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     *
     *\returns a row chain.
     */
    friend SparseChain get_row(const SparseMatrix &matrix, size_t index);

    /**
     * \brief Gets a constant reference over the column of  index`i` from a column matrix.
     * \pre `matrix.is_column()` must be `true`.
     *
     * \returns a constant reference to a column chain.
     */
    const SparseChain& cget_column(const SparseMatrix& matrix, size_t i);

    /**
     * \brief Gets a constant reference over the row of  index`i` from a row matrix.
     * \pre `matrix.is_row()` must be `true`.
     *
     * \returns a constant reference to a row chain.
     */
    const SparseChain& cget_row(const SparseMatrix& matrix, size_t i);

    /**
     * \brief Sets a column in the matrix (whatever the `StorageFormat` of the matrix).
     *
     * Set the `i`th column of `matrix` to `column`.
     * For column-matrices, it should be equivalent to an assignment, however, for row-matrices, a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     * \pre `column.is_column()` must be `true`
     */
    void set_column(SparseMatrix &matrix, size_t i, const SparseChain &column);

    /**
     * \brief Sets a row in `matrix` (whatever the `StorageFormat` of the matrix).
     *
     * Set the `i`th row of `matrix` to `chain`.
     * For row-matrices, it should be equivalent to an assignment, however, for column-matrices, a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     * \pre `row.is_row()` must be `true`
     */
    void set_row(SparseMatrix &matrix, size_t i, const Sparse_chain &row);

    /**
     * \brief Gets a submatrix from the matrix.
     *
     * Nullifies the chain of index `i` along the major direction of a copy of the matrix amd returns it.
     */
    friend SparseMatrix operator/(const SparseMatrix &matrix, size_t i);

    /**
     * \brief Gets a submatrix from the matrix and assign.
     *
     * Removes (along the major dimension) all indices provided in the vector `indices` from the matrix and returns it.
     */
    SparseMatrix& operator/=(const std::vector<size_t> &indices);

    /**
     * \brief Gets a submatrix from the matrix and assign.
     *
     * Removes (along the major dimension) the chain of index `i` from the matrix and returns it.
     */
    SparseMatrix& operator/=(size_t i);

    /**
     * \brief Nullifies a column from the matrix.
     *
     * Removes column of index `i` whatever the `StorageFormat` of the matrix. For column matrices, it just comes to the `\=` operator and for row matrices, it entails a traversal of the matrix.
     */
    friend SparseMatrix& remove_column(SparseMatrix& matrix, size_t index);

    /**
     * \brief Nullifies a row from the matrix.
     *
     * Removes row of index `i` whatever the `StorageFormat` of the matrix. For row matrices, it just comes to the `\=` operator and for column matrices, it entails a traversal of the matrix.
     */
    friend SparseMatrix& remove_row(SparseMatrix& matrix, size_t index);

    /**
     * \brief Nullifies a coefficient of the matrix.
     *
     * Removes coefficient on row `i` / column `j` of the matrix.
     */
    friend SparseMatrix& remove_coefficient(SparseMatrix& matrix, size_t i, size_t j);

    /// @}

};
