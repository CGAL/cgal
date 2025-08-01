/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `SparseMatrix` describes the requirements for sparse matrices optimized for topological computations.
Traditionally, sparse matrices data structures encode non zero coefficients of (sparse) matrices in order to optimize either matrices memory footprint, or linear algebra operations (which usually comes to optimize iterators over non zero coefficients and access to coefficients). However, topological operations require slightly different features:
 - fast access to row or columns of matrices (which are actually the images under the application encoded by the matrix)
 - fast block operations (especially along row or columns)

 The `SparseMatrix` concept describes requirements for such sparse matrix. It relies on the model of `SparseChain` which encodes sparse row or column vectors. Matrices are either column major or row major (hence they either store column sparse chains or row sparse chains).

 The following constants, called `ChainTypeFlag`, encode the major direction of both sparse chains and sparse matrices.
 - `OSM::COLUMN` for column-major chains and matrices (which is the default),
 - `OSM::ROW` for row-major chains and matrices.

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
\cgalHasModelsBare{`CGAL::OSM::Sparse_matrix<Ring, ChainTypeFlag>`}
\cgalHasModelsEnd

 \sa `Ring`
 \sa `SparseChain`
*/

class SparseMatrix {

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
     * \brief A data structure storing indices of non empty chains.
     */
    typedef unspecified_type NonZeroChainIndices;

    /*!
     * \brief Type of the chains stored in the matrix.
     */
    typedef unspecified_type MatrixChain;

    /// @}


    /// \name Creation, filling
    /// @{

    /*!
     * \brief Create an empty new SparseMatrix object.
     *
     * Default constructor, initialize an empty Matrix of type `ChainTypeFlag` with coefficients of type `CoefficientType`.
     * The default matrix size is 0x0.
     */
    SparseMatrix() ;

    /*!
     * \brief Create a new SparseMatrix object with given rows/columns sizes.
     *
     * Constructor with sizes, initialize an empty Matrix of type `ChainTypeFlag` with coefficients of type `CoefficientType` and a given size along rows/columns.
     */
    SparseMatrix(const size_t rowCount, const size_t columnCount) ;

    /**
     * \brief Create a new SparseMatrix from another SparseMatrix object (with possibly a different `ChainTypeFlag`).
     *
     * Copy constructor, initialize a SparseMatrix of same sizes, containing the same coefficients (but not necessarly of the same `ChainTypeFlag`).
     * If types are different, the constructor performs conversion.
     */
    SparseMatrix(const SparseMatrix& otherToCopy);

    /*!
     * \brief Assign to other matrix.
     *
     * Assign to other matrix coefficient-wise, equivalent to copying it.
     *
     * Matrices must have the same type.
     */
    SparseMatrix& operator=(const SparseMatrix& _otherToCopy);

    /**
     * \brief Clean a SparseMatrix (set all coefficients to zero).
     *
     * Empty all structures of the sparse matrix.
     *
     */
    void nullify();

    /// @}

    /// \name Matrix informations and iterators
    /// @{

    /**
     * \brief Test is a SparseMatrix is null.
     *
     * The function return `true` is the SparseMatrix is null (that is, empty) and `false` otherwise.
    */
    bool is_null();

    /**
     * \brief Gets the matrix size.
     *
     * The matrix size as a row/column pair.
     */
    std::pair<size_t, size_t> dimensions() const;

    /**
     * \brief Iterator to the beginning of the chains indices (visited by increasing indices).
     *
     * The function returns an iterator to the beginning of the (non zero) chains indices.
     */
    inline NonZeroChainIndices::iterator begin() const noexcept;

    /**
     * \brief Iterator to the ending of the chains indices (visited by increasing indices).
     *
     * The function returns an iterator to the ending of the chains indices.
     */
    inline NonZeroChainIndices::iterator end() const noexcept;

    /**
     * \brief Reverse iterator to the beginning of the chains indices (visited by decreasing indices).
     *
     * The function returns a reverse iterator to the beginning of the (non zero) chains indices.
     */
    inline NonZeroChainIndices::iterator reverse_begin() const noexcept;

    /**
     * \brief Reverse iterator to the ending of the chains indices (visited by decreasing indices).
     *
     * The function returns a reverse iterator to the ending of the chains indices.
     */
    inline NonZeroChainIndices::iterator reverse_end() const noexcept;

    /// @}

    /// \name Output
    /// @{

    /**
     * \brief Displays a SparseMatrix in the output stream.
     */
    friend std::ostream& operator<<(std::ostream &_stream, const SparseMatrix &_matrix);

    /**
     * \brief Write a sparse matrix to an output stream.
     */
    friend std::ostream& write_matrix (const Sparse_matrix& M, std::ostream& out);

    /**
     * \brief Read a sparse matrix from an input stream.
     */

    template <typename _CT>
    friend std::istream& read_matrix (Sparse_matrix& M, std::istream& in);
    /// @}

    /// \name Linear algebra operators
    /// @{

    /**
     * \brief Add a matrix and assign.
     *
     * Adds each coefficient of the matrix together and stores the result in `this`.
     * Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`.
     */
    friend SparseMatrix& operator+=(SparseMatrix &matrix, const SparseMatrix &other);

    /**
     * \brief Adds two matrices together.
     *
     * Adds each coefficient of the matrices together and returns a new matrix (of the same type as `first`) representing the result (when possible, prefer `+=` for efficiency).
     * Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`.
     */
    friend SparseMatrix operator+(const SparseMatrix &first, const SparseMatrix &second);


    /**
     * \brief Substract a matrix and assign.
     *
     * Substracts each coefficient of the matrix together and stores the result in `matrix`.
     * Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`.
     */
    SparseMatrix& operator-=(SparseMatrix &matrix, const SparseMatrix &other);

    /**
     * \brief Substracts two matrices together.
     *
     * Substracts each coefficient of the matrix together and returns a new matrix (of the same type as `first`) representing the result (when possible, prefer `-=` for efficiency).
     * Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`.
     */
    friend SparseMatrix operator-(const SparseMatrix &first, const SparseMatrix &second);

    /**
     * \brief Compute the negative of a matrix (unary operator).
     *
     * \return The resulting matrix.
     */
    friend SparseMatrix operator-(const SparseMatrix& matrix);

    /**
     * \brief Apply factor on each coefficients into a new matrix.
     *
     * This method creates a new matrix obtained by multiplying the matrix by a scalar factor `lambda`. If `lambda` is zero, the function comes to nullify the matrix (when possible, prefer `*=` for efficiency).
     */
    friend Sparse_matrix operator*(const CoefficientType& lambda, const Sparse_matrix &matrix);

    /**
     * \brief Apply factor on each coefficients into a new matrix.
     *
     * This method creates a new matrix obtained by multiplying the matrix by a scalar factor `lambda`. If `lambda` is zero, the function comes to nullify the matrix (when possible, prefer `*=` for efficiency).
     */
    friend Sparse_matrix operator*(const Sparse_matrix &matrix, const CoefficientType& lambda);

    /**
     * \brief Apply factor on each coefficient and assign.
     *
     * This method multiplies the matrix by a scalar factor `lambda`.
     */
    SparseMatrix& operator*=(const CoefficientType& lambda);

    /**
     * \brief Multiply a matrix and assign.
     *
     * Multiply each coefficient of the matrix together and stores the result in `matrix`.
     * Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`.
     */
    friend SparseMatrix& operator*=(SparseMatrix &matrix, const SparseMatrix &other);

    /**
     * \brief Perform multiplication between matrices and returns a new column-major matrix.
     *
     * Perform standard linear algebra product between matrices and returns a new column-major matrix (when possible, prefer `*=` for efficiency). Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`. Efficiency of the product depends of `ChainTypeFlag` (when possible, prefer row-major by column-major products).
     */
    friend SparseMatrix operator*(const SparseMatrix &first, const SparseMatrix &second);

    /**
     * \brief Perform multiplication between matrices and returns a new row-major matrix.
     *
     * Perform standard linear algebra product between matrices and returns a new row-major matrix (when possible, prefer `*=` for efficiency). Matrices must have the same `CoefficientType` but can have different `ChainTypeFlag`. Efficiency of the product depends of `ChainTypeFlag`.
     */
    friend SparseMatrix operator%(const SparseMatrix &first, const SparseMatrix &second);

    /**
     * \brief Perform multiplication between a matrix and a column-based chain.
     *
     * Generates a column-based chain from the matrix multiplication and returns it.
     *
     * The matrix and the chain must have the same coefficient type.
     */
    friend Sparse_chain<CoefficientType, COLUMN> operator*(const SparseMatrix& first, const Sparse_chain<CoefficientType, COLUMN>& second);

    /**
     * \brief Perform multiplication between a row-based chain and a matrix.
     *
     * Generates a row-based chain from the matrix multiplication and returns it.
     *
     * The matrix and the chain must have the same coefficient type.
     */
    friend Sparse_chain<CoefficientType, ROW> operator*(const Sparse_chain<CoefficientType, ROW>& first, const SparseMatrix& second);

    /**
     * \brief Transpose a matrix.
     *
     * The function returns a new matrix where the chain type flag is changed.
     */
    SparseMatrix transpose();


    /// @}

    /// \name Access and blocks operations
    /// @{

    /**
     * \brief Get a chain from a const matrix.
     *
     * If the matrix is column-major, returns the `i`th column, and if the matrix is row-major, returns the `i`th row.
     */
    MatrixChain operator[](size_t index) const;

    /**
     * \brief Set a given coefficient.
     *
     * Assign the scalar `d` to the coefficient on row `i` and column `j`.
     */
    friend void set_coef(Sparse_matrix& matrix, size_t i, size_t j, const CoefficientType d);

    /**
     * \brief Get a given coefficient.
     *
     * Returns the coefficient on row `i` and column `j` of the matrix.
     */
    friend CoefficientType get_coef(const Sparse_matrix& matrix, size_t i, size_t j);

    /**
     * \brief Get the value of the column at a given `index` from the matrix (whatever the `ChainTypeFlag` of the matrix).
     *
     * For column-matrices, it is equivalent to `operator[]`, for row-matrices a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     */
    friend Sparse_chain<CoefficientType, COLUMN> get_column(const SparseMatrix &matrix, size_t index);

    /**
     * \brief Get the value of the row at a given `index` from the matrix (whatever the `ChainTypeFlag` of the matrix).
     *
     * For row-matrices, it is equivalent to `operator[]`, for column-matrices a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     */
    friend Sparse_chain<CoefficientType, ROW> get_row(const SparseMatrix &matrix, size_t index);

    /**
     * \brief Get a constant reference over the column of  index`i` from a column matrix.
     */
    const Sparse_chain<CoefficientType, COLUMN> & cget_column(const SparseMatrix<CoefficientType, COLUMN> &_matrix, size_t i);

    /**
     * \brief Get a constant reference over the row of  index`i` from a row matrix.
     */
    const Sparse_chain<CoefficientType, ROW> & cget_row(const SparseMatrix<CoefficientType, ROW> &_matrix, size_t i);

    /**
     * \brief Set a column in the matrix (whatever the `ChainTypeFlag` of the matrix).
     *
     * Set the `i`th column of `matrix` to `chain`.
     * For column-matrices, it should be equivalent to an affectation, however, for row-matrices, a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     */
    void set_column(SparseMatrix &matrix, size_t i, const Sparse_chain<CoefficientType, COLUMN> &chain);

    /**
     * \brief Set a row in the matrix (whatever the `ChainTypeFlag` of the matrix).
     *
     * Set the `i`th row of `matrix` to `chain`.
     * For row-matrices, it should be equivalent to an affectation, however, for column-matrices, a traversal of the matrix is required (in \f$\mathcal O(n)\f$).
     */
    void set_row(SparseMatrix &matrix, size_t i, const Sparse_chain<CoefficientType, ROW> &chain);

    /**
     * \brief Get a submatrix from the matrix.
     *
     * Nullifies the chain of index `i` along the major direction of a copy of the matrix returns it.
     */
    friend SparseMatrix operator/(const SparseMatrix &_matrix, size_t i);

    /**
     * \brief Get a submatrix from the matrix and assign.
     *
     * Removes (along the major dimension) all indexes provided in the vector `indexes` from the matrix and returns it.
     */
    SparseMatrix& operator/=(const std::vector<size_t> &indexes);

    /**
     * \brief Get a submatrix from the matrix and assign.
     *
     * Removes (along the major dimension) the chain of index `i` from the matrix and returns it.
     */
    SparseMatrix& operator/=(size_t i);

    /**
     * \brief Nullifies a column from the matrix.
     *
     * Removes column of index `i` whatever the `ChainTypeFlag` of the matrix. For column matrices, it just comes to the `\=` operator and for row matrices, it entails a traversal of the matrix.
     */
    friend Sparse_matrix& del_column(Sparse_matrix& matrix, size_t index);

    /**
     * \brief Nullifies a row from the matrix.
     *
     * Removes row of index `i` whatever the `ChainTypeFlag` of the matrix. For row matrices, it just comes to the `\=` operator and for column matrices, it entails a traversal of the matrix.
     */
    friend Sparse_matrix& del_row(Sparse_matrix& matrix, size_t index);

    /**
     * \brief Nullifies a coefficient of the matrix.
     *
     * Removes coefficient on row `i` / column `j` of the matrix.
     */
    friend Sparse_matrix& del_coef(Sparse_matrix& matrix, size_t i, size_t j);

    /// @}

};
