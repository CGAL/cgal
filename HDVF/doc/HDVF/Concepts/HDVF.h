/*!
 \ingroup PkgHDVFConcepts
 \cgalConcept

 The concept `HDVF` describes the requirements for Homological Discrete Vector Fields (%HDVF for short) , a theory of computational homology unifying discrete Morse theory and effective homology. HDVFs were introduced by Aldo Gonzalez-Lorenzo in his PhD (see [AGL,2017], [AGL,2016]).

 A HDVF is a combinatorial tool computing the homology of a chain complex (described by the `AbstractChainComplex` concept) including "geometric" information such as homology and cohomology generators.

 In order to build a perfect HDVF, a pairing operation called `A()` associates (under conditions) two critical cells (becoming PRIMARY and SECONDARY) and updates the reduction in time \f$\mathscr O(n^2)\f$. Intuitively, this operation bears some similitudes  with discrete Morse theory's DGVF arrows (see \ref secHDVFdiscrete_Morse "introduction to discrete Morse theory" in user manual). Actually, an HDVF can always be represented by a vector field - possibly with cycles.
 The *validity condition for A* is tested by various `find_pair_A()` and `find_pairs_A()` methods, which search for such valid pairs.


 Starting from this operation, `HDVF` provides two methods to generate perfect HDVFs (and hence compute homology):
 - `compute_perfect_hdvf()` which computes a perfect HDVF by chosing at each step the first A-pairings available (depending on cells ordering in the chain complex)
 - `compute_rand_perfect_hdvf()` which computes a perfect random HDVF by chosing at each step a rand A-pairing among all possible ones  (this option is thus slower)

 For efficiency, the HDVF reduction can be computed partially and we define constants to define the HDVF computation option required:
 - compute only the reduced boundary  (option `OPT_BND`): to get Betti numbers
 - compute \f$\partial'\f$ together with \f$h\f$ and \f$g\f$ (option `OPT_G`): to get homology generators
 - compute \f$\partial'\f$ together with \f$h\f$ and \f$f\f$ (option `OPT_F`): to get co-homology generators
 - compute the full reduction (option `OPT_FULL`): to get the full homological / cohomological information


 \cgalHasModelsBegin
 \cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Hdvf_core<IntegralDomainWithoutDivision, AbstractChainComplex, SparseChain, Sparse_matrix_structure>`}
 \cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Hdvf<IntegralDomainWithoutDivision, AbstractChainComplex>`}
 \cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Hdvf_persistence<IntegralDomainWithoutDivision, AbstractChainComplex>`}
 \cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Hdvf_duality<IntegralDomainWithoutDivision, AbstractChainComplex>`}
 \cgalHasModelsEnd

 \sa `AbstractChainComplex`
 \sa `IntegralDomainWithoutDivision`
 \sa `SparseChain`
 \sa `SparseMatrix`

 \warning How to describe constants declared in the namespace Homological_discrete_vector_field and used everywhere? PSC_flag, options, exporttype*

 \warning How mention more clearly the export functions from Geometric_chain_complex_tools which provide write_VTK functions ?

 [AGL, 2017] Aldo Gonzalez-Lorenzo, Alexandra Bac, Jean-Luc Mari, Pedro Real. Allowing cycles in discrete Morse theory, Topology and its Applications, Volume 228, 2017, Pages 1-35.

 [AGL, 2016] Aldo Gonzalez-Lorenzo. Computational Homology Applied to Discrete Objects. Discrete Mathematics [cs.DM]. Aix-Marseille Université; Universidad de Sevilla, 2016.

 */

class HDVF
{
public:
    /// \name Constants
    /// @{

    /** \brief HDVF computation option: compute only the reduced boundary. */
    const int OPT_BND = 0b0001;
    /** \brief HDVF computation option: compute only the reduced boundary and \f$f\f$. */
    const int OPT_F = 0b0010;
    /** \brief HDVF computation option: compute only the reduced boundary and \f$g\f$. */
    const int OPT_G = 0b0100;
    /** \brief HDVF computation option: compute the full reduction. */
    const int OPT_FULL = 0b1000;

    /** \brief HDVF Enum for the label of cells.
     *
     * Values can be `PRIMARY`, `SECONDARY`, `CRITICAL` or `NONE` (for discarded cells in persistent homology or duality).
     */
    enum PSC_flag {
        PRIMARY,
        SECONDARY,
        CRITICAL,
        NONE
    };

    /// @}

    /// \name Types
    /// @{

    /*!
     * \brief Type of underlying chain complex (a model of `AbstractChainComplex`).
     */
    typedef unspecified_type Complex;

    /*! \brief Type of coefficients used to compute homology (a model of `IntegralDomainWithoutDivision`). */
    typedef Complex::Coefficient_ring Coefficient_ring;

    /*!
     * \brief Type of sparse chains (a model of `SparseChain`).
     */
    typedef unspecified_type Sparse_chain;

    /*!
     * \brief Type of sparse matrices (a model of `SparseMatrix`).
     */
    typedef unspecified_type Sparse_matrix;


    /*! \brief Structure to represent data for HDVF operations (between pairs of cells).
     *
     * Cells are always sorted so that the dimension of `sigma` is lesser than the dimension of `tau`.
     */
    struct Cell_pair {
        /// Index of the first cell
        size_t sigma;
        /// Index of the second cell
        size_t tau;
        /// Lower dimension of the cells
        int dim;
        /// Comparison operator
        inline bool operator==(const Cell_pair& pair2);
        /// Output operator
        friend std::ostream& operator<<(std::ostream &out, const std::vector<Cell_pair>& pairs);
    };

    /*!
     Type of column-major chains
     */
    typedef Sparse_chain<CoefficientType, CGAL::OSM::COLUMN> Column_chain;

    /*!
     Type of row-major chains
     */
    typedef Sparse_chain<CoefficientType, CGAL::OSM::ROW> Row_chain ;

    /*!
     Type of column-major sparse matrices
     */
    typedef Sparse_matrix<CoefficientType, CGAL::OSM::COLUMN, Sparse_chain> Column_matrix;

    /*!
     Type of row-major sparse matrices
     */
    typedef Sparse_matrix<CoefficientType, CGAL::OSM::ROW, Sparse_chain> Row_matrix;

    /// @}

    /// \name Functions to find valid pairs for A() operations
    /// @{

    /*!
     * \brief Finds a valid Cell_pair of dimension q / q+1 for A.
     *
     * The function searches a pair of critical cells \f$(\gamma_1, \gamma2)\f$ of dimension q / q+1, valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma_2), \gamma_1 \rangle\f$ invertible). It returns the first valid pair found by iterators.
     *
     * If the dimension is incorrect (negative or larger than the dimension of the underlying complex), `found` is set to `false`.
     *
     * \param q Lower dimension of the pair.
     * \param found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    virtual Cell_pair find_pair_A(int q, bool &found) const;

    /*!
     * \brief Finds a valid Cell_pair for A containing `gamma` (a cell of dimension `q`)
     *
     * The function searches a cell \f$\gamma'\f$ such that one of the following conditions holds:
     * - \f$\gamma'\f$ has dimension q+1 and \f$(\gamma, \gamma')\f$ is valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma'), \gamma \rangle\f$ invertible),
     * - \f$\gamma'\f$ has dimension q-1 and \f$(\gamma', \gamma)\f$ is valid for A (ie.\ such that \f$\langle \partial_{q}(\gamma), \gamma' \rangle\f$ invertible).
     *
     * If the dimension is incorrect (negative or larger than the dimension of the underlying complex), `found` is set to `false`. If the cell index is incorrect (larger than the number of cells in the dimension), `found` is set to `false`.
     *
     * \param q Dimension of the cell `gamma`.
     * \param gamma Index of a cell to pair.
     * \param found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    virtual Cell_pair find_pair_A(int q, size_t gamma, bool &found) const;

    /**
     * \brief Finds *all* valid Cell_pair of dimension q / q+1 for A.
     *
     * The function searches all pairs of critical cells \f$(\gamma_1, \gamma2)\f$ of dimension q / q+1, valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma_2), \gamma_1 \rangle\f$ invertible).
     * It returns a vector of such pairs.
     *
     * If the dimension is incorrect (negative or larger than the dimension of the underlying complex), `found` is set to `false` and the function returns an empty vector.
     *
     * \param q Lower dimension of the pair.
     * \param found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    virtual std::vector<Cell_pair> find_pairs_A(int q, bool &found) const;

    /**
     * \brief Finds *all* valid Cell_pair for A containing `gamma` (a cell of dimension `q`)
     *
     * The function searches all `CRITICAL` cells \f$\gamma'\f$ such that one of the following conditions holds:
     * - \f$\gamma'\f$ has dimension q+1 and \f$(\gamma, \gamma')\f$ is valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma'), \gamma \rangle\f$ invertible),
     * - \f$\gamma'\f$ has dimension q-1 and \f$(\gamma', \gamma)\f$ is valid for A (ie.\ such that \f$\langle \partial_{q}(\gamma), \gamma' \rangle\f$ invertible).
     * It returns a vector of such pairs.
     *
     * If the dimension is incorrect (negative or larger than the dimension of the underlying complex), `found` is set to `false`. If the cell index is incorrect (larger than the number of cells in the dimension), `found` is set to `false`. In both cases, the function returns an empty vector.
     *
     * \param q Dimension of the cell `gamma`.
     * \param gamma Index of a cell to pair.
     * \param found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
     */
    virtual std::vector<Cell_pair> find_pairs_A(int q, size_t gamma, bool &found) const;

    /** \brief Checks if the pair of cells \f$(\gamma, \gamma')\f$, of dimensions q / q+1, is valid for A.
     *
     * The pair is valid if \f$\langle \mathrm d(\gamma),\gamma' \rangle\f$ is invertible.
     *
     * \param gamma Index of the first cell (dimension `q`).
     * \param gamma_prime Index of the second cell (dimension `q+1`).
     * \param q Lower dimension of cells.
     */

    bool is_valid_pair_for_A(size_t gamma, size_t gamma_prime, int q);


    /// @}


    /// \name HDVF functions and operations
    /// @{

    /*!
     * \brief A operation: pairs critical cells.
     *
     * A pair of critical cells \f$(\gamma_1, \gamma_2)\f$ of respective dimension q and q+1 is valid for A if \f$\langle \partial_{q+1}(\gamma_2), \gamma_1 \rangle\f$ is invertible. After the `A()` operation, \f$\gamma_1\f$ becomes `PRIMARY`, \f$\gamma_2\f$ becomes `SECONDARY`. The A method updates the reduction accordingly (in time \f$\mathcal O(n^2)\f$).
     *
     * \param gamma1 First cell of the pair (dimension `q`)
     * \param gamma2 Second cell of the pair (dimension `q+1`)
     * \param q Dimension of the pair
     */
    void A(size_t gamma1, size_t gamma2, int q);

    /**
     * \brief A operation (pairs critical cells).
     *
     * A pair of critical cells \f$(\gamma_1, \gamma_2)\f$ of respective dimension q and q+1 is valid for A if \f$\langle \partial_{q+1}(\gamma_2), \gamma_1 \rangle\f$ is invertible. After the `A()` operation, \f$\gamma_1\f$ becomes `PRIMARY`, \f$\gamma_2\f$ becomes `SECONDARY`. The A method updates the reduction accordingly (in time \f$\mathcal O(n^2)\f$).
     *
     * \param p Argument `Cell_pair` containing \f$\gamma_1\f$, \f$\gamma_2\f$ and the dimension of the pair.
     */
    void A(const Cell_pair& p);

    /*!
     * \brief Computes a perfect HDVF.
     *
     * As long as valid pairs for A exist, the function selects the first available pair (returned by `find_pair_A`()) and applies the corresponding `A()` operation.
     * If the `IntegralDomainWithoutDivision` of coefficients is a field, this operation always produces a perfect HDVF (ie.\ the reduced boundary is null and the reduction provides homology and cohomology information).
     * Otherwise the operation produces a maximal HDVF with a residual boundary matrix over critical cells.
     *
     * If the HDVF is initially not trivial (some cells have already been paired), the function completes it into a perfect HDVF.
     *
     * \param verbose If this parameter is `true`, all intermediate reductions are printed out.
     *
     * \return The vector of all `Cell_pair` paired with A.
     *
     */
    std::vector<Cell_pair> compute_perfect_hdvf(bool verbose = false);

    /*!
     * \brief Computes a random perfect HDVF.
     *
     * As long as valid pairs for A exist, the function selects a random pair (among pairs returned by `find_pairs_A()`) and applies the corresponding `A()` operation.
     * If the `IntegralDomainWithoutDivision` of coefficients is a field, this operation always produces a perfect HDVF (that  is the reduced boundary is null and the reduction provides homology and cohomology information).
     *
     * If the HDVF is initially not trivial (some cells have already been paired), the function randomly completes it into a perfect HDVF.
     *
     * \warning This method is slower that `compute_perfect_hdvf()` (finding out all possible valid pairs requires additional time).
     *
     * \param verbose If this  parameter is `true`, all intermediate reductions are printed out.
     *
     * \return The vector of all pairs of cells used for apply A.
     */
    std::vector<Cell_pair> compute_rand_perfect_hdvf(bool verbose = false);

    /*!
     * \brief Tests if a HDVF is perfect.
     *
     * The function returns `true` if the reduced boundary matrix is null and `false` otherwise.
     */
    bool is_perfect_hdvf();

    /** \brief Comparison operator (with a HDVF over the same underlying complex).
     *
     * Compares the HDVF with `other` HDVF (fast comparison, testing only HDVF flags).
     *
     * \param other Other HDVF to compare.
     */
    bool operator==(const Hdvf_core& other);

    /** \brief Compares the HDVF with another HDVF over the same underlying complex.
     *
     * \param other Other HDVF to compare.
     * \param full_compare Turns on "in depth" HDVF comparison (reduction matrices).
     */

    bool compare(const Hdvf_core& other, bool full_compare = false) const;

    /// @}

    /// \name Getters
    /// @{

    /*!
     * \brief Gets cells with a given `PSC_flag` in any dimension.
     *
     * The function returns a vector containing, for each dimension, the vector of cells with a given `PSC_flag`. Indices are sorted by increasing value.
     *
     * \param flag PSC_flag to select.
     */
    std::vector<std::vector<size_t> > psc_flag (PSC_flag flag) const;

    /*!
     * \brief Gets cells with a given `PSC_flag` in dimension `q`.
     *
     * The function returns the vector of cells of dimension `q` with a given `PSC_flag`. Indices are sorted by increasing value.
     *
     * \param flag PSC_flag to select.
     * \param q Dimension visited.
     */
    std::vector<size_t> psc_flag (PSC_flag flag, int q) const;

    /*!
     * \brief Gets the PSC_flag of the cell `tau` in dimension `q`.
     *
     * \param tau Index of the cell.
     * \param q Dimension of the cell.
     */
    PSC_flag psc_flag (size_t tau, int q) const;

    /**
     * \brief Gets the`PSC_flag` of all cells in dimension `q`.
     *
     * The function returns the vector of `PSC_flag` of cells of dimension `q`.
     *
     * \param q Dimension visited.
     */
    virtual const std::vector<PSC_flag>& psc_flags (int q) const;

    /**
     * \brief Gets the`PSC_flag` of all cells.
     *
     * The function returns a vector which `q`th element contains the vector of `PSC_flag` of cells of dimension `q`.
     */
    virtual const std::vector<std::vector<PSC_flag> >& psc_flags () const;

    /** \brief Returns the number of cells with a given flag in dimension `q`.
     *
     * \param flag `PSC_flag` tested.
     * \param q Dimension along which the number of cells is returned.
     */
    size_t number_of_cells_by_flag(PSC_flag flag, int q);

    /*!
     * \brief Gets HDVF option.
     */
    int hdvf_opts ();

    /*!
     * \brief Gets the row-major matrix of \f$f\f$ (from the reduction associated to the HDVF).
     */
    const Row_matrix& matrix_f (int q) const;

    /*!
     * \brief Gets the column-major matrix of \f$g\f$ (from the reduction associated to the HDVF).
     */
    const Column_matrix& matrix_g (int q) const;

    /*!
     * \brief Gets the column-major matrix of \f$h\f$ (from the reduction associated to the HDVF).
     */
    const Column_matrix& matrix_h (int q) const;

    /*!
     * \brief Gets the column-major matrix of \f$\partial'\f$, reduced boundary operator (from the reduction associated to the HDVF).
     */
    const Column_matrix& matrix_dd (int q) const;

    /**
     \brief Gets a constant reference over the underlying chain complex.
     */
    const Chain_complex& complex () const;

    /*!
     * \brief Returns a chain containing the homology generator associated to `cell_index` (critical cell) of dimension  `q`.
     *
     * The method exports the chain \f$g(\sigma)\f$ for \f$\sigma\f$ the cell of index `cell_index` and dimension `q`.
     *
     * The complexity is in the number of (non zero) coefficients of the generator.
     *
     * \return A column-major chain.
     */
    Column_chain homology_chain (size_t cell_index, int q) const;

    /*!
     * \brief Returns a chain containing the cohomology generators associated to `cell_index` (critical cell) of dimension  `q`.
     *
     * The method exports:
     * - the chain \f$f^\star(\sigma)\f$ for \f$\sigma\f$ the cell of index `cell_index` and dimension `q`,
     * - or the co-faces of this chain if `co_faces` is true (sometimes easier to visualize cohomology generators)
     *
     * The complexity is in the number of (non zero) coefficients of the generator.
     *
     * Below, a sample mesh with, (left) homology generators, (right) two examples of cohomology generators (corresponding generators/co-generators bear similar colours):
     *
     * <img src="HDVF_dtorus_homs.png" align="center" width=25%/>
     * <img src="HDVF_dtorus_cohom1.png" align="center" width=25%/>
     * <img src="HDVF_dtorus_cohom2.png" align="center" width=25%/>
     *
     * The same generators displayed through their co-faces:
     *
     * <img src="HDVF_dtorus_cohom1_co.png" align="center" width=25%/>
     * <img src="HDVF_dtorus_cohom2_co.png" align="center" width=25%/>
     *
     * All homology / cohomology generators:
     *
     *<img src="HDVF_dtorus_all.png" align="center" width=30%/>
     *
     * \return A column-major chain.
     */
    virtual Column_chain cohomology_chain (size_t cell_index, int dim, bool co_faces = false) const;

    /// @}

    /// \name Dimension restriction functions
    /// @{

    /** \brief Returns the dimension of the Hdvf computation.
     *
     * If HDVF computation is full (computation of homology in any dimension), returns -1; otherwise, returns the dimension along which the HDVF is computed.
     */
    int dimension_restriction() const;

    /** \brief Changes the dimension of Hdvf computation.
     *
     * \param dimension Dimension along which the computation is restricted.  Set -1 for full computation (computation of homology in any dimension), and a positive value lower than the complex dimension for a computation restricted to dimension `dimension`.
     */
    inline void dimension_restriction(int dimension);

    /// @}

    /// \name Input, output
    /// @{

    /*!
     * \brief Writes the matrices of the reduction.
     *
     * Writes the matrices of the reduction (that is \f$f\f$, \f$g\f$, \f$h\f$, \f$\partial'\f$ the reduced boundary).
     *
     * If `dimension_restriction` is -1, writes matrices in any dimensions. If `dimension_restriction` is positive, restricts output to `dimension_restriction`. If `dimension_restriction` is -2 (default value), uses the HDVF dimension restriction value.
     *
     * \param out Output stream (by default, writes the complex to `std::cout`).
     * \param dimension_restriction If positive, restricts the output to a single dimension, if -1, outputs matrices in any dimensions, if -2 (default value) uses the HDVF dimension restriction value.
     */
    std::ostream& write_matrices(std::ostream &out = std::cout, int dimension_restriction = -2) const;

    /*! \brief Writes the homology and cohomology reduction information.
     *
     * Writes \f$f^*\f$, \f$g\f$ \f$\partial'\f$ the reduced boundary over each critical cell.
     *
     * If `dimension_restriction` is -1, writes the reduction in any dimensions. If `dimension_restriction` is positive, restricts output to `dimension_restriction`. If `dimension_restriction` is -2 (default value), uses the HDVF dimension restriction value.
     *
     * \param out Output stream (by default, writes the complex to `std::cout`).
     * \param dimension_restriction If positive, restricts the output to a single dimension, if -1, outputs reductions in any dimensions, if -2 (default value) uses the HDVF dimension restriction value.
     */
    std::ostream& write_reduction(std::ostream &out = std::cout, int dimension_restriction = -2) const;

    /**
     * \brief Writes the PSC flags of cells along computed dimensions.
     *
     * If `dimension_restriction` is -1, writes flags in any dimensions. If `dimension_restriction` is positive, restricts output to `dimension_restriction`. If `dimension_restriction` is -2 (default value), uses the HDVF dimension restriction value.
     *
     * \param out Output stream (by default, writes the complex to `std::cout`).
     * \param dimension_restriction If positive, restricts the output to a single dimension, if -1, outputs reductions in any dimensions, if -2 (default value) uses the HDVF dimension restriction value.
     */
    std::ostream& write_flags(std::ostream &out = std::cout, int dimension_restriction = -2) const;

    /*! \brief Writes a HDVF together with the associated reduction (f, g, h, d matrices)
     *
     * Exports a HDVF to a stream in %HDVF file format (a simple text file format, see %HDVF_file_format "user manual" for a specification).
     *
     * \param out Output stream.
    */
   std::ostream& write_hdvf_reduction(std::ostream& out) const ;

    /** \brief Writes a HDVF together with the associated reduction to a file (f, g, h, d matrices).
     *
     * Exports a HDVF to a file in `hdvf` file format (a simple text file format, see %HDVF_file_format "user manual" for a specification).
     *
     * \param filename Output file name.
     *
     * \exception File_not_found If the file `filename` cannot be opened, raises a `%std::runtime_error`.
     */
    void write_hdvf_reduction(std::string filename) const;

    /*!
     * \brief Reads a HDVF together with the associated reduction (f, g, h, d matrices)
     *
     * Loads a HDVF and its reduction from a stream in `hdvf` file format, a simple text file format (see %HDVF_file_format "user manual" for a specification).
     * \warning The underlying complex is not stored in the file!
     *
     * \param in_stream Input stream.
     *
     * \exception Incompatible_prefix If the prefix of the stored HDVF is not 'HDVF_core' raises a `%std::runtime_error`.
     *
     * \exception Incompatible_type If the type of the stored HDVF is not 0 (HDVF and reduction storage) raises a `%std::runtime_error`.
     *
     * \exception Incompatible_dimension If the dimension of the HDVF stored does not match the dimension of the underlying complex, raises a `%std::runtime_error`.
     *
     * \exception Incoherent_number_of_cells If the number of cells loaded in a given dimension does not match that of the underlying complex, raises a `%std::runtime_error`.
     */
    std::istream& read_hdvf_reduction(std::istream& in_stream) ;

    /**
     * \brief Reads a HDVF together with the associated reduction from a file (f, g, h, d matrices)
     *
     * Loads a HDVF and its reduction from a file in `hdvf` file format, a simple text file format (see %HDVF_file_format "user manual" for a specification).
     * \warning The underlying complex is not stored in the file!
     *
     * \param filename Input file name.
     *
     * \exception File_not_found If the file `filename` cannot be opened, raises a `%std::runtime_error`.
     */
    void read_hdvf_reduction(std::string filename);

    /*!
     * \brief Exports an integer code for cells labels (in particular for vtk export).
     *
     * The method exports integers encoding cells labels for every cell in each dimension. Labels are encoded as follows:
     * - `PRIMARY` is encoded as -1
     * - `SECONDARY` is encoded as 1
     * - `CRITICAL` is encoded as 0
     * - `NONE` is encoded as -2
     *
     * \return A vector containing, for each dimension, the vector of encoded labels by cell index.
     */
    std::vector<std::vector<int> > psc_labels () const;

    /// @}

};
