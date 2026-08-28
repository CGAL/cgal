/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `HDVF` describes the requirements for Homological Discrete Vector Fields (%HDVF for short) , a theory of computational homology unifying discrete Morse theory and effective homology. HDVFs were introduced by Aldo Gonzalez-Lorenzo in his PhD (see [AGL,2017], [AGL,2016]).

A HDVF is a combinatorial tool computing the homology of a chain complex (described by the `AbstractChainComplex` concept) including "geometric" information such as homology and cohomology generators.

 In order to build a perfect HDVF, a pairing operation called `A()` associates (under conditions) two critical cells (becoming PRIMARY and SECONDARY) and updates the reduction in time \f$\mathscr O(n^2)\f$. Intuitively, this operation bears some similitudes (see users guide) with discrete Morse theory's DGVF arrows (and an HDVF can always be represented by such a vector field - possibly with cycles).
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
\cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Hdvf_core<IntegralDomainWithoutDivision, AbstractChainComplex, SparseChain, Sparse_matrix>`}
\cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Hdvf<IntegralDomainWithoutDivision, AbstractChainComplex>`}
\cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Hdvf_persistence<IntegralDomainWithoutDivision, AbstractChainComplex>`}
\cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Hdvf_duality<IntegralDomainWithoutDivision, AbstractChainComplex>`}
\cgalHasModelsEnd

 \sa `AbstractChainComplex`
 \sa `IntegralDomainWithoutDivision`
 \sa `SparseChain`
 \sa `SparseMatrix`

 *How to describe constants declared in the namespace Homological_discrete_vector_field and used everywhere? PSC_flag, options, exporttype*

 [AGL, 2017] Aldo Gonzalez-Lorenzo, Alexandra Bac, Jean-Luc Mari, Pedro Real. Allowing cycles in discrete Morse theory, Topology and its Applications, Volume 228, 2017, Pages 1-35.

 [AGL, 2016] Aldo Gonzalez-Lorenzo. Computational Homology Applied to Discrete Objects. Discrete Mathematics [cs.DM]. Aix-Marseille Université; Universidad de Sevilla, 2016.

*/

class HDVF
{
public:
/// \name Types
/// @{

/*!
 * \brief Type of underlying chain complex (a model of `AbstractChainComplex`).
 */
typedef unspecified_type Complex;

/*! \brief Type of coefficients used to compute homology. */
typedef Complex::Coefficient_ring Coefficient_ring;

    /** \brief HDVF option (compute only reduced boundary). */
    const int OPT_BND = 0b0001;
    /** \brief HDVF option (compute only reduced boundary and f). */
    const int OPT_F = 0b0010;
    /** \brief HDVF option (compute only reduced boundary and g). */
    const int OPT_G = 0b0100;
    /** \brief HDVF option (compute full reduction). */
    const int OPT_FULL = 0b1000;

/*!
 * \brief Type of sparse chains (a model of `SparseChain`).
 */
    typedef unspecified_type Sparse_chain;

/*!
 * \brief Type of sparse matrices (a model of `SparseMatrix`).
 */
    typedef unspecified_type Sparse_matrix;

/*! \brief Structure to represent data for HDVF operations (pairs of cells).
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
typedef Sparse_matrix<CoefficientType, CGAL::OSM::COLUMN> Column_matrix;

/*!
 Type of row-major sparse matrices
 */
typedef Sparse_matrix<CoefficientType, CGAL::OSM::ROW> Row_matrix;

/// @}

/// \name Functions to find valid pairs for `A()` operations
/// @{

/*!
 * \brief Finds a valid Cell_pair of dimension q / q+1 for A.
 *
 * The function searches a pair of critical cells \f$(\gamma_1, \gamma2)\f$ of dimension q / q+1, valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma_2), \gamma_1 \rangle\f$ invertible). It returns the first valid pair found by iterators.
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
 * \param q Dimension of the cell `gamma`.
 * \param gamma Index of a cell to pair.
 * \param found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
 */
virtual Cell_pair find_pair_A(int q, size_t gamma, bool &found) const;

/**
 * \brief Finds *all* valid Cell_pair of dimension q / q+1 for A.
 *
 * The function searches all pairs of critical cells \f$(\gamma_1, \gamma_2)\f$ of dimension q / q+1, valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma_2), \gamma_1 \rangle\f$ invertible).
 * It returns a vector of such pairs.
 *
 * \param q Lower dimension of the pair.
 * \param found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
 */
virtual std::vector<Cell_pair> find_pairs_A(int q, bool &found) const;

/**
 * \brief Finds *all* valid Cell_pair for A containing `gamma` (a cell of dimension `q`)
 *
 * The function searches all critical cells \f$\gamma'\f$ such that one of the following conditions holds:
 * - \f$\gamma'\f$ has dimension q+1 and \f$(\gamma, \gamma')\f$ is valid for A (ie.\ such that \f$\langle \partial_{q+1}(\gamma'), \gamma \rangle\f$ invertible),
 * - \f$\gamma'\f$ has dimension q-1 and \f$(\gamma', \gamma)\f$ is valid for A (ie.\ such that \f$\langle \partial_{q}(\gamma), \gamma' \rangle\f$ invertible).
 * It returns a vector of such pairs.
 *
 * \param q Dimension of the cell `gamma`.
 * \param gamma Index of a cell to pair.
 * \param found Reference to a %Boolean variable. The method sets `found` to `true` if a valid pair is found, `false` otherwise.
 */
virtual std::vector<Cell_pair> find_pairs_A(int q, size_t gamma, bool &found) const;


/// @}


/// \name HDVF functions and operations
/// @{

/*!
 * \brief A operation (pairs two critical cells).
 *
 * A pair of critical cells \f$(\gamma_1, \gamma_2)\f$ of respective dimension q and q+1 is valid for A if \f$\langle \partial_{q+1}(\gamma_2), \gamma_1 \rangle\f$ is invertible. After the `A()` operation, \f$\gamma_1\f$ becomes PRIMARY, \f$\gamma_2\f$ becomes SECONDARY. The `A()` functions updates the reduction accordingly (in time \f$\mathscr O(n^2)\f$).
 */
void A(size_t gamma1, size_t gamma2, int q);

/*!
 * \brief Computes a perfect HDVF.
 *
 * As long as valid pairs for A exist, the function selects the first available pair (returned by `find_pair_A()`) and applies the corresponding `A()` operation.
 * If the `IntegralDomainWithoutDivision` of coefficients is a field, this operation always produces a perfect HDVF (ie.\ the reduced boundary is null and the reduction provides homology and cohomology information).
 *
 * If the HDVF is initially not trivial (some cells have already been paired), the function completes it into a perfect HDVF.
 *
 * \param verbose If `true`, all intermediate reductions are printed out.
 *
 */
std::vector<Cell_pair> compute_perfect_hdvf(bool verbose = false);

/*!
 * \brief Computes a random perfect HDVF.
 *
 *As long as valid pairs for A exist, the function selects a random pair (among pairs returned by `find_pairs_A()`) and applies the corresponding `A()` operation.
 * If the `IntegralDomainWithoutDivision` of coefficients is a field, this operation always produces a perfect HDVF (ie.\ the reduced boundary is null and the reduction provides homology and cohomology information).
 *
 * If the HDVF is initially not trivial (some cells have already been paired), the function randomly completes it into a perfect HDVF.
 *
 * \param verbose If `true`, all intermediate reductions are printed out.
 */
std::vector<Cell_pair> compute_rand_perfect_hdvf(bool verbose = false);
/// @}

/// \name Getters
/// @{

/*!
 * \brief Gets cells with a given `PSC_flag` in any dimension.
 *
 * The function returns in each dimension the vector of cells with a given `PSC_flag`.
 */
std::vector<std::vector<size_t> > flag (PSC_flag flag) const;

/*!
 * \brief Gets cells with a given `PSC_flag` in dimension `q`.
 *
 * The function returns the vector of cells of dimension `q` with a given `PSC_flag`.
 */
std::vector<size_t> psc_flags (PSC_flag flag, int q) const;

/*!
 * \brief Gets the `PSC_flag` of the cell `tau` in dimension `q`.
 */
PSC_flag psc_flag (size_t tau, int q) const;

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

/*!
 * \brief Tests if a HDVF is perfect.
 *
 * The function returns `true` if the reduced boundary matrix is null and `false` otherwise.
 */
    bool is_perfect_hdvf();

/// @}

/// \name Output functions
/// @{

/*!
 * \brief Writes the matrices of the reduction.
 *
 * Writes the matrices of the reduction (that is \f$f\f$, \f$g\f$, \f$h\f$, \f$\partial'\f$ the reduced boundary).
 *
 * By default, writes the complex to `std::cout`.
*/
std::ostream& write_matrices(std::ostream &out = std::cout) const;

/*! \brief Writes the homology and cohomology reduction information.
 *
 * Writes the homology and cohomology reduction information (that is \f$f^*\f$, \f$g\f$ \f$\partial'\f$ the reduced boundary over each critical cell).
 *
 * By default, writes the complex to `std::cout`.
*/
std::ostream& write_reduction(std::ostream &out = std::cout) const;

/*!
 *\brief Writes a HDVF and its reduction to a stream.
 *
 * Writes them to a stream  in  `.hdvf` file format (see xxx for a specification).
 */

std::ostream& write_hdvf_reduction(std::ostream& out) ;

/*!
 * \brief Reads a HDVF and its reduction from a stream.
 *
 * Reads them from a `.hdvf` file format (see xxx for a specification).
 */
std::istream& read_hdvf_reduction(std::istream& in) ;

/*!
 * \brief Exports primary/secondary/critical labels (in particular for vtk export).
 *
 * The method exports the labels of every cell in each dimension.
 *
 * \return A vector containing, for each dimension, the vector of labels by cell index.
 */
std::vector<std::vector<int> > psc_labels () const;

/*!
 * \brief Returns a chain containing the homology generator associated to `cell_index` (critical cell) of dimension  `q` (in particular for vtk export).
 *
 * The method exports the chain \f$g(\sigma)\f$ for \f$\sigma\f$ the cell of index `cell_index` and dimension `q`.
 *
 * \return A column-major chain.
 */
Column_chain homology_chain (size_t cell_index, int q) const;

/*!
 * \brief Returns a chain containing the cohomology generators associated to `cell_index` (critical cell) of dimension  `q` (in particular for vtk export).
 *
 * The method exports:
 * - the chain \f$f^\star(\sigma)\f$ for \f$\sigma\f$ the cell of index `cell_index` and dimension `q`,
 * - or the co-faces of this chain if `co_faces` is true (sometimes easier to view cohomology generators)
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

};
