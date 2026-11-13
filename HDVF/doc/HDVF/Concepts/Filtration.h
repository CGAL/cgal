/*!
 \ingroup PkgHDVFConcepts
 \cgalConcept

 The concept `Filtration` describes the requirements for persistent filtrations associated to persistent homology computation.

 By definition, a filtration over a chain complex `K` is a sequence of sub complexes \f$\{K_t\,;\, t\in[d_1, \ldots, d_N]\}\f$, where \f$d_i\f$ are scalars, such that:

 - \f$\forall d\leqslant d'\f$, \f$K_d\subseteq K_{d'}\f$.
 - \f$|K_{d_{i+1}}| = |K_{d_{i}}|+1\f$.
 - As a consequence, any cell \f$\sigma\in K\f$, is "added" at a given *time* \f$i\f$ and the corresponding \f$d_i\f$ is called its *degree* (denoted by \f$\mathrm{deg}(\sigma)\f$).
 - \f$K_{d_1}\f$ is an empty complex.
 - \f$K_{d_N} = K\f$.
 - For any cell \f$\sigma\in K\f$, if \f$\tau\f$ is a proper face of \f$\sigma\f$, then \f$\mathrm{deg}(\tau)\leqslant \mathrm{deg}(\sigma)\f$.

 Hence a filtration is associated to an `AbstractChainComplex`. Each cell is equipped with a scalar value (called its degree) and the filtration is an enumeration of cells in any dimension by increasing degrees.

 In consequence, a filtration associated to a complex `K`:
 - maps each cell to a scalar called its degree
 - orders the cells of `K` (each cell has a unique index along the filtration called *time* along the filtration) such that:
 - the index of a cell is larger that the indices of its faces,
 - the degree map is increasing along the filtration.

 Let us remind that in all `AbstractChainComplex`, cells are indexed along each dimension and thus identified by their index together with their dimension.

 \cgalHasModelsBegin
 \cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Filtration_core<ChainComplex,Degree>`}
 \cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Filtration_lower_star<ChainComplex,Degree>`}
 \cgalHasModelsEnd

 */

class Filtration
{
public:

    /// \name Types
    /// @{

    /*! \brief Type of the underlying complex (a model of `AbstractChainComplex`). */
    typedef unspecified_type Chain_complex;

    /*! \brief (Scalar) type of degrees (model of `RealEmbeddable`).
     */
    typedef unspecified_type Degree ;

    /*! \brief Type of coefficients used to compute homology . */
    typedef Complex::Coefficient_ring Coefficient_ring;


    /*! \brief Type for indexing uniquely a cell.
     *
     * As stated in `AbstractChainComplex`, cells are identified by their dimension, together with their index along this dimension. The type `Cell` stores this pair:
     *
     * - First element of the pair: index of the cell.
     * - Second element of the pair: dimension of the cell.
     */
    typedef std::pair<size_t, int> Cell ;

    /*! \brief Value type of the filtration iterator.
     * Contains a cell (identified by its index and dimension in a `Cell`) and its associated degree.
     */
    typedef struct {
        Cell cell_dim ;
        Degree degree ;
    } Filtration_iter_value ;

    /// @}

public:

    /// \name Iterator
    /// @{
    /**
     * \brief The iterator over filtration.
     *
     * Iterate the filtration by increasing degrees.
     */
    struct iterator
    {
        // Iterator tags
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = Filtration_iter_value;

        /*! \brief Iterator constructor
         *
         * \param f Constant reference on the `Filtration_core` iterated.
         * \param i The initial index.
         */
        iterator(const Filtration& f, size_t i=0);

        /*! \brief Iterator dereference
         *
         * \returns A `Filtration_iter_value` structure containing the information of the current cell and its degree.
         */
        value_type operator*() const;

        /*!
         * \brief Prefix incrementation.
         *
         * Moves to next cell in the filtration.
         *
         * \returns The reference to the current iterator.
         */
        iterator& operator++();

        /*! \brief Postfix incrementation. Moves to next cell in the filtration.
         * \returns The pre-incremented iterator.
         */
        iterator operator++(int);

        /*!
         * \brief Equality check.
         * \returns `true` if the indices are equal.
         */
        friend bool operator== (const iterator& a, const iterator& b);

        /*!
         * \brief Inequality check.
         * \returns `true` if the indices are different.
         */
        friend bool operator!= (const iterator& a, const iterator& b);
    };

    /*!
     * \brief Iterator to the beginning of the filtration.
     *
     * \return The iterator to the beginning of the filtration.
     */
    iterator begin();

    /**
     * \brief Iterator past-the-end of the filtration.
     *
     * \return The iterator past-the-end of the filtration.
     */
    iterator end();

    /// @}

    /// \name Getters and tests
    /// @{

    /*! \brief Gets the filtration size.
     */
    size_t size () const;

    /*! \brief Gets the cell (that is cell index and dimension) at the index `i` of the filtration.
     */
    Cell cell (size_t i) const;

    /*! \brief Gets the degree of the `i`th element of the filtration.
     */
    Degree degree (size_t i) const;


    /*! \brief Checks that a filtration is valid.
     * Checks that cells are ordered in increasing degrees and all cells have indices larger than their faces.
     */
    bool is_valid() const;

    /** \brief Comparison operator. */
    bool operator==(const Filtration_core& f) const;


    /// @}

    // Input/output filtration

    /// \name Input, output
    /// @{

    /**
     * \brief Exports the filtration time indices.
     *
     * The method exports the time index of every cells in each dimension.
     *
     * \returns A vector containing, for each dimension, the vector of labels by cell index.
     */
    vector<vector<size_t> > export_filtration () const;

    /*!
     \brief Extracts a filtration from a stream.

     The filtration is encoded using the simple %FILTR file format (see secHDVFfiltr_file_format "user manual" for details).
     */
    friend std::istream & operator>>(std::istream & in, Filtration &f);

    /*!
     \brief Inserts a filtration into a stream.

     The filtration is encoded using the simple %FILTR file format (see secHDVFfiltr_file_format "user manual" for details).
     */
    friend std::ostream & operator<<(std::ostream & out, const Filtration &f);
    /// @}
};


