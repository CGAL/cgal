/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `Filtration` describes the requirements for persistent filtrations associated to persistent homology computation.

A filtration is associated to an `AbstractChainComplex`. Each cell is equiped with a scalar value (called its degree) and the filtration is an enumeration of cells in any dimension by increasing degrees.

A filtration class provides:

 - an iterator to visit all cells by increasing degrees.
 - getters to get the degree of a cell, and the cell of a given index along the filtration.
 - an overload of `operator<<()` to output filtrations.

Cells are indexed along each dimension and thus identified by their index together with their dimension.

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Filtration_core<ChainComplex,Degree>`}
\cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Filtration_lower_star<ChainComplex,Degree>`}
\cgalHasModelsEnd

*/

class Filtration
{
public:
    /*! \brief Type of the underlying complex. */
    typedef unspecified_type Complex;

    /*! \brief (Scalar) type of degrees (model of `RealEmbeddable`).
     */
    typedef unspecified_type Degree ;

    /*! \brief Type of coefficients used to compute homology. */
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

protected:
    /*!
     Type of column-major sparse matrices.
     */
    typedef SparseMatrix<Coefficient_ring,CGAL::OSM::COLUMN> Column_matrix ;

    /*!
     Type of row-major sparse matrices.
     */
    typedef SparseMatrix<Coefficient_ring,CGAL::OSM::ROW> Row_matrix ;

    /*!
     Type of column-major chains.
     */
    typedef SparseChain<Coefficient_ring,CGAL::OSM::COLUMN> Column_chain ;

    /*!
     Type of row-major chains.
     */
    typedef SparseChain<Coefficient_ring,CGAL::OSM::ROW> Row_chain ;

public:
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

        /*! Iterator constructor
         */
        iterator(const Filtration& f, size_t i=0);

        /*! Dereference operator.
         */
        value_type operator*() const;

        /*! Gets the index (time) associated to the iterator.
         */
        size_t time () const ;

        /*! Gets the `Cell` (cell and its dimension) associated to the iterator.
         */
        Cell cell_dim () const ;

        /*! Gets the degree associated to the iterator.
         */
        Degree degree () const ;

        /*!
         * \brief Prefix incrementation.
         * Move to next index in the filtration.
         */
        iterator& operator++();

        /*! \brief Postfix incrementation.
         */
        iterator operator++(int);

        /*!
         * \brief Equality check.
         * \returns True if the indices are equal.
         */
        friend bool operator== (const iterator& a, const iterator& b);

        /*!
         * \brief Inequality check.
         * \returns True if the indices are different.
         */
        friend bool operator!= (const iterator& a, const iterator& b);
    };

    /*!
     * \brief Iterator to the beginning of the filtration.
     */
    iterator begin();

    /*!
     * \brief Returns a past-the-end iterator.
     */
    iterator end();

    // getters
    /*! \brief Gets the filtration size.
     */
    size_t size () const;

    /*! \brief Gets the cell (that is cell index and dimension) at the index `i` of the filtration.
     */
    Cell cell_index_dimension (size_t i) const;

    /*! \brief Gets the degree of the `i`th element of the filtration.
     */
    Degree degree (size_t i) const;

    // Filtration verification
    /*! \brief Checks that a filtration is valid.
     * Checks that cells are ordered in increasing degrees and all cells have indices larger than their faces.
     */
    bool is_valid() const;

    // Input/output filtration


    /**
     * \brief Exports the filtration time indices.
     *
     * The method exports the time index of every cell in each dimension.
     */
    vector<vector<size_t> > export_filtration () const
};

    /*!
    \relates Filtration
    \brief extracts a filtration from a stream.
     */
    istream & operator>>(istream & in, Filtration &f);

    /*!
    \relates Filtration
    \brief Inserts a filtration into a stream.
     */
    ostream & operator<<(ostream & out, const Filtration &f);