/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `Filtration` describes the requirements for persistent filtrations associated to persistent homology computation.

A filtration is associated to an `AbstractChainComplex`. Each cells is equiped with a scalar value (called its degree) and the filtration is an enumeration of cells in any dimension by increasing degrees.

A filtration class provides:

 - an iterator to visit all cells by increasing degrees.
 - getters to get the degree of a cell, and the cell of a given index along the filtration.
 - an overload of `<<` to output filtrations.

Cells are indexed along each dimension and thus identified by their index together with their dimension.

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::HDVF::Filtration_core<Degree_type>`}
\cgalHasModelsBare{`CGAL::HDVF::Filtration_lower_star<Degree_type>`}
\cgalHasModelsEnd

*/

class Filtration
{
public:
    /*! \brief (Scalar) type of degrees.
     */
    typedef unspecified_type Degree_type ;

    /*! \brief Type for indexing uniquely a cell.
     *
     * As stated in `AbstractChainComplex`, cells are identified by their dimension, together with their index along this dimension. The type `Cell_index_dimension` stores this pair:
     *
     * - First element of the pair: index of the cell.
     * - Second element of the pair: dimension of the cell.
     */
    typedef std::pair<size_t, int> Cell_index_dimension ;

    /*! \brief Value type of the filtration iterator.
     * Contains a cell (identified by its index and dimension in a `Cell_index_dimension`) and its associated degree.
     */
    typedef struct {
        Cell_index_dimension cell_dim ;
        Degree_type degree ;
    } Filtration_iter_value ;

protected:
    /*!
     Type of column-major sparse matrices.
     */
    typedef CGAL::OSM::Sparse_matrix<_CoefType,CGAL::OSM::COLUMN> Column_matrix ;

    /*!
     Type of row-major sparse matrices.
     */
    typedef CGAL::OSM::Sparse_matrix<_CoefType,CGAL::OSM::ROW> Row_matrix ;

    /*!
     Type of column-major chains.
     */
    typedef CGAL::OSM::Sparse_chain<_CoefType,CGAL::OSM::COLUMN> Column_chain ;

    /*!
     Type of row-major chains.
     */
    typedef CGAL::OSM::Sparse_chain<_CoefType,CGAL::OSM::ROW> Row_chain ;

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

        /*! Gets the `Cell_index_dimension` (cell and its dimension) associated to the iterator.
         */
        Cell_index_dimension cell_dim () const ;

        /*! Gets the degree associated to the iterator.
         */
        Degree_type degree () const ;

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
    Cell_index_dimension cell_index_dimension (size_t i) const;

    /*! \brief Gets the degree of the `i`th element of the filtration.
     */
    Degree_type degree (size_t i) const;

    // Filtration verification
    /*! \brief Checks that a filtration is valid
     * Checks that cells are ordered in increasing degrees and all cells have indices larger than their faces.
     */
    bool is_valid() const;

    // Input/output filtration
    /*! \brief Overload of the `>>`operator for filtrations.
     */
    friend istream & operator>>(istream & in, Filtration &f);
    
    /*! \brief Overload of the `operator<<()`operator for filtrations.
     */
    friend ostream & operator<<(ostream & out, const Filtration &f);

    /**
     * \brief Exports the filtration time indices.
     *
     * The method exports the time index of every cells in each dimension.
     */
    vector<vector<size_t> > export_filtration () const
};
