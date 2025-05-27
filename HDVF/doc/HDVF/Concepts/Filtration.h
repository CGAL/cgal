/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `Filtration` describes the requirements for persistent filtrations associated to persistent homology computation.
 
A filtration is associated to an `AbstractSimplicialComplex`. Each cells is equiped with a scalar value (called its degree) and the filtration is an enumeration of cells in any dimension by increasing degrees.
 
A filtration class provides:
 
 - an iterator to visit all cells by increasing degrees.
 - getters to get the degree of a cell, and the cell of a given index along the filtration.
 - an overload of `<<` to output filtrations.
 
Cells are indexed along each dimension and thus identified by their index together with their dimension.
 
\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Filtration_core<DegreeType>`}
\cgalHasModelsBare{`CGAL::Filtration_lower_star<DegreeType>`}
\cgalHasModelsEnd

*/

class Filtration
{
public:
    /*! \brief (Scalar) type of degrees.
     */
    typedef unspecified_type DegreeType ;
    
    /*! \brief Type for indexing uniquely a cell.
     * - First element of the pair: index of the cell.
     * - Second element of the pair: dimension of the cell.
     */
    typedef std::pair<int, int> CellDim ;
    
    /*! \brief Value returned by the filtration iterator.
     * Contains a cell (identified by its index and dimenion in a `CellDim`) and its associated degree.
     */
    typedef struct {
        CellDim cell_dim ;
        _DegType degree ;
    } FiltrationIterValue ;
    
protected:
    /*!
     Type of column-major sparse matrices.
     */
    typedef OSM::Sparse_matrix<_CoefType,OSM::COLUMN> CMatrix ;
    
    /*!
     Type of row-major sparse matrices.
     */
    typedef OSM::Sparse_matrix<_CoefType,OSM::ROW> RMatrix ;
    
    /*!
     Type of column-major chains.
     */
    typedef OSM::Sparse_chain<_CoefType,OSM::COLUMN> CChain ;
    
    /*!
     Type of row-major chains.
     */
    typedef OSM::Sparse_chain<_CoefType,OSM::ROW> RChain ;
    
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
        using value_type        = FiltrationIterValue;
        
        /*! Iterator constructor
         */
        iterator(const Filtration& f, int i=0);
        
        /*! Dereference operator.
         */
        value_type operator*() const;
        
        /*! Get the index (time) associated to the iterator.
         */
        int time () const ;
        
        /*! Get the `CellDim` (cell and its dimension) associated to the iterator.
         */
        CellDim cell_dim () const ;
        
        /*! Get the degree associated to the iterator.
         */
        DegreeType degree () const ;
        
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
     * \brief Iterator to the ending of the filtration.
     */
    iterator end();
    
    // getters
    /*! \brief Get the filtration size.
     */
    int get_filtration_size () const;
    
    /*! \brief Get the cell (that is cell index and dimension) at the index `i` of the filtration.
     */
    CellDim get_cell_dim (int i) const;
    
    /*! \brief Get the degree of the `i`th element of the filtration.
     */
    DegreeType get_degree (int i) const;
    
    // Filtration verification
    /*! \brief Check that a filtration is valid
     * Checks that cells are ordered in increasing degrees and all cells have indices larger than their faces.
     */
    bool is_valid_filtration() const;
    
    // Output filtration
    /*! \brief Overload of the `<<`operator for filtrations.
     */
    friend ostream & operator<<(ostream & out, const Filtration &f);
    
    /**
     * \brief Export the filtration time indices.
     *
     * The method exports the time index of every cells in each dimension.
     */
    vector<vector<int> > export_filtration () const
};
