/*!
 \ingroup PkgHDVFConcepts
 \cgalConcept

 The concept `GeometricChainComplex` refines the concept `AbstractChainComplex` and describes the requirements for (topological) chain complexes associated to geometric complexes used in the concept `HDVF`.
 It adds to `AbstractChainComplex` methods to get vertex coordinates.

 \cgalRefines{AbstractChainComplex}

 \cgalHasModelsBegin
 \cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Simplicial_chain_complex<CoefficientRing,Traits>`}
 \cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Cubical_chain_complex<CoefficientRing,Traits>`}
 \cgalHasModelsEnd

 */

class GeometricChainComplex
{
public:
    /// \name Types
    /// @{

    /*!
     Type of points.
     */
    typedef unspecified_type Point ;
    /// @}


    /// \name Access functions
    /// @{

    /*!
     \brief Returns the vector of vertex coordinates.
     */
    const std::vector<Point>& points() const;

    /*!
     \brief Returns the coordinates of the cell of index `i` and dimension 0.
     */
    Point point(size_t i) const;

    /// @}

    /// \name VTK export
    /// @{

    /**
     * \brief Exports a column chain `chain` over a `GeometricChainComplex` to a %VTK file.
     *
     * Exports cells of `chain` with a non null coefficient to a %VTK file, with a distinguished cell (for instance the underlying critical cell to a generator). Two properties are exported:
     * - the `Label` property (2 for the distinguished cell, 0 for other cells)
     * - the `CellId` property containing the "basis" index of the cell.
     *
     * \param K Underlying `GeometricChainComplex`.
     * \param filename Output file name.
     * \param chain Column chain exported to %VTK.
     * \param q Dimension of the chain.
     * \param cellId Index of the distinguished cell.
     *
     * \exception File_not_found If the file `filename` cannot be created and opened, throw a `%std::runtime_error` exception.
     */
    friend void write_VTK(const GeometricChainComplex &K, const std::string &filename, const Column_chain& chain, int q, size_t cellId);

    /**
     * \brief Exports a `GeometricChainComplex`  (plus, optionally, labels) to a %VTK file.
     *
     * The method generates legacy text %VTK files. Labels are exported as such in a %VTK property, together with CellID property, containing the index of each cell.
     *
     * \tparam CoefficientRing Coefficient ring of the underlying complex `K`.
     * \tparam Traits Traits class of the underlying complex `K`.
     * \tparam LabelType Type of labels provided (default: int).
     *
     * \param K `GeometricChainComplex` exported.
     * \param filename Output file root (output filenames will be built from this root).
     * \param labels Pointer to a vector of labels in each dimension. `labels->at(q)` is the set of integer labels of cells of dimension q. If labels is NULL, only CellID property is exported.
     * \param label_type_name Typename used in %VTK export (e.g. "int" or "unsigned_long", see <a href = "https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html">VTK manual </a>).
     *
     * \exception File_not_found If the file `filename` cannot be created and opened, throw a `%std::runtime_error` exception.
     */
    template <typename LabelType = int>
    void write_VTK(const GeometricChainComplex &K, const std::string &filename, const std::vector<std::vector<LabelType> > *labels=NULL, std::string label_type_name = "int");

    /// @}

};
