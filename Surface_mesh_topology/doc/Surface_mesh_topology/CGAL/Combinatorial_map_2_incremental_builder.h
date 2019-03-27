
namespace CGAL {
  /*!
    \ingroup PkgSurfaceMeshTopologyClasses

    The class `Combinatorial_map_2_incremental_builder` is a tool enabling to build incrementally a 2D combinatorial map.

    \tparam CMap a model of `CombinatorialMap`
  */
  template<typename CMap>
  class Combinatorial_map_2_incremental_builder
  {
  public:
    /// Dart_handle type.
    typedef typename CMap::Dart_handle Dart_handle;
    
    /*! Constructor. Creates a Combinatorial_map_2_incremental_builder object, to create amap.
     */
    Combinatorial_map_2_incremental_builder(CMap& amap);
    
    /// Start a new surface
    void begin_surface();
      
    /// End of the surface. Return one dart of the created surface.
    /// @pre A surface is under creation.
    Dart_handle end_surface();

    /// Start a new facet.
    void begin_facet();
    
    /// End of the facet. Return the first dart of this facet.
    /// @pre A facet is under creation.
    Dart_handle end_facet();

    /// Add one edge to the current facet, given by its label (any string, using minus sign for orientation)
    /// @pre A facet is under creation.
    void add_edge_to_facet(const std::string& s);
    
    /// Add the given edges to the current facet
    /// s is a sequence of labels, add all the corresponding edges into the current facet.
    /// @pre A facet is under creation.
    void add_edges_to_facet(const std::string& s);
    
    /// Add directly one facet giving the sequence of labels 's' of all its edges.
    /// @pre A surface is under creation.
    void add_facet(const std::string& s);
      
   /// Start a path on the surface
   void begin_path();
    
    /// End the current path.
    /// @return The path created.
    /// @pre A path is under creation.
    CGAL::Path_on_surface<CMap> end_path();

    /// Add edge labeled e at the end of the current path
    /// @pre A path is under creation.
    void add_edge_to_path(const std::string& e);
    
    /// A shortcut enabling to create a path directly with a sequence of edge labels
    /// @return The path created.
    CGAL::Path_on_surface<CMap> create_path(const std::string& s);
    
  };
}
