
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
    
    /*! creates a `Combinatorial_map_2_incremental_builder` object, to create amap.
     */
    Combinatorial_map_2_incremental_builder(CMap& amap);
    
    /// starts a new surface
    void begin_surface();
      
    /// finishes the current surface. Returns one dart of the created surface.
    /// @pre A surface is under creation.
    Dart_handle end_surface();

    /// starts a new facet.
    void begin_facet();
    
    /// finishes the current facet. Returns the first dart of this facet.
    /// @pre A facet is under creation.
    Dart_handle end_facet();

    /// adds one edge to the current facet, given by its label `l` (any string containing no space, using minus sign for orientation)
    /// @pre A facet is under creation.
    void add_edge_to_facet(const std::string& l);
    
    /// adds the given edges to the current facet.
    /// `s` is a sequence of labels, separated by spaces. All the corresponding edges are added into the current facet.
    /// @pre A facet is under creation.
    void add_edges_to_facet(const std::string& s);
    
    /// adds directly one facet giving the sequence of labels `s` of all its edges (labels are separated by spaces).
    /// @pre A surface is under creation.
    void add_facet(const std::string& s);
      
   /// starts a path on the surface.
   void begin_path();
    
    /// finishes the current path. Returns the path created.
    /// @pre A path is under creation.
    CGAL::Path_on_surface<CMap> end_path();

    /// adds edge labeled `l` at the end of the current path.
    /// @pre A path is under creation.
    void add_edge_to_path(const std::string& l);
    
    /// create a path directly from a sequence of edge labels `s` (labels are separated by spaces). Returns the path created.
    CGAL::Path_on_surface<CMap> create_path(const std::string& s);
  };
}
