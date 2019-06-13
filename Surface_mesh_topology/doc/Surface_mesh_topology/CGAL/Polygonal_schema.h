
namespace CGAL {
namespace Surface_mesh_topology {
  /*!
    \ingroup PkgSurfaceMeshTopologyClasses

    The class `Polygonal_schema` is a model of `CombinatorialMap` with labeled edges. A  Polygonal_schema is created incrementally by adding facets one at a time.

    \tparam CMap a model of `CombinatorialMap`
  */
  template<typename CMap>
  class Polygonal_schema
  {
  public:
    /// Dart_handle type.
    typedef typename CMap::Dart_handle Dart_handle;
    
    /*! creates an empty `Polygonal_schema` object.
     */
    Polygonal_schema();
    
    /// starts a new surface
    void init_surface();
      
    /// finishes the current surface. Returns one dart of the created surface.
    /// @pre A surface is under creation.
    Dart_handle finish_surface();

    /// starts a new facet.
    void init_facet();
    
    /// finishes the current facet. Returns the first dart of this facet.
    /// @pre A facet is under creation.
    Dart_handle finish_facet();

    /// adds one edge to the current facet, given by its label `l` (any string containing no space, using minus sign for orientation).
    /// Since the surface is oriented, each label can be used only twice with opposite signs. If this method is called with a label already used, with same sign, an error message is given and this label is ignored. 
    /// @pre A facet is under creation.
    void add_edge_to_facet(const std::string& l);
    
    /// adds the given edges to the current facet.
    /// `s` is a sequence of labels, separated by spaces. All the corresponding edges are added into the current facet.
    /// @pre A facet is under creation.
    void add_edges_to_facet(const std::string& s);
    
    /// adds directly one facet giving the sequence of labels `s` of all its edges (labels are separated by spaces).
    /// @pre A surface is under creation.
    void add_facet(const std::string& s);
  };
}
}
