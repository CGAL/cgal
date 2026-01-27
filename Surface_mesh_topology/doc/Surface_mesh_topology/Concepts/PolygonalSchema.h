/*!
  \ingroup PkgSurfaceMeshTopologyConcepts
  \cgalConcept

  The concept `PolygonalSchema` defines a 2D polygonal schema, i.e. a combinatorial surface with labeled edges. A PolygonalSchema is created incrementally by adding facets one at a time. A label is any word, that does not contain a space.

  PolygonalSchema::Dart_info should be a class having a public data member std::string m_label.
  PolygonalSchema::dimension should be equal to 2.

  \cgalRefines{GenericMap}

  \cgalHasModelsBegin
  \cgalHasModelsBare{\link CGAL::Surface_mesh_topology::Polygonal_schema_with_combinatorial_map `CGAL::Surface_mesh_topology::Polygonal_schema_with_combinatorial_map<Items,Alloc>`\endlink}
  \cgalHasModelsBare{\link CGAL::Surface_mesh_topology::Polygonal_schema_with_generalized_map `CGAL::Surface_mesh_topology::Polygonal_schema_with_generalized_map<Items,Alloc>`\endlink}
  \cgalHasModelsEnd
  */

class PolygonalSchema
{
public:
  /// creates an empty `PolygonalSchema` object.
  PolygonalSchema();

  /// starts a new facet.
  void init_facet();

  /// finishes the current facet. Returns the first dart of this facet.
  /// @pre A facet is under creation.
  Dart_descriptor finish_facet();

  /// adds one edge to the current facet, given by its label `l` (any string containing no space, using minus sign for orientation).
  /// Since the surface is oriented, each label can be used only twice with opposite signs. If this method is called with a label already used, with same sign, an error message is given and this label is ignored.
  /// @pre A facet is under creation.
  void add_edge_to_facet(const std::string& l);

  /// adds the given edges to the current facet.
  /// `s` is a sequence of labels, separated by spaces. All the corresponding edges are added into the current facet.
  /// @pre A facet is under creation.
  void add_edges_to_facet(const std::string& s);

  /// adds directly one facet giving the sequence of labels `s` of all its edges (labels are separated by spaces).
  void add_facet(const std::string& s);

  /// returns the label of dart `d`.
  std::string get_label(Dart_descriptor d) const;

  /// returns dart with label `s`, NULL if this label is not used.
  Dart_descriptor get_dart_labeled(const std::string & s) const;

  /// returns true iff the facet containing `d` is perforated.
  bool is_perforated(Dart_const_descriptor d) const;

  /// Shortcut for `is_perforated(get_dart_labeled(s))`.
  bool is_perforated(const std::string & s) const;

  /// perforates the facet containing `d`. Returns the number of darts of the face; 0 if the facet was already perforated.
  size_type perforate_facet(Dart_descriptor d);

  /// Shortcut for perforate_facet(get_dart_labeled(s)).
  size_type perforate_facet(const std::string & s);

  /// fills the facet containing `d`. Returns the number of darts of the face; 0 if the facet was already filled.
  size_type fill_facet(Dart_descriptor d);

  /// Shortcut for `fill_facet(get_dart_labeled(s))`.
  size_type fill_facet(const std::string & s);
};
