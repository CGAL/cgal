namespace CGAL {
namespace Surface_mesh_topology {

/*!
  \ingroup PkgSurfaceMeshTopologyClasses
  
  The class `Shortest_noncontractible_cycle` provides methods to 
  
  \tparam Mesh_ a model of `GenericMap` or of `FaceGraph`
  \tparam Weight_ a model of `WeightFunctor`
*/

template <class Mesh_, class Weight_ = void>
class Shortest_noncontractible_cycle {

public:

  /// Default weight functor is used in the unweighted case.
  /// Every edge has weight of 1.
  struct Default_weight_functor {
    /// Number type of the weight
    using Weight_t = unsigned int;

    /// Return 1 for any edge's weight
    template <class T>
    Weight_t operator()(T);
  };

  /// The weight functor used.
  /// If the template parameter `Weight_` is void, set as `Default_weight_functor`; otherwise set as `Weight_`
  using Weight = unspecified_type;
  
  /// Number type of the weights.
  using Distance_type = Weight::Weight_t;

  /// Dart_handle of the input mesh.
  /// If the template parameter `Mesh_` is a model of `GenericMap`, set as \link GenericMap::Dart_handle `Mesh::Dart_handle`\endlink ; otherwise set as `boost::graph_traits<Mesh_>::halfedge_descriptor`
  using Dart_handle_orig = unspecified_type;

  /// The cycle type.
  using Path = std::vector<Dart_handle_orig>;
  
  /// Constructor
  Shortest_noncontractible_cycle(Mesh_& mesh, const Weight& wf = Weight());
  
  /// Find the shortest non-contractible cycle through a vertex.
  /// The vertex is depicted by `root_vertex`, a \link GenericMap::Dart_handle `Dart_handle` \endlink of the input mesh. The cycle found is returned in `cycle`. The total length of the cycle is calculated by adding up the weights of all edges in the cycle and is returned through the pointer `length`.
  bool find_cycle(Dart_handle_orig root_vertex, Path& cycle, Distance_type* length = NULL);

  /// Find the edge width of the whole mesh
  /// The edge width is returned in `cycle`. The total length of the cycle is calculated by adding up the weights of all edges in the cycle and is returned through the pointer `length`.
  void edge_width(Path& cycle, Distance_type* length = NULL);

};

}
}
