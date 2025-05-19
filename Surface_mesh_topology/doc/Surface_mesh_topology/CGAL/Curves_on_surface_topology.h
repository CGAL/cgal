
namespace CGAL {
namespace Surface_mesh_topology {
  /*!
    \ingroup PkgSurfaceMeshTopologyClasses

    The class `Curves_on_surface_topology` provides methods to compute shortest non contractible cycles and to test homotopy on paths. Each object of this class is constructed from an external mesh, either a \ref CombinatorialMap "2D combinatorial map" or a model of a FaceGraph. It maintains a correspondence between this mesh and an internal representation, computed the first time an homotopy test is called. The user must not modify the input surface as long as homotopy tests are performed with this `Curves_on_surface_topology`.

    \tparam Mesh a model of `CombinatorialMap` or of `FaceGraph`
  */
  template<typename Mesh>
  class Curves_on_surface_topology
  {
  public:

    /*!
      A descriptor to `Dart` for combinatorial/generalized maps, or a halfedge descriptor for models of the `FaceGraph` concept.
    */
    typedef unspecified_type halfedge_descriptor;

    /*! creates a `Curves_on_surface_topology` object using `amesh` as input.
     */
    Curves_on_surface_topology(const Mesh& amesh);

    /*! returns `true` if the closed paths `p1` and `p2` are freely homotopic.
     *  @pre `p1` and `p2` must be two paths on `amesh`.
     */
    bool are_freely_homotopic(const Path_on_surface<Mesh>& p1,
                              const Path_on_surface<Mesh>& p2) const;

    /*! returns `true` if the paths `p1` and `p2` are homotopic with fixed endpoints. The paths `p1` and `p2` must have the same endpoints but must not be closed. Equivalent to \link Surface_mesh_topology::Curves_on_surface_topology::is_contractible `is_contractible(q)`\endlink where `q` is the concatenation of `p1` and the reverse of `p2`.
     *  @pre `p1` and `p2` must be two paths on `amesh`.
     */
    bool are_homotopic_with_fixed_endpoints(const Path_on_surface<Mesh>& p1,
                                  const Path_on_surface<Mesh>& p2) const;

    /*! returns `true` if the closed path `p` is contractible.
     *  @pre `p` must be a closed path on `amesh`.
     */
    bool is_contractible(const Path_on_surface<Mesh>& p) const;

    /*! returns `true` if the closed path `p` is homotopic to some simple cycle.
     *  @pre `p` must be a closed path on `amesh`.
     */
    bool is_homotopic_to_simple_cycle(const Path_on_surface<Mesh>& p) const;

    /*!  returns a non-contractible cycle of type `Path_on_surface` with minimal number of edges. This number of edges is the edge width of the mesh.
     */
    Path_on_surface<Mesh> compute_edge_width() const;

    /*! returns a non-contractible cycle of type `Path_on_surface` with minimal length, where the length of a cycle is the sum of the weights of its edges computed thanks to the WeightFunctor `wf`. By default, all the edge weights are set to 1 (thanks to the `Unit_weight_functor` functor).
     */
    template <class WeightFunctor=Unit_weight_functor>
    Path_on_surface<Mesh> compute_shortest_non_contractible_cycle(const WeightFunctor& wf=WeightFunctor()) const;

    /*! returns a non-contractible cycle of type `Path_on_surface` with minimal length going through the source vertex of `d`, where the length of a cycle is the sum of the weights of its edges computed thanks to the WeightFunctor `wf`. By default, all the edge weights are set to 1 (thanks to the `Unit_weight_functor` functor).
     */
    template <class WeightFunctor=Unit_weight_functor>
    Path_on_surface<Mesh> compute_shortest_non_contractible_cycle_with_base_point(halfedge_descriptor d, const WeightFunctor& wf=WeightFunctor()) const;

    /*! returns a vector of darts representing a non-contractible curve with a minimal number of intersection with the graph of the mesh. This curve can be described by the alternating sequence of faces and vertices it goes through, so that each dart in the returned vector belongs to both a face and the next vertex in the alternating sequence. (Here, faces and vertices are viewed as subsets of darts.) The size of the returned vector is the face width of the mesh.
     */
    std::vector<halfedge_descriptor> compute_face_width() const;

    /*! set whether the function should output error message to `std::cerr` when the prerequisite of the argument(s) is not met.
     * Affects \link Surface_mesh_topology::Curves_on_surface_topology::are_freely_homotopic `are_freely_homotopic(p1, p2)`\endlink, \link Surface_mesh_topology::Curves_on_surface_topology::are_homotopic_with_fixed_endpoints `are_homotopic_with_fixed_endpoints(p1, p2)`\endlink, \link Surface_mesh_topology::Curves_on_surface_topology::is_contractible `is_contractible(p)`\endlink, and \link Surface_mesh_topology::Curves_on_surface_topology::is_homotopic_to_simple_cycle `is_homotopic_to_simple_cycle(p)`\endlink
     */
    void set_verbose(bool is_verbose);
  };

  /*!
    \ingroup PkgSurfaceMeshTopologyClasses

    A model of `WeightFunctor` assigning unit weight to every edge.
   */
  struct Unit_weight_functor
  {
    /// Number type of the weights.
    using Weight_t=unsigned int;
  };

  /*!
    \ingroup PkgSurfaceMeshTopologyClasses

    A model of `WeightFunctor` assigning its Euclidean length to every edge.
   *  @pre `amesh` should have points associated with their vertices.
   */
  template<typename Mesh>
  struct Euclidean_length_weight_functor
  {
    /// Number type of the weights.
    using Weight_t=double;

    /// creates a Euclidean_length_weight_functor given a mesh.
    Euclidean_length_weight_functor(const Mesh& m);
  };

 }
}
