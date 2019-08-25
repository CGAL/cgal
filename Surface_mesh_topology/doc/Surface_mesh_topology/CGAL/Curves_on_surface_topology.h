
namespace CGAL {
namespace Surface_mesh_topology {
  /*!
    \ingroup PkgSurfaceMeshTopologyClasses
    
    The class `Curves_on_surface_topology` provides methods to test homotopy on paths. Each object of this class is constructed from an external mesh, either a \ref CombinatorialMap "2D combinatorial map" or a model of a FaceGraph. It maintains a correspondence between this mesh and an internal representation, computed the first time an homotopy test is called. The user must not modify the input surface as long as homotopy tests are performed with this `Curves_on_surface_topology`.
    
    \tparam Mesh a model of `CombinatorialMap` or of `FaceGraph`
  */
  template<typename Mesh>
  class Curves_on_surface_topology
  {
  public:

    /*!
      %halfedge_descriptor type. A handle to Dart for combinatorial/generalized maps, or a halfedge descriptor for models of FaceGraph.
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

    /*! returns `true` if the paths `p1` and `p2` are homotopic with fixed endpoints. The paths `p1` and `p2` must have the same endpoints but must not be closed. Equivalent to `Curves_on_surface_topology::is_contractible(q)` where `q` is the concatenation of `p1` and the reverse of `p2`.
     *  @pre `p1` and `p2` must be two paths on `amesh`.
     */
    bool are_homotopic_with_fixed_endpoints(const Path_on_surface<Mesh>& p1,
                                  const Path_on_surface<Mesh>& p2) const;

    /*! returns `true` if the closed path `p` is contractible.
     *  @pre `p` must be a closed path on `amesh`.
     */
    bool is_contractible(const Path_on_surface<Mesh>& p) const;

    /*! returns the edgewidth of the unweighted mesh
     */
    Path_on_surface<Mesh> compute_edgewidth() const;

    /*! returns the edgewidth of the weighted mesh
     */
    template <class WeightFunctor>
    Path_on_surface<Mesh> compute_edgewidth(const WeightFunctor& wf) const;

    /*! returns the shortest non-contractible cycle going through the source vertex of `dh` of the unweighted mesh
     */
    Path_on_surface<Mesh> compute_shortest_noncontractible_cycle_with_basepoint(halfedge_descriptor dh) const;

    /*! returns the shortest non-contractible cycle going through the source vertex of `dh` of the weighted mesh
     */
    template <class WeightFunctor>
    Path_on_surface<Mesh> compute_shortest_noncontractible_cycle_with_basepoint(halfedge_descriptor dh, const WeightFunctor& wf) const;

    /*! returns a list of one dart per face of the facewidth
     */
    std::vector<halfedge_descriptor> compute_facewidth() const;
  };

}
}
