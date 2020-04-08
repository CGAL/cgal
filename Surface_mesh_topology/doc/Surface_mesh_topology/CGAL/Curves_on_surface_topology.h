
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

    /*!  returns a non-contractible cycle of type `Path_on_surface` with minimal number of edges. This edge number is the edgewidth of the mesh.
     */
    Path_on_surface<Mesh> compute_edgewidth() const;

    /*! returns a non-contractible cycle of type `Path_on_surface` with minimal length, where the length of a cycle is the sum of the weights of its edges computed thanks to the WeightFunctor `wf`.
     */
    template <class WeightFunctor>
    Path_on_surface<Mesh> compute_shortest_noncontractible_cycle(const WeightFunctor& wf) const;

    /*! returns a non-contractible cycle of type `Path_on_surface` with minimal length going through the source vertex of `dh`, where the length of a cycle is the sum of the weights of its edges computed thanks to the WeightFunctor `wf`. When omitted, `wf` defaults to  'Unit_weight_functor'
     */
    template <class WeightFunctor>
    Path_on_surface<Mesh> compute_shortest_noncontractible_cycle_with_basepoint(halfedge_descriptor dh, const WeightFunctor& wf) const;

    /*! returns a vector of darts representing a non-contractible curve with a minimal number of intersection with the graph of the mesh. This curve can be decribed by the alternating sequence of faces and vertices it goes through, so that each dart in the returned vector belongs to both a face and the next vertex in the alternating sequence. (Here, faces and vertices are viewes as subsets of darts.) The size of the returned vector is the 'facewidth' of the mesh. 
     */
    std::vector<halfedge_descriptor> compute_facewidth() const;

    /*! A model of WeightFunctor assigning unit weight to every edge.
     */
    struct Unit_weight_functor
    {
      using Weight_t=unsigned int;
      
      Weight_t operator() (halfedge_descriptor dh)
	};

     /*! A model of WeightFunctor assigning its Euclidean length to every edge.
      *  @pre The vertices of `amesh` should have coordinates in their attributes.
     */
      struct Euclidean_length_weight_functor
      {
	using Weight_t=double;
	
	Euclidean_length_weight_functor(const Mesh& m)
	{}
	
	Weight_t operator() (halfedge_descriptor hd) const
	  };

  };
 }
}
