
namespace CGAL {
  /*!
    \ingroup PkgSurfaceMeshTopologyClasses
    
    The class `Surface_mesh_curve_topology` represents a combinatorial map internally represented as a quadrangulation. Each object in this class is constructed from an external mesh, either a 2D combinatorial map or a model of a FaceGraph. It maintains a correspondence between this mesh and its internal quadrangulation.
    
    \tparam Mesh a model of `CombinatorialMap` or of `FaceGraph`
  */
  template<typename Mesh>
  class Surface_mesh_curve_topology
  {
  public:

    /*! Constructor. Creates a Surface_mesh_curve_topology object using amesh as input.
     */
    Surface_mesh_curve_topology(const Mesh& amesh);

    /*! Returns `true` if the closed paths p1 and p2 are freely homotopic.
     *  @pre p1 and p2 must be two paths on amesh.
     */
    bool are_freely_homotopic(const Path_on_surface<Mesh>& p1,
                              const Path_on_surface<Mesh>& p2) const;

    /*! Returns `true` if the paths p1 and p2 are homotopic with fixed endpoints. The paths p1 and p2 must have the same endpoints but must not be closed. Equivalent to `Surface_mesh_curve_topology::is_contractible` (q) where q is the concatenation of p1 and the reverse of p2.
     *  @pre p1 and p2 must be two paths on amesh.
     */
    bool are_base_point_homotopic(const Path_on_surface<Mesh>& p1,
                                  const Path_on_surface<Mesh>& p2) const;

    /*! Returns `true` if the closed path p is contractible.
     *  @pre p must be two paths on amesh.
     */
    bool is_contractible(const Path_on_surface<Mesh>& p) const;    
  };

}
