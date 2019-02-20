
namespace CGAL {
  /*!
    \ingroup PkgSurfaceMeshTopologyClasses
    
    The class `Surface_mesh_curve_topology` represents a combinatorial map internally represented as a quadrangulation. Each object in this class is constructed from an external combinatorial map. It maintains a correspondence between this combinatorial map and its internal quadrangulation.
    
    \tparam Map a model of `CombinatorialMap`
  */
  template<typename Map>
  class Surface_mesh_curve_topology
  {
  public:

    /*! Constructor. Creates a Surface_mesh_curve_topology object using amap as input.
     */
    Surface_mesh_curve_topology(Map& amap);

    /*! Returns true if the closed paths p1 and p2 are freely homotopic.
     */
    bool are_freely_homotopic(const Path_on_surface<Map>& p1,
                              const Path_on_surface<Map>& p2) const;

    /*! Returns true if the paths p1 and p2 are homotopic with fixed endpoints. The paths p1 and p2 should have the same endpoints but need not be closed. Equivalent to `is_contractible` (q) where q is the concatenation of p1 and the reverse of p2.
     */
    bool are_base_point_homotopic(const Path_on_surface<Map>& p1,
                                  const Path_on_surface<Map>& p2) const;

    /*! Returns true if the closed path p is contractible.
     */
    bool is_contractible(const Path_on_surface<Map>& p) const;    
  };

}
