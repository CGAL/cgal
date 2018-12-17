
namespace CGAL {
  /*!
    \ingroup PkgSurfaceMeshTopologyClasses
    
    The class `Surface_mesh_curve_topology` represents XXX
    
    \tparam Map XXX
  */
  template<typename Map>
  class Surface_mesh_curve_topology
  {
  public:

    /*! XXX
     */
    Surface_mesh_curve_topology(Map& amap);

    /*! XXX
     */
    bool are_freely_homotopic(const Path_on_surface<Map>& p1,
                              const Path_on_surface<Map>& p2) const;

    /*! XXX
     */
    bool are_base_point_homotopic(const Path_on_surface<Map>& p1,
                                  const Path_on_surface<Map>& p2) const;

    /*! XXX
     */
    bool is_contractible(const Path_on_surface<Map>& p) const;    
  };

}
