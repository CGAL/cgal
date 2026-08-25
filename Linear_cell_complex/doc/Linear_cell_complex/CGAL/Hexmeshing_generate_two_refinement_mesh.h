namespace CGAL {

/*!
\ingroup PkgLinearCellComplexConstructions

@brief generates a pure hexahedral mesh from a surface triangle mesh using the two-refinement algorithm described in \cgalCite{cgal:owen2017template-based}.

Starts to create a regular grid of `cube_cells_per_dim`\f$ ^3\f$ voxels. Then refines voxels intersected by the surface `nb_levels` times, while creating transitions between refined and non-refined voxels.

@tparam LCC the target Linear Cell Complex type, or `CGAL::Default`.
@tparam TriangleMesh a model of `FaceListGraph`.
@tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".

@param tmesh the input triangle mesh.
@param cube_cells_per_dim the initial grid resolution per dimension.
@param nb_levels the number of two-refinement levels.
@param np n optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:

\cgalNamedParamsBegin
  \cgalParamNBegin{trim_mesh}
    \cgalParamDescription{trims exterior mesh elements, i.e., remove volumes that are entirely outside of `tmesh`}
    \cgalParamType{`bool`}
    \cgalParamDefault{`true`}
  \cgalParamNEnd

  \cgalParamNBegin{smooth_mesh}
   \cgalParamDescription{applies Laplacian smoothing to project boundary vertices onto `tmesh`}
    \cgalParamType{`bool`}
    \cgalParamDefault{`true`}
  \cgalParamNEnd

 \cgalNamedParamsEnd

@return the resulting linear cell complex representing the hexahedral mesh.

\sa `CGAL::Linear_cell_complex_for_combinatorial_map`

*/
  template<typename LCC=Default, typename TriangleMesh,
           typename NamedParameters=parameters::Default_named_parameters>
  auto generate_hexahedral_mesh_using_two_refinement
  (const TriangleMesh& tmesh, int cube_cells_per_dim, int nb_levels, const NamedParameters& np = parameters::default_values());

} // end namespace CGAL
