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
\cgalParamNBegin{use_trimming}
  \cgalParamDescription{trims exterior mesh elements, i.e., remove volumes that are entirely outside of `tmesh`}
  \cgalParamType{`bool`}
  \cgalParamDefault{`true`}
\cgalParamNEnd

\cgalParamNBegin{use_smoothing}
  \cgalParamDescription{applies Laplacian smoothing to project boundary vertices onto `tmesh`}
  \cgalParamType{`bool`}
  \cgalParamDefault{`true`}
\cgalParamNEnd

\cgalParamNBegin{vertex_point_map}
  \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
  \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::vertex_descriptor`
                 as key type and `GeomTraits::Point_3` as value type, `GeomTraits` being the type of the parameter `geom_traits`}
  \cgalParamDefault{`get(CGAL::vertex_point, tmesh)`}
\cgalParamNEnd

\cgalParamNBegin{geom_traits}
  \cgalParamDescription{an instance of a geometric traits class}
  \cgalParamType{a class model of `Kernel`}
  \cgalParamDefault{a CGAL Kernel deduced from the point type using `CGAL::Kernel_traits`}
  \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
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
