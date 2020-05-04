
namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

/*!
 * \ingroup PkgBGLIOFct
 *
 * \brief reads a PolyData in the VTP format into a triangulated surface mesh.
 *
 * \tparam FaceGraph a model of `FaceListGraph`.
 *
 * \param fname the path to the file that will be read.
 * \param g the output mesh.
 *
 * \pre \cgal needs to be configured with the VTK Libraries for this function to be available.
 * \cgalNamedParamsBegin
 *  \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
 *    If this parameter is omitted, an internal property map for
 *    `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
 * \cgalNamedParamsEnd
 * \pre The data must represent a 2-manifold
 * \see \ref IOStreamVTK
*/
template<typename FaceGraph, typename NamedParameter>
bool read_VTP(const char* fname, FaceGraph& g, const NamedParameter& np);

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

/*! \ingroup PkgBGLIOFct
 *
 * \brief  writes a triangulated surface mesh in the `PolyData` XML format.
 *
 * \tparam FaceGraph a model of `FaceListGraph` with only triangle faces.
 * \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * \param os the stream used for writing.
 * \param g the triangle mesh to be written.
 * \param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the
 * ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{use_binary_mode} a Boolean indicating if the
 *    data should be written in binary (`true`, the default) or in ASCII (`false`).
 *     \cgalParamEnd
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to
 * the vertices of `g`. If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_Point` must be available in `FaceGraph`.
 *     \cgalParamEnd
 *    \cgalParamBegin{vertex_index_map} the property map with the indices associated to
 * the vertices of `g`.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 * \see \ref IOStreamVTK
 */
template<typename FaceGraph, typename NamedParameter>
void write_VTP(std::ostream& os,
               const FaceGraph& g,
               const NamedParameter& np);

/*! \ingroup PkgBGLIOFct
 *
 * \brief  writes a triangulated surface mesh the file `fname`, in the `PolyData` XML format.
 *
 * \tparam FaceGraph a model of `FaceListGraph` with only triangle faces.
 * \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * \param os the stream used for writing.
 * \param g the triangle mesh to be written.
 * \param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the
 * ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{use_binary_mode} a Boolean indicating if the
 *    data should be written in binary (`true`, the default) or in ASCII (`false`).
 *     \cgalParamEnd
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to
 * the vertices of `g`. If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_Point` must be available in `FaceGraph`.
 *     \cgalParamEnd
 *    \cgalParamBegin{vertex_index_map} the property map with the indices associated to
 * the vertices of `g`.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 * \see \ref IOStreamVTK
 */
template<typename FaceGraph, typename NamedParameter>
void write_VTP(const char* fname,
               const FaceGraph& g,
               const NamedParameter& np);
} // namespace CGAL
