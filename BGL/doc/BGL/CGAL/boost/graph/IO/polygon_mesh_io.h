namespace CGAL{
/*!
 * \ingroup PkgBGLIOFct
 * \brief reads a polygon mesh from a file.
 * \tparam FaceGraph a model of `FaceGraph`
 * \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * \param fname the name of the file. Its extension must be one of the following :
 * `.off` (\ref IOStreamOFF "OFF file format") , `.obj` (\ref IOStreamOBJ "OBJ file format"),
 * `.stl` (\ref IOStreamSTL "STL file format"), `.ply` (\ref IOStreamPLY "PLY file format"),
 * `.vtp`(\ref IOStreamVTK "VTP file format")  or `.ts`(\ref IOStreamGocad "GOCAD file format").
 * \param g the mesh
 * \param np optional \ref pmp_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 * \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
 * \cgalNamedParamsEnd
 * Other named parameters may be used according to the file extension.
 * See `PkgBGLIOFct` for an exhaustive list.

 * \return `true` if the reading worked, `false` otherwise.
 *
 * \pre The data must represent a 2-manifold
 * \see \ref IOStreamOFF
 *
 */
template <class FaceGraph, typename NamedParameters>
bool read_polygon_mesh(const std::string& fname,
                       FaceGraph& g,
                       const NamedParameters& np);

/*!
 * \ingroup PkgBGLIOFct
 * \brief writes a polygon mesh in a file.
 * \tparam FaceGraph a model of `FaceGraph`
 * \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * \param fname the name of the file. Its extension must be one of the following :
 * `.off` (\ref IOStreamOFF "OFF file format") , `.obj` (\ref IOStreamOBJ "OBJ file format"),
 * `.stl` (\ref IOStreamSTL "STL file format"), `.ply` (\ref IOStreamPLY "PLY file format"),
 * `.vtp`(\ref IOStreamVTK "VTP file format")  or `.ts`(\ref IOStreamGocad "GOCAD file format").
 * \param g the mesh
 * \param np optional \ref pmp_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 * \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
 * \cgalNamedParamsEnd
 * Other named parameters may be used according to the file extension.
 * See `PkgBGLIOFct`  for an exhaustive list.
 * \return `true` if the writing worked, `false` otherwise.
 *
 * \see \ref IOStreamOFF
 */
template <class FaceGraph, typename NamedParameters>
bool write_polygon_mesh(const std::string& fname,
                       FaceGraph& g,
                       const NamedParameters& np);
}
