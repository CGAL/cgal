
namespace CGAL {

/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the wrl format (VRML 2.0).

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`
    \cgalParamEnd
  \cgalNamedParamsEnd

  \see \ref IOStreamWRL
*/
template <typename FaceGraph, typename NamedParameters>
bool write_WRL(std::ostream& os,
               const FaceGraph& g,
               const NamedParameters& np);
} // namespace CGAL
