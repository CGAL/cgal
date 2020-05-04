namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from data in the OFF format. Ignores comment lines which start with a hash,
  and lines with whitespace.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd

    \cgalParamBegin{vertex_normal_map} the property map with the normals associated to the vertices of `g`.\cgalParamEnd

    \cgalParamBegin{vertex_color_map} the property map with the colors associated to the vertices of `g`.\cgalParamEnd

    \cgalParamBegin{vertex_texture_map} the property map with the textures associated to the vertices of `g`.\cgalParamEnd

    \cgalParamBegin{face_color_map} the property map with the colors associated to the faces of `g`.\cgalParamEnd
  \cgalNamedParamsEnd

  \pre The data must represent a 2-manifold

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \see \ref IOStreamOFF
*/
template <typename FaceGraph, typename NamedParameters>
bool read_OFF(std::istream& in, FaceGraph& g, const NamedParameters& np);

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from `fname`, a file in the OFF format. Ignores comment lines which start with a hash,
  and lines with whitespace.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd

    \cgalParamBegin{vertex_normal_map} the property map with the normals associated to the vertices of `g`.\cgalParamEnd

    \cgalParamBegin{vertex_color_map} the property map with the colors associated to the vertices of `g`.\cgalParamEnd

    \cgalParamBegin{vertex_texture_map} the property map with the textures associated to the vertices of `g`.\cgalParamEnd

    \cgalParamBegin{face_color_map} the property map with the colors associated to the faces of `g`.\cgalParamEnd
    \cgalNamedParamsEnd

  \pre The data must represent a 2-manifold

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \see \ref IOStreamOFF
*/
template <typename FaceGraph, typename NamedParameters>
bool read_OFF(const char* fname, FaceGraph& g, const NamedParameters& np);

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the OFF format.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd

    \cgalParamBegin{vertex_normal_map} the property map with the normals associated to the vertices of `g`.\cgalParamEnd

    \cgalParamBegin{vertex_color_map} the property map with the colors associated to the vertices of `g`.\cgalParamEnd

    \cgalParamBegin{vertex_texture_map} the property map with the textures associated to the vertices of `g`.\cgalParamEnd

    \cgalParamBegin{face_color_map} the property map with the colors associated to the faces of `g`.\cgalParamEnd
    \cgalNamedParamsEnd

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \see \ref IOStreamOFF
*/
template <typename FaceGraph, typename NamedParameters>
bool write_OFF(std::ostream& os, const FaceGraph& g, const NamedParameters& np);

/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the file `fname`, in the OFF format.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd

    \cgalParamBegin{vertex_normal_map} the property map with the normals associated to the vertices of `g`.\cgalParamEnd

    \cgalParamBegin{vertex_color_map} the property map with the colors associated to the vertices of `g`.\cgalParamEnd

    \cgalParamBegin{vertex_texture_map} the property map with the textures associated to the vertices of `g`.\cgalParamEnd

    \cgalParamBegin{face_color_map} the property map with the colors associated to the faces of `g`.\cgalParamEnd
    \cgalNamedParamsEnd

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \see \ref IOStreamOFF
*/
template <typename FaceGraph, typename NamedParameters>
bool write_OFF(const char* fname, const FaceGraph& g, const NamedParameters& np);

} // namespace CGAL
