namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from data in the PLY format.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map}
      the property map with the points associated to the vertices of `g` .
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`
    \cgalParamEnd
    \cgalParamBegin{vertex_index_map}
      is a property map containing for each vertex of `g` a unique index between `0` and `num_vertices(g)-1`.
    \cgalParamEnd
    \cgalParamBegin{face_color_map} the property map with the colors associated to the faces of `g`.\cgalParamEnd
    \cgalParamBegin{vertex_color_map} the property map with the colors associated to the vertices of `g`.\cgalParamEnd
  \cgalNamedParamsEnd

  \pre The data must represent a 2-manifold

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \see \ref IOStreamPLY
*/
template <typename FaceGraph, typename NamedParameters>
bool read_PLY(std::istream& in, FaceGraph& g, const NamedParameters& np);

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from a file named `fname`, in the PLY format.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map}
      the property map with the points associated to the vertices of `g` .
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`
    \cgalParamEnd
    \cgalParamBegin{vertex_index_map}
      is a property map containing for each vertex of `g` a unique index between `0` and `num_vertices(g)-1`.
    \cgalParamEnd
    \cgalParamBegin{face_color_map} the property map with the colors associated to the faces of `g`.\cgalParamEnd
    \cgalParamBegin{vertex_color_map} the property map with the colors associated to the vertices of `g`.\cgalParamEnd
  \cgalNamedParamsEnd

  \pre The data must represent a 2-manifold

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \see \ref IOStreamPLY
*/
template <typename FaceGraph, typename NamedParameters>
bool read_PLY(const char* fname, FaceGraph& g, const NamedParameters& np);

/*!
 \ingroup PkgBGLIOFct

 Inserts the graph in an output stream in PLY format.

 If provided, the `comments` string is included line by line in
 the header of the PLY stream (each line will be precedeed by
 "comment ").

 The `np` is an optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map}
      the property map with the points associated to the vertices of `g` .
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`
    \cgalParamEnd
    \cgalParamBegin{vertex_index_map}
      is a property map containing for each vertex of `g` a unique index between `0` and `num_vertices(g)-1`.
    \cgalParamEnd
    \cgalParamBegin{face_color_map} the property map with the colors associated to the faces of `g`.\cgalParamEnd
    \cgalParamBegin{vertex_color_map} the property map with the colors associated to the vertices of `g`.\cgalParamEnd
  \cgalNamedParamsEnd

  \see \ref IOStreamPLY
 */
template <class FaceGraph, class NamedParameters>
bool write_PLY(std::ostream& os,
               const FaceGraph& g,
               const std::string& comments,
               const NamedParameters& np);

/*!
  \ingroup PkgBGLIOFct

  Inserts the graph in the output file `fname` in PLY format.

  If provided, the `comments` string is included line by line in
  the header of the PLY stream (each line will be precedeed by
  "comment ").

  The `np` is an optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map}
      the property map with the points associated to the vertices of `g` .
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`
    \cgalParamEnd
    \cgalParamBegin{vertex_index_map}
      is a property map containing for each vertex of `g` a unique index between `0` and `num_vertices(g)-1`.
    \cgalParamEnd
    \cgalParamBegin{face_color_map} the property map with the colors associated to the faces of `g`.\cgalParamEnd
    \cgalParamBegin{vertex_color_map} the property map with the colors associated to the vertices of `g`.\cgalParamEnd
  \cgalNamedParamsEnd

  \see \ref IOStreamPLY
 */
template <class FaceGraph, class NamedParameters>
bool write_PLY(const char* fname,
               const FaceGraph& g,
               const std::string& comments,
               const NamedParameters& np);
} // namespace CGAL
