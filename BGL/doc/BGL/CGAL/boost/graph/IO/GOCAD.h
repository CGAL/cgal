namespace CGAL {

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from data in the TS format.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`
    \cgalParamEnd
  \cgalNamedParamsEnd

  \pre The data must represent a 2-manifold

  \attention The graph `g` is not cleared, and the data from the stream is added.

  \see \ref IOStreamGocad
*/
template <typename FaceGraph, typename NamedParameter>
bool read_GOCAD(std::istream& in,
                FaceGraph& g,
                const NamedParameter& np);

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from the file `fname` in the TS format.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`
    \cgalParamEnd
  \cgalNamedParamsEnd

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \pre The data must represent a 2-manifold

  \see \ref IOStreamGocad
*/
template <typename FaceGraph, typename NamedParameter>
bool read_GOCAD(const char* fname, FaceGraph& g, const NamedParameter& np);

/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the TS format into `os`. `fname` is the
  mandatory name that will be assigned to `g` in the file.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`
    \cgalParamEnd
  \cgalNamedParamsEnd

  \see \ref IOStreamGocad
*/
template <typename FaceGraph, typename NamedParameter>
bool write_GOCAD(std::ostream& os,
                 const char* fname,
                 const FaceGraph& g,
                 const NamedParameter& np);
/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the TS format into a file named `fname`.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`
    \cgalParamEnd
  \cgalNamedParamsEnd

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \see \ref IOStreamGocad
*/
template <typename FaceGraph, typename NamedParameter>
bool write_GOCAD(const char* fname,
                 const FaceGraph& g,
                 const NamedParameter& np);

/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the TS format into `os`. The name
 that will be assigned to `g` in the file is `anonymous`.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`
    \cgalParamEnd
  \cgalNamedParamsEnd

  \see \ref IOStreamGocad
*/
template <typename FaceGraph, typename NamedParameter>
bool write_GOCAD(std::ostream& os, const FaceGraph& g, const NamedParameter& np);
} // namespace CGAL
