namespace CGAL {


/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from the stream `in` in the OBJ format.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
  \cgalNamedParamsEnd
    `vertex_normal_map` the property map with the normals associated to the vertices of `g`.


  \returns `true` if the resulting mesh is valid.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \pre The data must represent a 2-manifold

  \see \ref IOStreamOBJ
*/
template <typename FaceGraph, typename NamedParameter>
bool read_OBJ(std::istream& in,
              FaceGraph& g,
              const NamedParameter& np);
/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from data in the OFF format. Ignores comment lines which start with a hash,
  and lines with whitespace.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
  \cgalNamedParamsEnd
  `vertex_normal_map` the property map with the normals associated to the vertices of `g`.

  \returns `true` if the resulting mesh is valid.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \pre The data must represent a 2-manifold

  \see \ref IOStreamOBJ
*/
template <typename FaceGraph, typename NamedParameter>
bool read_OBJ(const char* fname,
              FaceGraph& g,
              const NamedParameter& np);

/*!
 \ingroup PkgBGLIOFct

  writes the graph `g` in the OBJ format.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
  \cgalNamedParamsEnd
  `vertex_normal_map` the property map with the normals associated to the vertices of `g`.

  \returns `true` if writing was successful.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \see \ref IOStreamOBJ
*/
template <typename FaceGraph, typename NamedParameter>
bool write_OBJ(std::ostream& os,
               const FaceGraph& g,
               const NamedParameter& np);

/*!
\ingroup PkgBGLIOFct

 writes the graph `g` in the OFF format into a file named `fname`.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
  \cgalNamedParamsEnd
  `vertex_normal_map` the property map with the normals associated to the vertices of `g`.

  \returns `true` if writing was successful.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

 \see \ref IOStreamOBJ
*/
template <typename FaceGraph, typename NamedParameter>
bool write_OBJ(const char* fname,
               const FaceGraph& g,
               const NamedParameter& np);
} // namespace CGAL
