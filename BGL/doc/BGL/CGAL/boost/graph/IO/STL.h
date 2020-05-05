namespace CGAL {


/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from data in the STL format.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \pre The data must represent a 2-manifold

  \see \ref IOStreamSTL
*/
template <typename FaceGraph, typename NamedParameter>
bool read_STL(std::istream& in,
              FaceGraph& g,
              const NamedParameter& np);

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from the file `fname` in the STL format.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \pre The data must represent a 2-manifold

  \see \ref IOStreamSTL
*/
template <typename FaceGraph, typename NamedParameter>
bool read_STL(const char* fname, FaceGraph& g, const NamedParameter& np);


/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the stream `out` in the STL format.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
  \cgalNamedParamsEnd
  \pre The graph must contain only triangle faces.

  \see \ref IOStreamSTL
*/
template <typename FaceGraph, typename NamedParameters>
bool write_STL(std::ostream& out,
               const FaceGraph& g,
          const NamedParameters& np);
/*!
\ingroup PkgBGLIOFct

 writes the graph `g` in the STL format into a file named `fname`.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
  \cgalNamedParamsEnd

 \sa Overloads of this function for specific models of the concept `FaceGraph`.

 \see \ref IOStreamSTL
*/
template <typename FaceGraph, typename NamedParameter>
bool write_STL(const char* fname, const FaceGraph& g, const NamedParameter& np);
} // namespace CGAL
