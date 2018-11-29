//!\ingroup PkgBGLIOFct
//!
//! \brief  writes the content of a triangulated surface mesh in the .vtp
//! XML format.
//! 
//! \tparam TriangleMesh a model of `FaceListGraph` with triangle faces.
//! \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
//! 
//! \param os a `std::ostream`.
//! \param mesh an instance of `TriangleMesh` to be written.
//! \param binary decides if the data should be written in binary(`true`)
//!   or in ASCII(`false`).
//! \param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the 
//! ones listed below
//!
//! \cgalNamedParamsBegin
//!    \cgalParamBegin{vertex_point_map} the property map with the points associated to
//! the vertices of `mesh`. If this parameter is omitted, an internal property map for
//!       `CGAL::vertex_point_t` must be available in `TriangleMesh`.
//!     \cgalParamEnd
//!    \cgalParamBegin{vertex_index_map} the property map with the indices associated to 
//! the vertices of `mesh`. If this parameter is omitted, an internal property map for
//!       `CGAL::vertex_index_t` must be available in `TriangleMesh`.
//!     \cgalParamEnd
//! \cgalNamedParamsEnd
template<class TriangleMesh, 
         class NamedParameters>
void write_VTP(std::ostream& os,
                    const TriangleMesh& mesh,
                    bool binary,
                    const NamedParameters& np);
