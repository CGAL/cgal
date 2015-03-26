namespace CGAL {

/*!
\ingroup PkgArrExtensionTraits

The class `Arr_extension_default_traits_2` serves as a traits class for all visibility polygon calculation function.
This class extends the vertex, halfedges and face. User also may define their own traits class to choose which to extend.

\cgalModels `ArrExtensionTraits_2`

\sa `CGAL::Arr_extended_dcel`


*/

template< typename A_ >
class Arr_extension_default_traits_2 {
public:

    /// \name Types
    /// @{
    /*!
     *
     */
    typedef A_::Vertex_iterator Vertex_iterator;

    /*!
     *
     */
    typedef A_::Halfedge_iterator Halfedge_iterator;

    /*!
     *
     */
    typedef A_::Fayce_iterator Face_iterator;

    /// @}

    /// \name Functor classes
    /// @{

    /*!
     *
     */
    class Extend_vertex;

    /*!
     *
     */
    class Extend_halfedge;

    /*!
     *
     */
    class Extend_face;
    /// @}

    /// \name Operations
    /// @{

    /*!
     *
     */
    Extend_vertex extend_vertex_object();

    /*!
     *
     */
    Extend_halfedge extend_halfedge_object();

    /*!
     *
     */
    Extend_face extend_face_object();

    /// @}

};

}  /* end namespace CGAL */
