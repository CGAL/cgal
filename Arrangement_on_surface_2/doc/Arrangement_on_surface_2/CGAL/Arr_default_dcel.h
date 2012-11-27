
namespace CGAL {

/*!
\ingroup PkgArrangement2DCEL

The default <span class="textsc">Dcel</span> class used by the `Arrangement_2` class-template 
is parameterized by a traits class, which is a model of the 
`ArrangementBasicTraits_2` concept. It simply uses the nested 
`Traits::Point_2` and `Traits::X_monotone_curve_2` to instantiate 
the base vertex and halfedge types, respectively. Thus, the default 
<span class="textsc">Dcel</span> records store no other information, except for the topological 
incidence relations and the geometric data attached to vertices and edges. 

\cgalModels `ArrangementDcelWithRebind`

\sa `Arr_dcel_base<V,H,F>` 
*/
template< typename Traits >
class Arr_default_dcel : 
    public Arr_dcel_base< Arr_vertex_base<typename Traits_::Point_2>,
                          Arr_halfedge_base<typename Traits_::X_monotone_curve_2>,
                          Arr_face_base>
{
}; /* end Arr_default_dcel */
} /* end namespace CGAL */
