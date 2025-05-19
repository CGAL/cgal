namespace CGAL {

  /*!
\ingroup PkgDrawConstrainedTriangulation2

opens a new window and draws a constrained triangulation. If the triangulation has constraints they are drawn. The faces inside and outside of the domain, based on the property map, are drawn in different colors.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam CT2 which must be an instantiation of a `CGAL::Constrained_triangulation_2<...>`.
\tparam InDomainPMap a class model of `ReadablePropertyMap` with `CT2::Face_handle` as key type and `bool` as value type.

\param ct2 the constrained triangulation to draw.
\param ipm the property map defining the faces which are in the domain.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class Gt, class Tds, class Itag, class InDomainPMap>

 void CGAL::draw(const CGAL::Constrained_triangulation_2<Gt, Tds, Itag>& ct2, InDomainPMap ipm);
</code>
\cgalAdvancedEnd
*/
  template<class CT2, class InDomainPMap>
  void draw(const CT2& ct2, InDomainPMap ipm);

/*!
\ingroup PkgDrawConstrainedTriangulation2

adds the vertices, edges and faces of `ct2` into the given graphic scene `gs`. If the triangulation has constraints they are drawn. The faces inside and outside of the domain, based on the property map, are drawn in different colors. Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam CT2 which must be an instantiation of a `CGAL::Constrained_triangulation_2<...>`.
\tparam InDomainPMap a class model of `ReadablePropertyMap` with `CT2::Face_handle` as key type and `bool` as value type.

\param ct2 the constrained triangulation to draw.
\param ipm the property map defining the faces which are in the domain.
\param gs the graphic scene to fill.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class Gt, class Tds, class Itag, class InDomainPMap>

 void CGAL::add_to_graphics_scene(const CGAL::Constrained_triangulation_2<Gt, Tds, Itag>& ct2, InDomainPMap ipm, CGAL::Graphics_scene& gs);
</code>
\cgalAdvancedEnd
*/
template<class CT2, class InDomainPMap>
void add_to_graphics_scene(const CT2& ct2, InDomainPMap ipm,
                           CGAL::Graphics_scene& gs);

} /* namespace CGAL */
