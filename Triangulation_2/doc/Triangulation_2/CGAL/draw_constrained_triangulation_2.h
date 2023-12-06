namespace CGAL {

  /*!
\ingroup PkgDrawConstrainedTriangulation2

opens a new window and draws a constrained triangulation. If the triangulation has constraints they are drawn. The faces inside and outside of the domain, based on the property map, are drawn in different colors. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam CT2 which must be an instanciation of a `CGAL::Constrained_triangulation_2<...>`.
\tparam InDomainPMap a class model of `ReadablePropertyMap` with `CT2::Face_handle` as key type and `bool` as value type.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param ct2 the constrained triangulation to draw.
\param ipm the property map defining the faces which are in the domain.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class Gt, class Tds, class Itag, class InDomainPMap, class GSOptions>

 void CGAL::draw(const CGAL::Constrained_triangulation_2<Gt, Tds, Itag>& ct2, InDomainPMap ipm, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
  template<class CT2, class InDomainPMap, class GSOptions>
  void draw(const CT2& act2, InDomainPMap ipm, const GSOptions& gso);

/*!
\ingroup PkgDrawConstrainedTriangulation2

A shortcut to `CGAL::draw(at2, ipm, Graphics_scene_options{})`.
*/
  template<class T2, class InDomainPMap>
  void draw(const T2& at2, InDomainPMap ipm);

} /* namespace CGAL */
