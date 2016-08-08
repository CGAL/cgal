/*!
\ingroup PkgGeneralizedMapsConcepts
\cgalConcept

The concept `GeneralizedMapItems` allows to customize a <I>d</I>D generalized map by choosing the type of darts, and by enabling and disabling some attributes. For that, it defines an inner class template named \ref GeneralizedMapItems::Dart_wrapper "Dart_wrapper", with one template parameter, `GMap`, a model of the `GeneralizedMap` concept. This inner class must define two types: `%Dart` and `%Attributes`.

\cgalHasModel \ref CGAL::Generalized_map_min_items "CGAL::Generalized_map_min_items<d>"

\sa `GeneralizedMap`
\sa `GMapDart`

  \cgalHeading{Example}

  The following examples show two possible models of the `GeneralizedMapItems` concept: the first one for a 4D generalized map without enabled attributes, the second one for a 3D generalized map with edge attributes enabled, and associated with a \ref CGAL::Cell_attribute "Cell_attribute" containing an `int`.

  \code{.cpp}
  struct Exemple_Item_4
  {
    template < class GMap >
    struct Dart_wrapper
    {
      typedef CGAL::GMap_dart<4, GMap> Dart;
      typedef CGAL::cpp11::tuple<>     Attributes;
    };
  };

  struct Exemple_Item_3
  {
    template < class GMap >
    struct Dart_wrapper
    {
      typedef CGAL::GMap_dart<3, GMap>             Dart;
      typedef CGAL::Cell_attribute<GMap, int>      Edge_attrib;
      typedef CGAL::cpp11::tuple<void,Edge_attrib> Attributes;
    };
  };
  \endcode
*/
class GeneralizedMapItems {
public:

  /*!
    Wrapper class defining type of darts and types of attributes.
    The class `%Dart_wrapper<GMap>` must provide:

  - `%Dart_wrapper<GMap>::%Dart`, the type of dart, a model of the `GMapDart` concept.
  - `%Dart_wrapper<GMap>::%Attributes` The tuple of attributes, containing at most \ref GeneralizedMap::dimension "GMap::dimension+1" types (one for each possible cell of the generalized map). Each type of the tuple must be either a model of the `CellAttribute` concept or `void`. The first type corresponds to 0-attributes, the second to 1-attributes and so on. If the \f$ i^{\mbox{th}}\f$ type in the tuple is `void`, (<I>i</I>-1)-attributes are disabled. Otherwise, (<I>i</I>-1)-attributes are enabled and have the given type. If the size of the tuple is <I>k</I>, with <I>k</I><\ref GeneralizedMap::dimension "GMap::dimension+1", \f$ \forall\f$<I>i</I>: <I>k</I>\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref GeneralizedMap::dimension "GMap::dimension", <I>i</I>-attributes are disabled.

  \note It can be implemented using a nested template class.
  */
  template <typename GMap>
  using Dart_wrapper = unspecified_type;

}; /* end #GeneralizedMapItems */
