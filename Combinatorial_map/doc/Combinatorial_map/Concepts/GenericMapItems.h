/*!
\ingroup PkgCombinatorialMapsConcepts
\cgalConcept

The concept `GenericMapItems` allows to customize a <I>d</I>D basic map by choosing the type of darts, and by enabling and disabling some attributes. For that, it defines an inner class template named
\link GenericMapItems::Dart_wrapper `Dart_wrapper`\endlink, with one template parameter, `Map`, a model of the `GenericMap` concept. This inner class must define two types: `%Dart` and `%Attributes`.

\cgalHasModel \link CGAL::Combinatorial_map_min_items `CGAL::Combinatorial_map_min_items<d>`\endlink
\cgalHasModel \link CGAL::Generalized_map_min_items `CGAL::Generalized_map_min_items<d>`\endlink

\sa `GenericMap`
\sa `Dart`

  \cgalHeading{Example}

  The following examples show two possible models of the `GenericMapItems` concept: the first one for a 4D combinatorial map without enabled attributes, the second one for a 3D generalized map with edge attributes enabled, and associated with a \link CGAL::Cell_attribute `Cell_attribute`\endlink containing an `int`.

  \code{.cpp}
  struct Exemple_Item_4
  {
    template < class CMap >
    struct Dart_wrapper
    {
      typedef CGAL::Combinatorial_map_dart<4, CMap> Dart;
      typedef CGAL::cpp11::tuple<> Attributes;
    };
  };

  struct Exemple_Item_3
  {
    template < class GMap >
    struct Dart_wrapper
    {
      typedef CGAL::Generalized_map_dart<3, GMap> Dart;
      typedef CGAL::Cell_attribute<GMap, int> Edge_attrib;
      typedef CGAL::cpp11::tuple<void,Edge_attrib> Attributes;
    };
  };
  \endcode
*/
class GenericMapItems {
public:

  /*!
    Wrapper class defining type of darts and types of attributes. The class `%Dart_wrapper<CMap>` must provide:

  - `%Dart_wrapper<CMap>::%Dart`, the type of dart, a model of the `BasicDart` concept.
  - `%Dart_wrapper<CMap>::%Attributes` The tuple of attributes, containing at most \link GenericMap::dimension `Map::dimension+1`\endlink types (one for each possible cell of the basic map). Each type of the tuple must be either a model of the `CellAttribute` concept or `void`. The first type corresponds to 0-attributes, the second to 1-attributes and so on. If the \f$ i^{\mbox{th}}\f$ type in the tuple is `void`, (<I>i</I>-1)-attributes are disabled. Otherwise, (<I>i</I>-1)-attributes are enabled and have the given type. If the size of the tuple is <I>k</I>, with <I>k</I><\link GenericMap::dimension `Map::dimension+1`\endlink, \f$ \forall\f$<I>i</I>: <I>k</I>\f$ \leq\f$<I>i</I>\f$ \leq\f$\link GenericMap::dimension `Map::dimension`\endlink, <I>i</I>-attributes are disabled.

  \note It can be implemented using a nested template class.
  */
  template <typename Map>
  using Dart_wrapper = unspecified_type;

}; /* end #GenericMapItems */
