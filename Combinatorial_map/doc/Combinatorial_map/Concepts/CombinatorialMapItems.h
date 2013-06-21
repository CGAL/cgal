/*!
\ingroup PkgCombinatorialMapsConcepts
\cgalConcept

The concept `CombinatorialMapItems` allows to customize a <I>d</I>D combinatorial map
by choosing the type of darts, and by enabling and disabling some
attributes. For that, it defines an inner class template named
\ref CombinatorialMapItems::Dart_wrapper "Dart_wrapper", with one template parameter, `CMap`, a model
of the `CombinatorialMap` concept. This inner class must define
two types: `%Dart` and `%Attributes`.

\cgalHasModel \ref CGAL::Combinatorial_map_min_items "CGAL::Combinatorial_map_min_items<d>"

\sa `CombinatorialMap`
\sa `Dart`

  \cgalHeading{Example}

  The following examples show two possible models of the
  `CombinatorialMapItems` concept: the first one for a 4D
  combinatorial map without enabled attributes, the second one for a 3D
  combinatorial map with edge attributes enabled, and associated with a
  \ref CGAL::Cell_attribute "Cell_attribute" containing an `int`.

  \code{.cpp}
  struct Exemple_Item_4
  {
    template < class CMap >
    struct Dart_wrapper
    {
      typedef CGAL::Dart<4, CMap> Dart;
      typedef CGAL::cpp11::tuple<> Attributes;
    };
  };

  struct Exemple_Item_3
  {
    template < class CMap >
    struct Dart_wrapper
    {
      typedef CGAL::Dart<3, CMap> Dart;
      typedef CGAL::Cell_attribute<CMap, int> Edge_attrib;
      typedef CGAL::cpp11::tuple<void,Edge_attrib> Attributes;
    };
  };
  \endcode
*/
class CombinatorialMapItems {
public:

  /*!
    Wrapper class defining type of darts and types of attributes.
    The class `%Dart_wrapper<CMap>` must provide:

  - `%Dart_wrapper<CMap>::%Dart`, the type of dart, a model of the `Dart` concept.
  - `%Dart_wrapper<CMap>::%Attributes` The tuple of attributes, containing at most
    \ref CombinatorialMap::dimension "CMap::dimension+1" types (one for each possible cell of the combinatorial
    map). Each type of the tuple must be either a model of the
    `CellAttribute` concept or `void`.
    The first type corresponds to 0-attributes,
    the second to 1-attributes and so on.
    If the \f$ i^{\mbox{th}}\f$ type in the tuple is `void`,
    (<I>i</I>-1)-attributes are disabled. Otherwise, (<I>i</I>-1)-attributes are enabled and
    have the given type. If the size of the tuple is <I>k</I>,
    with <I>k</I><\ref CombinatorialMap::dimension "CMap::dimension+1",
    \f$ \forall\f$<I>i</I>: <I>k</I>\f$ \leq\f$<I>i</I>\f$ \leq\f$\ref CombinatorialMap::dimension "CMap::dimension",
    <I>i</I>-attributes are disabled.

  \note It can be implemented using a nested template class.
  */
  template <typename CMap>
  using Dart_wrapper = unspecified_type;

}; /* end #CombinatorialMapItems */
