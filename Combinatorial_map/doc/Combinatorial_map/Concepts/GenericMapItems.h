/*!
\ingroup PkgCombinatorialMapsConcepts
\cgalConcept

The concept `GenericMapItems` allows to customize a <I>d</I>D generic map by choosing the information associated with darts, by enabling and disabling some attributes, and by choosing to use indices or handles. For that, it defines an inner class template named \link GenericMapItems::Dart_wrapper `Dart_wrapper`\endlink, with one template parameter, `Map`, a model of the `GenericMap` concept. This inner class can define the two types `%Dart_info` and `%Attributes`.

\cgalHasModelsBegin
\cgalHasModelsBare{\link CGAL::Generic_map_min_items `CGAL::Generic_map_min_items`\endlink}
\cgalHasModelsEnd

\sa `GenericMap`

\cgalHeading{Example}

The following examples show two possible models of the `GenericMapItems` concept: the first one for a generic map without dart information, nor enabled attributes, the second one for a generic map that uses std::uint16_t as indices, with a `double` associated with each dart, and edge attributes enabled, and associated with a \link CGAL::Cell_attribute `Cell_attribute`\endlink containing an `int`.

  \code{.cpp}
  struct Exemple_Item_1
  {
    template < class CMap >
    struct Dart_wrapper
    {};
  };

  struct Exemple_Item_2
  {
    typedef CGAL::Tag_true Use_index; // CGAL::Tag_false by default
    typedef std::uint16_t Index_type; // std::uint32_t by default
    template < class GMap >
    struct Dart_wrapper
    {
      typedef double Dart_info;
      typedef CGAL::Cell_attribute<GMap, int> Edge_attrib;
      typedef std::tuple<void,Edge_attrib> Attributes;
    };
  };
  \endcode
*/
class GenericMapItems {
public:
  /*!
    If `Use_index` is  `CGAL::Tag_true`, uses indices as descriptors. Otherwise (if this type is not defined or different from `CGAL::Tag_true`), use handles.
  */
  typedef unspecified_type Use_index;

  /*!
    Number type used for indices when Use_index is CGAL::Tag_true. By default, use std::uint32_t.
  */
  typedef unspecified_type Index_type;

  /*!
    Wrapper class defining type of information associated with darts and types of attributes. The class `%Dart_wrapper<Map>` must provide:

  - `%Dart_wrapper<Map>::%Dart_info`, the type of information associated with darts. Not defined or equal to `void` to have no information.
  - `%Dart_wrapper<Map>::%Attributes` The tuple of attributes, containing at most \link GenericMap::dimension `Map::dimension+1`\endlink types (one for each possible cell of the generic map). Each type of the tuple must be either a model of the `CellAttribute` concept or `void`. The first type corresponds to 0-attributes, the second to 1-attributes and so on. If the \f$ i^{\mbox{th}}\f$ type in the tuple is `void`, (<I>i</I>-1)-attributes are disabled. Otherwise, (<I>i</I>-1)-attributes are enabled and have the given type. If the size of the tuple is <I>k</I>, with <I>k</I><\link GenericMap::dimension `Map::dimension`\endlink+1, \f$ \forall\f$<I>i</I>: <I>k</I>\f$ \leq\f$<I>i</I>\f$ \leq\f$\link GenericMap::dimension `Map::dimension`\endlink, <I>i</I>-attributes are disabled. If this type is not defined, all attributes are disabled.

  \note It can be implemented using a nested template class.
  */
  template <typename Map>
  using Dart_wrapper = unspecified_type;

}; /* end #GenericMapItems */
