#ifndef AABB_TREE_TESTS_AABB_TRAITS_CONSTRUCT_BY_SORTING_H
#define AABB_TREE_TESTS_AABB_TRAITS_CONSTRUCT_BY_SORTING_H

#include <CGAL/AABB_traits.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/property_map.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>

#include <boost/property_map/function_property_map.hpp>

namespace CGAL {

  // forward declaration
  template<typename AABBTraits>
  class AABB_tree;


  /// This traits class provides the same functionality of the default \ref AABB_traits class,
  /// it can be used to enable faster construction of trees, at the cost of lower traversal speed.
  ///
  /// This is done by sorting the primitives along the Hilbert curve,
  /// rather than repeatedly partitioning them along the longest axis of their bounding box.
  /// The result is a tree that may contain nodes that have very high aspect ratios,
  /// which means slower traversals on average.
  ///
  /// In practice, construction can be up to 50% faster, and traversal tends to be around 20% slower.
  /// The breakeven point at which building a tree using the more optimal approach becomes worthwhile
  /// varies depending on the number of primitives.
  /// Generally, the sort-based construction is useful until the tree is used for
  /// thousands to hundreds of thousands of traversals.
  ///
  /// \cgalModels AABBTraits
  /// \cgalModels AABBRayIntersectionTraits
  ///
  /// \tparam GeomTraits
  /// must  be a model of the concept \ref AABBGeomTraits,
  /// and provide the geometric types as well as the intersection tests and computations.
  ///
  /// \tparam AABBPrimitive
  /// provide the type of primitives stored in the AABB_tree.
  /// It is a model of the concept `AABBPrimitive` or `AABBPrimitiveWithSharedData`.
  ///
  /// \tparam BboxMap
  /// must be a model of `ReadablePropertyMap` that has as key type a primitive id,
  /// and as value type a `Bounding_box`.
  /// If the type is `Default` the `Datum` must have the
  /// member function `bbox()` that returns the bounding box  of the primitive.
  ///
  /// \tparam ConcurrencyTag
  /// Must be one of CGAL::Sequential_tag, CGAL::Parallel_tag, or CGAL::Parallel_if_available_tag.
  /// It is used to determine the algorithm used by the underlying hilbert sort.
  ///
  /// If the argument `GeomTraits` is a model of the concept \ref
  /// AABBRayIntersectionGeomTraits, this class is also a model of \ref
  /// AABBRayIntersectionTraits.
  ///
  template<typename GeomTraits, typename AABBPrimitive, typename BboxMap = Default, class ConcurrencyTag = Sequential_tag>
  class AABB_traits_construct_by_sorting : public AABB_traits<GeomTraits, AABBPrimitive, BboxMap> {

  public:
    typedef GeomTraits Geom_traits;
    typedef AABBPrimitive Primitive;
    typedef BboxMap Bbox_map;
    typedef ConcurrencyTag Concurrency_tag;

    typedef AABB_traits_construct_by_sorting<GeomTraits, AABBPrimitive, BboxMap, Concurrency_tag> Traits;
    typedef internal::Primitive_helper <Traits> Helper;
    typedef typename GeomTraits::Point_3 Point_3;

  public:

    /// \internal
    ///
    /// A splitting algorithm that sorts the range provided along the hilbert curve exactly once,
    /// further invocations of this functor will not attempt to sort the list
    ///
    class Split_primitives {
      const Traits &m_traits;
      mutable bool has_been_sorted = false;

    public:

      struct Get_reference_point : public std::unary_function<const Primitive &, typename Traits::Point_3> {
        const Traits &m_traits;

        explicit Get_reference_point(const Traits &traits)
                : m_traits(traits) {}

        typename Traits::Point_3 operator()(const Primitive &p) const {
          return Helper::get_reference_point(p, m_traits);
        }
      };

      explicit Split_primitives(const Traits &traits)
              : m_traits(traits) {}

      template<typename PrimitiveIterator>
      void operator()(PrimitiveIterator first,
                      PrimitiveIterator beyond,
                      const typename Traits::Bounding_box &bbox) const {

        // If this is our first time splitting the primitives, sort them along the hilbert curve
        // This should generally put nearby primitives close together in the list
        if (!has_been_sorted) {

          // Create a property map using our Get_reference_point functor
          auto property_map = boost::make_function_property_map<Primitive, Point_3, Get_reference_point>(
                  Get_reference_point(m_traits)
          );

          // Our search traits will use that property map
          typedef CGAL::Spatial_sort_traits_adapter_3<GeomTraits, decltype(property_map)> Search_traits_3;

          // Perform our hilbert sort using the search traits type with our custom property map
          CGAL::hilbert_sort<Concurrency_tag>(first, beyond, Search_traits_3(property_map));

          // In the future, it's not necessary to re-sort the primitives (we can blindly partition them in the middle)
          has_been_sorted = true;
        }
      }
    };

    /// Factory function that produces a splitting functor which is associated with this traits class.
    ///
    /// \return a new splitting function
    ///
    Traits::Split_primitives split_primitives_object() const { return Traits::Split_primitives(*this); }
  };
}

#endif //AABB_TREE_TESTS_AABB_TRAITS_CONSTRUCT_BY_SORTING_H
