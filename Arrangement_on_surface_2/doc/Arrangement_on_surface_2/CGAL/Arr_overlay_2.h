namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2Funcs
\brief Computes the overlay of two arrangements `arr1` and `arr2`, and sets
the output arrangement `res` to represent the overlaid arrangement.

\details
Computes the overlay of two input arrangement
objects, and returns the overlaid arrangement. All three arrangements
can be instantiated with different geometric traits classes and different
\dcel classes (encapsulated in the various topology-traits classes).
The geometry traits of the resulting arrangement is used to construct the
resulting arrangement. This means that all the types (e.g.,
`Traits::Point_2`, `Traits::Curve_2`, and `Traits::Point_2`)
of both input arrangements have to be convertible to the types in the
resulting arrangement. A given overlay-traits object is used to properly
construct the overlaid \dcel that represents the resulting arrangement.

\pre `res` does not refer to either `arr1` or `arr2` (that is, "self overlay" is not supported).

\pre The overlay-traits object `ovl_tr` must model the `OverlayTraits`
  concept, which is able to construct records of the `ResDcel` class on
  the basis of the `Dcel1` and `Dcel2` records that induce them.

\sa `OverlayTraits`
*/
template <class GeomTraitsA, class GeomTraitsB,
          class GeomTraitsRes, class TopTraitsA,
          class TopTraitsB, class TopTraitsRes,
          class OverlayTraits>
void overlay (const Arrangement_2<GeomTraitsA, TopTraitsA>& arr1,
              const Arrangement_2<GeomTraitsB, TopTraitsB>& arr2,
              Arrangement_2<GeomTraitsRes, TopTraitsRes>& arr_res,
              OverlayTraits& ovl_tr);

/*!
\ingroup PkgArrangementOnSurface2Funcs
\brief Computes the overlay of two arrangements with history `arr1` and
`arr2`, and sets the output arrangement with history `res` to
represent the overlaid arrangement. The function also constructs a
consolidated set of curves that induce `res`.

\details
Computes the overlay of two input arrangement
objects, and returns the overlaid arrangement. All three arrangements
can be instantiated with different geometric traits classes and different
\dcel classes (encapsulated in the various topology-traits classes).
The geometry traits of the resulting arrangement is used to construct the
resulting arrangement. This means that all the types (e.g.,
`Traits::Point_2`, `Traits::Curve_2`, and `Traits::Point_2`)
of both input arrangements have to be convertible to the types in the
resulting arrangement. A given overlay-traits object is used to properly
construct the overlaid \dcel that represents the resulting arrangement.

\pre `res` does not refer to either `arr1` or `arr2` (that is, "self overlay" is not supported).

\pre The overlay-traits object `ovl_tr` must model the `OverlayTraits`
  concept, which is able to construct records of the `ResDcel` class on
  the basis of the `Dcel1` and `Dcel2` records that induce them.

\sa `OverlayTraits`

*/
template<typename Traits, typename Dcel1, typename Dcel2,
         typename ResDcel, typename OverlayTraits>
void overlay (const Arrangement_with_history_2<Traits,Dcel1>& arr1,
              const Arrangement_with_history_2<Traits,Dcel2>& arr2,
              Arrangement_with_history_2<Traits,ResDcel>& res,
              OverlayTraits& ovl_tr);




} /* end namesapce CGAL */
