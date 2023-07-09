
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2TraitsClasses
\ingroup PkgArrangementOnSurface2Overlay

An instance of `Arr_default_overlay_traits` should be used for overlaying two arrangements
of type `Arrangement` that store no auxiliary data with their \dcel records, where the resulting overlaid arrangement stores no auxiliary
\dcel data as well. This class simply gives empty implementation for all
traits-class functions.

\cgalModels{OverlayTraits}

\sa `overlay`

*/
template< typename Arrangement >
class Arr_default_overlay_traits {
public:

}; /* end Arr_default_overlay_traits */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2TraitsClasses
\ingroup PkgArrangementOnSurface2Overlay

An instance of `Arr_face_overlay_traits` should be used for overlaying two arrangements
of types `Arr_A` and `Arr_B`, which are instantiated using the same
geometric traits-class and with the \dcel classes `Dcel_A` and
`Dcel_B` respectively, in order to store their overlay in an arrangement
of type `Arr_R`, which is instantiated using a third \dcel class
`Dcel_R`. All three \dcel classes are assumed to be instantiations of the
`Arr_face_extended_dcel` template with types `FaceData_A`,
`FaceData_B` and `FaceData_R`, respectively.

This class gives empty implementation for all overlay traits-class functions,
except the function that computes the overlay of two faces. In this case,
it uses the functor `OvlFaceData`, which accepts a `FaceData_A` object
and a `FaceData_B` object and computes a corresponding `FaceData_R`
object, in order to set the auxiliary data of the overlay face.

\cgalModels{OverlayTraits}

\sa `overlay`
\sa `CGAL::Arr_face_extended_dcel<Traits,FData,V,H,F>`

*/
template< typename Arr_A, typename Arr_B, typename Arr_R, typename OvlFaceData >
class Arr_face_overlay_traits {
public:

}; /* end Arr_face_overlay_traits */
} /* end namespace CGAL */
