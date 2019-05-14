namespace CGAL {

/*!
\ingroup PkgEnvelope3

The traits class `Env_sphere_traits_3` models the `EnvelopeTraits_3` 
concept, and is used for the construction of lower and upper envelopes 
of spheres. Note that when projecting the intersection curve of two 
spheres (a circle in 3D) onto the \f$ xy\f$-plane, the resulting curve is an 
ellipse. The traits class is therefore parameterized by an 
arrangement-traits class that is capable of handling conic 
curves - namely an instantiation of the `Arr_conic_traits_2` 
class-template - and inherits from it. 

The conic-traits class defines a nested type named `Rat_kernel`, 
which is a geometric kernel parameterized by an exact rational type. 
`Env_sphere_traits_3` defines its `Surface_3` type to be constructible from 
`Rat_kernel::Sphere_3`. Namely, it can handle spheres whose center 
points have rational coordinates (i.e., of the type `Rat_kernel::FT`), 
and whose squared radius is also rational. The `Surface_3` type is 
also convertible to a `Rat_kernel::Sphere_3` object. 

The `Xy_monotone_surface_3` type is the same as the nested 
`Surface_3` type. The traits-class simply ignores the upper 
hemisphere when it computes lower envelopes, and ignores the lower 
hemisphere when it computes upper envelopes. 

\cgalModels `EnvelopeTraits_3`

*/
template< typename ConicTraits >
class Env_sphere_traits_3 : public ConicTraits {
public:
}; /* end Env_sphere_traits_3 */
} /* end namespace CGAL */
