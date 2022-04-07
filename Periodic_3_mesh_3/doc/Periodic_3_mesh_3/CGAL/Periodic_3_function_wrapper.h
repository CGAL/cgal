namespace CGAL {

/*!
\ingroup PkgPeriodic3Mesh3Domains

The class `Periodic_3_function_wrapper` is a helper class designed to wrap
an (a priori non-periodic) implicit function describing a domain through the relationship
[`p` is inside if `f(p)<0`] and defined over the canonical cube to a function
defined over the whole Euclidean space and periodic, with the same period as the canonical cube.

More precisely, if `f` is the real function defined either over \f$ \mathbb R^3\f$
or over the canonical cube, we construct the periodic real function \f$ f^{\ast} \f$ defined over
\f$ \mathbb R^3\f$ as follows:
- For any point \f$(x,y,z)\f$ in the canonical cube, \f$ f^{\ast}(x,y,z) = f(x,y,z)\f$
- For any point \f$(x,y,z)\f$ outside the canonical cube, there is a unique canonical representative
\f$(x_0,y_0,z_0)\f$ of \f$(x,y,z)\f$ in the canonical cube, i.e.,
\f$(x,y,z)=(x_0 + k_x c, y_0 + k_y c, z_0 + k_z c)\f$ with \f$(k_x,k_y,k_z)\f$ in \f$ \mathbb Z^3\f$,
and \f$ f^{\ast}(x,y,z) = f(x_0,y_0,z_0) \f$.

For example, if considering the unit cube as canonical cube, an oracle answering a
query such as <I>"what is the value of the implicit function at this point?"</I>
at the point `(2.5, 2.5, 2.5)` will be in fact evaluated at the canonical representative, that is
`(0.5, 0.5, 0.5)`.
Consequently, it is then not required to provide an input domain that is defined over the whole
space or periodic, but only defined over the canonical cube.

\cgalFigureBegin{Periodic_3_mesh_3FromCanonicalToWhole, periodicity_base.svg}
Illustration in 2D (cut view) of a domain defined by an implicit function that is transformed
into a periodic implicit function.
Only the values of the implicit function that are in the canonical cube are used:
the values of the implicit function at \f$ P \f$ and \f$ Q \f$ are obtained by evaluating
instead at \f$ P' \f$ and \f$ Q' \f$, as shown on the right.
\cgalFigureEnd

In practice, the implicit function provided by the user is likely defined
over a larger domain than the canonical cube (in general, it is \f$ \mathbb R^3\f$).
Note that -- when constructing artificially periodic functions -- all the values of the implicit function
for points outside this canonical cube are unused since queries are always answered by looking at the canonical representative.
\cgalFigureRef{Periodic_3_mesh_3FromCanonicalToWholeDiscard} gives an example of such domain where some information is discarded.

\cgalFigureBegin{Periodic_3_mesh_3FromCanonicalToWholeDiscard, periodicity.svg}
Illustration in 2D (cut view) of a domain defined by an implicit function artificially made periodic.
Any value of the function outside of the canonical cube is ignored.
\cgalFigureEnd

Note also that when constructing artificially periodic functions, it is the responsability of the user
to provide an input function that is compatible with the canonical cube (that is, whose isovalues
are <em>periodically</em> continuous and without intersections).
\cgalFigureRef{Periodic_3_mesh_3ContinuityIssue} is an example of a bad choice
of input function and canonical cube: there is no continuity of the isovalues
at the border of the canonical cube. In such configuration, the mesher might
or might not finish and the result is likely to be non-manifold and to contain self-intersections.

\cgalFigureBegin{Periodic_3_mesh_3ContinuityIssue, periodicity_issue.svg}
Illustration in 2D (cut view) of a domain defined by an implicit function artificially made periodic.
The zero isovalue of the implicit function does not form a continuous curve.
\cgalFigureEnd

\tparam Function provides the definition of the function.
        This parameter stands for a model of the concept `ImplicitFunction`
        described in the surface mesh generation package.
        The number types `Function::FT` and `BGT::FT` are required to match.
\tparam BGT is a geometric traits class that provides
        the basic operations to implement
        intersection tests and intersection computations
        through a bisection method. This parameter must be instantiated
        with a model of the concept `BisectionGeometricTraits_3`.

*/
template<class Function, class BGT>
class Periodic_3_function_wrapper
{
public:
  /// \name Types
  /// @{
  //!
  typedef typename BGT::FT                                     FT;
  typedef typename BGT::Point_3                                Point_3;
  typedef typename BGT::Iso_cuboid_3                           Iso_cuboid_3;
  /// @}

  /// \name Creation
  /// @{
  /*!
   * \brief Construction from an implicit function and the canonical cube.
   */
  Periodic_3_function_wrapper(Function f, const Iso_cuboid_3& domain);
  /// @}

  /// \name Operations
  /// @{

  /*!
   * Evaluates the function \f$ f \f$ passed in input at the canonical representative of \f$ p \f$.
   */
  FT operator()(const Point_3& p) const;

  /// @}
}; /* end Implicit_to_labeled_subdomains_function_wrapper */

} /* end namespace CGAL */
