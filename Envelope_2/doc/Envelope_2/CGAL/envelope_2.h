namespace CGAL {

/*!
\ingroup PkgEnvelope2

Computes the lower envelope of a set of curves in \f$ \mathbb{R}^2\f$,
as given by the range `[begin, end)`. The lower envelope is
represented using the output minimization diagram `diag`,
which must be a model of the `EnvelopeDiagram_1` concept.
\pre The value-type of `InputIterator` is `EnvelopeDiagram::Traits_2::Curve_2`.
*/
template<class InputIterator, class EnvelopeDiagram>
void lower_envelope_2 (InputIterator begin, InputIterator end,
EnvelopeDiagram& diag);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgEnvelope2

Computes the lower envelope of a set of \f$ x\f$-monotone curves in 
\f$ \mathbb{R}^2\f$, as given by the range `[begin, end)`. The lower 
envelope is represented using the output minimization diagram `diag`,
which must be a model of the `EnvelopeDiagram_1` concept.
\pre The value-type of `InputIterator` is `EnvelopeDiagram::X_monotone_curve_2`.
*/
template<class InputIterator, class EnvelopeDiagram>
void lower_envelope_x_monotone_2 
(InputIterator begin, InputIterator end,
EnvelopeDiagram& diag);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgEnvelope2

Computes the upper envelope of a set of curves in \f$ \mathbb{R}^2\f$,
as given by the range `[begin, end)`. The upper envelope is
represented using the output maximization diagram `diag`,
which must be a model of the `EnvelopeDiagram_1` concept.
\pre The value-type of `InputIterator` is `EnvelopeDiagram::Traits_2::Curve_2`.
*/
template<class InputIterator, class EnvelopeDiagram>
void upper_envelope_2 (InputIterator begin, InputIterator end,
EnvelopeDiagram& diag);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgEnvelope2

Computes the upper envelope of a set of \f$ x\f$-monotone curves in 
\f$ \mathbb{R}^2\f$, as given by the range `[begin, end)`. The upper 
envelope is represented using the output maximization diagram `diag`,
which must be a model of the `EnvelopeDiagram_1` concept.
\pre The value-type of `InputIterator` is `EnvelopeDiagram::X_monotone_curve_2`.
*/
template<class InputIterator, class EnvelopeDiagram>
void upper_envelope_x_monotone_2 
(InputIterator begin, InputIterator end,
EnvelopeDiagram& diag);

} /* namespace CGAL */
