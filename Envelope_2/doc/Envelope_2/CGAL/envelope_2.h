namespace CGAL {

/*!
\ingroup PkgEnvelope2Ref

Computes the lower envelope of a set of curves in \f$ \mathbb{R}^2\f$,
as given by the range `[begin, end)`. The lower envelope is
represented using the output minimization diagram `diag`.

\tparam InputIterator must be an input iterator with value type `EnvelopeDiagram::Traits_2::Curve_2`.
\tparam EnvelopeDiagram must be a model of the concept `EnvelopeDiagram_1`.
*/
template<class InputIterator, class EnvelopeDiagram>
void lower_envelope_2 (InputIterator begin, InputIterator end,
EnvelopeDiagram& diag);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgEnvelope2Ref

Computes the lower envelope of a set of \f$ x\f$-monotone curves in
\f$ \mathbb{R}^2\f$, as given by the range `[begin, end)`. The lower
envelope is represented using the output minimization diagram `diag`.

\tparam InputIterator must be an input iterator with value type `EnvelopeDiagram::X_monotone_curve_2`.
\tparam EnvelopeDiagram must be a model of the concept `EnvelopeDiagram_1`.
*/
template<class InputIterator, class EnvelopeDiagram>
void lower_envelope_x_monotone_2
(InputIterator begin, InputIterator end,
EnvelopeDiagram& diag);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgEnvelope2Ref

Computes the lower envelope of a set of \f$ x\f$-monotone curves in
\f$ \mathbb{R}^2\f$, as given by the range `[begin, end)` with the help
of the arrangement traits object `traits` responsible for their creation.
Reusing the same traits object improves speed if the traits class caches data.
The lower envelope is represented using the output minimization diagram `diag`.

\tparam InputIterator must be an input iterator with value type `EnvelopeDiagram::X_monotone_curve_2`.
\tparam EnvelopeDiagram must be a model of the concept `EnvelopeDiagram_1`.
\tparam Traits must be a model of the concept `ArrangementXMonotoneTraits_2`.
*/
template<class InputIterator, class EnvelopeDiagram, class Traits>
void lower_envelope_x_monotone_2
(InputIterator begin, InputIterator end,
EnvelopeDiagram& diag, const Traits& traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgEnvelope2Ref

Computes the upper envelope of a set of curves in \f$ \mathbb{R}^2\f$,
as given by the range `[begin, end)`. The upper envelope is
represented using the output maximization diagram `diag`.

\tparam InputIterator must be an input iterator with value type `EnvelopeDiagram::Traits_2::Curve_2`.
\tparam EnvelopeDiagram must be a model of the concept `EnvelopeDiagram_1`.
*/
template<class InputIterator, class EnvelopeDiagram>
void upper_envelope_2 (InputIterator begin, InputIterator end,
EnvelopeDiagram& diag);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgEnvelope2Ref

Computes the upper envelope of a set of \f$ x\f$-monotone curves in
\f$ \mathbb{R}^2\f$, as given by the range `[begin, end)`. The upper
envelope is represented using the output maximization diagram `diag`.

\tparam InputIterator must be an input iterator with value type `EnvelopeDiagram::X_monotone_curve_2`.
\tparam EnvelopeDiagram must be a model of the concept `EnvelopeDiagram_1`.
*/
template<class InputIterator, class EnvelopeDiagram>
void upper_envelope_x_monotone_2
(InputIterator begin, InputIterator end,
EnvelopeDiagram& diag);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgEnvelope2Ref

Computes the upper envelope of a set of \f$ x\f$-monotone curves in
\f$ \mathbb{R}^2\f$, as given by the range `[begin, end)` with the help
of the arrangement traits object `traits` responsbile for their creation.
Reusing the same traits object improves speed if the traits class caches data.
The upper envelope is represented using the output maximization diagram `diag`.

\tparam InputIterator must be an input iterator with value type `EnvelopeDiagram::X_monotone_curve_2`.
\tparam EnvelopeDiagram must be a model of the concept `EnvelopeDiagram_1`.
\tparam Traits must be a model of the concept `ArrangementXMonotoneTraits_2`.
*/
template<class InputIterator, class EnvelopeDiagram, class Traits>
void upper_envelope_x_monotone_2
(InputIterator begin, InputIterator end,
EnvelopeDiagram& diag, const Traits& traits);

} /* namespace CGAL */
