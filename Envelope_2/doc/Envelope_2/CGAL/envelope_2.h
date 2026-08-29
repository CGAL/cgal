namespace CGAL {
namespace Envelope_2 {

/*! \ingroup PkgEnvelope2Ref
 *
 * computes the lower envelope of a set of curves in \f$\mathbb{R}^2\f$, given
 * by the range `[begin, end)`. The lower envelope is represented using the
 * output minimization diagram `diag`.
 *
 * \tparam InputIterator must be an input iterator with value type `EnvelopeDiagram::Traits_2::Curve_2`.
 * \tparam EnvelopeDiagram must be a model of the concept `EnvelopeDiagram_1`.
 */
template <typename InputIterator, typename EnvelopeDiagram>
void lower_envelope_2(InputIterator begin, InputIterator end, EnvelopeDiagram& diag);

} // namespace Envelope_2
} // namespace CGAL

namespace CGAL {
namespace Envelope_2 {

/*! \ingroup PkgEnvelope2Ref
 *
 * computes the lower envelope of a set of \f$x\f$-monotone curves in
 * \f$\mathbb{R}^2\f$, given by the range `[begin, end)`. The lower envelope is
 * represented using the output minimization diagram `diag`.
 *
 * \tparam InputIterator must be an input iterator with value type `EnvelopeDiagram::X_monotone_curve_2`.
 * \tparam EnvelopeDiagram must be a model of the concept `EnvelopeDiagram_1`.
 */
template <typename InputIterator, typename EnvelopeDiagram>
void lower_envelope_x_monotone_2(InputIterator begin, InputIterator end, EnvelopeDiagram& diag);

} // namespace Envelope_2
} // namespace CGAL

namespace CGAL {
namespace Envelope_2 {

/*! \ingroup PkgEnvelope2Ref
 *
 * computes the lower envelope of a set of \f$x\f$-monotone curves in
 * \f$\mathbb{R}^2\f$, given by the range `[begin, end)` with the help of the
 * arrangement traits object `traits` responsible for their creation.  Reusing
 * the same traits object improves speed if the traits class caches data.  The
 * lower envelope is represented using the output minimization diagram `diag`.
 *
 * \tparam InputIterator must be an input iterator with value type `EnvelopeDiagram::X_monotone_curve_2`.
 * \tparam EnvelopeDiagram must be a model of the concept `EnvelopeDiagram_1`.
 * \tparam Traits must be a model of the concept `AosXMonotoneTraits_2`.
 */
template <typename InputIterator, typename EnvelopeDiagram, typename Traits>
void lower_envelope_x_monotone_2(InputIterator begin, InputIterator end, EnvelopeDiagram& diag, const Traits& traits);

} // namespace Envelope_2
} // namespace CGAL

namespace CGAL {
namespace Envelope_2 {

/*! \ingroup PkgEnvelope2Ref
 *
 * computes the upper envelope of a set of curves in \f$\mathbb{R}^2\f$, given
 * by the range `[begin, end)`. The upper envelope is represented using the
 * output maximization diagram `diag`.
 *
 * \tparam InputIterator must be an input iterator with value type `EnvelopeDiagram::Traits_2::Curve_2`.
 * \tparam EnvelopeDiagram must be a model of the concept `EnvelopeDiagram_1`.
 */
template <typename InputIterator, typename EnvelopeDiagram>
void upper_envelope_2(InputIterator begin, InputIterator end, EnvelopeDiagram& diag);

} // namespace Envelope_2
} // namespace CGAL

namespace CGAL {
namespace Envelope_2 {

/*! \ingroup PkgEnvelope2Ref
 *
 * computes the upper envelope of a set of\f$x\f$-monotone curves in
 * \f$\mathbb{R}^2\f$, given by the range `[begin, end)`. The upper envelope is
 * represented using the output maximization diagram `diag`.
 *
 * \tparam InputIterator must be an input iterator with value type `EnvelopeDiagram::X_monotone_curve_2`.
 * \tparam EnvelopeDiagram must be a model of the concept `EnvelopeDiagram_1`.
 */
template <typename InputIterator, typename EnvelopeDiagram>
void upper_envelope_x_monotone_2(InputIterator begin, InputIterator end, EnvelopeDiagram& diag);

} // namespace Envelope_2
} // namespace CGAL

namespace CGAL {
namespace Envelope_2 {

/*! \ingroup PkgEnvelope2Ref
 *
 * computes the upper envelope of a set of \f$x\f$-monotone curves in
 * \f$\mathbb{R}^2\f$, given by the range `[begin, end)` with the help of the
 * arrangement traits object `traits` responsible for their creation.  Reusing
 * the same traits object improves speed if the traits class caches data.  The
 * upper envelope is represented using the output maximization diagram `diag`.
 *
 * \tparam InputIterator must be an input iterator with value type `EnvelopeDiagram::X_monotone_curve_2`.
 * \tparam EnvelopeDiagram must be a model of the concept `EnvelopeDiagram_1`.
 * \tparam Traits must be a model of the concept `AosXMonotoneTraits_2`.
 */
template <typename InputIterator, typename EnvelopeDiagram, typename Traits>
void upper_envelope_x_monotone_2(InputIterator begin, InputIterator end, EnvelopeDiagram& diag, const Traits& traits);

} // namespace Envelope_2
} // namespace CGAL
