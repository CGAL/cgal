
namespace CGAL {

/*!
\ingroup PkgSpatialSortingUtils

`Median` is a tag class. It can be used to parameterize a strategy policy
in order to specify the strategy for spatial sorting.
`Hilbert_policy<Median>` can be passed to
as parameter to `hilbert_sort()` to choose the sorting policy.

\cgalModels{DefaultConstructible,CopyConstructible}

\sa `Middle`
\sa `Hilbert_policy`
\sa `Hilbert_sort_median_policy`
\sa `Hilbert_sort_middle_policy`
*/
struct Median { };

/*!
\ingroup PkgSpatialSortingUtils

`Middle` is a tag class. It can be used to parameterize a strategy policy
in order to specify the strategy for spatial sorting.
`Hilbert_policy<Middle>` can be passed to
as parameter to `hilbert_sort()` to choose the sorting policy.

\cgalModels{DefaultConstructible,CopyConstructible}

\sa `Median`
\sa `Hilbert_policy`
\sa `Hilbert_sort_median_policy`
\sa `Hilbert_sort_middle_policy`
*/
struct Middle { };

/*!
\ingroup PkgSpatialSortingUtils

`Hilbert_policy` is a policy class which can be used to parameterize a strategy policy
in order to specify the strategy for spatial sorting.
`Hilbert_policy<Median>` or `Hilbert_policy<Middle>`
can be passed  as parameter to `hilbert_sort()` to choose the sorting policy.

\tparam Tag must be either `Median` or `Middle`.

\cgalModels{DefaultConstructible,CopyConstructible}

\sa `Median`
\sa `Middle`
\sa `Hilbert_sort_median_policy`
\sa `Hilbert_sort_middle_policy`
*/
template< typename Tag >
struct Hilbert_policy { };

/*!
\ingroup PkgSpatialSortingUtils

A typedef to `Hilbert_policy<Median>`.
*/
typedef Hilbert_policy<Median>  Hilbert_sort_median_policy;

/*!
\ingroup PkgSpatialSortingUtils

A typedef to `Hilbert_policy<Middle>`.
*/
typedef Hilbert_policy<Middle>  Hilbert_sort_middle_policy;

} /* end namespace CGAL */
