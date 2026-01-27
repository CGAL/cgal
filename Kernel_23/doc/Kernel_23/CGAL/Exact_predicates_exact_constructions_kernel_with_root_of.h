
namespace CGAL {

/*!
\ingroup kernel_predef

A typedef to a kernel which has the following properties:

<UL>
<LI>It uses %Cartesian representation.
<LI>It supports constructions of points from `double` %Cartesian
coordinates.
<LI>It provides both exact geometric predicates and exact geometric
constructions.
<LI>Its `FT` nested type is model of the `FieldWithRootOf` concept.
</UL>

Note that it requires CORE or LEDA installed.

\cgalModels{Kernel}

\sa `CGAL::Exact_predicates_exact_constructions_kernel`
\sa `CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt`
\sa `CGAL::Exact_predicates_exact_constructions_kernel_with_kth_root`
\sa `CGAL::Exact_predicates_inexact_constructions_kernel`
\sa `CGAL::Cartesian`

*/

class Exact_predicates_exact_constructions_kernel_with_root_of {
public:

}; /* end Exact_predicates_exact_constructions_kernel_with_root_of */
} /* end namespace CGAL */
