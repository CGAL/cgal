#ifndef CGAL_HDVF_TRAITS_3_H
#define CGAL_HDVF_TRAITS_3_H

#include <CGAL/license/HDVF.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Dimension.h>

namespace CGAL {
namespace Homological_discrete_vector_field {

  /*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `HDVF_traits_3` implements the `HDVFTraits` concept for 3D data, using a geometric kernel `K`.

 @tparam K a geometric kernel model of the `Kernel` concept.

 \cgalModels{HDVFTraits}

 */

template <typename K>
struct Hdvf_traits_3 {
    using Dimension = Dimension_tag< 3 >;
    typedef K Kernel;
    typedef typename K::Point_3 Point;
    typedef typename K::Vector_3 Vector;
    typedef typename K::FT FT;
    typedef CGAL::Bbox_3 Bbox;
};

} /* end namespace Homological_discrete_vector_field */
} /* end namespace CGAL */

#endif // CGAL_HDVF_TRAITS_3_H
