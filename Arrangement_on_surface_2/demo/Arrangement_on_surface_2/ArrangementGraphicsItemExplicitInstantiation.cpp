
#include "ArrangementGraphicsItem.cpp"
#include "ArrangementTypes.h"

namespace CGAL {
namespace Qt {


template class ArrangementGraphicsItem<Seg_arr>;

template class ArrangementGraphicsItem<Pol_arr>;

template class ArrangementGraphicsItem<Conic_arr>;

template class ArrangementGraphicsItem<Lin_arr>;

template class ArrangementGraphicsItem<Arc_arr>;

template class ArrangementGraphicsItem<Alg_seg_arr>;


} // namespace Qt
} // namespace CGAL