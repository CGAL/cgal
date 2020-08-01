#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_simple_point_location.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_walk_along_line_point_location.h>

#include <QPointF>

#include "PointLocationFunctions.h"
#include "Utils.h"

BOOST_MPL_HAS_XXX_TRAIT_DEF(Approximate_2)

template <typename Arr_, bool b = has_Approximate_2<Arr_>::value>
struct Supports_landmarks
{
  typedef CGAL::Boolean_tag<b> Tag;
};

template <typename Arr_>
struct Supports_landmarks<Arr_, true>
{
  typedef CGAL::Tag_true Tag;
};

template <typename Arrangement, typename SupportsLandmarks>
struct StrategyHelper
{
  using type = CGAL::Arr_walk_along_line_point_location<Arrangement>;
};

template <typename Arrangement>
struct StrategyHelper<Arrangement, CGAL::Tag_true>
{
  using type = CGAL::Arr_landmarks_point_location<Arrangement>;
};

template <typename Traits>
static auto toKernelPoint(const QPointF& pt)
{
  using Kernel = typename ArrTraitsAdaptor<Traits>::Kernel;
  return CGAL::Qt::Converter<Kernel>{}(pt);
}

template <typename Arr_>
CGAL::Object
PointLocationFunctions<Arr_>::locate(const Arrangement* arr, const QPointF& pt)
{
  using SupportsLandmarks = typename Supports_landmarks<Arrangement>::Tag;
  using PointLocationStrategy =
    typename StrategyHelper<Arrangement, SupportsLandmarks>::type;

  Arr_construct_point_2<Traits> toArrPoint;
  auto arr_point = toArrPoint(toKernelPoint<Traits>(pt));

  PointLocationStrategy pointLocationStrategy{*arr};
  return pointLocationStrategy.locate(arr_point);
}

template <typename Arr_>
auto PointLocationFunctions<Arr_>::getFace(
  const Arrangement* arr, const QPointF& pt) -> Face_const_handle
{
  CGAL::Object obj = locate(arr, pt);

  Face_const_handle f;
  if (CGAL::assign(f, obj)) return f;

  Halfedge_const_handle he;
  if (CGAL::assign(he, obj)) return (he->face());

  Vertex_const_handle v;
  CGAL_assertion(CGAL::assign(v, obj));

  CGAL::assign(v, obj);
  if (v->is_isolated()) return v->face();
  Halfedge_around_vertex_const_circulator eit = v->incident_halfedges();
  return (eit->face());
}

template <typename Arr_>
CGAL::Object PointLocationFunctions<Arr_>::rayShootUp(
  const Arrangement* arr, const QPointF& pt)
{
  using Walk_pl_strategy =
    typename CGAL::Arr_walk_along_line_point_location<Arrangement>;

  Arr_construct_point_2<Traits> toArrPoint;
  Walk_pl_strategy pointLocationStrategy{*arr};
  return pointLocationStrategy.ray_shoot_up(
    toArrPoint(toKernelPoint<Traits>(pt)));
}

template <typename Arr_>
CGAL::Object PointLocationFunctions<Arr_>::rayShootDown(
  const Arrangement* arr, const QPointF& pt)
{
  using Walk_pl_strategy =
    typename CGAL::Arr_walk_along_line_point_location<Arrangement>;

  Arr_construct_point_2<Traits> toArrPoint;
  Walk_pl_strategy pointLocationStrategy{*arr};
  return pointLocationStrategy.ray_shoot_down(
    toArrPoint(toKernelPoint<Traits>(pt)));
}

template class PointLocationFunctions<Seg_arr>;
template class PointLocationFunctions<Pol_arr>;
template class PointLocationFunctions<Conic_arr>;
template class PointLocationFunctions<Lin_arr>;
template class PointLocationFunctions<Alg_seg_arr>;
template class PointLocationFunctions<Bezier_arr>;
