#ifndef CGAL_DRAW_AOS_ARR_RENDER_CONTEXT_H
#define CGAL_DRAW_AOS_ARR_RENDER_CONTEXT_H
#include <cstdlib>
#include <limits>
#include <memory>
#include <algorithm>
#include <atomic>
#include <chrono>

#include <CGAL/Bbox_2.h>
#include <CGAL/Arr_point_location_result.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Draw_aos/Arr_approximate_point_2.h>
#include <CGAL/Draw_aos/Arr_approximation_cache.h>
#include <CGAL/Draw_aos/Arr_construct_curve_end.h>
#include <CGAL/Draw_aos/Arr_construct_segments.h>
#include <CGAL/Draw_aos/Arr_portals.h>
#include <CGAL/Draw_aos/type_utils.h>

#if defined(CGAL_DRAW_AOS_DEBUG)
#include <fstream>
#endif

namespace CGAL {
namespace draw_aos {

class Arr_cancellable_context_mixin
{
  using Clock = std::chrono::steady_clock;
  using Duration = Clock::duration;
  using Time_point = std::chrono::time_point<Clock, Duration>;

protected:
  Arr_cancellable_context_mixin()
      : m_start_time(Clock::now())
      , m_done(std::make_shared<std::atomic<bool>>(false)) {}

public:
  Time_point start_time() const { return m_start_time; }

  Time_point end_time() const { return m_end_time; }

  Duration elapsed_time() const { return Clock::now() - m_start_time; }

  bool is_cancelled() const { return m_done->load(); }

  void cancel() {
    m_done->store(true, std::memory_order_relaxed);
    m_end_time = Clock::now();
  }

private:
  Time_point m_start_time, m_end_time;
  std::shared_ptr<std::atomic<bool>> m_done;
};

template <typename GeomTraits>
class Arr_bounds_context_mixin
{
  using Approx_point = typename Arr_approximation_geometry_traits<GeomTraits>::Approx_point;
  using Point_2 = typename Traits_adaptor<GeomTraits>::Point_2;
  using FT = typename Traits_adaptor<GeomTraits>::FT;

protected:
  Arr_bounds_context_mixin(const Bbox_2& bbox)
      : m_bbox(bbox) {}

public:
  double xmin() const { return m_bbox.xmin(); }
  double xmax() const { return m_bbox.xmax(); }
  double ymin() const { return m_bbox.ymin(); }
  double ymax() const { return m_bbox.ymax(); }

  const Bbox_2& bbox() const { return m_bbox; }

  bool strictly_contains_x(double x) const { return xmin() < x && x <= xmax(); }
  bool strictly_contains_x(FT x) const { return to_ft(xmin()) < x && x <= to_ft(xmax()); }

  bool strictly_contains_y(double y) const { return ymin() < y && y <= ymax(); }
  bool strictly_contains_y(FT y) const { return to_ft(ymin()) < y && y <= to_ft(ymax()); }

  bool strictly_contains(Point_2 pt) const { return strictly_contains_x(pt.x()) && strictly_contains_y(pt.y()); }
  bool strictly_contains(Approx_point pt) const { return strictly_contains_x(pt.x()) && strictly_contains_y(pt.y()); }

  bool contains_x(double x) const { return xmin() <= x && x <= xmax(); }
  bool contains_x(FT x) const { return to_ft(xmin()) <= x && x <= to_ft(xmax()); }

  bool contains_y(double y) const { return ymin() <= y && y <= ymax(); }
  bool contains_y(FT y) const { return to_ft(ymin()) <= y && y <= to_ft(ymax()); }

  bool contains(Approx_point pt) const { return contains_x(pt.x()) && contains_y(pt.y()); }
  bool contains(Point_2 pt) const { return contains_x(pt.x()) && contains_y(pt.y()); }

  bool is_on_boundary(Approx_point pt) const {
    return (pt.x() == xmin() || pt.x() == xmax()) && contains_y(pt.y()) ||
           (pt.y() == ymin() || pt.y() == ymax()) && contains_x(pt.x());
  }
  bool is_on_boundary(Point_2 pt) const {
    return (pt.x() == to_ft(xmin()) || pt.x() == to_ft(xmax())) && contains_y(pt.y()) ||
           (pt.y() == to_ft(ymin()) || pt.y() == to_ft(ymax())) && contains_x(pt.x());
  }

private:
  const Bbox_2 m_bbox;
  const Construct_coordinate<GeomTraits> to_ft;
};

template <typename GeomTraits>
class Arr_geom_traits_context_mixin
{
public:
  Arr_geom_traits_context_mixin(const GeomTraits& _traits)
      : traits(_traits)
      , cst_curve_end(_traits)
      , cst_horizontal_segment(_traits)
      , cst_vertical_segment(_traits)
      , intersect_2(_traits.intersect_2_object())
      , compare_xy_2(_traits.compare_xy_2_object())
      , is_vertical_2(_traits.is_vertical_2_object())
      , approx_pt(_traits) {}

  const GeomTraits& traits;
  const Arr_construct_curve_end<GeomTraits> cst_curve_end;
  const Arr_construct_vertical_segment<GeomTraits> cst_vertical_segment;
  const Arr_construct_horizontal_segment<GeomTraits> cst_horizontal_segment;
  const typename Traits_adaptor<GeomTraits>::Intersect_2 intersect_2;
  const typename Traits_adaptor<GeomTraits>::Compare_xy_2 compare_xy_2;
  const typename Traits_adaptor<GeomTraits>::Is_vertical_2 is_vertical_2;
  const Arr_approximate_point_2<GeomTraits> approx_pt;
};

template <typename Arrangement>
class Arr_render_context : public Arr_cancellable_context_mixin,
                           public Arr_geom_traits_context_mixin<typename Arrangement::Geometry_traits_2>
{
  using Point_location = Arr_trapezoid_ric_point_location<Arrangement>;
  using Feature_portals_map = typename Arr_portals<Arrangement>::Feature_portals_map;
  using Cancellable_context_mixin = Arr_cancellable_context_mixin;
  using Geom_traits_context_mixin = Arr_geom_traits_context_mixin<typename Arrangement::Geometry_traits_2>;

public:
  Arr_render_context(const Arrangement& arr,
                     const Point_location& pl,
                     const Feature_portals_map& feature_portals,
                     double approx_error)
      : Cancellable_context_mixin()
      , Geom_traits_context_mixin(*arr.geometry_traits())
      , arr(arr)
      , point_location(pl)
      , feature_portals(feature_portals)
      , approx_error(approx_error) {
#if defined(CGAL_DRAW_AOS_DEBUG) && defined(CGAL_DRAW_AOS_TRIANGULATOR_DEBUG_FILE_DIR)
    std::filesystem::path debug_file_dir(CGAL_DRAW_AOS_TRIANGULATOR_DEBUG_FILE_DIR);
    // clear the index file.
    std::filesystem::remove(debug_file_dir / "index.txt");
#endif
  }

public:
  const double approx_error;
  const Arrangement& arr;
  const Point_location& point_location;
  const Feature_portals_map& feature_portals;

#if defined(CGAL_DRAW_AOS_DEBUG)
  std::shared_ptr<int> debug_counter = std::make_shared<int>(0);
#endif
};

template <typename Arrangement>
class Arr_bounded_render_context : public Arr_render_context<Arrangement>,
                                   public Arr_bounds_context_mixin<typename Arrangement::Geometry_traits_2>
{
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_point = typename Arr_approximation_geometry_traits<Geom_traits>::Approx_point;
  using Point_2 = typename Traits_adaptor<Geom_traits>::Point_2;
  using Render_context = Arr_render_context<Arrangement>;
  using Bounds_context_mixin = Arr_bounds_context_mixin<Geom_traits>;
  using Approx_cache = Arr_approximation_cache<Arrangement>;

  constexpr static double ep_base = std::numeric_limits<double>::epsilon();

public:
  Arr_bounded_render_context(const Render_context& ctx, const Bbox_2& bbox, Approx_cache& cache)
      : Render_context(ctx)
      , Bounds_context_mixin(bbox)
      , ep_xmin(std::max(std::abs(ep_base * bbox.xmin()), ep_base))
      , ep_xmax(std::max(std::abs(ep_base * bbox.xmax()), ep_base))
      , ep_ymin(std::max(std::abs(ep_base * bbox.ymin()), ep_base))
      , ep_ymax(std::max(std::abs(ep_base * bbox.ymax()), ep_base))
      , cache(cache) {}

  Approx_point make_on_boundary(const Approx_point& pt) const {
    double x = pt.x(), y = pt.y();

    if(std::abs(x - this->xmin()) < ep_xmin) {
      x = this->xmin();
    } else if(std::abs(x - this->xmax()) < ep_xmax) {
      x = this->xmax();
    } else if(std::abs(y - this->ymin()) < ep_ymin) {
      y = this->ymin();
    } else if(std::abs(y - this->ymax()) < ep_ymax) {
      y = this->ymax();
    } else {
      // We shall not call this function if not approximated from a boundary point.
      CGAL_assertion(false && "Failed to match point to the boundary");
    }

    return Approx_point(x, y);
  }

public:
  Approx_cache& cache;

private:
  // floating point epsilon at boundary coordinates
  double ep_xmin, ep_xmax, ep_ymin, ep_ymax;
};

template <typename Context>
class Arr_context_delegator
{
public:
  Arr_context_delegator(const Context& ctx)
      : m_ctx(ctx) {}

  const Context* operator->() const { return &m_ctx; }
  const Context& operator*() const { return m_ctx; }

  // implicit conversion to the context type
  operator const Context&() const { return m_ctx; }

private:
  const Context& m_ctx;
};

} // namespace draw_aos
} // namespace CGAL
#endif // CGAL_DRAW_AOS_ARR_RENDER_CONTEXT_H