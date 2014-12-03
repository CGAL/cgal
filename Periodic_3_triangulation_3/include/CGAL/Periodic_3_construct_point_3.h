#ifndef CGAL_PERIODIC_3_CONSTRUCT_POINT_3_H
#define CGAL_PERIODIC_3_CONSTRUCT_POINT_3_H

namespace CGAL
{
template < typename K, typename Construct_point_3_base>
class Periodic_3_construct_point_3 : public Construct_point_3_base
{
  typedef K Kernel;

public:
  typedef typename Kernel::Point_3       Point;
  typedef typename Kernel::Offset        Offset;
  typedef typename Kernel::Iso_cuboid_3  Iso_cuboid_3;

  typedef Point       result_type;

  Periodic_3_construct_point_3(const Iso_cuboid_3 & dom) : _dom(dom) { }

  Point operator() ( const Point& p, const Offset& o ) const {
    return Point(p.x()+(_dom.xmax()-_dom.xmin())*o.x(),
  p.y()+(_dom.ymax()-_dom.ymin())*o.y(),
  p.z()+(_dom.zmax()-_dom.zmin())*o.z());
  }

private:
  Iso_cuboid_3 _dom;
};
}

#endif
