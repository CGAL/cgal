#ifndef GENERATORS_H
#define GENERATORS_H

#include <CGAL/Random.h>
#include <CGAL/Point_with_normal_3.h>
#include <string>

template <typename K>
typename K::FT random_float(const double &min, const double &max) {
  static CGAL::Random rand;
  return K::FT(rand.get_double(min, max));
}

template <typename K>
typename CGAL::Vector_3<K> random_normal() {
  CGAL::Vector_3<K> n;
  do {
    n = CGAL::Vector_3<K>(random_float<K>(-1.0, 1.0),
               random_float<K>(-1.0, 1.0),
               random_float<K>(-1.0, 1.0));
  } while (n.squared_length() < 0.001);

  n = n * 1.0 / (CGAL::sqrt(n.squared_length()));

  return n;
}

template <class K>
typename CGAL::Point_3<K> random_point_in(const CGAL::Bbox_3& bbox)
{
  typedef typename K::FT FT;

  FT x = random_float<K>(bbox.xmin(), bbox.xmax());
  FT y = random_float<K>(bbox.ymin(), bbox.ymax());
  FT z = random_float<K>(bbox.zmin(), bbox.zmax());

  return typename CGAL::Point_3<K>(x, y, z);
}

template <class K>
typename CGAL::Point_with_normal_3<K> random_pwn_in(const CGAL::Bbox_3 &bbox) {
  return typename CGAL::Point_with_normal_3<K>(random_point_in<K>(bbox), random_normal<K>());
}

template <class K>
typename CGAL::Vector_3<K> normalize(typename const CGAL::Vector_3<K> &v) {
  K::FT l = CGAL::Sqrt(v.squared_length());
  if (l < 0.00001)
    return typename CGAL::Vector_3<K>(0, 0, 0);
  else return v * l;
}

template <typename K, typename OutputIterator>
void generatePointsOnSphere(const std::size_t num_points,
  typename const CGAL::Point_3<K> &center, typename K::FT radius,
  OutputIterator points)
{
  for (std::size_t i = 0;i < num_points;++i)
  {
    CGAL::Vector_3<K> direction(random_float(-1.f, 1.f), 
      random_float(-1.f, 1.f), 
      random_float(-1.f, 1.f));
    direction = direction * (1.0 / sqrt(direction.squared_length()));

    CGAL::Point_3<K> p = center + direction * radius;

    *points = CGAL::Point_with_normal_3<K>(p, direction);
    ++points;
  }
}

template <typename K, typename OutputIterator>
void sampleCylinderInBox(const std::size_t num_points,
  const CGAL::Bbox_3 &bbox, const typename CGAL::Point_3<K> &center,
  const typename CGAL::Vector_3<K> &axis, const typename K::FT radius,
  const typename K::FT length, OutputIterator points) {
     // Sample shape
    for (size_t i = 0 ; i < num_points ; ++i) {
      CGAL::Vector_3<K> normal(random_float<K>(-1.f, 1.f), 
        random_float<K>(-1.f, 1.f), 
        random_float<K>(-1.f, 1.f));
      normal = normal - ((normal * axis) * axis);
      if (normal.squared_length() < 0.0001) {
        i--;
        continue;
      }
      normal = normal * (1.0 / sqrt(normal.squared_length()));

      K::FT l = random_float<K>(-length, length);

      CGAL::Point_3<K> p = center + axis * l + radius * normal;
      *points = CGAL::Point_with_normal_3<K>(p, normal);
      ++points;
    }
}

template <typename K, typename OutputIterator>
void sampleRandomCylinderInBox(const std::size_t num_points,
  const CGAL::Bbox_3 &bbox, typename CGAL::Point_3<K> &center,
  typename CGAL::Vector_3<K> &axis, typename K::FT &radius,
  OutputIterator points) {
    // Generate random parameters
    axis = random_normal<K>();
    radius = random_float<K>(0.5, 5.0);
    K::FT length = random_float<K>(0.2, 10);

    // Find random center point placed on the plane through
    // the origin with 'axis' as normal.
    CGAL::Vector_3<K> u = random_normal<K>();
    CGAL::Vector_3<K> v = CGAL::cross_product(u, axis);
    v = v * 1.0 / CGAL::sqrt(v.squared_length());
    u = CGAL::cross_product(v, axis);

    center = CGAL::ORIGIN + random_float<K>(-5.0, 5.0) * u + random_float<K>(-5.0, 5.0) * v;

    sampleCylinderInBox(num_points, bbox, center, axis, radius, length, points);
}

template <typename K, typename OutputIterator>
void generatePointsOnCone(const std::size_t num_points, typename const CGAL::Point_3<K> &apex, typename const CGAL::Vector_3<K> &axis, typename K::FT angle, typename K::FT s, typename K::FT e, OutputIterator points)
{
  assert(s < e);
  K::FT radiusGrow = CGAL::tan(angle);
  axis = axis * 1.0 / (CGAL::sqrt(axis.squared_length()));

  K::FT cosAng = CGAL::cos(angle);
  K::FT sinAng = CGAL::sin(angle);

  for (size_t i = 0 ; i < num_points ; ++i)
  {
    Vector normal(random_float(-1.f, 1.f), 
      random_float(-1.f, 1.f), 
      random_float(-1.f, 1.f));
    normal = normal - ((normal * axis) * axis);
    if (normal.squared_length() < 0.0001) {
      i--;
      continue;
    }
    normal = normal * (1.0 / CGAL::sqrt(normal.squared_length()));

    K::FT l = random_float<K>(s, e);

    CGAL::Point_3<K> p = apex + axis * l + (l * radiusGrow) * normal;

    // axis is pointing down from apex
    normal = normal * 1.0 / (CGAL::sqrt(normal.squared_length()));
    // normal is pointing from axis to surface point
    normal = normal * cosAng - axis * sinAng;
    l = CGAL::sqrt(normal.squared_length());
    if (0.95 > l || l > 1.05)
      std::cout << "normal not normalized" << std::endl;

    *points = CGAL::Point_with_normal_3<K>(p, normal);
    ++points;
  }
}

template <typename K, typename OutputIterator>
void sampleRandomParallelogramInBox(const std::size_t num_points, const CGAL::Bbox_3 &bbox, typename CGAL::Vector_3<K> &normal, typename K::FT &dist, OutputIterator points) {
  // Generate random plane from 3 non collinear points.
  CGAL::Vector_3<K> u, v;
  CGAL::Point_3<K> p[3];
  do {
    p[0] = random_point_in<K>(bbox);
    p[1] = random_point_in<K>(bbox);
    p[2] = random_point_in<K>(bbox);

    CGAL::Vector_3<K> a = p[1] - p[0];
    CGAL::Vector_3<K> b = p[2] - p[0];

    if (a.squared_length() < 4.0 || b.squared_length() < 4.0) {
      normal = CGAL::Vector_3<K>(0, 0, 0);
      continue;
    }

    a = a * 1.0 / CGAL::sqrt(a.squared_length());
    b = b * 1.0 / CGAL::sqrt(b.squared_length());

    normal = CGAL::cross_product(a, b);
    // repeat if angle between a and b is small
  } while (normal.squared_length() < 0.2);

  normal = normal * 1.0 / CGAL::sqrt(normal.squared_length());
  dist = -((p[0] - CGAL::ORIGIN) * normal);

  // sample
  size_t i = 0;
  u = p[1] - p[0];
  v = p[2] - p[0];
  
  for (std::size_t i = 0;i < num_points; ++i) {
    double s = random_float<K>(0, 1.0);
    double t = random_float<K>(0, 1.0);

    *points = CGAL::Point_with_normal_3<K>(p[0] + s * u + t * v, normal);
    ++points;
  }
}

template <typename K, typename OutputIterator>
void generatePointsOnTorus(const std::size_t num_points, typename const CGAL::Point_3<K> &center, typename const CGAL::Vector_3<K> &axis, typename const K::FT majorRadius, typename const K::FT minorRadius, OutputIterator points)
{
  size_t i = 0;
  axis = axis / CGAL::sqrt(axis.squared_length());

  // calculate basis
  Vector b1, b2;
  b1 = CGAL::cross_product(axis, Vector(0, 0, 1));
  if (b1.squared_length() < 0.1)
    b1 = CGAL::cross_product(axis, Vector(0, 1, 0));
  b1 = b1 / sqrt(b1.squared_length());

  b2 = CGAL::cross_product(axis, b1);
  b2 = b2 / sqrt(b2.squared_length());

  for (size_t i = 0 ; i < num_points ; ++i)
  {
    FT tau = random_float<K>(0, 2 * M_PI);
    FT phi = random_float<K>(0, 2 * M_PI);

    Vector normal = sin(tau) * b1 + cos(tau) * b2;
    normal = normal / sqrt(normal.squared_length());
    Point p = center + majorRadius * normal;
    normal = sin(phi) * normal + cos(phi) * axis;
    p = p + minorRadius * normal;

    *points = Point_with_normal(p, normal);
    ++points;
  }
}

template <typename K>
void saveScene(const std::string &fn, const std::vector<CGAL::Point_with_normal_3<K>> &pts) {
  std::ofstream plyFile(fn);

  plyFile << "ply" << std::endl;
  plyFile << "format ascii 1.0" << std::endl;
  plyFile << "element vertex " << pts.size() << std::endl;
  plyFile << "property float x" << std::endl;
  plyFile << "property float y" << std::endl;
  plyFile << "property float z" << std::endl;
  plyFile << "property float nx" << std::endl;
  plyFile << "property float ny" << std::endl;
  plyFile << "property float nz" << std::endl;
  plyFile << "end_header" << std::endl;

  plyFile << std::setprecision(6);

  for (size_t i = 0;i<pts.size();i++) {
    plyFile << pts[i][0] << " " << pts[i][1] << " " << pts[i][2];
    CGAL::Vector_3<K> n = pts[i].normal();
    plyFile << " " << n[0] << " " << n[1] << " " << n[2];
    plyFile << std::endl;
  }

  plyFile.close();
}

#endif
