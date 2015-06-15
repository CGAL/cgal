#ifndef GENERATORS_H
#define GENERATORS_H

#include <CGAL/Random.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Kernel/global_functions.h>
#include <string>
#include <fstream>


template <typename fl_t>
fl_t random_float(fl_t min, fl_t max) {
  static CGAL::Random rand;
  return fl_t(rand.get_double(min, max));
}

template <typename K>
typename CGAL::Vector_3<K> random_normal() {
  CGAL::Vector_3<K> n;
  do {
    n = CGAL::Vector_3<K>(random_float((K::FT) -1.0, (K::FT) 1.0),
               random_float((K::FT) -1.0, (K::FT) 1.0),
               random_float((K::FT) -1.0, (K::FT) 1.0));
  } while (n.squared_length() < (K::FT) 0.001);

  n = n * (K::FT) 1.0 / (CGAL::sqrt(n.squared_length()));

  return n;
}

template <class K>
typename CGAL::Point_3<K> random_point_in(const CGAL::Bbox_3& bbox) {
  typedef typename K::FT FT;

  FT x = random_float((K::FT) bbox.xmin(), (K::FT) bbox.xmax());
  FT y = random_float((K::FT) bbox.ymin(), (K::FT) bbox.ymax());
  FT z = random_float((K::FT) bbox.zmin(), (K::FT) bbox.zmax());

  return typename CGAL::Point_3<K>(x, y, z);
}

template <class K>
typename CGAL::Point_with_normal_3<K> random_pwn_in(const CGAL::Bbox_3 &bbox) {
  return typename CGAL::Point_with_normal_3<K>(random_point_in<K>(bbox),
    random_normal<K>());
}

template <class K>
typename CGAL::Vector_3<K> normalize(typename CGAL::Vector_3<K> const& v) {
  typename K::FT l = CGAL::sqrt(v.squared_length());
  if (l < (K::FT) 0.00001)
    return typename CGAL::Vector_3<K>((K::FT) 0, (K::FT) 0, (K::FT) 0);
  else return v * l;
}

template <typename K, typename OutputIterator>
void sample_sphere(const std::size_t num_points,
  typename CGAL::Point_3<K> const& center, typename K::FT radius,
  OutputIterator points) {
  for (std::size_t i = 0;i < num_points;++i) {
    CGAL::Vector_3<K> direction(random_float((K::FT) -1, (K::FT) 1), 
      random_float((K::FT) -1, (K::FT) 1), 
      random_float((K::FT) -1,(K::FT)  1));
    direction = direction * ((K::FT) 1.0 / CGAL::sqrt(direction.squared_length()));

    CGAL::Point_3<K> p = center + direction * radius;

    *points = CGAL::Point_with_normal_3<K>(p, direction);
    ++points;
  }
}

template <typename K, typename OutputIterator>
void sample_random_sphere_in_box(const std::size_t num_points,
  const CGAL::Bbox_3 &bbox, typename CGAL::Point_3<K> &center,
  typename K::FT &radius, OutputIterator points) {
    // Generate random parameters
    radius = random_float((K::FT) 0.01, (K::FT) 5);
    center = random_point_in<K>(bbox);

    sample_sphere(num_points, center, radius, points);
}

template <typename K, typename OutputIterator>
void sample_cylinder(const std::size_t num_points,
  const typename CGAL::Point_3<K> &center,
  const typename CGAL::Vector_3<K> &axis, const typename K::FT radius,
  const typename K::FT length, OutputIterator points) {
     // Sample shape
    for (size_t i = 0 ; i < num_points ; ++i) {
      CGAL::Vector_3<K> normal(
        random_float((K::FT) -1, (K::FT) 1), 
        random_float((K::FT) -1, (K::FT) 1), 
        random_float((K::FT) -1, (K::FT) 1));
      normal = normal - ((normal * axis) * axis);
      if (normal.squared_length() < (K::FT)0.0001) {
        i--;
        continue;
      }
      normal = normal * ((K::FT)1.0 / CGAL::sqrt(normal.squared_length()));

      typename K::FT l = random_float(-length, length);

      CGAL::Point_3<K> p = center + axis * l + radius * normal;
      *points = CGAL::Point_with_normal_3<K>(p, normal);
      ++points;
    }
}

template <typename K, typename OutputIterator>
void sample_random_cylinder(const std::size_t num_points,
  typename CGAL::Point_3<K> &center,
  typename CGAL::Vector_3<K> &axis, typename K::FT &radius,
  OutputIterator points) {
    // Generate random parameters
    axis = random_normal<K>();
    radius = random_float((K::FT) 0.5, (K::FT) 5);
    typename K::FT length = random_float((K::FT) 0.2, (K::FT) 10);

    // Find random center point placed on the plane through
    // the origin with 'axis' as normal.
    CGAL::Vector_3<K> u = random_normal<K>();
    CGAL::Vector_3<K> v = CGAL::cross_product(u, axis);
    v = v * 1.0 / CGAL::sqrt(v.squared_length());
    u = CGAL::cross_product(v, axis);

    center = CGAL::ORIGIN + random_float((K::FT) -5, (K::FT) 5) * u + random_float((K::FT) -5, (K::FT) 5) * v;

    sample_cylinder(num_points, center, axis, radius, length, points);
}

template <typename K, typename OutputIterator>
void sample_cone(const std::size_t num_points,
  typename CGAL::Point_3<K> const& apex,
  typename CGAL::Vector_3<K> const& axis, typename K::FT angle,
  typename K::FT s, typename K::FT e, OutputIterator points) {
  typedef typename K::FT FT;
  typedef typename CGAL::Vector_3<K> Vector;
  FT radiusGrow = std::tan(angle);

  FT cosAng = std::cos(angle);
  FT sinAng = std::sin(angle);

  for (size_t i = 0 ; i < num_points ; ++i)
  {
    Vector normal(random_float((K::FT) -1, (K::FT) 1), 
      random_float((K::FT)-1, (K::FT) 1), 
      random_float((K::FT)-1, (K::FT) 1));
    normal = normal - ((normal * axis) * axis);
    if (normal.squared_length() < (K::FT) 0.0001) {
      i--;
      continue;
    }
    normal = normal * ((K::FT) 1.0 / CGAL::sqrt(normal.squared_length()));

    FT l = random_float(s, e);

    CGAL::Point_3<K> p = apex + axis * l + (l * radiusGrow) * normal;

    // axis is pointing down from apex
    normal = normal * (K::FT) 1.0 / (CGAL::sqrt(normal.squared_length()));
    // normal is pointing from axis to surface point
    normal = normal * cosAng - axis * sinAng;
    l = CGAL::sqrt(normal.squared_length());
    if ((K::FT) 0.95 > l || l > (K::FT) 1.05)
      std::cout << "normal not normalized" << std::endl;

    *points = CGAL::Point_with_normal_3<K>(p, normal);
    ++points;
  }
}

template <typename K, typename OutputIterator>
void sample_random_cone(const std::size_t num_points,
  typename CGAL::Point_3<K> &apex,  typename CGAL::Vector_3<K> &axis,
  typename K::FT &angle, OutputIterator points) {
    // Generate random parameters
    apex = random_point_in<K>(CGAL::Bbox_3(-5, -5, -5, 5, 5, 5));
    axis = random_normal<K>();
    angle = random_float((K::FT) 0.2, (K::FT) 1.4);
    K::FT start  = random_float((K::FT) 0, (K::FT) 2.5);
    K::FT end = start + random_float((K::FT) 0.5, (K::FT) 2.5);

    sample_cone(num_points, apex, axis, angle, start, end, points);
}

template <typename K, typename OutputIterator>
void sample_random_parallelogram_in_box(const std::size_t num_points,
  const CGAL::Bbox_3 &bbox, typename CGAL::Vector_3<K> &normal,
  typename K::FT &dist, OutputIterator points) {
  // Generate random plane from 3 non collinear points.
  CGAL::Vector_3<K> u, v;
  CGAL::Point_3<K> p[3];
  do {
    p[0] = random_point_in<K>(bbox);
    p[1] = random_point_in<K>(bbox);
    p[2] = random_point_in<K>(bbox);

    CGAL::Vector_3<K> a = p[1] - p[0];
    CGAL::Vector_3<K> b = p[2] - p[0];

    if (a.squared_length() < (K::FT) 4.0 || b.squared_length() < (K::FT) 4.0) {
      normal = CGAL::Vector_3<K>((K::FT) 0, (K::FT) 0, (K::FT) 0);
      continue;
    }

    a = a * 1.0 / CGAL::sqrt(a.squared_length());
    b = b * 1.0 / CGAL::sqrt(b.squared_length());

    normal = CGAL::cross_product(a, b);
    // repeat if angle between a and b is small
  } while (normal.squared_length() < (K::FT) 0.2);

  normal = normal * (K::FT) 1.0 / CGAL::sqrt(normal.squared_length());
  dist = -((p[0] - CGAL::ORIGIN) * normal);

  // sample
  size_t i = 0;
  u = p[1] - p[0];
  v = p[2] - p[0];
  
  for (std::size_t i = 0;i < num_points; ++i) {
    K::FT s = random_float((K::FT) 0, (K::FT) 1);
    K::FT t = random_float((K::FT) 0, (K::FT) 1);

    *points = CGAL::Point_with_normal_3<K>(p[0] + s * u + t * v, normal);
    ++points;
  }
}

template <typename K, typename OutputIterator>
void generate_points_on_torus(const std::size_t num_points,
  typename const CGAL::Point_3<K> &center,
  typename const CGAL::Vector_3<K> &axis,
  typename const K::FT major_radius,
  typename const K::FT minor_radius,
  OutputIterator points) {
  size_t i = 0;

  // calculate basis  
  CGAL::Vector_3<K> b1, b2;
  b1 = CGAL::cross_product(axis, CGAL::Vector_3<K>((K::FT) 0, (K::FT) 0, (K::FT) 1));
  if (b1.squared_length() < (K::FT) 0.1)
    b1 = CGAL::cross_product(axis, CGAL::Vector_3<K>((K::FT) 0, (K::FT) 1, (K::FT) 0));
  b1 = b1 / CGAL::sqrt(b1.squared_length());

  b2 = CGAL::cross_product(axis, b1);
  b2 = b2 / CGAL::sqrt(b2.squared_length());

  for (size_t i = 0 ; i < num_points ; ++i)
  {
    K::FT tau = random_float((K::FT) 0, (K::FT) (2 * 3.141592656));
    K::FT phi = random_float((K::FT) 0, (K::FT) (2 * 3.141592656));

    CGAL::Vector_3<K> normal = sin(tau) * b1 + cos(tau) * b2;
    normal = normal / CGAL::sqrt(normal.squared_length());
    CGAL::Point_3<K> p = center + major_radius * normal;
    normal = sin(phi) * normal + cos(phi) * axis;
    p = p + minor_radius * normal;

    *points = CGAL::Point_with_normal_3<K>(p, normal);
    ++points;
  }
}

template <typename K, typename OutputIterator>
void sample_random_torus(const std::size_t num_points, typename CGAL::Point_3<K> &center,
  typename CGAL::Vector_3<K> &axis, typename K::FT &major_radius,
  typename K::FT &minor_radius, OutputIterator points) {
    // Generate random parameters
    center = random_point_in<K>(CGAL::Bbox_3(-5, -5, -5, 5, 5, 5));
    axis = random_normal<K>();
    major_radius = random_float((K::FT) 1.0, (K::FT) 5.0);
    minor_radius = random_float((K::FT) 0.1, (K::FT) 1.0);

    generate_points_on_torus(num_points, center, axis, major_radius, minor_radius, points);
}

template <typename K, typename P>
void filter_by_distance(typename const CGAL::Plane_3<K> &plane, typename K::FT dist, std::vector<typename P> &points) {
  K::FT d2 = dist * dist;

  std::vector<P>::iterator it = points.begin();
  while (it != points.end()) {
    if (CGAL::squared_distance(plane, *it) < d2)
      it = points.erase(it);
    else it++;
  }
}

template <typename K>
void save_scene(const std::string &fn, const std::vector<CGAL::Point_with_normal_3<K> > &pts) {
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
