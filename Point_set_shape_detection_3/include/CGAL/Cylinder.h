#ifndef CGAL_EFFICIENT_RANSAC_CYLINDER_H
#define CGAL_EFFICIENT_RANSAC_CYLINDER_H

#include "Primitive.h"
#include <set>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

namespace CGAL {

  namespace Efficient_ransac {

    template <typename Kernel, class inputDataType>
    class Cylinder : public Primitive_ab<Kernel, inputDataType>
    {
    public:
      typedef typename Kernel::FT FT;
      typedef typename Kernel::Line_3 Line;
      typedef typename Kernel::Point_3 Point;
      typedef typename Kernel::Vector_3 Vector;
      typedef typename Kernel::Plane_3 Plane_3;
      typedef typename Kernel::Sphere_3 Sphere_3;
      FT m_radius;
      Line m_axis;
      Point m_point_on_axis;

    public:
      Cylinder() :  Primitive_ab<Kernel, inputDataType>(0.1, 0.9) {m_type = CYLINDER; m_type_name ="Cylinder";}
      Cylinder(FT _a, FT _b)  :  Primitive_ab<Kernel, inputDataType>(_a, _b)  {m_type = CYLINDER;m_type_name ="Cylinder";}

      void compute(std::set<int> &l_list_index_selected, InputConstIterator &m_it_Point_Normal) {
        if ( l_list_index_selected.size() < 3) {
          m_isValid = false;
          return;
        }

        std::vector<int> output(l_list_index_selected.begin(), l_list_index_selected.end());
        Point p1 = (m_it_Point_Normal + output[0])->first;
        Vector n1 = (m_it_Point_Normal + output[0])->second;
        Point p2 = (m_it_Point_Normal + output[1])->first;
        Vector n2 = (m_it_Point_Normal + output[1])->second;
        Point p3 = (m_it_Point_Normal + output[2])->first;
        Vector n3 = (m_it_Point_Normal + output[2])->second;

        Vector axis = CGAL::cross_product((m_it_Point_Normal + output[0])->second, (m_it_Point_Normal + output[1])->second);
        FT axisL = sqrt(axis.squared_length());
        if (axisL < 0.001) {
          m_isValid = false;
          return;
        }
        axis = axis * (1.0 / axisL);

        // establish two directions in the plane axis * x = 0, whereas xDir is the projected n1
        Vector xDir = (m_it_Point_Normal + output[0])->second - ((m_it_Point_Normal + output[0])->second * axis) * axis;
        xDir = xDir * (1.0 / sqrt(xDir.squared_length()));
        Vector yDir = CGAL::cross_product(axis, xDir);
        yDir = yDir * (1.0 / sqrt(yDir.squared_length()));

        FT n2x = (m_it_Point_Normal + output[1])->second * yDir;
        FT n2y = -((m_it_Point_Normal + output[1])->second * xDir);

        Vector dist = (m_it_Point_Normal + output[1])->first - (m_it_Point_Normal + output[0])->first;

        FT Ox = xDir * dist;
        FT Oy = yDir * dist;

        FT lineDist = n2x * Ox + n2y * Oy;

        m_radius = lineDist / n2x;
        m_point_on_axis = (m_it_Point_Normal + output[0])->first + m_radius * xDir;
        m_radius = abs(m_radius);

        m_axis = Line(m_point_on_axis, axis);

        if (squared_distance((m_it_Point_Normal + output[0])->first) > m_epsilon || (cos_to_normal((m_it_Point_Normal + output[0])->first, (m_it_Point_Normal + output[0])->second) < m_normalThresh)) {
          m_isValid = false;
          return;
        }

      /*       FT xDist = dist * n2y;

        Point c = p1 + xDist * xDir;

        Point proj1 = (m_it_Point_Normal + output[0])->first - (((m_it_Point_Normal + output[0])->first - CGAL::ORIGIN) * axis) * axis;
        Point proj2 = (m_it_Point_Normal + output[1])->first - (((m_it_Point_Normal + output[1])->first - CGAL::ORIGIN) * axis) * axis;

   */
      }

      std::string info() {
        std::stringstream sstr;
        Point c = m_axis.point();
        Vector a = m_axis.to_vector();

        sstr << "Type: " << m_type_name << " c: (" << c.x() << ", " << c.y() << ", " << c.z() << ") a: (" << a.x() << ", " << a.y() << ", " << a.z() << ") r:" << m_radius
          << " ev: " << ExpectedValue() << " s: " << m_nb_subset_used << " #Pts: " <<  m_indices.size()	<< std::endl;

        return sstr.str();
      }

      std::string type_str() const {return m_type_name;}

      void parameters(InputConstIterator first, std::vector<std::pair<FT, FT>> &parameterSpace, const std::vector<int> &indices, FT min[2], FT max[2]) const {
        Vector d1 = Vector(0, 0, 1);
        Vector a = m_axis.to_vector();
        a = a * (1.0 / sqrt(a.squared_length()));

        Vector d2 = CGAL::cross_product(a, d1);
        FT l = d2.squared_length();
        if (l < 0.0001) {
          d1 = Vector(1, 0, 0);
          d2 = CGAL::cross_product(m_axis.to_vector(), d1);
          l = d2.squared_length();
          if (l < 0.0001) {
            std::cout << "Cylinder::pointOnPrimitive() construction failed!" << std::endl;
          }
        }
        d2 = d2 / sqrt(l);

        d1 = CGAL::cross_product(m_axis.to_vector(), d2);
        d1 = d1 * (1.0 / sqrt(d1.squared_length()));

        // 1.0 / circumfence
        FT c = 1.0 / 2 * M_PI * m_radius;

        // first one separate for initializing min/max
        Vector vec = first[indices[0]].first - m_point_on_axis;
        FT v = vec * a;
        vec = vec - ((vec * a) * a);
        vec = vec * (1.0 / sqrt(vec.squared_length()));

        FT a1 = acos(vec * d1);
        FT a2 = acos(vec * d2);

        FT u = ((a2 < M_PI_2) ? 2 * M_PI - a1 : a1) * c;

        parameterSpace[0] = std::pair<FT, FT>(u, v);

        min[0] = max[0] = u;
        min[1] = max[1] = v;

        for (unsigned int i = 0;i<indices.size();i++) {
          Vector vec = first[indices[i]].first - m_point_on_axis;
          FT v = vec * a;
          vec = vec - ((vec * a) * a);
          vec = vec * (1.0 / sqrt(vec.squared_length()));

          FT a1 = acos(vec * d1);
          FT a2 = acos(vec * d2);

          FT u = ((a2 < M_PI_2) ? 2 * M_PI - a1 : a1) * c;

          min[0] = std::min<FT>(min[0], u);
          max[0] = std::max<FT>(max[0], u);
          min[1] = std::min<FT>(min[1], v);
          max[1] = std::max<FT>(max[1], v);

          parameterSpace[i] = std::pair<FT, FT>(u, v);
        }


      }

      void parameterExtend(const Point &center, FT width, FT min[2], FT max[2]) const {
        //V length of axis in box? not enough
        FT maxLambda = std::numeric_limits<double>::max(), minLambda = -std::numeric_limits<double>::max();
        Vector a = m_axis.to_vector();
        Point p = m_point_on_axis;

        for (unsigned int i = 0;i<3;i++) {
          if (abs(a[i]) > 0.001) {
            FT l1 = (center[i] + width + m_radius - p[i]) / a[i];
            FT l2 = (center[i] - width - m_radius - p[i]) / a[i];
            if (l1 * l2 > 0) {
              std::cout << "Cylinder::parameterExtend(): dim 0, l1*l2 > 0" << std::endl;
            }
            minLambda = std::max<FT>(minLambda, std::min<FT>(l1, l2));
            maxLambda = std::min<FT>(maxLambda, std::max<FT>(l1, l2));
          }
        }

        min[1] = minLambda;
        max[1] = maxLambda;
        //U circumfence
        min[0] = 0;
        max[0] = (2 * M_PI * m_radius);
      }

      Point pointOnPrimitive() const {
        Vector d1 = Vector(0, 0, 1);
        Vector d2 = Vector::cross_product(m_axis.to_vector(), d1);
        FT l = d2.squared_length();
        if (l < 0.0001) {
          d1 = Vector(1, 0, 0);
          d2 = Vector::cross_product(m_axis.to_vector(), d1);
          l = d2.squared_length();
          if (l < 0.0001) {
            std::cout << "Cylinder::pointOnPrimitive() construction failed!" << std::endl;
          }
        }
        return m_point_on_axis + d2 * m_radius / sqrt(l);
      }

      FT squared_distance(const Point &_p) const {
        Vector a = m_axis.to_vector();
        a = a * (1.0 / sqrt(a.squared_length()));
        Vector v = _p - m_point_on_axis;
        v = v - ((v * a) * a);
        FT d = sqrt(v.squared_length()) - m_radius;
        return d * d;
      }

      void squared_distance(InputConstIterator first, std::vector<FT> &dists, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) {
        Vector a = m_axis.to_vector();
        a = a * (1.0 / sqrt(a.squared_length()));
        for (unsigned int i = 0;i<indices.size();i++) {
          if (shapeIndex[indices[i]] == -1) {
            Vector v = first[indices[i]].first - m_point_on_axis;
            v = v - ((v * a) * a);
            dists[i] = sqrt(v.squared_length()) - m_radius;
            dists[i] = dists[i] * dists[i];
          }
        }
      }

      void cos_to_normal(InputConstIterator first, std::vector<FT> &angles, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) const {
        Vector a = m_axis.to_vector();
        a = a * (1.0 / sqrt(a.squared_length()));
        for (unsigned int i = 0;i<indices.size();i++) {
          if (shapeIndex[indices[i]] == -1) {
            Vector v = first[indices[i]].first - m_point_on_axis;
            v = v - ((v * a) * a);
            v = v * (1.0 / sqrt(v.squared_length()));
            angles[i] = abs(v * first[indices[i]].second);
          }
        }
      }

      FT cos_to_normal(const Point &_p, const Vector &_n) const {
        Vector a = m_axis.to_vector();
        a = a * (1.0 / sqrt(a.squared_length()));
        Vector v = _p - m_point_on_axis;
        v = v - ((v * a) * a);
        v = v * (1.0 / sqrt(v.squared_length()));
        return abs(v * _n);
      }

      virtual bool supportsConnectedComponent() {return true;}
      // U is longitude
      virtual bool wrapsU() const {return true;}
      // V is between caps
      virtual bool wrapsV() const {return false;}
    };
  }
}
#endif
