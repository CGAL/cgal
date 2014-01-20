#ifndef CGAL_EFFICIENT_RANSAC_CONE_H
#define CGAL_EFFICIENT_RANSAC_CONE_H

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
    class Cone : public Primitive_ab<Kernel, inputDataType>
    {
    public:
      typedef typename std::vector<inputDataType>::iterator inputIterator;
      typedef typename std::vector<inputDataType>::const_iterator InputConstIterator;
      typedef typename Kernel::FT FT;
      typedef typename Kernel::Line_3 Line;
      typedef typename Kernel::Point_3 Point;
      typedef typename Kernel::Point_2 Point_2d;
      typedef typename Kernel::Vector_3 Vector;
      typedef typename Kernel::Plane_3 Plane_3;
      typedef typename Kernel::Sphere_3 Sphere_3;
      FT m_angle;
      Point m_apex;
      Vector m_axis;
      Point m_pointOnPrimitive;
      FT m_nSinAng, m_cosAng;

    public:
      Cone() :  Primitive_ab<Kernel, inputDataType>(0.1, 0.9) {m_type = CONE; m_type_name ="Cone";}
      Cone(FT _a, FT _b)  :  Primitive_ab<Kernel, inputDataType>(_a, _b)  {m_type = CONE;m_type_name ="Cone";}

      void compute(std::set<int> &l_list_index_selected, InputConstIterator &m_it_Point_Normal) {
        m_isValid = false;
        if ( l_list_index_selected.size() < 3)
          return;

        std::vector<int> output(l_list_index_selected.begin(), l_list_index_selected.end());
        Point p1 = (m_it_Point_Normal + output[0])->first;
        Vector n1 = (m_it_Point_Normal + output[0])->second;
        Point p2 = (m_it_Point_Normal + output[1])->first;
        Vector n2 = (m_it_Point_Normal + output[1])->second;
        Point p3 = (m_it_Point_Normal + output[2])->first;
        Vector n3 = (m_it_Point_Normal + output[2])->second;

        // first calculate intersection of three planes -> apex

        Vector lineDir = CGAL::cross_product(n1, n2);
        lineDir = lineDir * 1.0 / (sqrt(lineDir.squared_length()));

        // lineDir not normalized direction of intersection lines of two planes (p1, n1) and (p2, n2)
        // get point on line by moving point p1 onto line
        Vector orthLineInPlane = CGAL::cross_product(n1, lineDir);
        orthLineInPlane = orthLineInPlane * 1.0 / (sqrt(orthLineInPlane.squared_length()));

        // distance of p1 to (p2, n2)
        FT d = (p1 - CGAL::ORIGIN) * n2 - (p2 - CGAL::ORIGIN) * n2;
        // projection of orthLineInPlane onto p2
        FT l = orthLineInPlane * n2;
        Point pointOnLine = p1 - (d/l) * orthLineInPlane;


        // checking
        d = (pointOnLine - p1) * n1; // should be 0
        d = (pointOnLine - p2) * n2; // should be 0


        // distance of pLineDir to (p3, n3)
        d = (pointOnLine - CGAL::ORIGIN) * n3 - (p3 - CGAL::ORIGIN) * n3;
        l = lineDir * n3;
        m_apex = pointOnLine - (d/l) * lineDir;


        // checking
        d = (m_apex - p1) * n1; // should be 0
        d = (m_apex - p2) * n2; // should be 0
        d = (m_apex - p3) * n3; // should be 0

        // 2. find axis
        Vector v1 = p1 - m_apex;
        v1 = v1 * 1.0 / (sqrt(v1.squared_length()));
        Point c1 = m_apex + v1;

        Vector v2 = p2 - m_apex;
        v2 = v2 * 1.0 / (sqrt(v2.squared_length()));
        Point c2 = m_apex + v2;

        Vector v3 = p3 - m_apex;
        v3 = v3 * 1.0 / (sqrt(v3.squared_length()));
        Point c3 = m_apex + v3;

        m_axis = CGAL::cross_product(c1 - c2, c1 - c3);
        m_axis = (orthLineInPlane * m_axis < 0) ? -m_axis : m_axis;
        m_axis = m_axis * 1.0 / sqrt(m_axis.squared_length());

        m_angle = acos(v1 * m_axis) + acos(v2 * m_axis) + acos(v3 * m_axis);
        m_angle /= 3;
        if (m_angle < 0 || m_angle > M_PI / 2.12)
          return;

        m_nSinAng = -sin(m_angle);
        m_cosAng = cos(m_angle);

        m_isValid = true;
      }

      std::string info() {
        std::stringstream sstr;

        sstr << "Type: " << m_type_name << " apex: (" << m_apex.x() << ", " << m_apex.y() << ", " << m_apex.z() << ") axis: (" << m_axis.x() << ", " << m_axis.y() << ", " << m_axis.z() << ") angle:" << m_angle
          << " ev: " << ExpectedValue() << " s: " << m_nb_subset_used << " #Pts: " <<  m_indices.size()	<< std::endl;

        return sstr.str();
      }

      std::string type_str() const {return m_type_name;}

      void parameters(InputConstIterator first, std::vector<std::pair<FT, FT>> &parameterSpace, const std::vector<int> &indices, FT min[2], FT max[2]) const {

      }

      void parameterExtend(const Point &center, FT width, FT min[2], FT max[2]) const {
      }

      Point pointOnPrimitive() const {
        return m_pointOnPrimitive;
      }

      FT squared_distance(const Point &_p) const {
        Vector toApex = _p - m_apex;
        FT a = toApex.squared_length();
        // projection on axis
        FT b = toApex * m_axis;
        // distance to axis
        FT l = sqrt(a - b * b); // should never be negative as m_axis is normalized
        FT c = m_cosAng * l;
        FT d = m_nSinAng * b;

        // far on other side?
        return (b < 0 && c - d < 0) ? a : abs(c + d) * abs(c + d);
      }

      void squared_distance(InputConstIterator first, std::vector<FT> &dists, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) {
        for (unsigned int i = 0;i<indices.size();i++) {
          if (shapeIndex[indices[i]] == -1) {
            Vector toApex = first[indices[i]].first - m_apex;
            FT a = toApex.squared_length();
            // projection on axis
            FT b = toApex * m_axis;
            // distance to axis
            FT l = sqrt(a - b * b); // should never be negative as m_axis is normalized
            FT c = m_cosAng * l;
            FT d = m_nSinAng * b;

            // far on other side?
            dists[i] = (b < 0 && c - d < 0) ? a : abs(c + d) * abs(c + d);
            if (dists[i] > 0.01) {
              int asd;
              asd = 2;
            }
          }
        }
      }

      void cos_to_normal(InputConstIterator first, std::vector<FT> &angles, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) const {
        for (unsigned int i = 0;i<indices.size();i++) {
          if (shapeIndex[indices[i]] == -1) {
            // construct vector orthogonal to axis in direction of the point
            Vector a = first[indices[i]].first - m_apex;
            Vector b = CGAL::cross_product(m_axis, CGAL::cross_product(m_axis, a));
            b = (a * b < 0) ? -b : b;
            b = b * 1.0 / sqrt(b.squared_length());
            b = m_cosAng * b + m_nSinAng * m_axis;

            angles[i] = abs(first[indices[i]].second * b);
          }
        }
      }

      FT cos_to_normal(const Point &_p, const Vector &_n) const {
        // construct vector orthogonal to axis in direction of the point
        Vector a = _p - m_apex;
        Vector b = CGAL::cross_product(m_axis, CGAL::cross_product(m_axis, a));
        b = (a * b < 0) ? -b : b;
        b = b * 1.0 / sqrt(b.squared_length());
        b = m_cosAng * b + m_nSinAng * m_axis;

        return abs(_n * b);
      }

      virtual bool supportsConnectedComponent() {return false;}
      // U is longitude
      virtual bool wrapsU() const {return true;}
      // V is between caps
      virtual bool wrapsV() const {return false;}
    };
  }
}
#endif
