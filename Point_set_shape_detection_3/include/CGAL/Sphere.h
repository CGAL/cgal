#ifndef CGAL_EFFICIENT_RANSAC_SPHERE_H
#define CGAL_EFFICIENT_RANSAC_SPHERE_H

#include "Primitive.h"
#include <set>

namespace CGAL {

  namespace Efficient_ransac {

    template <typename Kernel, class inputDataType>
    class Sphere : public Primitive_ab<Kernel, inputDataType>
    {
    public:
      typedef typename Kernel::FT FT;
      typedef typename Kernel::Line_3 Line;
      typedef typename Kernel::Point_3 Point;
      typedef typename Kernel::Vector_3 Vector;
      typedef typename Kernel::Plane_3 Plane_3;
      typedef typename Kernel::Sphere_3 Sphere_3;
      Sphere_3 m_sphere;
      FT m_radius;

    public:
      Sphere() :  Primitive_ab<Kernel, inputDataType>(0.1, 0.9) {m_type = SPHERE; m_type_name ="Sphere";}
      Sphere(FT _a, FT _b)  :  Primitive_ab<Kernel, inputDataType>(_a, _b)  {m_type = SPHERE;m_type_name ="Sphere";}

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


        // Determine center: select midpoint of shortest line segment between p1 and p2
        // implemented from "3D game engine design" by Eberly 2001
        Vector diff = (m_it_Point_Normal + output[0])->first - (m_it_Point_Normal + output[1])->first;
        FT a = (m_it_Point_Normal + output[0])->second * (m_it_Point_Normal + output[0])->second;
        FT b = -((m_it_Point_Normal + output[0])->second * (m_it_Point_Normal + output[1])->second);
        FT c = (m_it_Point_Normal + output[1])->second * (m_it_Point_Normal + output[1])->second;
        FT d = (m_it_Point_Normal + output[0])->second * diff;
        FT f = diff * diff;

        FT det = abs(a*c - b*b);
        
        // parallel?
        if (det < 0.00001) {
          m_isValid = false;
          return;
        }

        FT e = -(m_it_Point_Normal + output[1])->second * diff;
        FT invDet = 1.0 / det;
        FT s = (b*e - c*d) * invDet;
        FT t = (d*b - a*e) * invDet;
        
        Point center = CGAL::ORIGIN + 0.5 * ((((m_it_Point_Normal + output[0])->first + s * (m_it_Point_Normal + output[0])->second) - CGAL::ORIGIN) + (((m_it_Point_Normal + output[1])->first + t * (m_it_Point_Normal + output[1])->second) - CGAL::ORIGIN));

        Vector v1 = ((m_it_Point_Normal + output[0])->first - center);
        Vector v2 = ((m_it_Point_Normal + output[1])->first - center);
        FT d1 = sqrt(v1.squared_length());
        FT d2 = sqrt(v2.squared_length());

        if (abs(d1-d2) > 2 * m_epsilon) {
          m_isValid = false;
          return;
        }

        v1 = v1 * (1.0 / d1);
        v2 = v2 * (1.0 / d2);

        if ((m_it_Point_Normal + output[0])->second * v1 < m_normalThresh || (m_it_Point_Normal + output[1])->second * v2 < m_normalThresh) {
          m_isValid = false;
          return;
        }

        Vector v3 = ((m_it_Point_Normal + output[2])->first - center);
        FT d3 = sqrt(v3.squared_length());
        v3 = v3 * (1.0 / d3);

        m_radius = (d1 + d2) * 0.5;

        if (abs(d3 - m_radius) > m_epsilon || (m_it_Point_Normal + output[2])->second * v3 < m_normalThresh) {
          m_isValid = false;
          return;
        }

        m_sphere = Sphere_3(center, m_radius * m_radius);
      }

      operator Sphere_3 () {
        return m_sphere;
      }

      std::string info() {
        std::stringstream sstr;
        Point c = m_sphere.center();
        FT r = sqrt(m_sphere.squared_radius());

        sstr << "Type: " << m_type_name << " c: (" << c.x() << ", " << c.y() << ", " << c.z() << ") r:" << r 
          << " ev: " << ExpectedValue() << " s: " << m_nb_subset_used << " #Pts: " <<  m_indices.size()	<< std::endl;

        return sstr.str();
      }

      std::string type_str() const {return m_type_name;}

      Point pointOnPrimitive() const {
        return m_sphere.center() + Vector(0, 0, sqrt(m_sphere.squared_radius()));
      }

      void parameters(InputConstIterator first, std::vector<std::pair<FT, FT>> &parameterSpace, const std::vector<int> &indices, FT min[2], FT max[2]) const {
      }

      void parameterExtend(const Point &center, FT width, FT min[2], FT max[2]) const {
      }

      FT squared_distance(const Point &_p) const {
        FT d = sqrt((m_sphere.center() - _p).squared_length()) - m_radius;
        return d*d;
      }

      void squared_distance(InputConstIterator first, std::vector<FT> &dists, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) {
        for (unsigned int i = 0;i<indices.size();i++) {
          if (shapeIndex[indices[i]] == -1) {
            dists[i] = sqrt((m_sphere.center() - first[indices[i]].first).squared_length()) - m_radius;
            dists[i] = dists[i] * dists[i];
          }
        }
      }

      void cos_to_normal(InputConstIterator first, std::vector<FT> &angles, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) const {
        for (unsigned int i = 0;i<indices.size();i++) {
          if (shapeIndex[indices[i]] == -1) {
            Vector n = m_sphere.center() - first[indices[i]].first;
            n = n * (1.0 / (sqrt(n.squared_length())));
            angles[i] = abs(first[indices[i]].second * n);
          }
        }
      }

      FT cos_to_normal(const Point &_p, const Vector &_n) const {
        Vector n = m_sphere.center() - _p;
        n = n * (1.0 / (sqrt(n.squared_length())));
        return abs(_n * n);
      }

      // U is longitude
      virtual bool supportsConnectedComponent() {return false;}
      virtual bool wrapsU() const {return true;}
      virtual bool wrapsV() const {return false;}
    };
  }
}
#endif
