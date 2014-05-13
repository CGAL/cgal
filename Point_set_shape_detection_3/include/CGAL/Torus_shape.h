#ifndef CGAL_EFFICIENT_RANSAC_TORUS_H
#define CGAL_EFFICIENT_RANSAC_TORUS_H

#include "Shape_base.h"
#include <set>
#include <math.h>
#include <cmath>

namespace CGAL {

  namespace Efficient_ransac {

    template <typename Kernel, class inputDataType>
    class Torus : public Primitive_ab<Kernel, inputDataType> {
    public:
      Point m_center;
      Vector m_axis;
      FT m_majorRad;
      FT m_minorRad;

      Torus() :  Primitive_ab<Kernel, inputDataType>(0.1, 0.9) {m_type = TORUS; m_type_name ="Torus";}
      Torus(FT _a, FT _b) : Primitive_ab<Kernel, inputDataType>(_a, _b) {m_type = TORUS;m_type_name ="Torus";}
      void compute(std::set<int> &l_list_index_selected, InputConstIterator &first) 
      {
        m_isValid = false;
        if ( l_list_index_selected.size() < 3)
          return;

        std::vector<int> output(l_list_index_selected.begin(), l_list_index_selected.end()); 

        Point p1 = (first + output[0])->first;
        Point p2 = (first + output[1])->first;
        Point p3 = (first + output[2])->first;
        Point p4 = (first + output[3])->first;

        Vector n1 = (first + output[0])->second;
        Vector n2 = (first + output[1])->second;
        Vector n3 = (first + output[2])->second;
        Vector n4 = (first + output[3])->second;

        // Implemented method from 'Geometric least-squares fitting of spheres, cylinders, cones and tori' by G. Lukacs,A.D. Marshall, R. R. Martin
        double a01 = CGAL::cross_product(n1, n2) * n3;
        double b01 = CGAL::cross_product(n1, n2) * n4;
        double a0 = CGAL::cross_product(p3 - p2, n1) * n3;
        double b0 = CGAL::cross_product(p4 - p2, n1) * n4;
        double a1 = CGAL::cross_product(p1 - p3, n2) * n3;
        double b1 = CGAL::cross_product(p1 - p4, n2) * n4;
        double a = CGAL::cross_product(p1 - p3, p2 - p1) * n3;
        double b = CGAL::cross_product(p1 - p4, p2 - p1) * n4;

        double div = 1.0 / (b1 * a01 - b01 * a1);
        double p = ((a01 * b + b1 * a0 - b0 * a1 - b01 * a)) * div * 0.5;
        double q = (b * a0 - b0 * a) * div;

        FT root = p * p - q;
        if (p * p - q < 0)
          root = 0;

        double y1 = -p - sqrt(root);
        double y2 = -p + sqrt(root);
        double x1 = -(a1 * y1 + a) / (a01 * y1 + a0);
        double x2 = -(a1 * y2 + a) / (a01 * y2 + a0);

        // 1. center + axis
        FT majorRad1, minorRad1, dist1;
        Point c1;
        Vector axis1;
        if (isfinite(x1) && isfinite(y1)) {
          c1 = p1 + n1 * x1;
          axis1 = c1 - (p2 + n2 * y1);
          axis1 = axis1 / sqrt(axis1.squared_length());
          dist1 = getCircle(first, c1, axis1, majorRad1, minorRad1);
        }
        else dist1 = FLT_MAX;

        // 2. center + axis
        FT majorRad2, minorRad2, dist2;
        Point c2;
        Vector axis2;
        if (isfinite(x2) && isfinite(y2)) {
          c2 = p1 + n1 * x2;
          axis2 = c2 - (p2 + n2 * y2);
          axis2 = axis2 / sqrt(axis2.squared_length());
          dist2 = getCircle(first, c2, axis2, majorRad2, minorRad2);
        }
        else dist2 = FLT_MAX;
        

        if (dist1 < dist2) {
          if (dist1 > m_epsilon)
            return;
          m_center = c1;
          m_axis = axis1;
          m_majorRad = majorRad1;
          m_minorRad = sqrt(minorRad1);
        }
        else {
          if (dist2 > m_epsilon)
            return;
          m_center = c2;
          m_axis = axis2;
          m_majorRad = majorRad2;
          m_minorRad = sqrt(minorRad2);
        }

        m_isValid = true;
        
        //create primitive
        //validate points and normals
      }

      std::string info()
      {
        std::stringstream sstr;
        sstr << "Type: " << m_type_name;
        /*sstr << "Type: " << m_type_name << "(" << m_normal.x() << ", " << m_normal.y() << ", " << m_normal.z() << ")x - " << m_d << "= 0"
          << " ev: " << ExpectedValue() << " s: " << m_nb_subset_used << " #Pts: " <<  m_indices.size();
        if (uExt)
          std::cout << " uE: " << uExt << " vE: " << vExt;*/
        std::cout << std::endl;

        return sstr.str();
      };
      std::string type_str() const {return m_type_name;}

      Point pointOnPrimitive() const {return Point();}//m_point_on_primitive;}

      void parameters(InputConstIterator first, std::vector<std::pair<FT, FT>> &parameterSpace, const std::vector<int> &indices, FT min[2], FT max[2]) const {
        return;
      }

      void parameterExtend(const Point &center, FT width, FT min[2], FT max[2]) const {
        return;
      }

      FT squared_distance(const Point &_p) const {
        Vector d = _p - m_center;
        // height over symmetry plane
        FT p = d * m_axis;
        // distance from axis in plane
        FT l = sqrt(d * d - p * p);

/*
        Vector inPlane = CGAL::cross_product(m_axis, CGAL::cross_product(m_axis, d));
        if (inPlane * d < 0)
          inPlane = -inPlane;*/

        // inPlane distance from circle
        FT l2 = m_majorRad - l;

        // distance from torus
        l = sqrt(p * p + l2 * l2) - m_minorRad;

        return l * l;
      }

      void squared_distance(InputConstIterator first, std::vector<FT> &dists, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) {
        for (unsigned int i = 0;i<indices.size();i++) {
          if (shapeIndex[indices[i]] == -1) {
            Point po = first[i].first;
            Vector n = first[i].second;
            Vector d = first[i].first - m_center;
            // height over symmetry plane
            FT p = d * m_axis;
            // distance from axis in plane
            FT l = sqrt(d * d - p * p);

            // inPlane distance from circle
            FT l2 = m_majorRad - l;

            // distance from torus
            l = sqrt(p * p + l2 * l2) - m_minorRad;
            dists[i] = l * l;

            if (dists[i] > m_epsilon) {
              int asd;
              asd = 3;
            }
          }
        }
      }

      void cos_to_normal(InputConstIterator first, std::vector<FT> &angles, const std::vector<int> &shapeIndex, const std::vector<unsigned int> &indices) const {
        for (unsigned int i = 0;i<indices.size();i++) {
          if (shapeIndex[indices[i]] == -1) {
            Vector d = first[i].first - m_center;
            // height over symmetry plane
            FT p = d * m_axis;
            // distance from axis in plane
            FT l = sqrt(d * d - p * p);

            Vector inPlane = CGAL::cross_product(m_axis, CGAL::cross_product(m_axis, d));
            if (inPlane * d < 0)
              inPlane = -inPlane;
            inPlane = inPlane / sqrt(inPlane.squared_length());

            d = first[i].first - (m_center + inPlane * m_majorRad);
            d = d / sqrt(d.squared_length());
            angles[i] = abs(d * first[i].second);
          }
        }
      }

      FT cos_to_normal(const Point &_p, const Vector &_n) const{
        Vector d = _p - m_center;
        // height over symmetry plane
        FT p = d * m_axis;
        // distance from axis in plane
        FT l = sqrt(d * d - p * p);

        Vector inPlane = CGAL::cross_product(m_axis, CGAL::cross_product(m_axis, d));
        if (inPlane * d < 0)
          inPlane = -inPlane;
        inPlane = inPlane / sqrt(inPlane.squared_length());

        d = _p - (m_center + inPlane * m_majorRad);
        d = d / sqrt(d.squared_length());
        return abs(d * _n);
      } 

      virtual bool supportsConnectedComponent() {return false;}
      virtual bool wrapsU() const {return false;}
      virtual bool wrapsV() const {return false;}

      private:
        FT getCircle(InputConstIterator first, Point &center, const Vector &axis, FT &majorRad, FT &minorRad) const {
          // create spin image
          Kernel::Point_2 pts[4];
          for (unsigned int i = 0;i<4;i++) {
            Vector d = first[i].first - center;
            FT p = d * axis;
            pts[i] = Kernel::Point_2(p, sqrt(d * d - p * p));
          }
          Kernel::Circle_2 c(pts[0], pts[1], pts[2]);
          minorRad = c.squared_radius();
          majorRad = c.center().y();
          center = center + c.center().x() * axis;

          return abs((pts[3] - c.center()).squared_length() - c.squared_radius());
        }
    };
  }
}
#endif