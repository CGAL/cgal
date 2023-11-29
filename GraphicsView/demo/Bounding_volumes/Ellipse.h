#ifndef ELLIPSE_H
#define ELLIPSE_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

template <class R>
 class Ellipse_2 {

    typedef typename R::Point_2 Point_2;
    typedef typename R::Vector_2 Vector_2;
    typedef typename R::Conic_2 Conic_2;

    Point_2 ce;
    Vector_2 va_, vb_;
    double a_,b_;

 public:

    template <class T>
    Ellipse_2(const CGAL::Min_ellipse_2<T>& me)
    {
      CGAL::Simple_cartesian<double>::Conic_2 co;
      me.ellipse().double_conic(co);

      // https://en.wikipedia.org/wiki/Matrix_representation_of_conic_sections
      double A = co.r();
      double B = co.t();
      double C = co.s();
      double D = co.u();
      double E = co.v();
      double F = co.w();

      assert(A * C - B*B/4 >= 0);

      Eigen::Matrix2d A33, A33i;
      A33 << A, B/2, B/2, C;
      A33i = A33.inverse();

      Eigen::Vector2d v;
      v << -D/2 , -E/2;

      v = A33i * v;

      ce = Point_2(v[0],v[1]);

      Eigen::EigenSolver<Eigen::Matrix2d> es(A33);

      Eigen::Matrix2cd ev;
      ev =  es.eigenvectors();
      double x0 = ev.col(0)[0].real();
      double y0 = ev.col(0)[1].real();
      double x1 = ev.col(1)[0].real();
      double y1 = ev.col(1)[1].real();

      double lambda1 = es.eigenvalues()[0].real();
      double lambda2 = es.eigenvalues()[1].real();

      Eigen::Matrix3d AQ;
      AQ << A, B/2, D/2,
            B/2, C, E/2,
            D/2, E/2, F;

      double K = AQ.determinant() / A33.determinant();

      a_ = sqrt(-K/lambda1);
      b_ = sqrt(-K/lambda2);

      va_ =  a_ * Vector_2(x0,y0);
      vb_ =  b_ * Vector_2(x1,y1);
    }



    Point_2 center() const
    {
      return ce;
    }

    double a() const
    {
      return a_;
    }

    double b() const
    {
      return b_;
    }

    double eccentricity() const
    {
      if(a_ == b_){
        return 0;
      }
      return std::sqrt(1.0 - (b_*b_)/(a_*a_));
    }

    const Vector_2& va() const
    {
      return va_;
    }
    const Vector_2& vb() const
    {
      return vb_;
    }

 };
#endif // ELLIPSE_H
