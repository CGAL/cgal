#ifndef ELLIPSE_H
#define ELLIPSE_H

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

      double a11 = co.r();
      double a12 = co.t()/2.0;
      double a22 = co.s();
      double b1 = co.u();
      double b2 = co.v();
      double c = co.w();

      assert(a11*a22-a12*a12 >= 0);

      double delta = (a11>a22)? a11-a22 : a22-a11;
      if( (delta < 0.00000000001) && (a12 == 0)){
        ce = Point_2(-b1/(2*a11), -b2/(2*a11));
        a_ = b_ = std::sqrt(b1*b1+b2*b2-4*a11*c)/(2*a11);
        va_ = Vector_2(a_, 0);
        vb_ = Vector_2(0, b_);
        return;
        }

      double kden = 2.0*(a12*a12 - a11*a22);
      double k1 = (a22*b1 - a12*b2)/kden;
      double k2 = (a11*b2 - a12*b1)/kden;
      ce = Point_2(k1,k2);
      double mu = 1/(a11*k1*k1 + 2.0*a12*k1*k2 + a22*k2*k2 - c);
      double m11 = mu*a11;
      double m12 = mu*a12;
      double m22 = mu*a22;

      double r = std::sqrt((m11-m22)*(m11-m22) + 4.0*m12*m12);
      double lambda1 = ((m11+m22) + r)/2.0;
      double lambda2 = ((m11+m22) - r)/2.0;

      b_ = 1.0/sqrt(lambda1);
      a_ = 1.0/sqrt(lambda2);
      double omega = 1.0/
 std::sqrt((lambda1-m22)*(lambda1-m22) + m12*m12);

      double u1x = - a_ * omega*m12;
      double u1y = a_ * omega * (lambda1-m22);

      double u2x = b_* omega * (lambda1-m22);
      double u2y = b_* omega*m12;

      va_ = Vector_2(u1x,u1y);
      vb_ = Vector_2(u2x,u2y);
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
