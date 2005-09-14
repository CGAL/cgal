#ifndef RUNGE_KUTTA_INTEGRATOR_2_H_
#define RUNGE_KUTTA_INTEGRATOR_2_H_

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

CGAL_BEGIN_NAMESPACE

// The class Runge_kutta_integrator_2 is a model of the concept Integrator
template <class EulerIntegrator_2,class VectorField_2>
class Runge_kutta_integrator_2{
public:
	typedef Runge_kutta_integrator_2<EulerIntegrator_2, VectorField_2> self;
  typedef typename VectorField_2::FT FT;
  typedef typename VectorField_2::Point_2 Point_2;
  typedef typename VectorField_2::Vector_2 Vector_2;
  typedef typename VectorField_2::Vector_field_2 Vector_field_2;
  typedef typename EulerIntegrator_2::self Euler_integrator_2;
protected:
  Euler_integrator_2 * euler_integrator_2;
  FT default_integration_step;
public:
  Runge_kutta_integrator_2(Euler_integrator_2 & euler_integrator_2);
  Runge_kutta_integrator_2(Euler_integrator_2 & euler_integrator_2, const FT & integration_step);
  Point_2 operator()(const Point_2 & p, const Vector_field_2 & vector_field_2, const bool & index) const;
  Point_2 operator()(const Point_2 & p, const Vector_field_2 & vector_field_2, const FT & integration_step, const bool & index) const;
  Point_2 operator()(const Point_2 & p, const Vector_field_2 & vector_field_2, const FT & integration_step, Vector_2 v, const bool & index) const;

  inline FT get_default_integration_step(){return default_integration_step;}
  Euler_integrator_2 * get_euler_integrator_2(){return euler_integrator_2;}
  // just for debugging
  inline FT distance(const Point_2 & p, const Point_2 & q)
    {
      return sqrt(
		  (
		   (
		    p.x() - q.x()
		    )
		   *
		   (
		    p.x() - q.x()
		    )
		   )
		  +
		  (
		   (
		    p.y() - q.y()
		    )
		   *
		   (
		    p.y() - q.y()
		    )
		   )
		  );
    };
};


// The default constructor required by the concept
template <class EulerIntegrator_2, class VectorField_2>
Runge_kutta_integrator_2<EulerIntegrator_2, VectorField_2>::
Runge_kutta_integrator_2(Euler_integrator_2& integrator)
{
  default_integration_step =
    integrator.get_default_integration_step();
  euler_integrator_2 = &integrator;
}


// An additional constructor to specify the default integration step
template <class EulerIntegrator_2, class VectorField_2>
Runge_kutta_integrator_2<EulerIntegrator_2, VectorField_2>::Runge_kutta_integrator_2(Euler_integrator_2&
							  integrator, const FT & integration_step)
{
  default_integration_step = integration_step;
  euler_integrator_2 = &integrator;
}

// The integrator functor
template <class EulerIntegrator_2, class VectorField_2>
inline
typename Runge_kutta_integrator_2<EulerIntegrator_2, VectorField_2>::
Point_2 Runge_kutta_integrator_2<EulerIntegrator_2, VectorField_2>::operator()
  (const Point_2 & p, const Vector_field_2& vector_field_2, const FT & integration_step, Vector_2 v, const bool & index) const
{
  Point_2 p1 = (*euler_integrator_2)(p, vector_field_2, 0.5*integration_step, v, index);
  if(!vector_field_2.is_in_domain(p1))
    return p1;
  v = vector_field_2.get_field(p1).first;
  Point_2 p2 = (*euler_integrator_2)(p, vector_field_2, integration_step,v, index);
  return p2;
};

// The integrator functor with additional argument: the integration step
template <class EulerIntegrator_2, class VectorField_2>
inline
typename Runge_kutta_integrator_2<EulerIntegrator_2, VectorField_2>::Point_2 Runge_kutta_integrator_2<EulerIntegrator_2, VectorField_2>::operator()
(const Point_2 & p, const Vector_field_2& vector_field_2, const FT & integration_step, const bool & index) const{
  Vector_2 v;
  v = vector_field_2.get_field(p).first;
  Runge_kutta_integrator_2<EulerIntegrator_2, VectorField_2>
    runge_kutta_integrator_2(euler_integrator_2, integration_step);
  return runge_kutta_integrator_2(p, vector_field_2, integration_step, v, index);
};

// The integrator functor with additional argument: the integration step
template <class EulerIntegrator_2, class VectorField_2>
inline
typename Runge_kutta_integrator_2<EulerIntegrator_2, VectorField_2>::Point_2 Runge_kutta_integrator_2<EulerIntegrator_2, VectorField_2>::operator()
(const Point_2 & p, const Vector_field_2& vector_field_2, const bool & index) const
{
  Vector_2 v;
  v = vector_field_2.get_field(p).first;
  Runge_kutta_integrator_2<EulerIntegrator_2, VectorField_2> runge_kutta_integrator_2(*euler_integrator_2, default_integration_step);
  return runge_kutta_integrator_2(p, vector_field_2, default_integration_step, v, index);
};

CGAL_END_NAMESPACE

#endif
