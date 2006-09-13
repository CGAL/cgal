#ifndef CGAL_ARR_OF_S_3_H_SIM_FK
#define CGAL_ARR_OF_S_3_H_SIM_FK

#include <CGAL/Arrangement_of_spheres_traits_3.h>
#include <CGAL/Polynomial/internal/Kernel/Rational_between_roots.h>
#include <CGAL/Arrangement_of_spheres_3/coordinates.h>

struct Function_kernel {
  typedef Function_kernel This;
  typedef Arrangement_of_spheres_traits_3::Event_point_3 Root;
  typedef Root::NT NT;

  struct To_isolating_interval {
    typedef std::pair<NT, NT> result_type;
    typedef Root argument_type;
    result_type operator()(const argument_type &a) const {
      std::pair<double, double> i= CGAL::to_interval(a);
      return result_type(NT(i.first), NT(i.second));
    }
  };

  To_isolating_interval to_isolating_interval_object() const {
    return To_isolating_interval();
  }

  typedef CGAL::POLYNOMIAL::internal::Rational_between_roots<This>
  Rational_between_roots;
  Rational_between_roots rational_between_roots_object() const {
    return Rational_between_roots(*this);
  }
  
  struct Is_rational {
    typedef bool result_type;
    typedef Root argument_type;
    result_type operator()(const argument_type &a) const {
      return a.has_simple_coordinate(sweep_coordinate());
    }
  };

  Is_rational is_rational_object() const {
    return Is_rational();
  }
  
  struct To_rational {
    typedef NT result_type;
    typedef Root argument_type;
    result_type operator()(const argument_type &a) const {
      return a.simple_coordinate(sweep_coordinate());
    }
  };

  To_rational to_rational_object() const {
    return To_rational();
  }

};

#endif
