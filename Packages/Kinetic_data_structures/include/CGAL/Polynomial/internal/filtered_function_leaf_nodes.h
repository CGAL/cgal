#ifndef CGAL_POLYNOMIAL_INTERNAL_VIRTUAL_FUNCTION_REPS_H
#define CGAL_POLYNOMIAL_INTERNAL_VIRTUAL_FUNCTION_REPS_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/filtered_function_node_bases.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template <class Traits>
class Filtered_function_node_constant: public Filtered_function_node<Traits> {
  typedef Filtered_function_node_constant<Traits> This;
  typedef Filtered_function_node<Traits> P;
public:
  Filtered_function_node_constant(const typename P::Exact_function::NT& val,
			    const typename P::Interval_function_converter &ifc): P(ifc){
    P::set_interval_function(typename P::Interval_function(ifc.nt_converter()(val)));
    P::set_exact_function(typename P::Exact_function(this->num_));
  }

  virtual ~Filtered_function_node_constant(){}

  virtual void write(std::ostream &out) const {
    out << P::exact_function();
  }
protected:

  virtual  void generate_exact_function() const {
    CGAL_Polynomial_assertion(0);
  }
};

template <class Traits>
class Filtered_function_node_double_constant: public Filtered_function_node<Traits> {
  typedef Filtered_function_node_double_constant<Traits> This;
  typedef Filtered_function_node<Traits> P;
public:
  Filtered_function_node_double_constant(double d, const typename P::Interval_function_converter &ifc): P(ifc){
    P::set_interval_function(typename P::Interval_function(to_interval(d)));
 
  }

  virtual ~Filtered_function_node_double_constant(){}

  virtual void write(std::ostream &out) const {
    out << P::interval_function()[0].sup();
  }
protected:

  virtual  void generate_exact_function() const {
    P::set_exact_function(typename P::Exact_function(typename P::Exact_function::NT(P::interval_function()[0].sup())));  
  }
};

template <class Traits>
class Filtered_function_node_explicit: public Filtered_function_node<Traits> {
  typedef Filtered_function_node_explicit<Traits> This;
  typedef Filtered_function_node<Traits> P;
public:
  Filtered_function_node_explicit(const typename P::Exact_function f, const typename P::Interval_function_converter &ifc): P(ifc){
    P::set_interval_function(ifc(f));
    P::set_exact_function(f);
  }

  virtual ~Filtered_function_node_explicit(){}

  virtual void write(std::ostream &out) const {
    if (P::exact_function().degree() >0){
      out << "("<< P::exact_function() << ")";
    } else {
      out << P::exact_function();
    }
  }
protected:

  virtual  void generate_exact_function() const {
    CGAL_Polynomial_assertion(0);
  }
};

template <class Traits>
class Filtered_function_node_double_explicit: public Filtered_function_node<Traits> {
  typedef Filtered_function_node_double_explicit<Traits> This;
  typedef Filtered_function_node<Traits> P;
public:
  template <class DF>
  Filtered_function_node_double_explicit(const DF& f, const typename P::Interval_function_converter &ifc): P(ifc){
    typename P::Interval_function fi;
    fi.set_nominal_degree(f.degree());
    for (int i=0; i<= f.degree(); ++i){
      fi.set_coef(i, to_interval(f[i]));
    }
    P::set_interval_function(fi);
  }
  virtual ~Filtered_function_node_double_explicit(){}

  virtual void write(std::ostream &out) const {
    if (P::interval_function().nominal_degree() >0){
      out << "("<< P::exact_function() << ")";
    } else {
      out << P::exact_function();
    }
  }
protected:

  virtual  void generate_exact_function() const {
    typename P::Exact_function fe;
    int deg= P::interval_function().nominal_degree();
    fe.set_nominal_degree(deg);
    for (int i=0; i<= this->f.degree(); ++i){
      fe.set_coef(i, typename P::Exact_function::NT(this->f[i].sup()));
    }
    P::set_exact_function(fe);
  }
};


CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif
