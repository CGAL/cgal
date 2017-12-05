#ifndef CGAL_KINETIC_INTERVAL_SIMULATOR_TRAITS_H
#define CGAL_KINETIC_INTERVAL_SIMULATOR_TRAITS_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Ref_counted.h>


namespace CGAL { namespace Kinetic {

template <class Refiner>
class Interval_root {
  typedef Interval_root< Refiner> This;
public:
  typedef typename Refiner::Exact_root Exact_root;
  /*struct Refiner: public Ref_counted<Refiner> {
    typedef Exact Exact_root;
    virtual const Exact &exact_root(const std::pair<double,double> &iv) const=0;
    virtual bool refine(std::pair<double,double> &iv) const=0;
    virtual void write(std::ostream&out) const {
      out << "";
    }
    virtual ~Refiner(){};
    };*/

  //typedef typename Refiner::Handle RHandle;

  template <class NT>
  Interval_root(NT a) {
    iv_= CGAL::to_interval(a);
    CGAL_precondition(iv_.first == iv_.second);
  }

  Interval_root(): iv_(1,-1){
  }

  Interval_root(double a, double b, Refiner d): data_(d) {
    iv_= std::make_pair(a,b);
  }


  const Refiner &refiner() const {
    return data_;
  }
  Refiner &refiner() {
    return data_;
  }

  /*void refine(std::pair<double,double> nv) const {
    iv_.first = std::max(iv_.first, nv.first);
    iv_.second = std::min(iv_.second, nv.second);
    }*/

  const std::pair<double,double>& interval_value() const {  
    CGAL_precondition(iv_.first <= iv_.second);
    return iv_;
  }

  bool is_point() const {
    return iv_.first == iv_.second;
  }

  double double_value() const {
    CGAL_precondition(iv_.first <= iv_.second);
    return (iv_.first+iv_.second)/2.0;
  }
  const Exact_root& exact_root() const {
#if 0
    std::pair<double,double> oiv= iv_;
    Exact_root er=data_.exact_root(iv_);
    if (er < oiv.first || er > oiv.second) {
      std::cerr << "Exact root in bounds. Root is " << er 
		<< " bounds are " << oiv.first << " " << oiv.second << std::endl;
    }
    CGAL_precondition(er >= oiv.first);
    CGAL_precondition(er <= oiv.second);
    //CGAL_precondition(!is_point());
#endif
    if (data_.ensure_exact_root(iv_) && iv_.first != iv_.second) {
      std::pair<double,double> niv= to_interval(data_.exact_root());
      iv_.first = std::max(niv.first, iv_.first);
      iv_.second = std::min(niv.second, iv_.second);
    }
    return data_.exact_root();
  }

  bool is_rational() const {
    return is_point();
  }

  bool operator==(const This &o) const {
    return iv_== o.iv_;
  }

  bool operator!=(const This &o) const {
    return iv_!= o.iv_;
  }

  CGAL::Comparison_result compare(const This &o) const {
    //std::cout << "Comparing " << *this << " and " << o << std::endl;
    CGAL_precondition(iv_.first <= iv_.second);
    CGAL_precondition(o.iv_.first <= o.iv_.second);
    CGAL::Comparison_result cr= compare_interval(o);
    if (cr != CGAL::EQUAL /*|| 
			    data_== RHandle() && o.data_==RHandle()*/){
      //std::cout << "Got " << cr << std::endl;
      return cr;
    }
    bool changed=false;
    do {
      changed=false;
      if (iv_.second-iv_.first > o.iv_.second- o.iv_.first) {
	// this is not a point since it is larger than 0
	changed=data_.refine(iv_);
      } else if (!o.is_point()) {
	changed=o.data_.refine(o.iv_);
      }
    
      if (!changed) break;
      cr= compare_interval(o);
      if (cr != CGAL::EQUAL) {
	//std::cout << "Got " << cr << std::endl;
	return cr;
      }
    } while (true);
    if (is_point() && o.is_point()){
      cr = CGAL::EQUAL;
    }
    else /*if (is_point()) {
      cr= CGAL::compare(Exact_root(iv_.first), o.exact_root());
      } else if (o.is_point()){
    cr= CGAL::compare(exact_root(), Exact_root(o.iv_.first));
      } else*/ {
      //CGAL_assertion(!data_.equal_description(o.data_));
      if (data_.equal_description(o.data_)) {
	//std::cout << "Equal descriptions in comparisons" << std::endl;
	return CGAL::EQUAL;
      }
      cr= CGAL::compare(exact_root(), o.exact_root());
    }
    //std::cout << "Got " << EQUAL<< std::endl;
    return cr;
  }

  This operator-() const {
    CGAL_error();
    //return This(-iv_.first, -iv_.second, data_->negate);
    return This();
  }

  void write(std::ostream &out) const {
    if (iv_.first > iv_.second) {
      out << "[XXX]";
    } else {
      out << "[" << iv_.first << "..." << iv_.second 
	  << "] ";
      data_.write(out);
    }
  }

  void set_interval(const std::pair<double,double>& iv) {
    iv_= iv;
  }

  void refine() const {
    data_.refine(iv_);
  }

private:
  CGAL::Comparison_result compare_interval(const This &o) const {
    if (iv_.second < o.iv_.first) {
      return SMALLER;
    } else if (iv_.first > o.iv_.second) {
      return LARGER;
    } else {
      return EQUAL;
    }
  }

  mutable std::pair<double,double> iv_;
  Refiner data_;
};

template <class R>
std::ostream &operator<<(std::ostream &out, const Interval_root<R>& i) {
  i.write(out);
  return out;
}


/*template <class D>
Comparison_result compare(const Interval_root<D>& a,
			  const Interval_root<D>& b) {
  if (a.interval_value().second < b.interval_value().first) {
    return SMALLER;
  } else if (a.interval_value().first > b.interval_value().second) {
    return LARGER;
  } else {
    return EQUAL;
  }
}

template <class D>
double to_double(const Interval_root<D>& a) {
  return a.double_value();
}

template <class D>
std::pair<double,double> to_interval(const Interval_root<D>& a) {
  return a.interval_value();
  }*/


} } //namespace CGAL::Kinetic

namespace CGAL { 

template <class D>
Comparison_result compare(const typename Kinetic::Interval_root<D>& a,
			  const typename Kinetic::Interval_root<D>& b) {
  return a.compare(b);
}

template <class D>
double to_double(const typename Kinetic::Interval_root<D>& a) {
  return a.double_value();
}

template <class D>
std::pair<double,double> to_interval(const typename Kinetic::Interval_root<D>& a) {
  return a.interval_value();
}

} //namespace CGAL

namespace std {
  template <class D>
  class numeric_limits<CGAL::Kinetic::Interval_root<D> >: public numeric_limits<double>
  {
    typedef CGAL::Kinetic::Interval_root<D> Rt;
  public:
    static const bool is_specialized = true;
    static const bool has_infinity=true;
    static Rt infinity() throw() {return Rt(std::numeric_limits<double>::infinity());}
    static Rt max() throw() {return Rt(std::numeric_limits<double>::max());}
  };
};

namespace CGAL { namespace Kinetic {

template <class Data>
struct Interval_simulator_traits {
  typedef Interval_root<Data> Root;
  typedef double FT;
  struct Is_rational {
    typedef bool result_type;
    typedef Root argument_type;
    bool operator()(Root rt) const {
      return rt.is_rational();
    }
  };

  Is_rational is_rational_object() const {
    return Is_rational();
  }

  struct To_rational {
    typedef double result_type;
    typedef Root argument_type;
    double operator()(Root rt) const {
      return CGAL::to_interval(rt).second;
    }
  };

  To_rational to_rational_object() const {
    return To_rational();
  }

  struct Rational_between_roots {
    typedef double result_type;
    typedef Root first_argument_type;
    typedef Root second_argument_type;
    result_type operator()(first_argument_type a, 
			   second_argument_type b) const {
      CGAL_precondition(CGAL::compare(a,b) == CGAL::SMALLER);
      return .5*(CGAL::to_interval(a).second + CGAL::to_interval(b).first);
    }
  };

  Rational_between_roots rational_between_roots_object() const {
    return Rational_between_roots();
  }

  struct To_isolating_interval {
    typedef std::pair<double,double> result_type;
    typedef Root argument_type;
    result_type operator()(argument_type a) const {
      return CGAL::to_interval(a);
    }
  };

  To_isolating_interval to_isolating_interval_object() const {
    return To_isolating_interval();
  }
};

} } //namespace CGAL::Kinetic

#endif
