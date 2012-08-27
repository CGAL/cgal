// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_LAZY_UPPER_BOUND_ROOT_STACK_H
#define CGAL_POLYNOMIAL_LAZY_UPPER_BOUND_ROOT_STACK_H
#include <CGAL/Polynomial/Upper_bound_root_stack.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Polynomial/internal/interval_arithmetic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

/*
  Root has pointer to Root_rep, solver_rep (ref counted), a double approximation;

  Root_rep is
  - an optional Root
  - a pointer to the solver (non-ref_counted)

  Solver_rep is
  - an normal solver
  - pointer to current root (ref counted)

Top does:
- creates a new root rep with a pointer to the current root

Pop does:
- isolates current root if necessary
- puts it in to Root_rep
- creates new Root_rep (if more roots)

Compare:
- compares double approximations
- if root is not isolated, isolates root,
- updates double approx
- compares double approximation
- compares isolated root

Optimizations:

- keep popping intervals off of the front of the solver to try to
perform the comparison, stop when isolated. Try this only initially with a nt or an isolated root.
*/

struct Lazy_stats
{
    Lazy_stats(): created_(0), isolated_(0),
        refine_attempted_(0),
        refine_succeeded_(0){}
    int created_;
    int isolated_;
    int refine_attempted_;
    int refine_succeeded_;
};

Lazy_stats lazy_stats;

template <class Rep, class Solver_rep>
class Lazy_upper_bound_root_stack_root
{
    typedef Lazy_upper_bound_root_stack_root<Rep, Solver_rep> This;

    Comparison_result compare_double(const This &o) const
    {
/*if (iv_.sup() == std::numeric_limits<double>::infinity()
|| o.iv_.sup() == std::numeric_limits<double>::infinity()){
  int i=-5;++i;
  }*/
        if (iv_.inf() > o.iv_.sup()) return LARGER;
        else if (iv_.sup() < o.iv_.inf()) return SMALLER;
        else if (iv_.is_point() && o.iv_.is_point()) return EQUAL;
        else return UNKNOWN;
    }

    Comparison_result compare(const This &o) const
    {
        int cd= compare_double(o);
        if (cd != UNKNOWN) {
            return cd;
        }
        else {
            if (rep_->has_root() && !o.rep_->has_root()) {
                bool ret= o.srep_->refine_top(iv_.sup());
                o.update_interval();
                if (!ret) {
                    return SMALLER;
                }
            }
            else if (o.rep_->has_root() && !rep_->has_root()) {
                bool ret= srep_->refine_top(o.iv_.sup());
//update_interval();
                if (!ret) {
                    return LARGER;
                }
            }
            ensure_exact();
            o.ensure_exact();
            Comparison_result ret= rep_->root().compare(o.rep_->root());
            return ret;
        }
        CGAL_Polynomial_postcondition(0);
        return  UNKNOWN;
    }

    public:
        template <class NT>
            Lazy_upper_bound_root_stack_root(const NT &nt): iv_(CGAL_POLYNOMIAL_TO_INTERVAL(nt)),
            rep_(new Rep(nt)) {}

        Lazy_upper_bound_root_stack_root(double d): iv_(d), rep_(new Rep(d)) {}

        Lazy_upper_bound_root_stack_root(typename Solver_rep::Pointer sp): srep_(sp),
        rep_(sp->current_root_rep()) {
            iv_=srep_->double_interval();
//++lazy_stats.created_;
        }
        Lazy_upper_bound_root_stack_root(): iv_(std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()){}

        ~Lazy_upper_bound_root_stack_root() {
            if (srep_ && rep_->has_root()) {
//++lazy_stats.isolated_;
            }
        }

/*typedef typename Rep::Root Root;
  Lazy_upper_bound_root_stack_root(const Root &rt): iv_(CGAL_POLYNOMIAL_TO_INTERVAL(rt), rep_(rt) {
  }*/

        bool operator==(const This &o) const
        {
            return compare(o)==EQUAL;
        }
        bool operator!=(const This &o) const
        {
            return compare(o)!=EQUAL;
        }
        bool operator<(const This &o) const
        {
            return compare(o)==SMALLER;
        }
        bool operator>(const This &o) const
        {
            return compare(o)==LARGER;
        }
        bool operator>=(const This &o) const
        {
            return compare(o)!=SMALLER;
        }
        bool operator<=(const This &o) const
        {
            return compare(o)!=LARGER;
        }
        const typename Rep::Root &troot() const
        {
            ensure_exact();
            return rep_->root();
        }
        bool is_even_multiplicity() const
        {
            ensure_exact();
            return rep_->root().is_even_multiplicity();
        }
        This operator-() const
        {
            ensure_exact();
            This ret;
            ret.iv_= -iv_;
            ret.rep_= new Rep(-rep_->root());
            return ret;
        }
        double to_double() const
        {
            ensure_exact();
            return CGAL_POLYNOMIAL_TO_DOUBLE(rep_->root());
        }
        std::pair<double,double> double_interval(double tol= std::numeric_limits<double>::infinity()) const
        {
            if (tol != std::numeric_limits<double>::infinity()) {
                ensure_exact();
//iv_= CGAL_POLYNOMIAL_TO_INTERVAL(rep_->root());
            }
            return std::pair<double,double>(iv_.inf(), iv_.sup());
        }
        bool is_rational() const
        {
            ensure_exact();
            return rep_->root().is_rational();
        }
        typename Solver_rep::NT to_rational() const
        {
            ensure_exact();
            return rep_->root().to_rational();
        }
        static This infinity() {
            This ret;
            ret.iv_= std::numeric_limits<double>::infinity();
        }

        void write(std::ostream &out) const
        {
            if (!rep_) out << "null";
            if (rep_->has_root()) out << rep_->root();
            else out << "Unisolated (" << iv_ << ")";
        }

    protected:
        void ensure_exact() const
        {
            if (!rep_->has_root()) {
                srep_->check_current(rep_);
                srep_->isolate();
//++lazy_stats.isolated_;
            }
            iv_= CGAL_POLYNOMIAL_TO_INTERVAL(rep_->root());
        }

        void update_interval() const
        {
            if (rep_->has_root()) {
                iv_= rep_->root().double_interval(std::numeric_limits<double>::infinity());
            }
            else {
                srep_->check_current(rep_);
                iv_= srep_->double_interval();
            }
        }

        mutable Interval_nt iv_;
        mutable typename Solver_rep::Pointer srep_;
        mutable typename Rep::Pointer rep_;
};

template <class Rep, class Solver_rep>
std::ostream& operator<<(std::ostream &out, const Lazy_upper_bound_root_stack_root<Rep, Solver_rep> &rt)
{
    rt.write(out);
    return out;
}


} } } //namespace CGAL::POLYNOMIAL::internal
namespace CGAL { namespace POLYNOMIAL {

/*template <class Rep, class Solver_rep>
  double to_double(const Lazy_upper_bound_root_stack_root<Rep, Solver_rep> &rt){
  return rt.to_double();
  }

  template <class Rep, class Solver_rep>
  std::pair<double,double> to_interval(const Lazy_upper_bound_root_stack_root<Rep, Solver_rep> &rt){
  return rt.double_interval();
  }*/

template <class Traits_t>
class Lazy_upper_bound_root_stack
{
    typedef typename Traits_t::Root TRoot;
    protected:
        typedef Upper_bound_root_stack<Traits_t> Solver;
        struct Solver_rep;
        struct Root_rep: public CGAL::Kinetic::Ref_counted<Root_rep>
	{
	  typedef TRoot Root;
	  typedef CGAL::Kinetic::Ref_counted<Root_rep> RC;
	  Root_rep(): has_root_(false) {
	    CGAL_Polynomial_postcondition(!has_root());
	  }
	  template <class NT>
	  Root_rep(const NT &nt): root_(nt), has_root_(true){}
	  
	  Root_rep(double d): root_(d), has_root_(true){ /*CGAL_Polynomial_precondition(d != -std::numeric_limits<double>::infinity());*/}
	  
	  Root_rep(const Root_rep &o): RC(), root_(o.root_), has_root_(o.has_root_){}
	  Root_rep(const TRoot &tr): root_(tr), has_root_(true) {
	    //CGAL_Polynomial_precondition(tr != -std::numeric_limits<TRoot>::infinity());
	  }
	  
	  void set_root(const TRoot &tr) {
	    CGAL_Polynomial_precondition(tr != -std::numeric_limits<TRoot>::infinity());
	    CGAL_Polynomial_precondition(!has_root());
	    root_=tr;
	    has_root_=true;
	  }
	  
	  const TRoot& root() const
	  {
	    CGAL_Polynomial_precondition(has_root());
	    return root_;
	  }
	  
	  bool has_root() const
	  {
	    CGAL_precondition(std::numeric_limits<Root>::has_quiet_NaN);
	    return has_root_;
	  }
	  
	  /*std::pair<double,double> current_interval(typename Solver_rep::Pointer srep) const {
	    if (is_isolated()){
	    return root_.double_interval();
	    } else {
	    srep->check_current(this);
	    return srep->double_interval();
	    }
	    }*/
	  
	  mutable TRoot root_;
	  bool has_root_;
	};
        typedef Lazy_upper_bound_root_stack<Traits_t> This;
    public:
        typedef Traits_t Traits;
        typedef internal::Lazy_upper_bound_root_stack_root<Root_rep, Solver_rep> Root;
        typedef typename Traits::Function Function;

        Lazy_upper_bound_root_stack() {
        };
        Lazy_upper_bound_root_stack(const Function &f,
            const Root &lb,
            const Root &ub,
            const Traits &tr): rep_(new Solver_rep(f, lb, ub, tr)){}
        Lazy_upper_bound_root_stack(const Lazy_upper_bound_root_stack &o) {
            copy_data(o);
        }
        void pop() {
//CGAL_precondition(!empty());
            if (!empty()) rep_->pop();
        }
        const Root& top() const
        {
//    CGAL_precondition(!empty());
            static Root ret;
            if (rep_) {
                Root retc(rep_);
                ret=retc;
            }
            else {
                ret= std::numeric_limits<Root>::infinity();
            }
//CGAL_postcondition(ret != std::numeric_limits<Root>::infinity());
            return ret;
        }

        bool empty() const
        {
            return !rep_ || rep_->empty();
        }

        This operator=(const This &o) {
            copy_data(o);
            return *this;
        }

    protected:
        void copy_data(const This &o) {
            if (o.rep_) {
                rep_= new Solver_rep(*o.rep_);
            }
            else {
                rep_=NULL;
            }
        }

        mutable typename Solver_rep::Pointer rep_;
};

template <class Traits_t>
struct Lazy_upper_bound_root_stack<Traits_t>::Solver_rep: public CGAL::Kinetic::Ref_counted<Solver_rep>
{
    typedef Solver_rep This;
    typedef CGAL::Kinetic::Ref_counted<Solver_rep> RC;
    Solver_rep(const Function &f,
        const Root &lb,
        const Root &ub,
    const Traits &tr): solver_(f, lb.troot(), ub.troot(), tr, false) {
        cur_= new Root_rep();
        ++internal::lazy_stats.created_;
        if (solver_.empty()) {
            cur_->set_root(std::numeric_limits<TRoot>::infinity());
        }
    }

    typedef typename Traits_t::NT NT;

    Solver_rep(const Solver_rep &o): RC() {
        copy_data(o);
    }

    const This& operator=(const This &o) {
        copy_data(o);
        return *this;
    }

    void check_current(typename Root_rep::Pointer cur) const
    {
        CGAL_assertion(cur== cur_);
    }

    typename Root_rep::Pointer current_root_rep() const
    {
        CGAL_precondition(cur_);
        return cur_;
    }

    typename Root_rep::Pointer current_root_rep() {
        CGAL_precondition(cur_);
        return cur_;
    }

    void pop() {
        if (!is_isolated()) {
            isolate();
        }

        if (!solver_.empty()) {
            cur_= new Root_rep();
        }
        else {
            cur_=new Root_rep();
            cur_->set_root(std::numeric_limits<TRoot>::infinity());
            CGAL_Polynomial_postcondition(empty());
        }
    }

    std::pair<double, double> double_interval() const
    {
        if (cur_->has_root()) {
            return CGAL_POLYNOMIAL_TO_INTERVAL(cur_->root());
        }
        else {
/*CGAL_Polynomial_precondition(solver_.double_interval().first
!= std::numeric_limits<double>::infinity());*/
            return std::pair<double,double>(solver_.double_interval().first,
                std::numeric_limits<double>::infinity());
        }
    }

    bool empty() const
    {
        return cur_->has_root() &&  cur_->root() == std::numeric_limits<TRoot>::infinity();
    }

    void isolate() {
        solver_.isolate_top();
        set_root();
        CGAL_Polynomial_postcondition(cur_->has_root());
    }

    bool is_isolated() const
    {
        return cur_->has_root();
    }

    bool refine_bottom(double ub) {
        bool ret= solver_.refine_bottom(ub);
        if (solver_.top_is_isolated()) {
            set_root();
        }
        return ret;
    }

    bool refine_top(double lb) {
        ++internal::lazy_stats.refine_attempted_;
        bool ret= solver_.refine_top(lb);
        if (solver_.top_is_isolated()) {
            set_root();
        }
        else {
            ++internal::lazy_stats.refine_succeeded_;
//++lazy_stats.isolated_;
        }
        return ret;
    }

    protected:
        void copy_data(const This &o) {
            solver_=o.solver_;

            cur_= new Root_rep(*o.cur_);

        }
        void set_root() {
            CGAL_precondition(!cur_->has_root());
            if (solver_.empty()) {
                cur_->set_root(std::numeric_limits<TRoot>::infinity());
            }
            else {
                ++internal::lazy_stats.isolated_;
                cur_->set_root(solver_.top());
                solver_.pop_no_isolate();
            }
        }

        Solver solver_;
        typename Root_rep::Pointer cur_;
};
/*template <class Traits_t>
struct Lazy_upper_bound_root_stack<Traits_t>::Root_rep: public CGAL::Kinetic::Ref_counted<Root_rep>
{
    typedef TRoot Root;
    typedef CGAL::Kinetic::Ref_counted<Root_rep> RC;
    Root_rep(): has_root_(false) {
        CGAL_Polynomial_postcondition(!has_root());
    }
    template <class NT>
        Root_rep(const NT &nt): root_(nt), has_root_(true){}

	Root_rep(double d): root_(d), has_root_(true){}

    Root_rep(const Root_rep &o): RC(), root_(o.root_), has_root_(o.has_root_){}
    Root_rep(const TRoot &tr): root_(tr), has_root_(true) {
//CGAL_Polynomial_precondition(tr != -std::numeric_limits<TRoot>::infinity());
    }

    void set_root(const TRoot &tr) {
        CGAL_Polynomial_precondition(tr != -std::numeric_limits<TRoot>::infinity());
        CGAL_Polynomial_precondition(!has_root());
        root_=tr;
        has_root_=true;
    }

    const TRoot& root() const
    {
        CGAL_Polynomial_precondition(has_root());
        return root_;
    }

    bool has_root() const
    {
        CGAL_precondition(std::numeric_limits<Root>::has_quiet_NaN);
        return has_root_;
    }


    mutable TRoot root_;
    bool has_root_;
};*/


} } //namespace CGAL::POLYNOMIAL

namespace CGAL
{
    template <class Rep, class Solver_rep>
    double to_double(const CGAL_POLYNOMIAL_NS::internal::Lazy_upper_bound_root_stack_root<Rep, Solver_rep> &rt) {
        return rt.to_double();
    }

    template <class Rep, class Solver_rep>
    std::pair<double,double> to_interval(const typename CGAL_POLYNOMIAL_NS::internal::Lazy_upper_bound_root_stack_root<Rep, Solver_rep> &rt) {
        return rt.double_interval();
    }
}


namespace std
{
    template <class R, class SR>
        struct numeric_limits<CGAL_POLYNOMIAL_NS::internal::Lazy_upper_bound_root_stack_root<R, SR> >
    {
        typedef numeric_limits<typename R::Root> Rnl;
        typedef CGAL_POLYNOMIAL_NS::internal::Lazy_upper_bound_root_stack_root<R, SR> T;
        static const bool is_specialized = true;
        static T min BOOST_PREVENT_MACRO_SUBSTITUTION () throw () {return -infinity();}
        static T max BOOST_PREVENT_MACRO_SUBSTITUTION () throw () {return infinity();}
        static const int digits =0;
        static const int digits10 =0;
        static const bool is_signed = true;
        static const bool is_integer = false;
        static const bool is_exact = true;
        static const int radix =0;
        static T epsilon() throw(){return T(0);}
        static T round_error() throw(){return T(0);}
        static const int min_exponent=0;
        static const int min_exponent10=0;
        static const int max_exponent=0;
        static const int max_exponent10=0;
        static const bool has_infinity=true;
        static const bool has_quiet_NaN = true;
        static const bool has_signaling_NaN= false;
        static const float_denorm_style has_denorm= denorm_absent;
        static const bool has_denorm_loss = false;
        static T infinity() throw() {T inf(std::numeric_limits<double>::infinity()); return inf;}
        static T quiet_NaN() throw(){return T();}
        static T denorm_min() throw() {return T(0);}
        static const bool is_iec559=false;
        static const bool is_bounded =false;
        static const bool is_modulo= false;
        static const bool traps = false;
        static const bool tinyness_before =false;
        static const float_round_style round_stype = round_toward_zero;
    };
};
#endif
