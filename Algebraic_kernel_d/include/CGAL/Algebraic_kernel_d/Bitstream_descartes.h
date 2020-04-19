// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_GENERIC_DESCARTES
#define CGAL_GENERIC_DESCARTES 1

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_type_generator.h>

// NOTE: If this flag is set, you need EXACUS!
#if CGAL_ACK_BITSTREAM_USES_E08_TREE
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_E08_tree.h>
#else
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree.h>
#endif

#include <CGAL/Algebraic_kernel_d/exceptions.h>

namespace CGAL {

namespace internal {

//! enum to distinguish between different descartes instances
enum Bitstream_descartes_type {
  GENERIC_DESCARTES = 0,
  SQUARE_FREE_DESCARTES = 1, //!< uses Square_free_descartes_tag constructor
  M_K_DESCARTES = 2, //!< uses M_k_descartes_tag constructor
  BACKSHEAR_DESCARTES = 3, //!< uses Backshear_descartes_tag constructor
  VERT_LINE_ADAPTER_DESCARTES = 4 // ! < uses Vert_line_adapter_descartes
};

// pre-declaration
template<typename BitstreamDescartesRndlTreeTraits>
class Bitstream_descartes;

//! Tag for the square free Bitstream Descartes method
struct Square_free_descartes_tag {};

//! Tag for  the Bitstream m-k-Descartes method
struct M_k_descartes_tag {};

//! Tag for the Backshear Descartes method
struct Backshear_descartes_tag {};

//! Tag for the Exchange-Descartes method
struct Vert_line_adapter_descartes_tag {};

//! forward declaration
template<typename BitstreamDescartesRndlTreeTraits>
class Bitstream_descartes;


/*
 * \brief Thrown whenever a non-specialised virtual member function is called
 */
class Virtual_method_exception {};


/*
 * \brief The base class for all variants of the Bitstream Descartes method.
 *
 */
template<typename BitstreamDescartesRndlTreeTraits,
      typename Policy=CGAL::Handle_policy_no_union>
class Generic_descartes_rep
  : public Policy::template Hierarchy_base<CGAL_ALLOCATOR(char) >::Type {

public:

  //! The traits class for approximations
  typedef BitstreamDescartesRndlTreeTraits
  Bitstream_descartes_rndl_tree_traits;

  //! The Handle class
  typedef Bitstream_descartes<Bitstream_descartes_rndl_tree_traits> Handle;

  //! The Coeeficient type of the input polynomial
  typedef typename Bitstream_descartes_rndl_tree_traits::Coefficient
  Coefficient;

  //! The polynomial type
  typedef typename Bitstream_descartes_rndl_tree_traits::POLY Polynomial;

  typedef Generic_descartes_rep<Bitstream_descartes_rndl_tree_traits> Self;

  //! The type of the used Bitstream Descartes tree
#if CGAL_ACK_BITSTREAM_USES_E08_TREE
  typedef CGAL::internal::Bitstream_descartes_E08_tree
  <Bitstream_descartes_rndl_tree_traits>
  Bitstream_tree;
#else
  typedef CGAL::internal::Bitstream_descartes_rndl_tree
  <Bitstream_descartes_rndl_tree_traits>
  Bitstream_tree;
#endif


  //! The used integer type
  typedef typename Bitstream_descartes_rndl_tree_traits::Integer Integer;

  //! The type for the iterator of the nodes of the bitstream tree
  typedef typename Bitstream_tree::Node_iterator Node_iterator;

  //! The same as constant iterator
  typedef typename Bitstream_tree::Node_const_iterator Node_const_iterator;

  //! How the boundaries of the isolating intervals are represented
  typedef typename Bitstream_descartes_rndl_tree_traits::Bound Bound;

  //! Default constructor (does nothing)
  Generic_descartes_rep(Bitstream_descartes_type type = GENERIC_DESCARTES) :
    type_(type) {
  };

  /*!
   * Constructor computing an interval containing all real roots of \c f,
   * and initialising the Bitstream Descartes tree
   */
  Generic_descartes_rep(Bitstream_descartes_type type,
                        Polynomial f,
                        Bitstream_descartes_rndl_tree_traits traits) :
    type_(type),
    f_(f),
    traits_(traits),
    is_isolated_(false) {

    Integer lower,upper;
    long log_div;
    this->get_interval(f,lower,upper,log_div,traits);
    //AcX_DSTREAM("f: " << f << std::endl);
    if (CGAL::degree(f) > 0) {
      bitstream_tree
#if CGAL_ACK_BITSTREAM_USES_E08_TREE
        = Bitstream_tree(-log_div,
                         f.begin(),
                         f.end(),
                         typename Bitstream_tree::Monomial_basis_tag(),
                         traits);
#else
      = Bitstream_tree(lower,upper,log_div,
                       f.begin(),
                       f.end(),
                       typename Bitstream_tree::Monomial_basis_tag(),
                       traits);
#endif

      if (bitstream_tree.begin() == bitstream_tree.end()) {
        number_of_intervals = 0;
      } else {
        number_of_intervals = 1;
      }
    } else {
      number_of_intervals=0;
    }
  }

  /*!
   * Constructor that copies the Bitstream tree given from outside
   * and initialising the Bitstream Descartes tree
   * The tree must "fit" to the polynomial
   */
  Generic_descartes_rep(Bitstream_descartes_type type,
                        Polynomial f,
                        Bitstream_tree tree,
                        Bitstream_descartes_rndl_tree_traits traits) :
    type_(type),
    f_(f),
    traits_(traits),
    bitstream_tree(tree),
    is_isolated_(false) {

    tree.set_traits(traits);

    number_of_intervals = 0;

    for(Node_iterator curr = bitstream_tree.begin();
        curr != bitstream_tree.end();
        curr++) {
      number_of_intervals++;
    }
  }

  //! Destructor (does nothing)
  virtual ~Generic_descartes_rep() {
  }

  //! Needed for the referencing counting mechanism
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Generic_descartes_rep(*this);
  }

  /*!
   * \brief Computes a better approximation of the \c i th root of the
   * polynomial
   */
  virtual void refine_interval(int i) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < number_of_intervals);
    Node_iterator curr = bitstream_tree.begin(), begin, end,
      new_begin, helper;
    std::advance(curr,i);
    int intervals = 1;
    end = curr;
    ++end;
    begin=curr;
    do {
      //std::cout << bitstream_tree.lower(begin) << " " << bitstream_tree.upper(begin) << std::endl;
      //std::cout << bitstream_tree.min_var(begin) << " " << bitstream_tree.max_var(begin) << std::endl;
      int new_intervals = bitstream_tree.subdivide(begin,new_begin,helper);
      intervals += new_intervals-1;
      begin = new_begin;
      curr = helper;

      // Fixes the bug when a interval splits, and the leftmost subintervals
      // has no children with sign variation >=1
      if (intervals == 1) {
        break;
      }
      if (new_intervals == 0) {
        continue;
      }

      while(curr != end) {
        intervals += bitstream_tree.subdivide(curr,new_begin,helper)-1;
        curr = helper;
      }

    }
    while (intervals != 1);
    //std::cout << "Refined " << left_bound(i) << " " << right_bound(i) << std::endl;

  }

  /*!
   * \brief isolates the root of \c f
   *
   * The mechanism is the following: The \c bitstream_tree member of the
   * object is transformed via subdivision until the
   * \c termination_condition routine of the object returns true. When this
   * happens, the \c process_nodes routine of the object is called.
   */
  virtual void isolate() {

    //AcX_DSTREAM("Starting isolation" << std::endl);

    Node_iterator curr = bitstream_tree.begin(),dummy,new_curr;

    if(curr == bitstream_tree.end()) {
      is_isolated_ = true;
      return;
    }

    int newly_created;

    while (!this->termination_condition()) {

      if (curr == bitstream_tree.end()) {
        curr = bitstream_tree.begin();
      }

      if (bitstream_tree.max_var(curr) == 1) {
        ++curr;
      }
      else {
        //AcX_DSTREAM("Subdivision at "
        //<< CGAL::to_double(bitstream_tree.lower(curr)) << " "
        //<< CGAL::to_double(bitstream_tree.upper(curr)) << std::flush);
        newly_created = bitstream_tree.subdivide(curr,dummy,new_curr);
        number_of_intervals += newly_created-1;
        curr = new_curr;
        //AcX_DSTREAM("done" << std::endl);
      }

    }
    this->process_nodes();
    is_isolated_ = true;
  }

  /*!
   * \brief Computes an interval containing all real roots of \c p,
   * using the Fujiwara root bound.
   *
   * So far, the \c log_div variable is always set to zero, this means
   * that <i>[lower,upper]</i> is the interval containing all real roots
   */
  virtual void get_interval(const Polynomial& p, Integer& lower,
                            Integer& upper, long& log_div,
                            Bitstream_descartes_rndl_tree_traits traits) {


    typename Bitstream_descartes_rndl_tree_traits::Lower_bound_log2_abs
      lower_bound_log2_abs = traits.lower_bound_log2_abs_object();
    typename
      Bitstream_descartes_rndl_tree_traits::Upper_bound_log2_abs_approximator
      upper_bound_log2_abs_approximator
      = traits.upper_bound_log2_abs_approximator_object();
    //AcX_DSTREAM("Fujiwara bound.." << p <<  std::endl);
#if CGAL_ACK_BITSTREAM_USES_E08_TREE
    log_div = -CGAL::internal::Fujiwara_root_bound_log
      (p.begin(),
       p.end(),
       lower_bound_log2_abs,
       upper_bound_log2_abs_approximator
      );
#else

    log_div = -CGAL::internal
      ::Fujiwara_root_bound_log
      (p.begin(),
       p.end(),
       lower_bound_log2_abs,
       upper_bound_log2_abs_approximator
      );
#endif

    //AcX_DSTREAM("Fujiwara returns " << log_div << std::endl);
    // To be sure
    log_div--;
    lower=Integer(-1);
    upper=Integer(1);
    return;

  }

  //! returns the number of detected isolating intervals
  virtual int number_of_real_roots() const {
    return number_of_intervals;
  }

  //! The lower bound of the \c i th root
  virtual Bound left_bound(int i) const  {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < number_of_intervals);
    Node_const_iterator curr = bitstream_tree.begin();
    std::advance(curr,i);
    return bitstream_tree.lower(curr);
  }

  //! The upper bound of the \c i th root
  virtual Bound right_bound(int i) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < number_of_intervals);
    Node_const_iterator curr = bitstream_tree.begin();
    std::advance(curr,i);
    return bitstream_tree.upper(curr);
  }

  //! Returns the polynomial which is isolated
  Polynomial polynomial() const {
    return f_;
  }

  /*!
   * \brief When does the isolation algorithm terminate?
   *
   * This method must be specialised by derived classes
   */
  virtual bool termination_condition() {
    throw Virtual_method_exception();
    return false;
  }

  /*!
   * \brief Gives an opportunity to process the nodes after
   * the subdivision steps are finished
   *
   * This method must be specialised by derived classes, but can
   * remain empty in many cases.
   */
  virtual void process_nodes() {
    throw Virtual_method_exception();
    return;
  }

  /*! \brief Returns whether the \c i th root is definitely a simple root
   * of the isolated polynomial
   *
   * Must be specialised by derived class
   */
  virtual bool is_certainly_simple_root(int) const {
    throw Virtual_method_exception();
    return false;
  }

  /*! \brief Returns whether the \c i th root is definitely a multiple root
   * of the isolated polynomial
   *
   * Must be specialised by derived class
   */
  virtual bool is_certainly_multiple_root(int) const {
    throw Virtual_method_exception();
    return false;
  }


  virtual int multiplicity_of_root(int CGAL_assertion_code(i)) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < number_of_intervals);
    return -1;
  }

  virtual int get_upper_bound_for_multiplicity(int i) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < number_of_intervals);
    Node_const_iterator curr = bitstream_tree.begin();
    std::advance(curr,i);
    return bitstream_tree.min_var(curr);
  }

  //! Must be specialized by the derived class
  virtual int degree_of_gcd() const {
    throw Virtual_method_exception();
    return -1;
  }

  //! Must be specialized by the derived class
  virtual Polynomial square_free_part() const {
    throw Virtual_method_exception();
    return Polynomial();
  }

  //! Must be specialized by the derived class
  virtual Handle inverse_transform_isolator() const {
    throw Virtual_method_exception();
    return Handle();
  }

  bool is_isolated() const {
    return is_isolated_;
  }

  Bitstream_descartes_rndl_tree_traits traits() const {
    return traits_;
  }

  Bitstream_tree get_tree() const {

    return bitstream_tree;

  }

  //! type to distinguish used constructor
  Bitstream_descartes_type type_;

protected:

  //! Polynomial which is isolated
  Polynomial f_;

  //! The traits class
  Bitstream_descartes_rndl_tree_traits traits_;

  //! The tree of the Bitstream Descartes method
  mutable Bitstream_tree bitstream_tree;

  //! The number of detected isolating intervals
  int number_of_intervals;

  //! Has isolation already taken place
  mutable bool is_isolated_;

};

/*
 * \brief Representation for square free polynomials
 */
template<typename BitstreamDescartesRndlTreeTraits,
      typename Policy=CGAL::Handle_policy_no_union>
class Square_free_descartes_rep
  : public Generic_descartes_rep<BitstreamDescartesRndlTreeTraits> {


public:

  //! Traits type
  typedef BitstreamDescartesRndlTreeTraits
  Bitstream_descartes_rndl_tree_traits;

  //! The generic representation
  typedef Generic_descartes_rep<BitstreamDescartesRndlTreeTraits,
        Policy> Base;

  //! Polynomial type
  typedef typename Base::Polynomial Polynomial;

  //! Iterator for the leaves in the bitstream tree
  typedef typename Base::Node_iterator Node_iterator;

  //! The type of the tree that controls the Bitstream instance
  typedef typename Base::Bitstream_tree Bitstream_tree;

  /*!
   * \brief Constructor with the square free polynomial <tt>f<tt>.
   */
  Square_free_descartes_rep(
      Polynomial f,
      Bitstream_descartes_rndl_tree_traits traits) :
    Base(SQUARE_FREE_DESCARTES, f,traits) {
  }

  /*!
   * \brief Constructor with the square free polynomial <tt>f<tt>.
   */
  Square_free_descartes_rep(
      Polynomial f,
      Bitstream_tree tree,
      Bitstream_descartes_rndl_tree_traits traits) :
    Base(SQUARE_FREE_DESCARTES, f, tree, traits) {
  }

  //! Needed for reference counting
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Square_free_descartes_rep(*this);
  }

  /*!
   * \brief Terminates when all detected roots are simple
   */
  virtual bool termination_condition() {
    for(Node_iterator curr=Base::bitstream_tree.begin();
        curr != Base::bitstream_tree.end();curr++) {
      if(Base::bitstream_tree.max_var(curr)!=1) {
        return false;
      }
    }
    return true;
  }

  //! nothing to do here
  virtual void process_nodes() {
    return;
  }

  //! Polynomial is square free, so gcd is 1
  virtual int degree_of_gcd() const {
    return 0;
  }

  //! Polynomial is square free
  virtual Polynomial square_free_part() const {
    return this->f_;
  }

  //! Always true
  virtual bool is_certainly_simple_root(int ) const {
    return true;
  }

  //! Always false
  virtual bool is_certainly_multiple_root(int ) const {
    return false;
  }

};

/*
 * \brief Representation for polynomials with at most one multiple root
 */
template<typename BitstreamDescartesRndlTreeTraits,
      typename Policy=CGAL::Handle_policy_no_union>
class M_k_descartes_rep
  : public Generic_descartes_rep<BitstreamDescartesRndlTreeTraits> {


public:

  //! Traits class
  typedef BitstreamDescartesRndlTreeTraits
  Bitstream_descartes_rndl_tree_traits;

  //! Generic representation
  typedef Generic_descartes_rep<BitstreamDescartesRndlTreeTraits,
                                        Policy> Base;

  //! Polynomial type
  typedef typename Base::Polynomial Polynomial;

  //! Iterator for the leaves of the Bitstream Descartes tree
  typedef typename Base::Node_iterator Node_iterator;

  //! Constant iterator for the leaves
  typedef typename Base::Node_const_iterator Node_const_iterator;

  //! The interval boundaries are represented in this type
  typedef typename  Bitstream_descartes_rndl_tree_traits::Bound
  Bound;

  //! The type of the tree that controls the Bitstream instance
  typedef typename Base::Bitstream_tree Bitstream_tree;

  /*!
   * \brief Constructor for a polynomial <tt>f<tt>, not necessarily square
   * free
   *
   * The values <tt>m</tt>
   * and <tt>k</tt> need to be the exact number of real roots of <tt>f</tt>
   * counted without multiplicity and the degree of the greatest common
   * divisor of <tt>f</tt> with its partial derivative, respectively.
   */
  M_k_descartes_rep(Polynomial f,int m, int k,
                    Bitstream_descartes_rndl_tree_traits traits) :
    Base(M_K_DESCARTES, f,traits),
    number_of_roots(m),
    gcd_degree(k),
    index_of_multiple(-1) {
  }

  M_k_descartes_rep(Polynomial f,int m, int k,
                    Bitstream_tree tree,
                    Bitstream_descartes_rndl_tree_traits traits) :
    Base(M_K_DESCARTES, f, tree, traits),
    number_of_roots(m),
    gcd_degree(k),
    index_of_multiple(-1) {
  }

  //! Default constructor
  M_k_descartes_rep() {
  }

  //! Needed for reference counting
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new M_k_descartes_rep(*this);
  }

  /*!
   * \brief Termination condition
   *
   * If <tt>m-1</tt> simple and one more leaf is detected, the Bitstream
   * Descartes method is stopped. If the minimal sign
   * variation drops under <tt>k</tt> in each leaf, a
   * \c Non_generic_position_exception is thrown.
   */
  virtual bool termination_condition() {
    int counted_simple_roots=0;
    int max_max_var = 0;
    for(Node_iterator curr=Base::bitstream_tree.begin();
        curr != Base::bitstream_tree.end(); curr++) {
      int max_var = Base::bitstream_tree.max_var(curr);
      if(max_var > max_max_var) {
        max_max_var = max_var;
      }
      if(max_var == 1) { // && Base::bitstream_tree.max_var(curr)==1) {
        ++counted_simple_roots;
      }
    }
    //AcX_DSTREAM("Situation: " << this->number_of_intervals << " intervals " << this->number_of_roots << " are expected" << std::endl);
    if (this->number_of_intervals == this->number_of_roots
        && counted_simple_roots >= number_of_roots-1) {
      return true;
    }
    if (max_max_var <= gcd_degree) {
      throw CGAL::internal::Non_generic_position_exception();
    }

    return false;

  }

  //! The index of the (possibly) multiple root is computed here.
  virtual void process_nodes() {
    int i = 0;
    for (Node_iterator curr=Base::bitstream_tree.begin();
        curr != Base::bitstream_tree.end(); curr++) {
      if(Base::bitstream_tree.max_var(curr) > 1 ) {
        index_of_multiple = i;
        return;
      } else {
        ++i;
      }
    }
    return;
  }

  //! Returns k
  virtual int degree_of_gcd() const {
    return gcd_degree;
  }

  //! True for all roots except for the candidate
  virtual bool is_certainly_simple_root(int i) const {
    return (i!=index_of_multiple);
  }

  //! Always false
  virtual bool is_certainly_multiple_root(int) const {
    return false;
  }


protected:

  //! The "m"
  int number_of_roots;

  //! The "k"
  int gcd_degree;

  //! The candidate's index
  int index_of_multiple;

};


template<typename BitstreamDescartesRndlTreeTraits,
             typename EventRefinement,
             typename Policy=CGAL::Handle_policy_no_union>
class Backshear_descartes_rep
  : public Generic_descartes_rep<BitstreamDescartesRndlTreeTraits> {


public:

  typedef EventRefinement Event_refinement;

  typedef BitstreamDescartesRndlTreeTraits
  Bitstream_descartes_rndl_tree_traits;

  typedef Generic_descartes_rep<BitstreamDescartesRndlTreeTraits,
                                        Policy> Base;

  typedef typename Base::Polynomial Polynomial;

  typedef typename Base::Node_iterator Node_iterator;

  typedef std::list<int>::iterator Marking_iterator;

  typedef std::list<int>::const_iterator Marking_const_iterator;

  typedef typename Base::Node_const_iterator Node_const_iterator;

  typedef typename Base::Bound Bound;

  Backshear_descartes_rep(
      Polynomial f,
      int number_of_non_event_points,
      int number_of_events,
      Event_refinement event_refinement,
      Bitstream_descartes_rndl_tree_traits traits) :
    Base(BACKSHEAR_DESCARTES,f,traits),
    number_of_non_event_points(number_of_non_event_points),
    number_of_events(number_of_events),
    event_refinement(event_refinement) {
  }

  Backshear_descartes_rep() {
  }

  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Backshear_descartes_rep(*this);
  }

  virtual void isolate() {

    Node_iterator curr = Base::bitstream_tree.begin(),sub_begin,new_curr;

    if(curr == Base::bitstream_tree.end()) {
      this->is_isolated_ = true;
      return;
    }
    markings.clear();
    markings.push_back(this->check_marking(curr));

    Marking_iterator curr_mark = markings.begin(),mark_helper;

    int newly_created;

    while(! this->termination_condition()) {
      //AcX_DSTREAM("Subdivision..." << number_of_intervals << std::endl);
      if (curr == Base::bitstream_tree.end()) {
        curr = Base::bitstream_tree.begin();
        CGAL_assertion(curr_mark == markings.end());
        curr_mark = markings.begin();
      }
      if(Base::bitstream_tree.max_var(curr) == 1) {
        ++curr;
        ++curr_mark;
        //AcX_DSTREAM("nothing happend" << std::endl);
      }
      else {
        newly_created =
          Base::bitstream_tree.subdivide(curr,sub_begin,new_curr);
        mark_helper = markings.erase(curr_mark);
        curr_mark = mark_helper;
        for(Node_iterator tmp_curr = sub_begin;
            tmp_curr != new_curr;
            tmp_curr++) {
          markings.insert(curr_mark,check_marking(tmp_curr));
        }
        Base::number_of_intervals += newly_created-1;
        curr = new_curr;
        //AcX_DSTREAM(newly_created << " new intervals, marking size: " << markings.size() << std::endl);

      }
    }
    this->process_nodes();
    this->is_isolated_ = true;
  }



  virtual bool termination_condition() {
    int marked_intervals = 0;
    int unmarked_odd_intervals = 0;
    Node_iterator curr = Base::bitstream_tree.begin();
    Marking_iterator curr_mark = markings.begin();
    for(;curr != Base::bitstream_tree.end(); curr++) {
      if((*curr_mark) >= 0) {
        ++marked_intervals;
      }
      else {
        if (Base::bitstream_tree.min_var(curr) % 2 == 1) { // odd
          ++unmarked_odd_intervals;
        }
      }
      ++curr_mark;
    }
    CGAL_assertion(curr_mark == markings.end());
    return ((marked_intervals == number_of_events)
            && (unmarked_odd_intervals == number_of_non_event_points));
  }

  virtual void process_nodes() {
    Node_iterator curr=Base::bitstream_tree.begin(),curr_helper;
    Marking_iterator curr_mark = markings.begin();
    while(curr!=Base::bitstream_tree.end()) {
      if(((*curr_mark) == -1) &&
         (Base::bitstream_tree.min_var(curr) % 2 == 0)) {
        ++curr;
        curr_helper = curr;
        curr_helper--;
        Base::bitstream_tree.erase(curr_helper);
        curr_mark = markings.erase(curr_mark);
        Base::number_of_intervals--;
      } else {
        ++curr_mark;
        ++curr;
      }
    }
    CGAL_assertion(curr_mark == markings.end());

    //AcX_DSTREAM(markings.size() << " " << number_of_non_event_points << " " << number_of_events << std::endl);
    CGAL_assertion(static_cast<int>(markings.size())
                   ==number_of_non_event_points + number_of_events);
    return;
  }

  virtual bool is_certainly_simple_root(int i) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < Base::number_of_intervals);
    Node_const_iterator curr=Base::bitstream_tree.begin();
    std::advance(curr,i);
    return (Base::bitstream_tree.max_var(curr) == 1);
  }

  virtual bool is_certainly_multiple_root(int i) const {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < Base::number_of_intervals);
    Marking_const_iterator curr = markings.begin();
    std::advance(curr,i);
    return (*curr>=0);
  }


protected:

  int number_of_non_event_points;

  int number_of_events;

  Event_refinement event_refinement;

  std::list<int> markings;

protected:

  int check_marking(Node_iterator node) {
    Bound lower = Base::bitstream_tree.lower(node),
      upper = Base::bitstream_tree.upper(node);
    for(int i = 0; i < number_of_events; i++) {
      while(true) {
        if(CGAL::compare(event_refinement.lower_bound(i),lower)
           !=CGAL::NEGATIVE
           &&
           CGAL::compare(event_refinement.upper_bound(i),upper)
           !=CGAL::POSITIVE) {
          //Event inside the interval
          return i;
        }
        if(CGAL::compare(event_refinement.lower_bound(i),upper)
           ==CGAL::POSITIVE
           ||
           CGAL::compare(event_refinement.upper_bound(i),lower)
           ==CGAL::NEGATIVE) {
          //This event is outside
          break;
        }
        event_refinement.refine(i);

      }
    }
    return -1;
  }

};

/*
 * \brief Adaptor for roots of a vert line
 * (needed as dummy in surface analysis)
 *
 */
template<typename BitstreamDescartesRndlTreeTraits,
      typename VertLine,
      typename Policy=CGAL::Handle_policy_no_union>
class Vert_line_adapter_descartes_rep
  : public Generic_descartes_rep<BitstreamDescartesRndlTreeTraits> {

public:

  //! The traits class for approximations
  typedef BitstreamDescartesRndlTreeTraits
  Bitstream_descartes_rndl_tree_traits;

  //! type of vert line
  typedef VertLine Vert_line;

  //! type of Curve_analysis_2
  typedef typename Vert_line::Curve_analysis_2 Curve_analysis_2;

  //! type of Curve_kernel_2;
  typedef typename Curve_analysis_2::Algebraic_kernel_with_analysis_2
  Curve_kernel_2;

  //! The Coeeficient type of the input polynomial
  typedef typename Bitstream_descartes_rndl_tree_traits::Coefficient
  Coefficient;

  //! The polynomial type
  typedef typename CGAL::Polynomial_type_generator<Coefficient,1>::Type
  Polynomial;

  typedef Vert_line_adapter_descartes_rep
  <Bitstream_descartes_rndl_tree_traits, Vert_line, Policy> Self;

  //! The used integer type
  typedef typename Bitstream_descartes_rndl_tree_traits::Integer Integer;

  //! How the boundaries of the isolating intervals are represented
  typedef typename Bitstream_descartes_rndl_tree_traits::Bound
  Bound;

  typedef Generic_descartes_rep<Bitstream_descartes_rndl_tree_traits>
  Base;

  //! The type for the inverse isolator
  typedef typename Base::Handle Handle;

  /*!
   * \brief Constructor
   */
  template<typename InputIterator>
  Vert_line_adapter_descartes_rep(InputIterator begin,
                                  InputIterator end,
                                  Bitstream_descartes_rndl_tree_traits traits)
    : Base(VERT_LINE_ADAPTER_DESCARTES)
  {
    for (InputIterator it = begin; it != end; it++) {
      root_vec.push_back(std::make_pair(*it, 4));
    }

    this->is_isolated_ = true;
    this->traits_ = traits;
    this->f_  = Polynomial(0);
    this->number_of_intervals
      = static_cast<int>(root_vec.size());
    // Isolate all real roots until intervals are disjoint:
    for (int i = 1; i < this->number_of_real_roots(); i++ ){
      while(left_bound(i) < right_bound(i-1) ) {
        if (right_bound(i)-left_bound(i) <
            right_bound(i-1) - left_bound(i-1) ) {
          refine_interval(i-1);
        } else {
          refine_interval(i);
        }
      }

    }
  }


  //! Destructor (does nothing)
  virtual ~Vert_line_adapter_descartes_rep() {
  }

  //! Needed for the referencing counting mechanism
  virtual CGAL::Reference_counted_hierarchy<>* clone() {
    return new Vert_line_adapter_descartes_rep(*this);
  }

  virtual void refine_interval(int i) const {
    root_vec[i] = std::make_pair(root_vec[i].first, root_vec[i].second * 2);
  }

  virtual void isolate() const {
  }


  //! The lower bound of the \c i th root
  virtual Bound left_bound(int i) const  {
    typename Curve_kernel_2::Approximate_absolute_y_2
      approx_y = Base::traits_.point().xy().kernel()
      ->approximate_absolute_y_2_object();
    return approx_y(root_vec[i].first.first.algebraic_real_2
                    (root_vec[i].first.second),
                    root_vec[i].second).first;
  }

  //! The upper bound of the \c i th root
  virtual Bound right_bound(int i) const {
    typename Curve_kernel_2::Approximate_absolute_y_2
      approx_y = Base::traits_.point().xy().kernel()
      ->approximate_absolute_y_2_object();
    return approx_y(root_vec[i].first.first.algebraic_real_2
                    (root_vec[i].first.second),
                    root_vec[i].second).second;
  }


  /*! \brief Returns whether the \c i th root is definitely a simple root
   * of the isolated polynomial
   *
   */
  virtual bool is_certainly_simple_root(int /* i */) const {
    return false;
  }

  /*! \brief Returns whether the \c i th root is definitely
   * a multiple root
   * of the isolated polynomial
   *
   */
  virtual bool is_certainly_multiple_root(int /* i */) const {
    return false;
  }

protected:

  //! Roots stored as pair of a AcX::Vert_line and an integer denoting the
  //! index. Also, current precision of each root is stored
  mutable std::vector<std::pair<std::pair<Vert_line, int>,int> > root_vec;


};

/*!
 * \brief Class for the Bitstream Descartes method
 *
 * Class for the real root isolation of polynomials, using the Bitstream
 * Descartes method. The polynomials coefficient type is arbitrary, the
 * approximations of the coefficient type are obtained with the
 * \c BitstreamDescartesRndlTreeTraits parameter. For the requirements
 * of this traits class, see the documentation of
 * CGAL::Bitstream_descartes_rndl_tree.
 *
 * Internally, an instance of CGAL::Bitstream_descartes_rndl_tree is explored
 * in a specific way. That exploration strategy depends on the constructor
 * that is used to create the object. A tag is passed that defines the
 * variant of the Bitstream Descartes method: The Square_free_descartes_tag
 * starts the usual Bitstream method for square free integer polynomials.
 * With the M_k_descartes tag, it is able to handle one multiple root in
 * favourable situations, the Backshear_descartes_tag allows to isolate
 * even more complicated polynomials, if the multiple roots with even
 * multiplicity can be refined from outside. See the corresponding
 * constructors for more information.
 *
 */
template<typename BitstreamDescartesRndlTreeTraits>
class Bitstream_descartes
  : ::CGAL::Handle_with_policy<
    CGAL::internal::Generic_descartes_rep<BitstreamDescartesRndlTreeTraits> > {

public:

  //! Traits class
  typedef BitstreamDescartesRndlTreeTraits
  Bitstream_descartes_rndl_tree_traits;

  // The generic representation class
  typedef
  CGAL::internal::Generic_descartes_rep<BitstreamDescartesRndlTreeTraits> Rep;

  // The Handle type
  typedef ::CGAL::Handle_with_policy<Rep> Base;

  //! The coefficients of the polynomial
  typedef typename Bitstream_descartes_rndl_tree_traits::Coefficient
  Coefficient;

  //! The polynomial's type
  typedef typename CGAL::Polynomial_type_generator<Coefficient,1>::Type
  Polynomial;

  typedef Bitstream_descartes<Bitstream_descartes_rndl_tree_traits>
  Self;

  // Type for the Bitstream Descartes tree
#if CGAL_ACK_BITSTREAM_USES_E08_TREE
  typedef CGAL::internal::Bitstream_descartes_E08_tree
  <Bitstream_descartes_rndl_tree_traits>
  Bitstream_tree;
#else
  typedef CGAL::internal::Bitstream_descartes_rndl_tree
  <Bitstream_descartes_rndl_tree_traits>
  Bitstream_tree;
#endif

  //! Type for Integers
  typedef typename Bitstream_descartes_rndl_tree_traits::Integer Integer;

  //! Iterator type for the leaves of the Descartes tree
  typedef typename Bitstream_tree::Node_iterator
  Node_iterator;

  //! Const iterator for the leaves
  typedef typename Bitstream_tree::Node_const_iterator
  Node_const_iterator;

  //! Type for the interval boundaries of the isolating intervals
  typedef typename Bitstream_descartes_rndl_tree_traits::Bound
  Bound;

  //! Default constructor
  Bitstream_descartes() : Base(new Rep()) {}

  //! Copy constructor
  Bitstream_descartes(const Self& other) : Base(static_cast<const Base&>(other))
  {}

  /*!
   * \brief Constructor for a polynomial \c f
   *
   * See the documentation of the constrctor
   * with \c Square_free_descartes_tag
   */
  Bitstream_descartes(Polynomial f,
                      Bitstream_descartes_rndl_tree_traits traits
                      = Bitstream_descartes_rndl_tree_traits(),
                      bool isolate=true)
    : Base(new CGAL::internal::Square_free_descartes_rep
           <Bitstream_descartes_rndl_tree_traits>(f,traits))
  {
    if (isolate) {
      this->isolate();
    }
  }

  /*!
   * \brief Constructor for the square free Descartes method
   *
   * The polynomial \c f must not have multiple real roots. The
   * Bitstream Descartes tree is traversed in a bfs manner until
   * all leaves have sign variation zero or one.
   */
  Bitstream_descartes(Square_free_descartes_tag ,
                      Polynomial f,
                      Bitstream_descartes_rndl_tree_traits traits
                      = Bitstream_descartes_rndl_tree_traits(),
                      bool isolate=true)
    : Base(new CGAL::internal::Square_free_descartes_rep
           <Bitstream_descartes_rndl_tree_traits>(f,traits))
  {
    if (isolate) {
      this->isolate();
    }
  }

  /*!
   * \brief Constructor for the square free Descartes method,
   * using a precomputed tree
   *
   * The polynomial \c f must not have multiple real roots. The
   * Bitstream Descartes tree is traversed in a bfs manner until
   * all leaves have sign variation zero or one.
   * The tree must be adequate for the polynomial.
   * Use that constructor only if you know what you're doing!
   */
  Bitstream_descartes(Square_free_descartes_tag ,
                      Polynomial f,
                      Bitstream_tree tree,
                      Bitstream_descartes_rndl_tree_traits traits
                      = Bitstream_descartes_rndl_tree_traits(),
                      bool isolate=true)
    : Base(new CGAL::internal::Square_free_descartes_rep
           <Bitstream_descartes_rndl_tree_traits>(f, tree, traits))
  {
    if(isolate) {
      this->isolate();
    }
  }

  /*!
   * \brief Constructor for the m-k-Descartes method
   *
   * The polynomial \c f must have exactly \c m real roots, counted without
   * multiplicity, and the degree of <tt>gcd(f,f')</tt> must be \c k. In this
   * case, the constructor either isolates the real roots of \c f sucessfully
   * or a Non_generic_position_exception is thrown. Such an exception
   * certainly occurs if \c f has more than one multiple real root. If \c f
   * has at most one multiple root over the complex numbers, the roots are
   * certainly isolated with success.
   */
  Bitstream_descartes(M_k_descartes_tag /* t */,
                      Polynomial f,int m,int k,
                      Bitstream_descartes_rndl_tree_traits traits
                      = Bitstream_descartes_rndl_tree_traits(),
                      bool isolate = true)
    : Base(new CGAL::internal::M_k_descartes_rep
           <Bitstream_descartes_rndl_tree_traits>(f,m,k,traits))
  {
    if (isolate) {
      this->isolate();
    }
  }


  Bitstream_descartes(M_k_descartes_tag /* t */,
                      Polynomial f,int m,int k,
                      Bitstream_tree tree,
                      Bitstream_descartes_rndl_tree_traits traits
                      = Bitstream_descartes_rndl_tree_traits(),
                      bool isolate = true)
    : Base(new CGAL::internal::M_k_descartes_rep
           <Bitstream_descartes_rndl_tree_traits>(f,m,k,tree,traits))
  {
    if (isolate) {
      this->isolate();
    }
  }


  /*!
   * \brief Constructor for the Backshear-Decartes method
   *
   * The polynomial \c f must have exactly \c number_of_real_roots
   * many real roots, counted without multiplicity. Additionally, a set of
   * \c number_of_events root can be refined to arbitrary precision with the
   * \c event_refinement object. This must support three operations
   * for each <tt>0<=i<number_of_events</tt>:
   * <ul><li>lower_bound(i), upper_bound(i) gives an interval (not
   * necessarily isolating) of some root of \c f</li>
   * <li>refine(i) refines the corresponding interval</li></ul>
   * Note that the roots in \c event_refinement need not be sorted. All roots
   * which are not covered by \c event_refinement must have odd multiplicity.
   */
  template<typename EventRefinement>
  Bitstream_descartes(Backshear_descartes_tag ,
                      Polynomial f,
                      int number_of_real_roots,
                      int number_of_events,
                      EventRefinement event_refinement,
                      Bitstream_descartes_rndl_tree_traits traits
                      = Bitstream_descartes_rndl_tree_traits(),
                      bool isolate = true)
    : Base(new
           CGAL::internal::Backshear_descartes_rep
           <Bitstream_descartes_rndl_tree_traits,EventRefinement>
           (f,number_of_real_roots-number_of_events,
            number_of_events,event_refinement,traits))
  {
    if (isolate) {
      this->isolate();
    }
  }


  /*!
   * \brief Constructor for the Vert-line-adapter-Descartes method
   *
   */
  template<typename InputIterator>
  Bitstream_descartes(Vert_line_adapter_descartes_tag /* t */,
                      InputIterator begin,
                      InputIterator end,
                      Bitstream_descartes_rndl_tree_traits traits)
    : Base(new CGAL::internal::Vert_line_adapter_descartes_rep
           <Bitstream_descartes_rndl_tree_traits,
           typename InputIterator::value_type::first_type>
           (begin, end, traits) )
  {
    // No isolation necessary
  }

  //! return the type of the used descartes method
  Bitstream_descartes_type type() const {
    return this->ptr()->type_;
  }

  //! Return the polynomial
  Polynomial polynomial() const {
    CGAL_assertion(is_isolated());
    return this->ptr()->polynomial();
  }

  //! Returns the traits class
  Bitstream_descartes_rndl_tree_traits traits() const {
    return this->ptr()->traits();
  }

  //! Number of real roots of the polynomial
  int number_of_real_roots() const {
    CGAL_assertion(is_isolated());
    return this->ptr()->number_of_real_roots();
  }

  //! Refine the <tt>i</tt>th isolating interval
  void refine_interval(int i) const {
    CGAL_assertion(is_isolated());
    this->ptr()->refine_interval(i);
  }

  //! The left bound of the <tt>i</tt>th isolating interval
  Bound left_bound(int i) const  {
    CGAL_assertion(is_isolated());
    return this->ptr()->left_bound(i);
  }

  //! The left bound of the <tt>i</tt>th isolating interval
  void left_bound(int i,
                  Integer& numerator,
                  Integer& denominator) const {
    typedef CGAL::Fraction_traits<Bound> Fraction_traits;
    typename Fraction_traits::Decompose decompose;
    decompose(left_bound(i),numerator,denominator);
  }

  //! The right bound of the <tt>i</tt>th isolating interval
  Bound right_bound(int i) const  {
    CGAL_assertion(is_isolated());
    return this->ptr()->right_bound(i);
  }

  //! The right bound of the <tt>i</tt>th isolating interval
  void right_bound(int i,
                   Integer& numerator,
                   Integer& denominator) const {
    typedef CGAL::Fraction_traits<Bound> Fraction_traits;
    typename Fraction_traits::Decompose decompose;
    decompose(right_bound(i),numerator,denominator);
  }

  //! The length of the <tt>i</tt>th isolating interval
  Bound length(int i) const {
    CGAL_assertion(is_isolated());
    return (this->ptr()->right_bound(i) -
            this->ptr()->left_bound(i));
  }

  bool is_exact_root(int) const { return false; }

  /*!
   * \brief Returns true if the <tt>i</tt>th root is known to be a simple
   * root of the curve.
   */
  bool is_certainly_simple_root(int i) const {
    CGAL_assertion(is_isolated());
    return this->ptr()->is_certainly_simple_root(i);
  }

  /*!
   * \brief Returns true if the <tt>i</tt>th root is known to be a multiple
   * root of the curve.
   */
  bool is_certainly_multiple_root(int i) const {
    CGAL_assertion(is_isolated());
    return this->ptr()->is_certainly_multiple_root(i);
  }


  /*!
   * \brief Returns the multiplicity of the root if know, otherwise -1
   */
  int multiplicity_of_root(int i) const {
    CGAL_assertion(is_isolated());
    return this->ptr()->multiplicity_of_root(i);
  }

  /*!
   * Returns an upper bound for the multiplicity of the ith root
   */
  int get_upper_bound_for_multiplicity(int i) const {
    CGAL_assertion(is_isolated());
    return this->ptr()->get_upper_bound_for_multiplicity(i);
  }

  /*!
   * \brief Returns the isolator of the polynomial f(1/x + q), if known
   */
  Self inverse_transform_isolator() const {
    return this->ptr()->inverse_transform_isolator();
  }


public:

  //! Starts the isolation of the real roots.
  void isolate() {
    CGAL_assertion(!is_isolated());
    this->ptr()->isolate();
  }

  bool is_isolated() const {
    return this->ptr()->is_isolated();
  }

  Bitstream_tree get_tree() const {
    return this->ptr()->get_tree();
  }

  //! returns the degree of the gcd of f and its derivative, if known
  int degree_of_gcd() const {
    return this->ptr()->degree_of_gcd();
  }

  //! returns the square free part of f, if known
  Polynomial square_free_part() const {
    return this->ptr()->square_free_part();
  }

};

} // namespace internal

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
