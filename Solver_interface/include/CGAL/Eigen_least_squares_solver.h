// Copyright (c) 2005-2008  Inria Loria (France).
/*
 * author:  Bruno Levy, INRIA, project ALICE
 * website: http://www.loria.fr/~levy/software
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Scientific work that use this software can reference the website and
 * the following publication:
 *
 * @INPROCEEDINGS {levy:NMDGP:05,
 *    AUTHOR = Bruno Levy,
 *    TITLE  = Numerical Methods for Digital Geometry Processing,
 *    BOOKTITLE =Israel Korea Bi-National Conference,
 *    YEAR=November 2005,
 *    URL=http://www.loria.fr/~levy/php/article.php?pub=../publications/papers/2005/Numerics
 * }
 *
 *  Laurent Saboret 2005-2006: Changes for CGAL:
 *      - Added OpenNL namespace
 *      - DefaultLinearSolverTraits is now a model of the SparseLinearAlgebraTraits_d concept
 *      - Added SymmetricLinearSolverTraits
 *      - copied Jacobi preconditioner from Graphite 1.9 code
 */


#ifndef __EIGEN_LEAST_SQUARES_SOLVER__
#define __EIGEN_LEAST_SQUARES_SOLVER__

#include <CGAL/Eigen_solver_traits.h>

#include <vector>
#include <iostream>
#include <cstdlib>


namespace CGAL {


/*
 * Solves a linear system or minimizes a quadratic form.
 */
class EigenLeastSquaresSolver
{
protected:
  
    enum State {
        INITIAL, IN_SYSTEM, IN_ROW, CONSTRUCTED, SOLVED
    } ;

public:
  // typedef Eigen_solver_traits<Eigen::ConjugateGradient<Eigen_sparse_symmetric_matrix<double>::EigenType,
  // 						       Eigen::Lower|Eigen::Upper> > Solver;
  typedef Eigen_solver_traits<Eigen::SimplicialLDLT<Eigen_sparse_symmetric_matrix<double>::EigenType> > Solver;
    typedef typename Solver::Matrix Matrix ;
    typedef typename Solver::Vector Vector ;
    typedef typename Solver::NT CoeffType ;

    class Variable {
    public:
        Variable() : x_(0), index_(-1), locked_(false) { }
        double value() const { return x_; }
        void set_value(double x_in) { x_ = x_in ; }
        void lock()   { locked_ = true ; }
        void unlock() { locked_ = false ; }
        bool is_locked() const { return locked_ ; }
        unsigned int index() const {
            CGAL_assertion(index_ != -1) ;
            return (unsigned int)(index_) ;
        }
        void set_index(unsigned int index_in) {
            index_ = index_in ;
        }
    private:
        CoeffType x_ ;
        int index_ ;
        bool locked_ ;
    }  ;


    EigenLeastSquaresSolver(unsigned int nb_variables) {
        state_ = INITIAL ;
        nb_variables_ = nb_variables ;
        variable_ = new Variable[nb_variables] ;
        A_ = NULL ;
        x_ = NULL ;
        b_ = NULL ;
    }

    ~EigenLeastSquaresSolver() {
        delete[] variable_ ;
        delete A_ ;
        delete x_ ;
        delete b_ ;
    }

    // __________________ Parameters ________________________

  void set_least_squares(bool x) { }

    // __________________ Access ____________________________

    int nb_variables() const { return nb_variables_ ; }

    Variable& variable(unsigned int idx) {
        CGAL_assertion(idx < nb_variables_) ;
        return variable_[idx] ;
    }

    const Variable& variable(unsigned int idx) const {
        CGAL_assertion(idx < nb_variables_) ;
        return variable_[idx] ;
    }

    // _________________ Construction _______________________

    void begin_system() {
        current_row_ = 0 ;
        transition(INITIAL, IN_SYSTEM) ;
        // Enumerate free variables.
        unsigned int index = 0 ;
        for(int ii=0; ii < nb_variables() ; ii++) {
            Variable& v = variable(ii) ;
            if(!v.is_locked()) {
                v.set_index(index) ;
                index++ ;
            }
        }
        unsigned int n = index ;
        A_ = new Matrix(static_cast<int>(n),static_cast<int>(n)) ;
        x_ = new Vector(n) ;
        b_ = new Vector(n) ;
        for(unsigned int i=0; i<n; i++) {
	  x_->set(i, 0.);
	  b_->set(i, 0.);
        }
        variables_to_vector() ;
    }

    void begin_row() {
        transition(IN_SYSTEM, IN_ROW) ;
        af_.clear() ;
        if_.clear() ;
        al_.clear() ;
        xl_.clear() ;
        bk_ = 0 ;
    }

    void set_right_hand_side(double b) {
        check_state(IN_ROW) ;
        bk_ = b ;
    }

    void add_coefficient(unsigned int iv, double a) {
        check_state(IN_ROW) ;
        Variable& v = variable(iv) ;
        if(v.is_locked()) {
            al_.push_back(a) ;
            xl_.push_back(v.value()) ;
        } else {
            af_.push_back(a) ;
            if_.push_back(v.index()) ;
        }
    }

    void normalize_row(CoeffType weight = 1) {
        check_state(IN_ROW) ;
        CoeffType norm = 0.0 ;
        unsigned int nf = af_.size() ;
        for(unsigned int i=0; i<nf; i++) {
            norm += af_[i] * af_[i] ;
        }
        unsigned int nl = al_.size() ;
        for(unsigned int i=0; i<nl; i++) {
            norm += al_[i] * al_[i] ;
        }
        norm = sqrt(norm) ;
        CGAL_assertion( fabs(norm)>1e-40 );
        scale_row(weight / norm) ;
    }

    void scale_row(CoeffType s) {
        check_state(IN_ROW) ;
        unsigned int nf = af_.size() ;
         for(unsigned int i=0; i<nf; i++) {
             af_[i] *= s ;
         }
         unsigned int nl = al_.size() ;
         for(unsigned int i=0; i<nl; i++) {
             al_[i] *= s ;
         }
         bk_ *= s ;
    }

    void end_row() {
      unsigned int nf = af_.size() ;
      unsigned int nl = al_.size() ;
      for(unsigned int i=0; i<nf; i++) {
	for(unsigned int j=0; j<nf; j++) {
	  A_->add_coef (if_[i], if_[j], af_[i] * af_[j]);
	}
      }
      CoeffType S = - bk_ ;
      for(unsigned int j=0; j<nl; j++) {
	S += al_[j] * xl_[j] ;
      }
      for(unsigned int i=0; i<nf; i++) {
	b_->set(if_[i], b_->get(if_[i]) - af_[i] * S);
      }
        current_row_++ ;
        transition(IN_ROW, IN_SYSTEM) ;
    }

    void end_system() {
      A_->assemble_matrix ();
        transition(IN_SYSTEM, CONSTRUCTED) ;
    }

    // ----------------------------- Solver -------------------------------

    // Solves a linear system or minimizes a quadratic form.
    // Return true on success.
    // (modified for SparseLinearAlgebraTraits_d concept)
    bool solve()
    {
        check_state(CONSTRUCTED) ;

        // Solve the sparse linear system "A*X = B". On success, the solution is (1/D) * X.	
        Solver solver_;
        CoeffType D;

        bool success = solver_.linear_solver(*A_, *b_, *x_, D) ;
        CGAL_assertion(D == 1.0);   // WARNING: this library does not support homogeneous coordinates!

        vector_to_variables() ;

        transition(CONSTRUCTED, SOLVED) ;

        delete A_ ; A_ = NULL ;
        delete b_ ; b_ = NULL ;
        delete x_ ; x_ = NULL ;

        return success;
    }

protected:

    // ----------- Converting between user representation and the internal representation -----

    void vector_to_variables() {
        for(int ii=0; ii < nb_variables(); ii++) {
            Variable& v = variable(ii) ;
            if(!v.is_locked()) {
	      v.set_value( x_->get(v.index()) ) ;
            }
        }
    }

    void variables_to_vector() {
        for(int ii=0; ii < nb_variables(); ii++) {
            Variable& v = variable(ii) ;
            if(!v.is_locked()) {
	      x_->set(v.index(), v.value());
            }
        }
    }

    // ----------- Finite state automaton (checks that calling sequence is respected) ---------

    std::string state_to_string(State s) {
            switch(s) {
            case INITIAL:
                return "initial" ;
            case IN_SYSTEM:
                return "in system" ;
            case IN_ROW:
                return "in row" ;
            case CONSTRUCTED:
                return "constructed" ;
            case SOLVED:
                return "solved" ;
            }
            // Should not go there.
            CGAL_error();
            return "undefined" ;
    }

    void check_state(State s) {
            CGAL_USE(s);
            CGAL_assertion(state_ == s) ;
    }

    void transition(State from, State to) {
            check_state(from) ;
            state_ = to ;
    }

private:

    // --------------- user representation --------------
    unsigned int nb_variables_ ;
    Variable* variable_ ;

    // --------------- construction -----------------------
    State state_ ;
    unsigned int current_row_ ;
    std::vector<CoeffType> af_ ;
    std::vector<unsigned int> if_ ;
    std::vector<CoeffType> al_ ;
    std::vector<CoeffType> xl_ ;
    double bk_ ;

    // --------------- internal representation ---------
    Matrix* A_ ;
    Vector* x_ ;
    Vector* b_ ;

} ;


} // namespace OpenNL

#endif
