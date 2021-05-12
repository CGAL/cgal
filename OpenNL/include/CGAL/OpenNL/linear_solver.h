// Copyright (c) 2005-2008  Inria Loria (France).
/*
 * author:  Bruno Levy, INRIA, project ALICE
 * website: http://www.loria.fr/~levy/software
 *
 * This file is part of CGAL (www.cgal.org)
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
 *
 * $URL$
 * $Id$
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */


#ifndef __OPENNL_LINEAR_SOLVER__
#define __OPENNL_LINEAR_SOLVER__

#include <CGAL/OpenNL/conjugate_gradient.h>
#include <CGAL/OpenNL/bicgstab.h>
#include <CGAL/OpenNL/preconditioner.h>
#include <CGAL/OpenNL/sparse_matrix.h>
#include <CGAL/OpenNL/full_vector.h>

#include <vector>
#include <iostream>
#include <cstdlib>

#include <CGAL/use.h>

namespace OpenNL {



// Class DefaultLinearSolverTraits
// is a traits class for solving general sparse linear systems.
// It uses BICGSTAB solver with Jacobi preconditioner.
//
// Concept: Model of the SparseLinearAlgebraTraits_d concept.

template
<
    class COEFFTYPE,                        // type of matrix and vector coefficients
    class MATRIX = SparseMatrix<COEFFTYPE>, // model of SparseLinearSolverTraits_d::Matrix
    class VECTOR = FullVector<COEFFTYPE>    // model of SparseLinearSolverTraits_d::Vector
>
class DefaultLinearSolverTraits
{
// Public types
public:
    typedef COEFFTYPE                       CoeffType ;
    typedef COEFFTYPE                       NT;
    typedef MATRIX                          Matrix ;
    typedef VECTOR                          Vector ;

// Private types
private:
    typedef Jacobi_Preconditioner<NT>       Preconditioner ;
    typedef Solver_preconditioned_BICGSTAB<Matrix, Preconditioner, Vector>
                                            Preconditioned_solver ;
    typedef Solver_BICGSTAB<Matrix, Vector> Solver ;

// Public operations
public:
    // Default contructor, copy constructor, operator=() and destructor are fine

    // Solve the sparse linear system "A*X = B"
    // Return true on success. The solution is then (1/D) * X.
    //
    // Preconditions:
    // - A.row_dimension()    == B.dimension()
    // - A.column_dimension() == X.dimension()
    bool linear_solver (const Matrix& A, const Vector& B, Vector& X, NT& D)
    {
        D = 1;              // OpenNL does not support homogeneous coordinates

        // Solve using BICGSTAB solver with preconditioner
        Preconditioned_solver preconditioned_solver ;
        NT omega = 1.5;
        Preconditioner C(A, omega);
        X = B;
        if (preconditioned_solver.solve(A, C, B, X))
            return true;

        // On error, solve using BICGSTAB solver without preconditioner
#ifdef DEBUG_TRACE
        std::cerr << "  Failure of BICGSTAB solver with Jacobi preconditioner. "
                  << "Trying BICGSTAB." << std::endl;
#endif
        Solver solver ;
        X = B;
        return solver.solve(A, B, X) ;
    }
} ;

// Class SymmetricLinearSolverTraits
// is a traits class for solving symmetric positive definite sparse linear systems.
// It uses Conjugate Gradient solver with Jacobi preconditioner.
//
// Concept: Model of the SparseLinearAlgebraTraits_d concept.

template
<
    class COEFFTYPE,                        // type of matrix and vector coefficients
    class MATRIX = SparseMatrix<COEFFTYPE>, // model of SparseLinearSolverTraits_d::Matrix
    class VECTOR = FullVector<COEFFTYPE>    // model of SparseLinearSolverTraits_d::Vector
>
class SymmetricLinearSolverTraits
{
// Public types
public:
    typedef COEFFTYPE                       CoeffType ;
    typedef COEFFTYPE                       NT;
    typedef MATRIX                          Matrix ;
    typedef VECTOR                          Vector ;

// Private types
private:
    typedef Jacobi_Preconditioner<NT>       Preconditioner ;
    typedef Solver_preconditioned_CG<Matrix, Preconditioner, Vector>
                                            Preconditioned_solver ;
    typedef Solver_CG<Matrix, Vector>       Solver ;

// Public operations
public:
    // Default contructor, copy constructor, operator=() and destructor are fine

    // Solve the sparse linear system "A*X = B"
    // Return true on success. The solution is then (1/D) * X.
    //
    // Preconditions:
    // - A.row_dimension()    == B.dimension()
    // - A.column_dimension() == X.dimension()
    bool linear_solver (const Matrix& A, const Vector& B, Vector& X, NT& D)
    {
        D = 1;              // OpenNL does not support homogeneous coordinates

        // Solve using Conjugate Gradient solver with preconditioner
        Preconditioned_solver preconditioned_solver ;
        NT omega = 1.5;
        Preconditioner C(A, omega);
        X = B;
        if (preconditioned_solver.solve(A, C, B, X))
            return true;

        // On error, solve using Conjugate Gradient solver without preconditioner
#ifdef DEBUG_TRACE
        std::cerr << "  Failure of Conjugate Gradient solver with Jacobi preconditioner. "
                  << "Trying Conjugate Gradient." << std::endl;
#endif
        Solver solver ;
        X = B;
        return solver.solve(A, B, X) ;
    }
};


/*
 * Solves a linear system or minimizes a quadratic form.
 *
 * Requirements for its traits class: must be a model of SparseLinearSolverTraits_d concept
 */
template <class TRAITS>
class LinearSolver
{
protected:
    enum State {
        INITIAL, IN_SYSTEM, IN_ROW, CONSTRUCTED, SOLVED
    } ;

public:
    typedef TRAITS Traits ;
    typedef typename Traits::Matrix Matrix ;
    typedef typename Traits::Vector Vector ;
    typedef typename Traits::NT CoeffType ;

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


    LinearSolver(unsigned int nb_variables) {
        state_ = INITIAL ;
        least_squares_ = false ;
        nb_variables_ = nb_variables ;
        variable_ = new Variable[nb_variables] ;
        A_ = nullptr ;
        x_ = nullptr ;
        b_ = nullptr ;
    }

    ~LinearSolver() {
        delete[] variable_ ;
        delete A_ ;
        delete x_ ;
        delete b_ ;
    }

    // __________________ Parameters ________________________

    void set_least_squares(bool x) { least_squares_ = x ; }

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
        A_ = new Matrix(static_cast<int>(n)) ;
        x_ = new Vector(n) ;
        b_ = new Vector(n) ;
        for(unsigned int i=0; i<n; i++) {
            (*x_)[i] = 0 ;
            (*b_)[i] = 0 ;
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
        if(least_squares_) {
            std::size_t nf = af_.size() ;
            std::size_t nl = al_.size() ;
            for(std::size_t i=0; i<nf; i++) {
                for(std::size_t j=0; j<nf; j++) {
                    A_->add_coef(if_[i], if_[j], af_[i] * af_[j]) ;
                }
            }
            CoeffType S = - bk_ ;
            for(std::size_t j=0; j<nl; j++) {
                S += al_[j] * xl_[j] ;
            }
            for(std::size_t i=0; i<nf; i++) {
                (*b_)[if_[i]] -= af_[i] * S ;
            }
        } else {
            std::size_t nf = af_.size() ;
            std::size_t nl = al_.size() ;
            for(std::size_t i=0; i<nf; i++) {
                A_->add_coef(current_row_, if_[i], af_[i]) ;
            }
            (*b_)[current_row_] = bk_ ;
            for(std::size_t i=0; i<nl; i++) {
                (*b_)[current_row_] -= al_[i] * xl_[i] ;
            }
        }
        current_row_++ ;
        transition(IN_ROW, IN_SYSTEM) ;
    }

    void end_system() {
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
        Traits solver_traits;
        CoeffType D;
        bool success = solver_traits.linear_solver(*A_, *b_, *x_, D) ;
        CGAL_assertion(D == 1.0);   // WARNING: this library does not support homogeneous coordinates!

        vector_to_variables() ;

        transition(CONSTRUCTED, SOLVED) ;

        delete A_ ; A_ = nullptr ;
        delete b_ ; b_ = nullptr ;
        delete x_ ; x_ = nullptr ;

        return success;
    }

protected:

    // ----------- Converting between user representation and the internal representation -----

    void vector_to_variables() {
        for(int ii=0; ii < nb_variables(); ii++) {
            Variable& v = variable(ii) ;
            if(!v.is_locked()) {
                v.set_value( (*x_)[v.index()] ) ;
            }
        }
    }

    void variables_to_vector() {
        for(int ii=0; ii < nb_variables(); ii++) {
            Variable& v = variable(ii) ;
            if(!v.is_locked()) {
                (*x_)[v.index()] = v.value() ;
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

    // --------------- parameters --------------------------
    bool least_squares_ ;

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
