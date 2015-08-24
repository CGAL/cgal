#ifndef __EIGEN_LEAST_SQUARES_SOLVER__
#define __EIGEN_LEAST_SQUARES_SOLVER__

#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/tuple.h>

#include <vector>
#include <iostream>
#include <cstdlib>


namespace CGAL {

  namespace internal
  {
    class Eigen_least_squares_solver
    {

    public:
      typedef Eigen_solver_traits<Eigen::SimplicialLDLT<Eigen_sparse_symmetric_matrix<double>::EigenType> > Solver;
      typedef typename Solver::Matrix Matrix ;
      typedef typename Solver::Vector Vector ;
      typedef typename Solver::NT Scalar ;
      typedef typename CGAL::cpp11::tuple<Scalar, std::size_t, bool> Indexed_scalar;


    private:

      std::vector<Indexed_scalar> _indexed_scalars;

      std::vector<Scalar>       _af;
      std::vector<unsigned int> _if;
      std::vector<Scalar>       _al;
      std::vector<Scalar>       _xl;

      Matrix* _A ;
      Vector* _x ;
      Vector* _b ;

    public:

      Eigen_least_squares_solver (unsigned int n_bvariables)
      {
        _indexed_scalars.resize (n_bvariables);
	for (std::size_t i = 0; i < _indexed_scalars.size (); ++ i)
	  {
	    _indexed_scalars[i].get<0>() = 0.;
	    _indexed_scalars[i].get<1>() = -1;
	    _indexed_scalars[i].get<2>() = false;
	  }
        _A = NULL ;
        _x = NULL ;
        _b = NULL ;
      }

      ~Eigen_least_squares_solver()
      {
        if (_A != NULL) delete _A ;
        if (_x != NULL) delete _x ;
        if (_b != NULL) delete _b ;
      }

      const Scalar& value (std::size_t index) const
      {
	return _indexed_scalars[index].get<0>();
      }
      
      void set_value (std::size_t index, Scalar value)
      {
	_indexed_scalars[index].get<0>() = value;
      }
      
      void lock (std::size_t index)
      {
	_indexed_scalars[index].get<2>() = true;
      }
      
      void store_scalars_in_x ()
      {
        for (std::size_t i = 0; i < _indexed_scalars.size (); ++ i)
	  if (!(_indexed_scalars[i].get<2>()))
	    _x->set (_indexed_scalars[i].get<1>(), _indexed_scalars[i].get<0>());
      }
      
      void get_result ()
      {
        for (std::size_t i = 0; i < _indexed_scalars.size (); ++ i)
	  if (!(_indexed_scalars[i].get<2>()))
	    _indexed_scalars[i].get<0>() = _x->get(_indexed_scalars[i].get<1>());
      }

      void begin_system ()
      {
	std::size_t index = 0;
        for (std::size_t i = 0; i < _indexed_scalars.size (); ++ i)
	  if(!(_indexed_scalars[i].get<2>()))
	    _indexed_scalars[i].get<1>() = index ++;
	
        _A = new Matrix (index, index);
        _x = new Vector (index);
        _b = new Vector (index);
	
        for (std::size_t i = 0; i < index; ++ i)
	  {
	    _x->set(i, 0.);
	    _b->set(i, 0.);
	  }
	
        store_scalars_in_x() ;
      }

      void begin_row ()
      {
        _af.clear() ;
        _if.clear() ;
        _al.clear() ;
        _xl.clear() ;
      }

      void add_coefficient (unsigned int iv, double a)
      {
        if(_indexed_scalars[iv].get<2>())
	  {
	    _al.push_back(a) ;
	    _xl.push_back(_indexed_scalars[iv].get<0>());
	  }
	else
	  {
	    _af.push_back(a) ;
	    _if.push_back(_indexed_scalars[iv].get<1>());
	  }
      }

      void end_row ()
      {
	std::size_t nf = _af.size() ;
	std::size_t nl = _al.size() ;
	
	for (std::size_t i = 0; i < nf; ++ i)
	  for (std::size_t j = 0; j < nf; ++ j)
	    _A->add_coef (_if[i], _if[j], _af[i] * _af[j]);

	Scalar S = 0.;
	
	for (std::size_t j = 0; j < nl; ++ j)
	  S += _al[j] * _xl[j] ;

	for (std::size_t i = 0; i < nf; ++ i)
	  _b->set(_if[i], _b->get(_if[i]) - _af[i] * S);
      }

      void end_system ()
      {
	_A->assemble_matrix ();
      }

      bool solve ()
      {
	Scalar D;
        Solver solver;
        bool success = solver.linear_solver(*_A, *_b, *_x, D) ;

        get_result() ;

        delete _A ; _A = NULL ;
        delete _b ; _b = NULL ;
        delete _x ; _x = NULL ;

        return success;
      }



    } ;


  } // namespace internal

} // namespace CGAL

#endif
