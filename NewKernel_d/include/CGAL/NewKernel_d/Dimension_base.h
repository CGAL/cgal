// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_KD_DIMENSION_BASE_h
#define CGAL_KD_DIMENSION_BASE_h
#include <CGAL/Dimension.h>
#include <CGAL/assertions.h>
#include <CGAL/NewKernel_d/utils.h>
#include <CGAL/use.h>
namespace CGAL {
struct Store_dimension_base {
        //TODO: add some assertions
        Store_dimension_base(int dim=UNKNOWN_DIMENSION):dim_(dim){}
        int dimension()const{return dim_;}
        void set_dimension(int dim){dim_=dim;}
        private:
        int dim_;
};
template<class=Dynamic_dimension_tag>
struct Dimension_base {
        Dimension_base(int = UNKNOWN_DIMENSION){}
        int dimension() const { return UNKNOWN_DIMENSION; }
        void set_dimension(int) {}
};
template<int dim_>
struct Dimension_base<Dimension_tag<dim_> > {
        Dimension_base(){}
        Dimension_base(int CGAL_assertion_code(dim)){CGAL_assertion(dim_==dim);}
        int dimension()const{return dim_;}
        void set_dimension(int dim){
          CGAL_assertion(dim_==dim);
          CGAL_USE(dim);
        }
};
}
#endif

