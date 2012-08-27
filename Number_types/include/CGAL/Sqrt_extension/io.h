// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>


#ifndef CGAL_SQRT_EXTENSION_IO_H
#define CGAL_SQRT_EXTENSION_IO_H

#include <sstream>

#include <CGAL/basic.h>
#include <CGAL/Sqrt_extension/Sqrt_extension_type.h>


namespace CGAL {


template<class NT, class ROOT, class ACDE_TAG, class FP_TAG>
void
input_ascii(std::istream& is , Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>& result){

  typedef Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG> EXT; 

    char c;
    NT a0;
    NT a1;
    ROOT root;

    swallow(is, 'E');
    swallow(is, 'X');
    swallow(is, 'T');
    swallow(is, '[');
    is >> iformat(a0);
    do c = is.get(); while (isspace(c));
    if (c != ',') CGAL_error_msg( "input error: , expected" );

    is >> iformat(a1);
    do c = is.get(); while (isspace(c));
    if (c != ',') CGAL_error_msg( "input error: , expected" );

    is >> iformat(root);
    do c = is.get(); while (isspace(c));
    if (c != ']') CGAL_error_msg( "input error: ] expected" );

    if ( root  < ROOT(0)) CGAL_error_msg("input error: non-negative root expected");

    if ( root == ROOT(0)) 
        result =  EXT(a0);
    else
        result = EXT(a0,a1,root);
}

template<class NT, class ROOT, class ACDE_TAG, class FP_TAG>
void
output_maple(std::ostream& os, const Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>& x){
    CGAL::IO::Mode o_mode=::CGAL::get_mode(os);
    ::CGAL::set_mode(os,CGAL::IO::PRETTY);
    
    if ( x.a0() != NT(0)){
        if ( x.a1() != NT(0)){
            os << x.a0()
               << "+" << CGAL::oformat(x.a1(),CGAL::Parens_as_product_tag())
               << "*sqrt(" << x.root() << ")";
        }else{
            os << x.a0();
        }
    }
    else{
        if (x.a1() != NT(0)){
            os << CGAL::oformat(x.a1(),CGAL::Parens_as_product_tag())
               << "*sqrt(" << x.root() << ")";
        }else{
            os << 0;
        }
    }
    ::CGAL::set_mode(os,o_mode);
    return;
}

template< class NT, class ROOT, class ACDE_TAG, class FP_TAG >
void
output_benchmark( std::ostream& os, const Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>& x ) {
    os << "Sqrt_extension( " << bmformat( x.a0() ) << ", " << bmformat( x.a1() )
       << ", " << bmformat( x.root()) << " )";
}

// Benchmark_rep specialization 
template < class NT, class ROOT, class ACDE_TAG, class FP_TAG >
class Benchmark_rep< CGAL::Sqrt_extension< NT,ROOT, ACDE_TAG, FP_TAG> > {
    const CGAL::Sqrt_extension< NT,ROOT,ACDE_TAG,FP_TAG>& t;
public:
    //! initialize with a const reference to \a t.
    Benchmark_rep( const CGAL::Sqrt_extension< NT,ROOT,ACDE_TAG,FP_TAG>& tt) : t(tt) {}
    //! perform the output, calls \c operator\<\< by default.
    std::ostream& operator()( std::ostream& out) const { 
        output_benchmark( out, t );
        return out;
    }
    
    static std::string get_benchmark_name() {
        std::stringstream ss;
        ss << "Sqrt_extension< " << Benchmark_rep< NT >::get_benchmark_name() 
           << ", " << Benchmark_rep< ROOT>::get_benchmark_name() << " >";
        return ss.str();
    }
};


template <class COEFF, class ROOT, class ACDE_TAG,class FP_TAG>
struct Needs_parens_as_product< Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG> >{
public:
    typedef Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG> NT;
    bool operator()(const NT& t){
        if( t.a0() != NT(0) && t.a1() != NT(0)){
            return true;
        }
        if( t.a1() == NT(0) ){
            Needs_parens_as_product<COEFF> npap;
            return npap(t.a0());
        }
        return false;
    }
};

template <class NT,class ROOT, class ACDE_TAG,class FP_TAG>
std::ostream& operator << (std::ostream& os,
        const Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>& ext){
    switch(CGAL::get_mode(os)) {
    case CGAL::IO::PRETTY:
        output_maple(os,ext); break; 
    default:
        os<<"EXT["<<ext.a0()<<","<<ext.a1()<<","<<ext.root()<<"]"; break;
    }
    return os;
}

/*! \relates CGAL::Sqrt_extension
 *  \brief try to read a CGAL::Sqrt_extension from \c is into \c ext
 *
 *  \c is must be in a mode that supports input of CGAL::Sqrt_extension
 *  (\c LiS::IO::ASCII or \c LiS::IO::BINARY) and the input from
 *  \c is must have the format of output to a stream of the same mode.
 */
template <class NT,class ROOT, class ACDE_TAG, class FP_TAG>
std::istream& operator >> (std::istream& is, Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>& ext) {
    CGAL_precondition(!CGAL::is_pretty(is));
    input_ascii(is,ext);
    return is;
}

} //namespace CGAL

#endif
