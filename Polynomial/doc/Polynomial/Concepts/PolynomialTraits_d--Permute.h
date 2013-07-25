
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableFunctor` permutes the variables of the given polynomial 
with respect to a permutation \f$ \sigma\f$, that is, each monomial 
\f$ \prod x_i^{e_i}\f$ will be mapped to the monomial \f$ \prod x_{\sigma(i)}^{e_i}\f$. 
The permutation \f$ \sigma\f$ is given by the iterator range of 
length `PolynomialTraits_d::d`, which is supposed to contain 
the second row of the permutation. 

For instance, let \f$ p\f$ be a polynomial in 4 variables and it is intended to 
change the order of the variables such that 
\f$ x_0 \mapsto x_2\f$, \f$ x_1 \mapsto x_0\f$, \f$ x_2 \mapsto x_1\f$ and \f$ x_3 \mapsto x_3\f$. 
In this case the iterator range should contain the sequence \f$ [2,0,1,3]\f$. 

\cgalRefines `AdaptableFunctor` 
\cgalRefines `CopyConstructible` 
\cgalRefines `DefaultConstructible` 

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::Permute {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef PolynomialTraits_d::Polynomial_d result_type; 

/// @} 

/// \name Operations 
/// @{

/*!

Returns \f$ p\f$ with interchanged variables as defined by the iterator range. 
\pre (end-begin == `PolynomialTraits_d::d`) 
\pre `std::iterator_traits< InputIterator >::%value_type` is convertible to int. 
\pre The iterator range contains each value in \f$ \{0,\dots,d-1\}\f$ exactly once. 

*/ 
template<class Input_iterator> 
result_type operator()(PolynomialTraits_d::Polynomial_d p, 
Input_iterator begin, Input_iterator end); 

/// @}

}; /* end PolynomialTraits_d::Permute */

