
namespace CGAL {

/*!
\ingroup nt_util

\brief The class `Protect_FPU_rounding` allows to reduce the number of rounding mode changes when evaluating sequences of interval arithmetic operations.

\anchor protect_fpu_rouding 

Floating-point arithmetic, as specified by the IEEE-754 standard, allows to use 
so-called directed rounding for the following arithmetic operations: addition, 
subtraction, multiplication, division and square root. The default behavior is 
that the result of such an arithmetic operation is the closest floating-point 
number to the exact real result of the operation (rounding to the nearest). 
The other rounding modes are: round towards plus infinity, round towards minus 
infinity, and round towards zero. 

Interval arithmetic uses such directed rounding modes to offer guaranteed 
enclosures for the evaluation of real functions, such as with \cgal's 
`Interval_nt` class. 

In order to efficiently evaluate sequences of interval arithmetic operations, 
such as a geometric predicate computing for example a determinant, it is 
advised to reduce the number of rounding mode changes, which otherwise are 
performed for each arithmetic operation. \cgal exploits the fact that it is 
possible to compute a sequence of interval arithmetic operations by doing only 
one rounding mode change around the whole function evaluation in order to 
benefit from this optimization. 

The class `Protect_FPU_rounding` allows to easily benefit from this. 
Its constructor saves the current rounding mode in the object, and then sets 
the current rounding mode to the value provided as argument to the constructor. 
The destructor sets the rounding mode back to the saved value. 
This allows to protect a block of code determined by a C++ scope, and have 
the destructor take care of restoring the value automatically. 

The related class `Set_ieee_double_precision` allows to similarly protect 
a block of code from excess precision on some machines (x86 typically with 
the traditional FPU, not the more recent SSE2). Note that 
`Protect_FPU_rounding`, when changing rounding modes, also sets the precision 
to the correct 64 bit precision, hence providing a similar effect to 
`Set_ieee_double_precision`. This notably affects the `Residue` class. 

Note for Visual C++ 64-bit users: due to a compiler bug, the stack unwinding 
process happenning when an exception is thrown does not correctly execute the 
rounding mode restoration when the `Protect_FPU_rounding` object is 
destroyed. Therefore, for this configuration, some explicit code has to be 
added. 

\cgalHeading{Parameters}

The template parameter `Protected` is a Boolean parameter, which defaults 
to `true`. It follows the same parameter of the `Interval_nt` class. 
When it is `false`, the constructor and the destructor of the class do 
nothing (this is meant to be used in a context where you know that the rounding 
mode change has been taken care of at a higher level in the call stack. 

What follows describes the behavior when the parameter has its default value, 
`true`. 

\sa `CGAL::Set_ieee_double_precision` 

*/
template< typename Protected >
class Protect_FPU_rounding {
public:

/// \name Creation 
/// @{

/*!
The current rounding mode is saved in the object, and rounding mode is set to `r` 
which can be any of `CGAL_FE_TONEAREST`, `CGAL_FE_TOWARDZERO`, 
`CGAL_FE_UPWARD` (the default) and `CGAL_FE_DOWNWARD`. 
*/ 
Protect_FPU_rounding(FPU_CW_t r = CGAL_FE_UPWARD); 



/*!
The rounding mode is restored to the saved value. 
*/ 
~Protect_FPU_rounding(); 

/// @}

}; // end class
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup nt_util

\brief The class `Set_ieee_double_precision` provides a mechanism to set 
the correct 53 bits precision for a block of code. It does so by having 
a default constructor that sets a particular mode on the FPU which corrects 
the problem, and have its destructor reset the mode to its previous state. 

\details 
\anchor set_ieee_double_precision 

The IEEE754 standard specifies that the precision of double precision 
floating-point numbers should be 53 bits, with 11 bits for the exponent range. 

Some processors violate this rule by providing excess precision during some 
computations (when values are in registers). This is the case of the x86 
Intel processor and compatible processors (note that the SSE2 more recent 
alternative FPU is fortunately not affected by this issue). The effect of such 
excess precision can be a problem for some computations, since it can produce 
so-called double rounding effects, where actually <I>less</I> precision is 
actually provided! It can also be the root of non-deterministic computations 
depending on compiler optimizations or not (since this affects how long 
variables are kept in registers), for example numerical floating-point 
values get computed with slightly different results. Finally, it affects code 
that carefully makes use of cancellation properties, like `Residue`. 



If the platform is not affected by the excess precision problem, this class becomes an empty class doing nothing. 

Note that nothing can be done for the excess range of the exponent, which 
affects underflow and overflow cases, fortunately less frequent. 

Note also that in the process of setting the correct precision, the rounding 
mode is also set to the nearest. 

Moreover, some compilers provide a flag that performs this setting at the 
time of program startup. For example, GCC provides the option <TT>-mpc64</TT> 
since release 4.3 which does exactly this. Other compilers may have similar 
options. 

Similarly, some third-party libraries may do the same thing as part of their 
startup process, and this is notably the case of LEDA (at least some versions 
of it). \cgal does not enforce this at startup as it would impact 
computations with long double performed by other codes in the same program. 

Note that this property is notably required for proper functionning of the 
`Residue` class that performs modular arithmetic using efficient 
floating-point operations. 

Note concerning Visual C++ 64-bit: due to a compiler bug, the stack unwinding 
process happenning when an exception is thrown does not correctly execute the 
restoring operation when the `Set_ieee_double_precision` object is 
destroyed. Therefore, for this configuration, some explicit code has to be 
added if you care about the state being restored. 

\sa `CGAL::Protect_FPU_rounding<Protected>`
\sa `CGAL::Residue` 

*/

class Set_ieee_double_precision {
public:

/// \name Creation 

/// @{

/*!
Sets the precision of operations on double to 53bits. 
Note that the rounding mode is set to the nearest in the same process. 
*/ 
Set_ieee_double_precision(); 

/*!
The precision and rounding modes are reset to the values they held before the 
constructor was called. 
*/ 
~Set_ieee_double_precision(); 

}; /* end Set_ieee_double_precision */


/*!

Sets the precision of operations on double to 53bits. 
Note that the rounding mode is set to the nearest in the same process. 

The function does the same thing as the default constructor of `Set_ieee_double_precision` 
except that it does not perform the save and restore of the previous state.

\relates Set_ieee_double_precision 
*/ 
void force_ieee_double_precision(); 

/// @}


} /* end namespace CGAL */
