// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef OPERATOR_ID_DWA2002531_HPP
# define OPERATOR_ID_DWA2002531_HPP

namespace boost { namespace python { namespace detail {

enum operator_id
{ 
    op_add, 
    op_sub, 
    op_mul, 
    op_div, 
    op_mod, 
    op_divmod,
    op_pow, 
    op_lshift, 
    op_rshift, 
    op_and, 
    op_xor, 
    op_or, 
    op_neg, 
    op_pos, 
    op_abs, 
    op_invert, 
    op_int, 
    op_long, 
    op_float, 
    op_str,
    op_cmp,
    op_gt,
    op_ge,
    op_lt,
    op_le,
    op_eq,
    op_ne,
    op_iadd,
    op_isub,
    op_imul,
    op_idiv,
    op_imod,
    op_ilshift,
    op_irshift,
    op_iand,
    op_ixor,
    op_ior,
    op_complex,
    op_nonzero
};

}}} // namespace boost::python::detail

#endif // OPERATOR_ID_DWA2002531_HPP
