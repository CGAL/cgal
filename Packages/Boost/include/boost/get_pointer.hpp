// Copyright Peter Dimov and David Abrahams 2002. Permission to copy,
// use, modify, sell and distribute this software is granted provided
// this copyright notice appears in all copies of the source. This
// software is provided "as is" without express or implied warranty,
// and with no claim as to its suitability for any purpose.
#ifndef GET_POINTER_DWA20021219_HPP
# define GET_POINTER_DWA20021219_HPP

# include <memory>

namespace boost { 

// get_pointer(p) extracts a ->* capable pointer from p

template<class T> T * get_pointer(T * p)
{
    return p;
}

// get_pointer(shared_ptr<T> const & p) has been moved to shared_ptr.hpp

template<class T> T * get_pointer(std::auto_ptr<T> const& p)
{
    return p.get();
}


} // namespace boost

#endif // GET_POINTER_DWA20021219_HPP
