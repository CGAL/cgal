// Copyright (C) 2001
// Mac Murrett
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee,
// provided that the above copyright notice appear in all copies and
// that both that copyright notice and this permission notice appear
// in supporting documentation.  Mac Murrett makes no representations
// about the suitability of this software for any purpose.  It is
// provided "as is" without express or implied warranty.
//
// See http://www.boost.org for most recent version including documentation.

#ifndef BOOST_SINGLETON_MJM012402_HPP
#define BOOST_SINGLETON_MJM012402_HPP

namespace boost {

namespace detail {

namespace thread {

// class singleton has the same goal as all singletons: create one instance of a
//  class on demand, then dish it out as requested.

template<class T>
class singleton: private T
{
private:
    singleton();
    ~singleton();

public:
    static T &instance();
};


template<class T>
inline singleton<T>::singleton()
{   /* no-op */ }

template<class T>
inline singleton<T>::~singleton()
{   /* no-op */ }


template<class T>
/*static*/ T &singleton<T>::instance()
{
    // function-local static to force this to work correctly at static
    // initialization time.
    static singleton<T> s_oT;
    return(s_oT);
}


} // namespace thread

} // namespace detail

} // namespace boost


#endif // BOOST_SINGLETON_MJM012402_HPP
