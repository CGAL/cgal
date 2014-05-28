
#ifndef REFERENCED_ARGUMENT_H
#define REFERENCED_ARGUMENT_H

template< typename T>
struct Referenced_argument
{
    Referenced_argument() : arg_() {}
    template<typename A>
    Referenced_argument(A a) : arg_(a) {}
    template<typename A, typename B>
    Referenced_argument(A a, B b) : arg_(a, b) {}
    T & arg() { return arg_; }

protected:
    T arg_;
};

#endif // REFERENCED_ARGUMENT_H
