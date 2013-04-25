#include <CGAL/Combination_enumerator.h>
#include <iostream>
#include <vector>
#include <cassert>


using namespace std;

long fac(long from, long to)
{
    long result(1);
    while( from <= to )
    {
        result *= from;
        ++from;
    }
    return result;
}

template< typename T >
void test(const int K, const T & first, const T & beyond)
{
    long n(0);

    CGAL::Combination_enumerator<T> combi(K, first, beyond);
    assert( first  == combi.min_element() );
    assert( beyond == combi.beyond_element() );
    assert( K      == combi.number_of_elements() );
    while( ! combi.finished() )
    {
        ++n;
        ++combi;
    }
    long nelem = static_cast<long>(beyond - first);
    long num = fac(nelem - K + 1, nelem) / fac(2, K);
    cout << endl << "Enumerated " << n << " combinations. Should be " << num;
    assert(n == num);
}

int main()
{
    test(3,  10, 21); // triples in [10,20]
    test(1, -10, 21); // singletons in [-10,20]
    test(4,   3,  7); // the unique set {3,4,5,6}
    char name[] = {'s', 'a', 'm', 'u', 'e', 'l', 0};
    test(2, name+0, name+6);
    vector<int> l;
    for( int i = 0; i < 10; ++i )
        l.push_back(rand());
    test(3, l.begin(), l.end());
    cout << endl;
    return EXIT_SUCCESS;
}
