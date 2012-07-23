#include "CGAL/Combination_enumerator.h"
#include <iostream>

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

void test(const int K, const int Min, const int Max)
{
    unsigned int n(0);
    cout << "Taking " << K << " distinct integers in the range [" <<
            Min << ", " << Max << "]:";

    CGAL::Combination_enumerator combi(K, Min, Max);
    while( ! combi.finished() )
    {
        cout << " (";
        for(int i = 0; i < K - 1; ++i)
            cout << combi[i] << ' ';
        cout << combi[K-1] << ')';
        ++n;
        ++combi;
    }
    long nelem = Max - Min + 1;
    long num = fac(nelem - K + 1, nelem) / fac(2, K);
    cout << endl << "Enumerated " << n << " combinations. Should be " << num << endl;
    CGAL_assertion(n == num);
}

int main()
{
    test(3,  10, 20);
    test(1, -10, 20);
    test(4,   3,  6);
    return EXIT_SUCCESS;
}
