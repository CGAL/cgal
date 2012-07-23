#include "CGAL/Combination_enumerator.h"
#include <iostream>
using namespace std;
int main()
{
    unsigned int n(0);
    const int K(3), Min(10), Max(20);
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
    cout << endl << "Enumerated " << n << " combinations." << endl;
    return EXIT_SUCCESS;
}
