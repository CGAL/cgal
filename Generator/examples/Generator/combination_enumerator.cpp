#include "CGAL/Combination_enumerator.h"
#include <iostream>
using namespace std;
int main()
{
    unsigned int n(0);
    const int K(3), First(10), Last(20);
    cout << "Taking " << K << " distinct integers in the range [" <<
            First << ", " << Last << "]:";

    CGAL::Combination_enumerator<int> combi(K, First, Last+1);
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
