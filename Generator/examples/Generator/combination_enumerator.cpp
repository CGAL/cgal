#include "CGAL/Combination_enumerator.h"
#include <iostream>
using namespace std;
int main()
{
    unsigned int n(0);
    const int k(3), first(10), last(15);
    cout << "Taking " << k << " distinct integers in the range [" <<
            first << ", " << last << "]:";

    CGAL::Combination_enumerator<int> combi(k, first, last + 1);
    while( ! combi.finished() ) {
        cout << " {";
        for(int i = 0; i < k; ++i) {
            cout << combi[i];
            if( i < k - 1 )
                cout << ' ';
        }
        cout << '}';
        ++n;
        ++combi;
    }
    cout << endl << "Enumerated " << n << " combinations." << endl;
    return EXIT_SUCCESS;
}
