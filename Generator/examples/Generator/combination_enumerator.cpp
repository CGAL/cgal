#include <CGAL/Combination_enumerator.h>
#include <iostream>

int main()
{
    unsigned int n(0);
    const int k(3), first(10), last(15);
    std::cout << "Taking " << k << " distinct integers in the range ["
              <<  first << ", " << last << "]:";

    CGAL::Combination_enumerator<int> combi(k, first, last + 1);
    while( ! combi.finished() ) {
        std::cout << " {";
        for(int i = 0; i < k; ++i) {
            std::cout << combi[i];
            if( i < k - 1 )
                std::cout << ' ';
        }
        std::cout << '}';
        ++n;
        ++combi;
    }
    std::cout << std::endl << "Enumerated " << n << " combinations." << std::endl;
    return EXIT_SUCCESS;
}
