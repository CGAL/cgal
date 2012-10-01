#include "CGAL/Combination_enumerator.h"
#include <iostream>
#include <vector>
using namespace std;
int main()
{
    vector<string> names;
    names.push_back("Sun");    names.push_back("Shannon");
    names.push_back("Hurley"); names.push_back("Sawyer");
    names.push_back("Kate");   names.push_back("Claire");
    names.push_back("John");   names.push_back("Jack");
    CGAL::Combination_enumerator<vector<string>::iterator>
        combi(2, names.begin(), names.end());
    while( ! combi.finished() ) {
        cout << " {" << *combi[0] << ' ' << *combi[1] << '}';
        ++combi;
    }
    cout << endl;
    return EXIT_SUCCESS;
}
