

#include "Bops_types.h"
#include "Bops_io.h"
#include <iostream>


#include <CGAL/boolean_operations_2.h>

using namespace CGAL;
using namespace std;

int main()
{
    TestPolygon pgn1, pgn2;
    if (!read_input(pgn1, pgn2)) {
	CGAL_STD::cerr << "Input error!\n";
	return 1;
    }

    list<Object> result;
    Union(pgn1, pgn2, back_inserter(result));
    
    write_result(result);
    return 0;
}
