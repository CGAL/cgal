#include <list>

int main()
{
    // Reading non necessarily convex polygons from file
    std::list<Polygon> P;
    std::ifstream ifs("input_polygons");
    std::istream_iterator<Polygon> ifs_it(ifs),ifs_end;
    std::copy(ifs_it,ifs_end,back_inserter(P));

    // Creating a scene and inserting the polygons in P
    Scene S;
    copy( P.begin() , P.end() , back_inserter(S) );

    // Computing the Visibility Complex
    Visibility_complex V(S);

    // Output bitangents
    for(
}
