int main()
{
    // Reading non necessarily convex polygons from file
    list<Polygon> D;

    // Creating a scene and inserting the polygons in D
    Scene S;
    copy( D.begin() , D.end() , back_inserter(S) );

    // Computing the Visibility Complex
    Visibility_complex V(S);
}
