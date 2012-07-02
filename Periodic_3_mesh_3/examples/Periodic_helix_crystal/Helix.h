
class Helix {

public:
    Helix(const Point_3& basepoint, const Vector_3& axisvec, const Vector_3& radiusvec,
             const Vector_3& voxelradius, t_chirality chirality, int ctrotations, int ctintervalsperrotation);

    Point_3 GetNearestPointOnPath(Point_3 p);
    Bbox_3 GetBoundingBox();
    FT ImplicitFunction(Point_3 p);

protected:
	double ImplicitFunction(Point_3 p, Point_3 path);
	double ImplicitFunction(Point_3 p, double t);
    Point_3 Path(double t);
	Vector_3 Velocity(double t);
	double TimeDerivative(Point_3 p, double t);
	double FindMaximum(Point_3 p, double tmin, double tmax);

    Point_3 m_basepoint;
    Vector_3 m_axisvec;
    Vector_3 m_r1, m_r2;
    Vector_3 m_voxelradius;
    double m_radius;
    double m_handedness;
    int m_ctrotations;
    int m_ctintervals;

    double m_dt;

	std::vector<Point_3> m_helixpathcache;
    Iso_cuboid_3 m_boundingbox;
    
    const double pi;
};
