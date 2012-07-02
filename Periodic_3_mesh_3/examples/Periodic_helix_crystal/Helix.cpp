
#include "stdafx.h"

#include "Helix.h"

Helix::Helix(const Point_3& basepoint, const Vector_3& axisvec, const Vector_3& radiusvec,
             const Vector_3& voxelradius, t_chirality chirality, int ctrotations, int ctintervalsperrotation)
             : pi(acos(-1.0))
{
    m_basepoint = basepoint;

    m_axisvec = axisvec;
    
    m_radius = sqrt(radiusvec.squared_length());
    m_r1 = radiusvec / m_radius;
    m_r2 = cross_product(axisvec / sqrt(axisvec.squared_length()), m_r1);

    m_voxelradius = voxelradius;
    
    m_handedness = (chirality == righthanded) ? +1.0 : -1.0;
    m_ctrotations = ctrotations;
    m_ctintervals = ctintervalsperrotation * ctrotations;

    m_dt = 1.0 / m_ctintervals;

	// populate helix path cache
	for (int k = 0; k <= m_ctintervals; k++)
	{
		m_helixpathcache.push_back(Path(k * m_dt));
	}

    Iso_cuboid_3 bbox = bounding_box(m_helixpathcache.begin(), m_helixpathcache.end());

    m_boundingbox = Iso_cuboid_3(bbox.min() - 2 * voxelradius, bbox.max() + 2 * voxelradius);   // safety margin
}

Point_3 Helix::GetNearestPointOnPath(Point_3 p)
{
    Point_3 pnearest = m_helixpathcache.front();

    for (std::vector<Point_3>::iterator iter = m_helixpathcache.begin(); iter != m_helixpathcache.end(); ++iter)
    {
        if (Vector_3(p, *iter).squared_length() < Vector_3(p, pnearest).squared_length())
        {
            pnearest = *iter;
        }
    }

    return pnearest;
}

Bbox_3 Helix::GetBoundingBox()
{
    return m_boundingbox.bbox();
}


FT Helix::ImplicitFunction(Point_3 p)
{
    if (m_boundingbox.has_on_unbounded_side(p))
    {
      return -1000.0;
    }

    std::vector<double> f;
    double fmax;
    int k, kmax;

    for (k = 0; k <= m_ctintervals; k++)
    {
      f.push_back(ImplicitFunction(p, m_helixpathcache[k]));
    }

    for (kmax = 0, fmax = f[kmax], k = 1; k <= m_ctintervals; k++)
    {
        if (f[k] > f[kmax])
        {
	  kmax = k;
        }
    }

    k = kmax;

    if (k == 0)
        k = 1;
    else if (k == m_ctintervals)
        k = m_ctintervals - 1;

	/*
	// interpolate maximum using Newton's polynomials
    double a = f[k] - f[k-1];
    double b = (f[k+1] - f[k]) - (f[k] - f[k-1]);
    double tmax = (k - 1) * m_dt + (b - 2*a)/(2*b) * m_dt;
	*/

	double tmax = FindMaximum(p, (k-1) * m_dt, (k+1) * m_dt);

    if (tmax < 0.0)
        tmax = 0.0;
    if (tmax > 1.0)
        tmax = 1.0;

    return ImplicitFunction(p, tmax);
}

double Helix::ImplicitFunction(Point_3 p, Point_3 path)
{
    const Vector_3 q = p - path;

	const FT qvx = q.x() / m_voxelradius.x();
	const FT qvy = q.y() / m_voxelradius.y();
	const FT qvz = q.z() / m_voxelradius.z();

    return 1.0 - (qvx * qvx + qvy * qvy + qvz * qvz);
}

double Helix::ImplicitFunction(Point_3 p, double t)
{
    return ImplicitFunction(p, Path(t));
}

Point_3 Helix::Path(double t)
{
    assert(0 <= t && t <= 1);

	double omega = 2 * pi * m_handedness * m_ctrotations;
    double phi = omega * t;

    return m_basepoint + t * m_axisvec + m_radius * (cos(phi) * m_r1 + sin(phi) * m_r2);
}

Vector_3 Helix::Velocity(double t)
{
    assert(0 <= t && t <= 1);

	// time derivative of Path(t)
	double omega = 2 * pi * m_handedness * m_ctrotations;
    double phi = omega * t;

	return m_axisvec + omega * m_radius * (-sin(phi) * m_r1 + cos(phi) * m_r2);
}

double Helix::TimeDerivative(Point_3 p, double t)
{
    assert(0 <= t && t <= 1);

    const Vector_3 q = p - Path(t);
	const Vector_3 dh = Velocity(t);

	const FT rx2 = m_voxelradius.x() * m_voxelradius.x();
	const FT ry2 = m_voxelradius.y() * m_voxelradius.y();
	const FT rz2 = m_voxelradius.z() * m_voxelradius.z();

	return 2.0 * ((q.x() * dh.x() / rx2) + (q.y() * dh.y() / ry2) + (q.z() * dh.z() / rz2));
}

double Helix::FindMaximum(Point_3 p, double tmin, double tmax)
{
	// find zero of TimeDerivative(p, [tmin tmax]) using Regula Falsi

    assert(0 <= tmin && tmin <= 1);
    assert(0 <= tmax && tmax <= 1);

	double dmin = TimeDerivative(p, tmin);
	double dmax = TimeDerivative(p, tmax);

	if (dmin * dmax > 0)	// tmin, tmax do not bracket an extremum
	{
		if (ImplicitFunction(p, tmin) > ImplicitFunction(p, tmax))
			return tmin;
		else
			return tmax;
	}

	const double eps = 5e-15;
	const int maxsteps = 100;
	double t, d;
	int side = 0;

	for (int k = 0; k < maxsteps; ++k)
	{
		// compute the root of the secant through (tmin,dmin) and (tmax,dmax)
		t = (dmin * tmax - dmax * tmin) / (dmin - dmax);

		// check error range
		if (fabs(dmax - dmin) < eps * fabs(tmax + tmin))
			break;

		// compute function value at root of secant
		d = TimeDerivative(p, t);

		if (d * dmax > 0)	// t > actual root
		{
			tmax = t;
			dmax = d;

			if (side == -1)	// same end-point retained twice in a row
				dmin /= 2;	// modify false position to improve convergence

			side = -1;
		}
		else if (dmin * d > 0)	// t < actual root
		{
			tmin = t;
			dmin = d;

			if (side == +1)	// same end-point retained twice in a row
				dmax /= 2;	// modify false position to improve convergence

			side = +1;
		}
		else
			break;
	}

	return t;
}
