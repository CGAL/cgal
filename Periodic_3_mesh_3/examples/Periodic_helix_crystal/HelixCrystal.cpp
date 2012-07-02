
#include "stdafx.h"

#include "HelixCrystal.h"

HelixCrystal::HelixCrystal(t_crystalsize cells, FT gratingconstant, FT helixradius, const Vector_3& voxelradius, t_chirality corner, t_chirality helix, bool unitcell)
{  
  m_unitcell=unitcell; 
  m_cells = cells;
  m_gratingconstant = gratingconstant;
  m_radius = helixradius;
  m_voxel = voxelradius;
  m_chiralcorner = corner;
  m_chiralhelix = helix;
  
  FT a = m_gratingconstant;
  FT da = m_radius / sqrt(2.0);
  FT s = (m_chiralcorner == righthanded) ? +1.0 : -1.0;

  Point_3 DP0(a/2 + s * da + a, a/2 - s * da + a, -3*m_voxel.z());
  Point_3 DP1(DP0.x() + a, DP0.y() + a, DP0.z() + cells.z*a +6*m_voxel.z());
  unitcell_bbox = Iso_cuboid_3(DP0, DP1); 			       


	for (int kz = 0; kz < m_cells.z; ++kz)
	{
		for (int ky = 0; ky < m_cells.y; ++ky)
		{
			m_helices.push_back(new Helix(
				Point_3(0, a/2 + s * da, a/2 + s * da) + a * Vector_3(0, ky, kz),
				a * m_cells.x * Vector_3(1,0,0),
				Vector_3(0, s * da, s * da),
				m_voxel,
				m_chiralhelix,
				m_cells.x,
				m_ctintervalsperrotation
				));
		}
	}


	for (int kz = 0; kz < m_cells.z; ++kz)
	{
		for (int kx = 0; kx < m_cells.x; ++kx)
		{
			m_helices.push_back(new Helix(
				Point_3(a/2 - s * da, 0, a/2 - s * da) + a * Vector_3(kx, 0, kz),
				a * m_cells.y * Vector_3(0,1,0),
				Vector_3(-s * da, 0, -s * da),
				m_voxel,
				m_chiralhelix,
				m_cells.y,
				m_ctintervalsperrotation
				));
		}
	}

	for (int kx = 0; kx < m_cells.x; ++kx)
	{
		for (int ky = 0; ky < m_cells.y; ++ky)
		{
			m_helices.push_back(new Helix(
				Point_3(a/2 + s * da, a/2 - s * da, 0) + a * Vector_3(kx, ky, 0),
				a * m_cells.z * Vector_3(0,0,1),
				Vector_3(s * da, -s * da, 0),
				m_voxel,
				m_chiralhelix,
				m_cells.z,
				m_ctintervalsperrotation
				));
		}
	}
}

HelixCrystal::~HelixCrystal()
{
	for (std::vector<Helix*>::iterator iter = m_helices.begin(); iter != m_helices.end(); ++iter)
	{
		delete(*iter);
	}
}

Sphere_3 HelixCrystal::GetBoundingSphere()
{
  Bbox_3 bbox=GetBoundingBox();  
  // find point on helix paths that is nearest to the center of the bounding box of all helix paths
  
  Point_3 bboxcenter = Point_3(0.5 * (bbox.xmin() + bbox.xmax()),
                                    0.5 * (bbox.ymin() + bbox.ymax()),
                                    0.5 * (bbox.zmin() + bbox.zmax()));

    std::vector<Helix*>::iterator iter = m_helices.begin();
    Point_3 center = (*iter)->GetNearestPointOnPath(bboxcenter);

    ++iter;

    while (iter != m_helices.end())
    {
        Point_3 p = (*iter)->GetNearestPointOnPath(bboxcenter);

        if (Vector_3(bboxcenter, p).squared_length() < Vector_3(bboxcenter, center).squared_length())
        {
            center = p;
        }

        ++iter;
    }

    assert(ImplicitFunction(center) > 0);

    Iso_cuboid_3 cuboid(bbox.xmin(), bbox.ymin(), bbox.zmin(), bbox.xmax(), bbox.ymax(), bbox.zmax());

    double squared_radius = 0.0;

    for (int k = 0; k < 8; ++k)
    {
        double dist = Vector_3(center, cuboid.vertex(k)).squared_length();

        if (dist > squared_radius)
        {
            squared_radius = dist;
        }
    }

    return Sphere_3(center, squared_radius);
}


Bbox_3 HelixCrystal::GetBoundingBox()
{
  if (m_unitcell)
    return unitcell_bbox.bbox();
  else {
    std::vector<Helix*>::iterator iter = m_helices.begin();
    Bbox_3 commonbox = (*iter)->GetBoundingBox();

    ++iter;

    while (iter != m_helices.end())
    {
      commonbox = commonbox + (*iter)->GetBoundingBox();
      
      ++iter;
    }      
    return commonbox;
  }
  
}

t_crystalsize HelixCrystal::GetCrystalSize()
{
    return m_cells;
}

double HelixCrystal::GetGratingConstant()
{
    return m_gratingconstant;
}

double HelixCrystal::GetRadius()
{
    return m_radius;
}

Vector_3 HelixCrystal::GetVoxel()
{
    return m_voxel;
}

FT HelixCrystal::ImplicitFunction(Point_3 p)
{
  if (m_unitcell)
    if (unitcell_bbox.has_on_unbounded_side(p) || (unitcell_bbox.has_on_boundary(p)))
      return -1000.0;
  
  FT f;
  FT fmax = -1000;

  for (std::vector<Helix*>::iterator iter = m_helices.begin(); iter != m_helices.end(); ++iter)
    {
      f = (*iter)->ImplicitFunction(p);
      
      if (f > fmax)
	{
	  fmax = f;
	}
    }
  
  return fmax;
}
