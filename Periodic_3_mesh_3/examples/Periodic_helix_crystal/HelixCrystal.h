
#include "Helix.h"


class HelixCrystal {

public:

  HelixCrystal(t_crystalsize cells, FT gratingconstant, FT helixradius, const Vector_3& voxelradius, t_chirality corner, t_chirality helix, bool unitcell=false);
	~HelixCrystal();
	
    Sphere_3 GetBoundingSphere();
    Bbox_3 GetBoundingBox();

    t_crystalsize GetCrystalSize();
    double GetGratingConstant();
    double GetRadius();
    Vector_3 GetVoxel();

	FT ImplicitFunction(Point_3 p);

protected:
	Iso_cuboid_3 unitcell_bbox;
	bool m_unitcell;
	t_crystalsize m_cells;
	FT m_gratingconstant;
	FT m_radius;
	Vector_3 m_voxel;
	t_chirality m_chiralcorner;
	t_chirality m_chiralhelix;
	
	std::vector<Helix*> m_helices;

	static const int m_ctintervalsperrotation = 100;
};
