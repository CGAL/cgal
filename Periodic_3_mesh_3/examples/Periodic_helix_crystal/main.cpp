
#include "stdafx.h"

#include "HelixCrystal.h"
#include "tclap/CmdLine.h" // in /opt/local/include/

typedef struct t_meshcriteria {
	t_meshcriteria() { angle = 30.0; radius = 0.03; distance = 0.03; }	// default values
	std::string str()
	{
		std::stringstream s;
		s << angle << ',' << radius << ',' << distance;
		return s.str();
	}
	double angle;
	double radius;
	double distance;
} t_meshcriteria;

typedef enum t_filetype {
	filetype_unknown,
	filetype_off,
	filetype_pov,
	filetype_stl,
	filetype_stl_ascii,
	filetype_stl_binary
};

typedef std::map<std::string, t_filetype> t_fextensionmap;

class BoxToCubeTranslator
{
public:
  BoxToCubeTranslator()
  {
  }
  
  void set_box(Bbox_3 box)
  {
    m_box = box;
    
    make_cube();
  }
  
  Bbox_3 cube() const
  {
    return m_cube;
  }
  
  // p - from a cube
  Point fromCubeToBox(Point p) const
  {
    double px = p.x() - m_cube.xmin();
    double py = p.y() - m_cube.ymin();
    double pz = p.z() - m_cube.zmin();
    
    return Point(m_box.xmin() + (px * wx)/minw, m_box.ymin() + (py * wy)/minw, m_box.zmin() + (pz * wz)/minw);
  }

private:
  void make_cube()
  {
    wx = m_box.xmax() - m_box.xmin();
    wy = m_box.ymax() - m_box.ymin();
    wz = m_box.zmax() - m_box.zmin();
    
    minw = wx;
    if(minw > wy) {
      minw = wy;
    }
    if(minw > wz) {
      minw = wz;
    }
    
    m_cube = Bbox_3(m_box.xmin(), m_box.ymin(), m_box.zmin(), m_box.xmin() + minw, m_box.ymin() + minw, m_box.zmin() + minw);
  }
  
  double wx, wy, wz, minw;
  
  Bbox_3 m_box, m_cube;
};

HelixCrystal* hc;

// extra data
Tr h_tr;

BoxToCubeTranslator boxToCube;
//

FT implicit_function(const Point& p)
{
  // p comes from the cube
  // q = stretched p
  h_tr.canonicalize_point(p);
  
  Point q = boxToCube.fromCubeToBox(p);
  
  return -1.0 * hc->ImplicitFunction(q);
}

struct Plane_equation {
	template <class Facet>
	typename Facet::Plane_3 operator() (Facet& f)
	{
		typename Facet::Halfedge_handle h = f.halfedge();
		typedef typename Facet::Plane_3 Plane;
		return Plane(h->vertex()->point(),
						h->next()->vertex()->point(),
						h->next()->next()->vertex()->point());
	}
};

std::string getFileExtension(const std::string& s)
{
	size_t ixdot = s.find_last_of('.');
  
	if (ixdot == std::string::npos)
	{
		return std::string();
	}
  
	return s.substr(ixdot+1);
}

template<class T, int dim>
class VectorConstraint
	: public TCLAP::Constraint<std::string>
{
public:
  VectorConstraint(const std::string& tname)
	{
		desc << tname << "[" << dim << "]";
	}

	std::string description() const
	{
		return desc.str();
	}

	std::string shortID() const
	{
		return desc.str();
	}

	bool check(const std::string& value) const
	{
		std::stringstream str(value);
		T val;
		char ch;
		
		for (int k = 0; k < dim; k++)
		{
			if (k > 0)
			{
				str >> ch;

				if (str.fail() || ch != ',')
					return false;
			}

			str >> val;

			if (str.fail() || val <= 0)
				return false;
		}

		return true;
	}

protected:
	std::stringstream desc;
};

class FileExtensionConstraint
	: public TCLAP::Constraint<std::string>
{
public:
	FileExtensionConstraint(t_fextensionmap fmap)
	{
		m_map = fmap;

		for (t_fextensionmap::iterator iter = fmap.begin(); iter != fmap.end(); ++iter)
		{
			if (iter != fmap.begin())
			{
				m_desc << '|';
			}

			m_desc << "*." << iter->first;
		}
	}

	std::string description() const
	{
		return m_desc.str();
	}

	std::string shortID() const
	{
		return m_desc.str();
	}

	bool check(const std::string& s) const
	{
		return m_map.find(getFileExtension(s)) != m_map.end();
	}

protected:
	t_fextensionmap m_map;
	std::stringstream m_desc;
};

Vector_3 parseVector(const std::string& s)
{
	double x, y, z;
	std::stringstream str(s);
	char ch;

	str >> x;
	str >> ch;	// comma
	str >> y;
	str >> ch;	// comma
	str >> z;

	return Vector_3(x, y, z);
}

t_crystalsize parseCrystalSize(const std::string& s)
{
	t_crystalsize cs;
	std::stringstream str(s);
	char(ch);

	str >> cs.x;
	str >> ch;	// comma
	str >> cs.y;
	str >> ch;	// comma
	str >> cs.z;

	return cs;
}

t_chirality parseHandedness(const std::string& s)
{
	return (s == "r" || s == "right") ? righthanded : lefthanded ; 
}

int main(int argc, char** argv) {

	bool flascii = false;
	bool flsmooth = false;
	t_meshcriteria mc;
	std::vector<std::string> fnames;

	t_fextensionmap fmap;
	fmap["inc"] = filetype_pov;
	fmap["off"] = filetype_off;
	fmap["stl"] = filetype_stl;
	fmap["stla"] = filetype_stl_ascii;
	fmap["stlb"] = filetype_stl_binary;

	// parse command line arguments
	try
	{
		TCLAP::CmdLine cmd("Compute a 3D surface mesh of a helix crystal.", ' ', "0.9a1", false);

		TCLAP::SwitchArg helpArg("", "help", "Displays usage information and exits.");
		cmd.add(helpArg);

		TCLAP::SwitchArg versionArg("", "version", "Displays version information and exits.");
		cmd.add(versionArg);

		TCLAP::SwitchArg asciiArg("", "ascii", "Switch from binary to ascii mode (for *.off and *.stl files only).");
		cmd.add(asciiArg);

		TCLAP::SwitchArg smoothArg("", "smooth", "Switch from flat to smooth triangles (for povray .inc files only).");
		cmd.add(smoothArg);

		VectorConstraint<double,3> meshConstraint("float");
		TCLAP::ValueArg<std::string> meshArg("", "criteria", "Sets the angular, radius, and distance bounds of the meshing algorithm (default: " + mc.str() + ").", false, mc.str(), &meshConstraint);
		cmd.add(meshArg);

		VectorConstraint<double,3> voxelConstraint("float");
		TCLAP::ValueArg<std::string> voxelArg("v", "voxel", "Sets the voxel diameter (default: 0.1,0.1,0.27).", false, "0.1,0.1,0.27", &voxelConstraint);
		cmd.add(voxelArg);

		std::vector<std::string> handedness;
		handedness.push_back("l");
		handedness.push_back("r");
		handedness.push_back("left");
		handedness.push_back("right");
		TCLAP::ValuesConstraint<std::string> allowedHandedness(handedness);

		TCLAP::ValueArg<std::string> helixArg("h", "helix", "Sets handedness of helices (default: left).", false, "l", &allowedHandedness);
		cmd.add(helixArg);

		TCLAP::ValueArg<std::string> cornerArg("c", "corner", "Sets handedness of corners (default: left).", false, "l", &allowedHandedness);
		cmd.add(cornerArg);

		TCLAP::ValueArg<double> radiusArg("r", "radius", "Sets the radius of the helices (default: 0.45 * a).", false, 0.0, "float");
		cmd.add(radiusArg);

		TCLAP::ValueArg<double> gratingArg("a", "gratingconstant", "Sets the grating constant and helix pitch (default: 1.0).", false, 1.0, "float");
		cmd.add(gratingArg);

		VectorConstraint<int,3> sizeConstraint("int");
		TCLAP::ValueArg<std::string> sizeArg("s", "size", "Sets the number of unit cells (default: 1,1,1).", false, "1,1,1", &sizeConstraint);
		cmd.add(sizeArg);

		FileExtensionConstraint fileConstraint(fmap);
		TCLAP::UnlabeledMultiArg<std::string> fileArg("file", "Output file(s). The file format is deduced from the extension.", false, &fileConstraint);
		cmd.add(fileArg);

		cmd.parse(argc, argv);

		if (helpArg.isSet())
		{
			cmd.getOutput()->usage(cmd);
			return 0;
		}

		if (versionArg.isSet())
		{
			cmd.getOutput()->version(cmd);
			return 0;
		}

		// output files
		flascii = asciiArg.getValue();
		flsmooth = smoothArg.getValue();
		fnames = fileArg.getValue();

		// helix crystal properties
		t_crystalsize cs = parseCrystalSize(sizeArg.getValue());
		double a = gratingArg.getValue();
		double r = radiusArg.getValue();
		Vector_3 voxel = 0.5 * parseVector(voxelArg.getValue());
		t_chirality corner = parseHandedness(cornerArg.getValue());
		t_chirality helix = parseHandedness(helixArg.getValue());

		if (r == 0.0)
		{
			r = 0.45 * a;	// actual default value
		}

		hc = new HelixCrystal(cs, a, r, voxel, corner, helix, /*true*/false);

		// meshing algorithm bounds
		char ch;
		std::stringstream critstring(meshArg.getValue());
		critstring >> mc.angle;
		critstring >> ch;
		critstring >> mc.radius;
		critstring >> ch;
		critstring >> mc.distance;
	}
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
		exit(1);
	}

  boxToCube.set_box(hc->GetBoundingBox());
  
  // dump info
  std::cout << hc->GetBoundingBox() << std::endl;
  std::cout << boxToCube.cube() << std::endl;
  
  std::cout << hc->GetBoundingSphere() << std::endl;
  
  Periodic_mesh_domain domain(implicit_function, boxToCube.cube(), 1e-9);
  h_tr.set_domain(domain.periodic_cuboid());
  
  Mesh_criteria criteria(domain, facet_angle=mc.angle, facet_size = mc.radius, facet_distance= mc.distance,
                         cell_radius_edge_ratio=2, cell_size = 0.05);  
  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_periodic_mesh_3<C3t3>(domain, criteria);
  
  // Output
  std::ofstream medit_file("out.mesh");
  
  // we set stretched values of the vertices of tr
  // Attantion: tr will become invalid
  typedef C3t3::Triangulation Triangulation;
  typedef Triangulation Tr;
  typedef Tr::Vertex_iterator Vertex_iterator;
  
  Tr& tr = c3t3.triangulation();
  
  for(Vertex_iterator vit = tr.vertices_begin(); vit != tr.vertices_end(); ++vit) {
    // point from a cube
    Point p = vit->point();
    
    // stretch
    Point q = boxToCube.fromCubeToBox(p);
    
    vit->set_point(q);
  }
  write_complex_to_medit(medit_file, c3t3);
  
  medit_file.close();

   
	delete hc;
}
