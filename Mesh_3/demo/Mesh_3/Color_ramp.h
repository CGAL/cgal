#ifndef _COLOR_RAMP_H
#define _COLOR_RAMP_H

class Color_ramp
{
public :
	Color_ramp();
	~Color_ramp() {}

public :
	unsigned char r(unsigned int index) const { return m_colors[0][index%256]; }
	unsigned char g(unsigned int index) const { return m_colors[1][index%256]; }
	unsigned char b(unsigned int index) const { return m_colors[2][index%256]; }
  
	void add_node(unsigned int index,
                unsigned char r,
                unsigned char g,
                unsigned char b);
  
	void build_red();
	void build_thermal();
	void build_blue();

private:
	void rebuild();
	void reset();

private :
	int m_nodes[256];
	unsigned char m_colors[4][256];
};

#endif // _COLOR_RAMP_H
