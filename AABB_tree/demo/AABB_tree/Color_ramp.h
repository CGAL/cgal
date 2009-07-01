#ifndef _COLOR_RAMP_H
#define _COLOR_RAMP_H

class Color_ramp
{
private :
  unsigned char m_colors[4][256];
  int m_nodes[256];

public :
  // life cycle
  Color_ramp() 
  { 
    build_thermal(); 
  } 
  ~Color_ramp() {}

public :
  unsigned char* color(unsigned int index) { return m_colors[index%256]; }
  const unsigned char r(unsigned int index) const { return m_colors[0][index%256]; }
  const unsigned char g(unsigned int index) const { return m_colors[1][index%256]; }
  const unsigned char b(unsigned int index) const { return m_colors[2][index%256]; }

private:

  void rebuild()
  {
    // build nodes
    m_colors[3][0] = 1;
    m_colors[3][255] = 1;
    unsigned int nb_nodes = 0;
    for(int i=0;i<256;i++)
      if(m_colors[3][i])
      {
        m_nodes[nb_nodes] = i;
        nb_nodes++;
      }
    // build Color_ramp
    for(int k=0;k<3;k++)
      for(unsigned int i=0;i<(nb_nodes-1);i++)
      {
        int x1 = m_nodes[i];
        int x2 = m_nodes[i+1];
        int y1 = m_colors[k][x1];
        int y2 = m_colors[k][x2];
        float a = (float)(y2-y1) / (float)(x2-x1);
        float b = (float)y1 - a*(float)x1;
        for(int j=x1;j<x2;j++)
          m_colors[k][j] = (unsigned char)(a*(float)j+b);
      }
  }

  void reset()
  {
    for(int i=1;i<=254;i++)
      m_colors[3][i] = 0;
    m_colors[3][0] = 1;
    m_colors[3][255] = 1;
  }

public:
  // predefined
  void build_default()
  {
    reset();
    add_node(0,0,0,0);
    add_node(255,255,255,255);
    rebuild();
  }
  void add_node(unsigned int index,
                unsigned char r,
                unsigned char g,
                unsigned char b)
  {
    m_colors[3][index] = 1;
    m_colors[0][index] = r;
    m_colors[1][index] = g;   
    m_colors[2][index] = b;
  }
  void build_menthol()
  {
    reset();
    add_node(0,0,0,128);
    add_node(48,0,0,254);
    add_node(176,0,255,255);
    add_node(208,128,255,255);
    add_node(255,255,255,255);
    rebuild();
  }

  void build_fire()
  {
    reset();
    m_colors[3][0] = 1;m_colors[0][0] = 0;m_colors[1][0] = 0;m_colors[2][0] = 0;
    m_colors[3][92] = 1;m_colors[0][92] = 255;m_colors[1][92] = 0;m_colors[2][92] = 0;
    m_colors[3][185] = 1;m_colors[0][185] = 255;m_colors[1][185] = 255;m_colors[2][185] = 0;
    m_colors[3][255] = 1;m_colors[0][255] = 255;m_colors[1][255] = 255;m_colors[2][255] = 0;
    rebuild();
  }
  void build_rainbow()
  {
    reset();
    add_node(0,0,0,0);
    add_node(32,128,0,255);
    add_node(80,0,0,255);
    add_node(112,0,255,255);
    add_node(144,0,255,0);
    add_node(176,255,255,0);
    add_node(208,255,128,0);
    add_node(240,255,0,0);
    add_node(255,255,255,255);
    rebuild();
  }

  void build_infrared()
  {
    reset();
    add_node(0,0,0,0);
    add_node(32,0,0,128);
    add_node(64,128,0,255);
    add_node(112,255,0,255);
    add_node(160,255,128,0);
    add_node(176,255,128,0);
    add_node(208,255,255,0);
    add_node(255,255,255,255);
    rebuild();
  }

	void build_light_rainbow()
  {
    reset();
    m_colors[3][0]   = 1; m_colors[0][0]   = 128; m_colors[1][0]   = 255; m_colors[2][0] = 255;
    m_colors[3][64]  = 1; m_colors[0][64]  = 168; m_colors[1][64] = 255;  m_colors[2][64] = 168;
    m_colors[3][128] = 1; m_colors[0][128] = 255; m_colors[1][128] = 255; m_colors[2][128] = 179;
    m_colors[3][192] = 1; m_colors[0][192] = 255; m_colors[1][192] = 175; m_colors[2][192] = 175;
    m_colors[3][255] = 1; m_colors[0][255] = 255; m_colors[1][255] = 255; m_colors[2][255] = 255;
    rebuild();
  }
  void build_thermal()
  {
    reset();
    m_colors[3][0] = 1;m_colors[0][0] = 0;m_colors[1][0] = 0;m_colors[2][0] = 0;
    m_colors[3][70] = 1;m_colors[0][70] = 128;m_colors[1][70] = 0;m_colors[2][70] = 0;
    m_colors[3][165] = 1;m_colors[0][165] = 255;m_colors[1][165] = 128;m_colors[2][165] = 0;
    m_colors[3][255] = 1;m_colors[0][255] = 255;m_colors[1][255] = 255;m_colors[2][255] = 255;
    rebuild();
  }
	void build_multicolor()
  {
    reset();
    m_colors[3][0]   = 1; m_colors[0][0]   = 255;   m_colors[1][0]  = 128;   m_colors[2][0] = 192;
    m_colors[3][32]  = 1; m_colors[0][32]  = 0;   m_colors[1][32]  = 0;   m_colors[2][32] = 255;
    m_colors[3][64]  = 1; m_colors[0][64]  = 0;   m_colors[1][64]  = 255; m_colors[2][64] = 255;
    m_colors[3][96] = 1; m_colors[0][96] = 0;     m_colors[1][96] = 252; m_colors[2][96] = 0;
    m_colors[3][128] = 1; m_colors[0][128] = 255; m_colors[1][128] = 254; m_colors[2][128] = 0;
    m_colors[3][176] = 1; m_colors[0][176] = 255; m_colors[1][176] = 128; m_colors[2][176] = 0;
    m_colors[3][208] = 1; m_colors[0][208] = 255; m_colors[1][208] = 0; m_colors[2][208] = 0;
    m_colors[3][240] = 1; m_colors[0][240] = 192; m_colors[1][240] = 192;   m_colors[2][240] = 192;
    m_colors[3][255] = 1; m_colors[0][255] = 255; m_colors[1][255] = 255; m_colors[2][255] = 255;
    rebuild();
  }

	void build_random_pastel()
  {
    reset();
    for(int i=0;i<256;i++)
    {
      m_colors[3][i] = 1;
      m_colors[1][i] = 128 + (unsigned char) ((double)rand() / (double)RAND_MAX * 127.0);
      m_colors[2][i] = 128 + (unsigned char) ((double)rand() / (double)RAND_MAX * 127.0);
      m_colors[3][i] = 128 + (unsigned char) ((double)rand() / (double)RAND_MAX * 127.0);
    }
  }
};

#endif // _COLOR_RAMP_H
