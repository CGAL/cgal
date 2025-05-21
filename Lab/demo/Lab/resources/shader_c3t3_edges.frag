#version 150

const int number_of_bitset = 4;

in vec4 color;
in float dist[6];
flat in vec2 subdomain_out;
uniform bool is_clipbox_on;
uniform bool is_surface;
uniform bool is_filterable;
uniform int is_visible_bitset[number_of_bitset];
out vec4 out_color;
void main(void) 
{ 
  if(is_clipbox_on)
    if(dist[0]>0.0 ||
      dist[1]>0.0 ||
      dist[2]>0.0 ||
      dist[3]>0.0 ||
      dist[4]>0.0 ||
      dist[5]>0.0)
        discard;
  if(is_filterable)
  {
    uint domain1 = uint(subdomain_out.x);
    uint domain2 = uint(subdomain_out.y);
    uint i1 = domain1/32u;
    uint j1 = domain1%32u;
    uint i2 = domain2/32u;
    uint j2 = domain2%32u;
    uint visible1 = uint(is_visible_bitset[i1]);
    uint visible2 = uint(is_visible_bitset[i2]);
    if(((visible1>>j1)&1u) == 0u && ((visible2>>j2)&1u) == 0u)
    {
      discard;
    }
  }

  if(color.w<0 || is_surface)
    out_color = vec4(0,0,0,1.0); 
  else
    discard;
}  
