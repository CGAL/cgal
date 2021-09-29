#version 150
in vec4 color;
flat in vec2 subdomain_out;
uniform bool is_surface;
uniform vec4 is_visible_bitset;
uniform bool is_filterable;
out vec4 out_color;
void main(void) 
{ 
  if(is_filterable)
  {
    uint domain1 = uint(subdomain_out.x);
    uint domain2 = uint(subdomain_out.y);
    uint i1 = domain1/32u;
    uint i2 = domain2/32u;
    uint visible1 = uint(is_visible_bitset[i1]);
    uint visible2 = uint(is_visible_bitset[i2]);
    if((visible1>>(domain1%32u))%2u == 0u && (visible2>>(domain2%32u))%2u == 0u)
    {
      discard;
    }
  }

  if(color.w<0 || is_surface)
    out_color = vec4(0,0,0,1.0); 
  else
    discard;
}  
