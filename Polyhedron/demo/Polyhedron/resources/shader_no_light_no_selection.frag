#version 150
in vec4 color;
in float dist[6];
uniform bool is_clipbox_on;
uniform float alpha;
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
  out_color = vec4(color.xyz, alpha); 
}  
