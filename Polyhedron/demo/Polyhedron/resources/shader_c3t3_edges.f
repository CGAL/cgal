#version 150
in vec4 color;
uniform bool is_surface;
out vec4 out_color;
void main(void) 
{ 
  if(color.w<0 || is_surface)
    out_color = vec4(0,0,0,1.0); 
  else
    discard;
}  
