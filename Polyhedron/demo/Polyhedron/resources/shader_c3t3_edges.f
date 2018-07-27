#version 150
in vec4 color;
out vec4 out_color;
void main(void) 
{ 
  if(color.w<0)
    out_color = vec4(0,0,0,1.0); 
  else
    discard;
}  
