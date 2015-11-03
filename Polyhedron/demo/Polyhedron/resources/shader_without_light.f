#version 120
varying highp vec4 color;
uniform bool is_selected;
void main(void) 
{ 
if(is_selected)
  gl_FragColor = vec4(0,0,0,1.0); 
else
  gl_FragColor = color; 
}  
