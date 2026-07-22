
varying highp vec4 color;
varying highp float dist[6];
uniform bool is_selected;
uniform bool is_clipbox_on;
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

if(is_selected)
  gl_FragColor = vec4(0,0,0,1.0);
else
  gl_FragColor = color; 
}  
