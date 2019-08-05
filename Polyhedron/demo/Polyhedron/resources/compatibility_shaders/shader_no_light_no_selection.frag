
varying highp vec4 color;
varying highp float dist[6];
uniform bool is_clipbox_on;
uniform highp float alpha;
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
  gl_FragColor = vec4(color.xyz, alpha);
}  
