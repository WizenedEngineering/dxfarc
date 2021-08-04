# dxfarc
Determines world coordinates (x, y, z) of start and end points of arcs from dxf-files

I need to extract x,y,z coordinates (World Coordinate System, WCS) of starting or end point of an arc that is somewhere in 3D space.
Lines can be extracted easily because the coordinates of start/end point are given directly in WCS. 
However arcs are given in a different way. It gives coordinates of center point, start angle, end angle, radius in Object (arbitrary) Coordinate system (OCS). Additional the 3 coordinates of an extrusion vector are given. 
Assuming that the conversion of OCS data to WCS data is something very simple. I tried some approaches but failed. 
Ok, next step searching internet. The only conversion I could find were arcs laying on x,y-plane. My arcs can have any position and orientation in 3D space. 
I'm a mechanical engineer. Complex (imaginary) numbers are something I came across before but I gave up to wrap my brain around imaginary spaces (quaternions). 
I tried different ways of Euler rotations and voilà I got results. 

This little program reads a .dxf file and extracts line data as well as arc data. The coordinates of start and end point of arc are converted to WCS and a line between those point drawn. Those lines are saved in another .dxf file. I used this new .dxf file to overlay to the original drawing as a visual check. So far it worked for me.  
