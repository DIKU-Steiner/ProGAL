package ProGAL.geom3d.tess;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.viewer.ClickListener;
import ProGAL.geom3d.viewer.J3DScene;

public class TPoint{
	double x,y,z;
	double sq;//lifted value = x*x+y*y+z*z-t*t
	double xv,yv,zv,sqv;// space for differences
	double p; // "pressure" sample from scalar field

	public String toString(){
		return String.format("TPoint[x=%.1f,y=%.1f,z=%.1f,p=%.1f]",x,y,z,p);
	}
}
