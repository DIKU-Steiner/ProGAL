package ProGAL.geom3d.surface;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Shape;
import ProGAL.geom3d.Vector;

public abstract class ParametricSurface implements Shape{
	
	public abstract Point getPoint(double u, double v);

	public Point getPoint(ProGAL.geom2d.Point p) {
		return getPoint(p.x(), p.y());
	}

	public Vector getNormal(double u, double v){
		Point p1 = getPoint(u,v);
		Point p2 = getPoint(u+0.001,v);
		Point p3 = getPoint(u,v+0.001);
		return p1.vectorTo(p2).crossThis(p1.vectorTo(p3)).normalizeThis();
		
	}

	public Vector getNormal(ProGAL.geom2d.Point p){
		return getNormal(p.x(), p.y());
	}
	
	public Point getCenter(){ return getPoint(0,0); }

	public abstract ParametricSurface clone();
}
