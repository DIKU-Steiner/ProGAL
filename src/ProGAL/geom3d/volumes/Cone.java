package ProGAL.geom3d.volumes;

import ProGAL.geom3d.Point;

/** 
 * A three-dimensional cone represented by two points on its central axis. The first point 
 * is the bottom and the second is the tip of the cone.   
 */
public class Cone implements Volume{
	public Point p1, p2;
	public float rad;

	public Cone(Point p1, Point p2, float r){
		this.p1 = p1;
		this.p2 = p2;
		this.rad = r;
	}

	@Override
	public double getVolume() {
		return Math.PI*rad*rad*p1.distance(p2)/3f;
	}


	public Point getCenter() {
		return Point.getMidpoint(p1,p2);
	}

	public Cone clone(){
		return new Cone(p1.clone(), p2.clone(), rad);
	}

	@Override
	public boolean overlaps(Volume vol) {
		throw new RuntimeException("overlaps not implemented");
	}


}

