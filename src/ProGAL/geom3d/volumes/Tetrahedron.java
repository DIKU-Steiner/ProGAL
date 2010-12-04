package ProGAL.geom3d.volumes;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;

/** 
 * A tetrahedron is a polyhedron with three triangular faces. It is defined using 
 * four corner-points.
 *  
 * @author R. Fonseca
 */
public class Tetrahedron implements Volume {
	protected Point[] corners = new Point[4];

	public Tetrahedron(Point p1, Point p2, Point p3, Point p4){
		corners[0] = p1;
		corners[1] = p2;
		corners[2] = p3;
		corners[3] = p4;
	}
	public Tetrahedron(Point[] corners){
		this.corners = corners; 
	}
	
	/** Return the specified corner. Throws an error if <code>c<0 || c>3</code>. */
	public Point getCorner(int c){
		if(c<0 || c>3) throw new IllegalArgumentException();
		return corners[c];
	}
	
	/** Return the specified corner-point. Throws an error if <code>c<0 || c>3</code>. */
	public Point getPoint(int c){
		if(c<0 || c>3) throw new IllegalArgumentException();
		return corners[c];
	}
	
	/** Get the volume of the tetrahedron. */
	public double getVolume() {
		Vector a = corners[3].vectorTo(corners[0]);
		Vector b = corners[3].vectorTo(corners[1]);
		Vector c = corners[3].vectorTo(corners[2]);
		return Math.abs(a.dot(b.crossThis(c)))/6.0;
	}

	/** Calculate the radius of the circumsphere. */
	public double getCircumradius(){
		Vector a = corners[3].vectorTo(corners[0]);
		Vector b = corners[3].vectorTo(corners[1]);
		Vector c = corners[3].vectorTo(corners[2]);
		Vector O = b.cross(c).multiplyThis(a.dot(a));
		O.addThis(c.cross(a).multiplyThis(b.dot(b)));
		O.addThis(a.cross(b).multiplyThis(c.dot(c)));
		O.multiplyThis(1.0/(2*a.dot(b.crossThis(c))));
		return O.length();
	}
	
	public boolean overlaps(Volume vol) {
		// TODO Auto-generated method stub
		return false;
	}

	/** Find the center of the circumscribing sphere. */
	public Point getCircumcenter(){
		Vector a = corners[3].vectorTo(corners[0]);
		Vector b = corners[3].vectorTo(corners[1]);
		Vector c = corners[3].vectorTo(corners[2]);
		Vector O = b.cross(c).multiplyThis(a.dot(a));
		O.addThis(c.cross(a).multiplyThis(b.dot(b)));
		O.addThis(a.cross(b).multiplyThis(c.dot(c)));
		O.multiplyThis(1.0/(2*a.dot(b.crossThis(c))));
		return corners[3].add(O);
	}
	
	public Point getCenter() {
		return getCircumcenter();
	}

	/** Returns true if the point p is inside this tetrahedron. */
	public boolean isInside(Point p) {
//		return isBehind(p,p1,p3,p2) && isBehind(p,p0,p2,p3) && isBehind(p,p0,p3,p1) && isBehind(p,p0,p1,p2);
		// TODO implement 
		return false;
	}

	
	
	public Volume clone(){
		return new Tetrahedron(corners[0].clone(), corners[1].clone(), corners[2].clone(), corners[3].clone());
	}
}
