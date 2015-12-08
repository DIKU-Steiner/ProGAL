package ProGAL.geom3d.volumes;

import java.awt.Color;

import ProGAL.geom3d.Plane;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Simplex;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.viewer.J3DScene;

/** 
 * A tetrahedron is a polyhedron with four triangular faces. It is defined using 
 * four corner-points.
 *  
 * @author R.Fonseca
 */
public class Tetrahedron implements Simplex, Volume {
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
	
	/* returns a tetrahedron with its circumscribing circle at (0,0,0) and four corners at unit distance */
	public static Tetrahedron regularTetrahedron() {
		return new Tetrahedron(new Point( 1.0,      0.0,               0.0),
						       new Point(-1.0/3.0,  Math.sqrt(8)/3.0,  0.0),
					     	   new Point(-1.0/3.0, -Math.sqrt(2)/3.0,  Math.sqrt(2.0/3.0)),
							   new Point(-1.0/3.0, -Math.sqrt(2)/3.0, -Math.sqrt(2.0/3.0)));
	}
	
	/** Return the specified corner. Throws an error if <code>c<0 || c>3</code>. */
	public Point getCorner(int c){
		if(c<0 || c>3) throw new IllegalArgumentException();
		return corners[c];
	}
	
	/** Return all four corners */
	public Point[] getCorners() {
		return corners;
	}
	
	
	public void translate(Vector v) {
		for (Point p : corners) p.translateThis(v.x(), v.y(), v.z());
	}
	public void translate(Point q) {
		for (Point p : corners) p.translateThis(q.x(), q.y(), q.z());
	}
	public void translate(double x, double y, double z) {
		for (Point p : corners) p.translateThis(x, y, z);
	}
	
	public void blowUp(double t) {
		Point center = circumCenter();
		translate(-center.x(), -center.y(), -center.z());
		for (Point p : corners) p.scaleThis(t);
		translate(center);
	}
	
	/** Return the specified corner-point. Throws an error if <code>c<0 || c>3</code>. */
	public Point getPoint(int c){
		if(c<0 || c>3) throw new IllegalArgumentException();
		return corners[c];
	}
	
	/** Return the 'dimension' of this object. Required by the interface Simplex. */
	public int getDimension() { return 3; }

	/** TODO: Comment */
	public void setPoint(int c, Point point) {
		if(c<0 || c>3) throw new IllegalArgumentException();
		corners[c] = point;
	}
	
	public boolean overlaps(Volume vol) {
		// TODO Auto-generated method stub
		return false;
	}
	
	/** Get the volume of the tetrahedron. */
	public double getVolume() {
		Vector a = corners[3].vectorTo(corners[0]);
		Vector b = corners[3].vectorTo(corners[1]);
		Vector c = corners[3].vectorTo(corners[2]);
		return Math.abs(a.dot(b.crossThis(c)))/6.0;
	}


	/** Return common triangle of 2 tetrahedra */
	public Triangle getCommonTriangle(Tetrahedron t) {
		Point[] common = new Point[4];
		int count = 0;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				if (corners[i].equals(t.corners[j])) {
					common[count++] = corners[i];
					j = 4;
				}
			}
		}
		if (count != 3) return null;
		return new Triangle(common[0], common[1], common[2]);	
	}
	
	/** Calculate the radius of the insphere. */
	public double getInradius(){
		Vector a = corners[3].vectorTo(corners[0]);
		Vector b = corners[3].vectorTo(corners[1]);
		Vector c = corners[3].vectorTo(corners[2]);
		Vector bXc = b.cross(c);
		double sixV = Math.abs(a.dot(bXc));
		Vector cXa = c.crossThis(a);
		Vector aXb = a.crossThis(b);
		double denom = bXc.length()+cXa.length()+aXb.length()+( bXc.addThis(cXa).addThis(aXb).length() );
		return sixV/denom;
	}
	
	/** Calculate the radius of the circumsphere. */
	public double circumRadius(){
		Vector a = corners[3].vectorTo(corners[0]);
		Vector b = corners[3].vectorTo(corners[1]);
		Vector c = corners[3].vectorTo(corners[2]);
		Vector O = b.cross(c).multiplyThis(a.dot(a));
		O.addThis(c.cross(a).multiplyThis(b.dot(b)));
		O.addThis(a.cross(b).multiplyThis(c.dot(c)));
		O.multiplyThis(1.0/(2*a.dot(b.crossThis(c))));
		return O.length();
	}
	

	/** Find the center of the circumscribing sphere. */
	public Point circumCenter(){
		Vector a = corners[3].vectorTo(corners[0]);
		Vector b = corners[3].vectorTo(corners[1]);
		Vector c = corners[3].vectorTo(corners[2]);
		Vector O = b.cross(c).multiplyThis(a.dot(a));
		O.addThis(c.cross(a).multiplyThis(b.dot(b)));
		O.addThis(a.cross(b).multiplyThis(c.dot(c)));
		O.multiplyThis(1.0/(2*a.dot(b.crossThis(c))));
		return corners[3].add(O);
	}
	
	
	/** Find the circumscribing sphere */
	public Sphere circumSphere(){
		Vector a = corners[3].vectorTo(corners[0]);
		Vector b = corners[3].vectorTo(corners[1]);
		Vector c = corners[3].vectorTo(corners[2]);
		Vector O = b.cross(c).multiplyThis(a.dot(a));
		O.addThis(c.cross(a).multiplyThis(b.dot(b)));
		O.addThis(a.cross(b).multiplyThis(c.dot(c)));
		O.multiplyThis(1.0/(2*a.dot(b.crossThis(c))));
		return new Sphere(corners[3].add(O), O.length());
	}

	/** Find the center of the inscribed sphere. */
	public Point incenter(){
		Vector a = corners[3].vectorTo(corners[0]);
		Vector b = corners[3].vectorTo(corners[1]);
		Vector c = corners[3].vectorTo(corners[2]);
		Vector bXc = b.cross(c);
		Vector cXa = c.cross(a);
		Vector aXb = a.cross(b);
		double bXcLength = bXc.length();
		double cXaLength = cXa.length();
		double aXbLength = aXb.length();
		double dLength = bXc.addThis(cXa).addThis(aXb).length();
		Vector O = a.multiplyThis(bXcLength);
		O.addThis(b.multiplyThis(cXaLength));
		O.addThis(c.multiplyThis(aXbLength));
		O.divideThis(bXcLength+cXaLength+aXbLength+dLength );
		return corners[3].add(O);
	}
	
	
	public Point getCenter() {
		Vector v = corners[0].vectorTo(corners[1]);
		v.addThis(corners[0].vectorTo(corners[2]));
		v.addThis(corners[0].vectorTo(corners[3]));
		return corners[0].add(v.multiplyThis(0.25));
	}

	/** Returns true if the point p is inside this tetrahedron. */
//	public boolean isInside(Point p) {
//		return isBehind(p,p1,p3,p2) && isBehind(p,p0,p2,p3) && isBehind(p,p0,p3,p1) && isBehind(p,p0,p1,p2);
		// TODO implement 
//		return false;
//	}
	
	public boolean isInside(Point p) {
		Plane pl012 = new Plane(getCorner(0), getCorner(1), getCorner(2));
		Plane pl013 = new Plane(getCorner(0), getCorner(1), getCorner(3));
		Plane pl023 = new Plane(getCorner(0), getCorner(2), getCorner(3));
		Plane pl123 = new Plane(getCorner(1), getCorner(2), getCorner(3));
		return (((pl012.above(p) == 1) && (pl013.above(p) == 1) && (pl023.above(p) == 1) && (pl123.above(p) == 1)) ||
				((pl012.below(p) == 1) && (pl013.below(p) == 1) && (pl023.below(p) == 1) && (pl123.below(p) == 1)));
	}

	/*
	 * returns TRUE if the tetrahedron is acute. Tetrahedron is acute if all its dihedral angles are acute (< 90�)
	 * added by pawel 12-11-2011
	 */

	public boolean isAcute() {
		return ((Point.getCosDihedralAngle(corners[0], corners[1], corners[2], corners[3]) > 0.0) &&
				(Point.getCosDihedralAngle(corners[0], corners[1], corners[3], corners[2]) > 0.0) &&
				(Point.getCosDihedralAngle(corners[0], corners[2], corners[3], corners[1]) > 0.0) &&
				(Point.getCosDihedralAngle(corners[2], corners[0], corners[1], corners[3]) > 0.0) &&
				(Point.getCosDihedralAngle(corners[1], corners[0], corners[2], corners[3]) > 0.0) &&
				(Point.getCosDihedralAngle(corners[1], corners[0], corners[3], corners[2]) > 0.0));
	}

	
	public Volume clone(){
		return new Tetrahedron(corners[0].clone(), corners[1].clone(), corners[2].clone(), corners[3].clone());
	}


	/** Return a string representation of this tetrahedron. */
	public String toString() {
		return toString(2);
	}
	
	/** Return a string representation of this tetrahedron with <code>dec</code> decimals precision */
	public String toString(int dec) {
		return String.format("Tetrahedron[%s,%s,%s,%s]",corners[0].toString(dec),corners[1].toString(dec),corners[2].toString(dec),corners[3].toString(dec));
	}
	
	/** Writes this tetrahedron to <code>System.out</code>. */
	public void toConsole() { System.out.println(toString()); }
	
	/** Writes this tetrahedron to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec) { System.out.println(toString(dec)); }
	
	public void toScene(J3DScene scene) {
		for (int i = 0; i < 4; i++)
			for (int j = i+1; j < 4; j++) scene.addShape(new LSS(corners[i], corners[j], 0.01), Color.BLUE, 3);	
	}
	
	public static void main(String[] args) {
		Point p1 = new Point(1,0,0);
		Point p2 = new Point(1,1,0);
		Point p3 = new Point(2,2,3);
		Point p4 = new Point(3,4,2);
		Point p5 = new Point(0,0,3);
		Tetrahedron tetr = new Tetrahedron(p2, p1, p3, p4);
		if (tetr.isInside(p5)) System.out.println("inside"); else System.out.println("outside");

	}
}
