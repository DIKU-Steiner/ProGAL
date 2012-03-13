package ProGAL.geom3d;


import java.awt.Color;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;

import ProGAL.math.Constants;

/** 
 *  A point in (x,y,z)-space represented with double precision. 
 */
public class Point extends ProGAL.geomNd.Point implements Simplex{
	private static final long serialVersionUID = -2120468832687547475L;

	/** Construct a point with the specified coordinates. */
	public Point(double x, double y, double z) { 
		super(new double[]{x,y,z});
	}
	/** Construct a point with the specified coordinates. */
	public Point(double[] coords) { 
		super(coords);
	}
	/** Construct a point that is a clone of p. */
	public Point(Point p) { 
		super(new double[]{p.coords[0],p.coords[1],p.coords[2]});
	}
	/** Construct a point at the coordinates of v. */
	public Point(Vector v) { 
		super(new double[]{v.x(),v.y(),v.z()});
	}

	
//	/** Get the i'th coordinate. */
//	public double getCoord(int i) {
//		switch(i){
//		case 0: return coords[0];
//		case 1: return coords[1];
//		case 2: return coords[2];
//		}
//		throw new Error("Trying to get invalid coordinate");
//	}
	
//	/** Get the i'th coordinate. */
//	public double get(int i) { return getCoord(i); }

	/** Get the first coordinate. */
	public double x() { return coords[0]; }
	/** Get the second coordinate. */
	public double y() { return coords[1]; }
	/** Get the third coordinate. */
	public double z() { return coords[2]; }
	
	/** Set the first coordinate */
	public void setX(double x) { this.coords[0] = x; }
	/** Set the second coordinate */
	public void setY(double y) { this.coords[1] = y; }
	/** Set the third coordinate */
	public void setZ(double z) { this.coords[2] = z; }
	
	/** 
	 * Return the 'dimension' of this object. Required by the interface Simplex.
	 * Beware not to confuse this method with getDimensions from geomNd.Point.  
	 */
	public int getDimension() { return 0; }
	
	/** Get the vector that points from this point to p */
	public Vector vectorTo(Point p){
		return new Vector(p.coords[0]-coords[0], p.coords[1]-coords[1], p.coords[2]-coords[2]);
	}

	/** 
	 * Returns true if three points are on the same line. This implies that
	 * overlapping points are considered collinear.
	 */
	public static boolean collinear(Point p0, Point p1, Point p2) {
		Vector v1v0 = p1.vectorTo(p0);
		Vector v1v2 = p1.vectorTo(p2);
		return v1v0.cross(v1v2).getLengthSquared()<Constants.EPSILON;
	}

	/** Returns true if four specified points are in the same plane	 */
	public static boolean coplanar(Point p0, Point p1, Point p2, Point p3) {
		double ax = p0.coords[0]; double ay = p0.coords[1]; double az = p0.coords[2];
		double bx = p1.coords[0]; double by = p1.coords[1]; double bz = p1.coords[2];
		double cx = p2.coords[0]; double cy = p2.coords[1]; double cz = p2.coords[2];
		double dx = p3.coords[0]; double dy = p3.coords[1]; double dz = p3.coords[2];
		return   Math.abs(
				-az*by*cx + ay*bz*cx + az*bx*cy - ax*bz*cy
				-ay*bx*cz + ax*by*cz + az*by*dx - ay*bz*dx
				-az*cy*dx + bz*cy*dx + ay*cz*dx - by*cz*dx
				-az*bx*dy + ax*bz*dy + az*cx*dy - bz*cx*dy
				-ax*cz*dy + bx*cz*dy + ay*bx*dz - ax*by*dz
				-ay*cx*dz + by*cx*dz + ax*cy*dz - bx*cy*dz	) < Constants.EPSILON; 
	}


	/** Translates this point by (x,y,z). */
	public void translateThis(double dx, double dy, double dz) {
		this.coords[0] += dx;
		this.coords[1] += dy;
		this.coords[2] += dz;
	}
	
	/** Scale this point by a factor s */
	public void scaleThis(double s){
		this.coords[0]*=s;
		this.coords[1]*=s;
		this.coords[2]*=s;
	}
	
	/** Returns p added to this (changing this object). */
	public Point addThis(Vector p) { translateThis(p.x(),p.y(),p.z()); return this; }

	public Point addThis(Point p) { translateThis(p.x(), p.y(), p.z()); return this; }
	
	/** Returns p added to this (without changing this object). */
	public Point add(Vector p) { return new Point(coords[0]+p.x(), coords[1]+p.y(), coords[2]+p.z()); }
	
	/** Returns p subtracted from this (changing this object). */
	public Point subtractThis(Vector p) { coords[0]-=p.x(); coords[1]-=p.y(); coords[2]-=p.z(); return this;	}
	
	/** Returns p subtracted from this (without changing this object). */
	public Point subtract(Vector p) {return new Point(coords[0]-p.x(),coords[1]-p.y(),coords[2]-p.z());	}
	
	/** Reflects this point through origo. */
	public Point reflectThroughOrigoThis() { coords[0]*=-1; coords[1]*=-1; coords[2]*=-1; return this; }

	/** Get the squared distance from this point to point q. */
	public double distanceSquared(Point q) {
		double dx = coords[0]-q.coords[0];
		double dy = coords[1]-q.coords[1];
		double dz = coords[2]-q.coords[2];
		return dx*dx+dy*dy+dz*dz;
	}

	/** Get the distance from this point to point q */
	public double distance(Point q) { 
		double dx = coords[0]-q.coords[0];
		double dy = coords[1]-q.coords[1];
		double dz = coords[2]-q.coords[2];
		return Math.sqrt(dx*dx+dy*dy+dz*dz); 
	}

	public double dot(Point p) { return x()*p.x() + y()*p.y() + z()*p.z(); }
	
	/** Creates a bisector between points p and q */
	public static Plane getBisector(Point p, Point q) {
		if (!p.equals(q)) 
			return new Plane(getMidpoint(p, q), p.vectorTo(q).normalizeThis()); 
		else return null;
	}
	
	/** Creates the midpoint of two points. */
	public static Point getMidpoint(Point p, Point q) {
		return new Point( (p.coords[0] + q.coords[0])/2,(p.coords[1] + q.coords[1])/2,(p.coords[2] + q.coords[2])/2 );		
	}


	// ANGLE METHODS

	/** Get the angle between the line segments p2->p1 and p2->p3. */
	public static double getAngle(Point p1, Point p2, Point p3) {
		return p2.vectorTo(p1).angle(p2.vectorTo(p3));
	}

	/** Get the dihedral angle defined by the 4 non-collinear points p1, p2, p3, p4. */
	public static double getDihedralAngle(Point p1, Point p2, Point p3, Point p4) {
		return Vector.getDihedralAngle(p1.vectorTo(p2), p2.vectorTo(p3), p3.vectorTo(p4));
	}

  	/*
  	 * returns cosinus of the dihedral angle between 4 non-collinear points p1, p2, p3, p4
  	 * added by pawel 12-11-2011
  	 */
  	public static double getCosDihedralAngle(Point p1, Point p2, Point p3, Point p4) {
 		return Vector.getCosDihedralAngle(new Vector(p1, p2), new Vector(p2, p3), new Vector(p3, p4));
  	}

	
	// COMPARISON METHODS

	/**
	 * Returns true if this point dominates point q. One point is said to dominate another 
	 * if it has a higher x-coordinate. If two points have identical x-coordinates, the 
	 * y-coordinate is considered and so forth.
	 */
	public boolean dominates(Point q) { 
		if (coords[0] > q.coords[0]) return true;
		if (coords[0] < q.coords[0]) return false;
		if (coords[1] > q.coords[1]) return true;
		if (coords[1] < q.coords[1]) return false; 
		return coords[2] > q.coords[2];
	}

	/**
	 * Returns true if this point dominates point q (i=0,1,2 is the most important coordinate,
	 * j=0,1,2 is the second most important coordinate and k=0,1,2 is the least important coordinate).
	 */
	public boolean dominates(Point q, int i, int j, int k) {
		if(i==j || i==k || j==k) 
			throw new Error(String.format("i, j and k must be distinct coordinate indices (%d,%d,%d)",i,j,k));
		if (this.getCoord(i) > q.getCoord(i)) return true;
		if (this.getCoord(i) < q.getCoord(i)) return false;
		if (this.getCoord(j) > q.getCoord(j)) return true;
		if (this.getCoord(j) < q.getCoord(j)) return false; 
		return this.getCoord(k) > q.getCoord(k);
	}

	/** 
	 * Returns a clone of this point. Since a point can be interpreted as a geometric shape 
	 * (a 0-simplex) the Shape interface requires the getCenter method to be implemented.
	 * TODO: Test
	 */
	public Point getCenter() {
		return clone();
	}
	
	public Point getPoint(int i){
		if(i!=0) throw new IllegalArgumentException("Invalid index ("+i+") 0-simplex has one point only");
		return this;
	}
	
	/** Returns true iff o is a point that equals this point. */
	public boolean equals(Object o){
		if(o instanceof Point) return equals((Point)o);
		else return false;
	}
	
	/** Returns true iff this point and point p are overlapping. */
	public boolean equals(Point p) {
		if(Math.abs(coords[0]-p.coords[0])>Constants.EPSILON) return false;
		if(Math.abs(coords[1]-p.coords[1])>Constants.EPSILON) return false;
		if(Math.abs(coords[2]-p.coords[2])>Constants.EPSILON) return false;
		return true;
	}

	/** Return a new object that equals this object. */
	public Point clone(){ return new Point(coords[0], coords[1], coords[2]); }
	
	/** Returns the vector from origo to this point. Converts this point to a vector. */
	public Vector toVector() { return new Vector(coords[0], coords[1], coords[2]); }

	/** Returns a string-representation of this point formatted with two decimals precision. */ 
	public String toString() {	return toString(2); }

	/** Returns a string-representation of this point formatted with <code>dec</code> decimals precision. */ 
	public String toString(int dec) {
		return String.format("Point[%."+dec+"f,%."+dec+"f,%."+dec+"f]", coords[0], coords[1], coords[2]); 
	}	
	
	public void toScene(J3DScene scene, double r, Color clr) { scene.addShape(new Sphere(this, r), clr);}
//	public void draw(J3DScene scene, double r) { draw(scene, r, Color.BLUE); }
//	public void draw(J3DScene scene) { draw(scene,0.1f, Color.BLUE); }


	/** Writes this point to <code>System.out</code>. */
	public void toConsole() { System.out.println(toString()); }
	
	/** Writes this point to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec) { System.out.println(toString(dec)); }

}



