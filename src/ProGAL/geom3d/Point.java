package ProGAL.geom3d;

import ProGAL.math.Constants;

/** 
 *  A point in (x,y,z)-space represented with double precision. 
 */
public class Point {
	protected double x,y,z;

	/** Construct a point with the specified coordinates. */
	public Point(double x, double y, double z) { this.x = x; this.y = y; this.z = z; }
	/** Construct a point that is a clone of p. */
	public Point(Point p) { x = p.x; y = p.y; z = p.z; }
	/** Construct a point at the coordinates of v. */
	public Point(Vector v) { x = v.x; y = v.y; z = v.z; }

	
	/** Get the i'th coordinate. */
	public double getCoord(int i) {
		switch(i){
		case 0: return x;
		case 1: return y;
		case 2: return z;
		}
		throw new Error("Trying to get invalid coordinate");
	}
	
	/** Get the i'th coordinate. */
	public double get(int i) { return getCoord(i); }

	/** Get the first coordinate. */
	public double getX() { return x; }
	/** Get the second coordinate. */
	public double getY() { return y; }
	/** Get the third coordinate. */
	public double getZ() { return z; }
	
	/** Set the i'th coordinate to v */
	public void setCoord(int i, double v){
		switch(i){
		case 0: this.x = v;return;
		case 1: this.y = v;return;
		case 2: this.z = v;return;
		}
	}
	
	/** Set the i'th coordinate to v */
	public void set(int i, double v){ setCoord(i,v); }

	/** Set the first coordinate */
	public void setX(double x) { this.x = x; }
	/** Set the second coordinate */
	public void setY(double y) { this.y = y; }
	/** Set the third coordinate */
	public void setZ(double z) { this.z = z; }
	
	/** Get the vector that points from this point to p */
	public Vector vectorTo(Point p){
		return new Vector(p.x-x, p.y-y, p.z-z);
	}

	/** 
	 * Returns true if three points are on the same line. This implies that
	 * overlapping points are considered collinear.
	 */
	public static boolean collinear(Point p0, Point p1, Point p2) {
		Vector v1v0 = p1.vectorTo(p0);
		Vector v1v2 = p1.vectorTo(p2);
		return v1v0.cross(v1v2).lengthSquared()<Constants.EPSILON;
	}

	/** Returns true if four specified points are in the same plane	 */
	public static boolean coplanar(Point p0, Point p1, Point p2, Point p3) {
		double ax = p0.x; double ay = p0.y; double az = p0.z;
		double bx = p1.x; double by = p1.y; double bz = p1.z;
		double cx = p2.x; double cy = p2.y; double cz = p2.z;
		double dx = p3.x; double dy = p3.y; double dz = p3.z;
		return   Math.abs(
				-az*by*cx + ay*bz*cx + az*bx*cy - ax*bz*cy
				-ay*bx*cz + ax*by*cz + az*by*dx - ay*bz*dx
				-az*cy*dx + bz*cy*dx + ay*cz*dx - by*cz*dx
				-az*bx*dy + ax*bz*dy + az*cx*dy - bz*cx*dy
				-ax*cz*dy + bx*cz*dy + ay*bx*dz - ax*by*dz
				-ay*cx*dz + by*cx*dz + ax*cy*dz - bx*cy*dz	) < Constants.EPSILON; 
	}


	/** Translates this point by (x,y,z). */
	public void translate(double dx, double dy, double dz) {
		this.x += dx;
		this.y += dy;
		this.z += dz;
	}
	
	/** Scale this point by a factor s */
	public void scale(double s){
		this.x*=s;
		this.y*=s;
		this.z*=s;
	}
	
	/** Returns p added to this (changing this object). */
	public Point addThis(Vector p) { translate(p.x,p.y,p.z); return this; }
	
	/** Returns p added to this (without changing this object). */
	public Point add(Vector p) { return new Point(x+p.x, y+p.y, z+p.z); }
	
	/** Returns p subtracted from this (changing this object). */
	public Point subtractThis(Vector p) { x-=p.x; y-=p.y; z-=p.z; return this;	}
	
	/** Returns p subtracted from this (without changing this object). */
	public Point subtract(Vector p) {return new Point(x-p.x,y-p.y,z-p.z);	}
	
	/** Reflects this point through origo. */
	public Point reflectThroughOrigoThis() { x*=-1; y*=-1; z*=-1; return this; }

	/** Get the squared distance from this point to point q. */
	public double getDistanceSquared(Point q) {
		double dx = x-q.x;
		double dy = y-q.y;
		double dz = z-q.z;
		return dx*dx+dy*dy+dz*dz;
	}

	/** Get the distance from this point to point q */
	public double getDistance(Point q) { 
		double dx = x-q.x;
		double dy = y-q.y;
		double dz = z-q.z;
		return Math.sqrt(dx*dx+dy*dy+dz*dz); 
	}

	/** Creates a bisector between points p and q */
	public static Plane getBisector(Point p, Point q) {
		if (!p.equals(q)) 
			return new Plane(getMidpoint(p, q), p.vectorTo(q).normalizeThis()); 
		else return null;
	}
	
	/** Creates the midpoint of two points. */
	public static Point getMidpoint(Point p, Point q) {
		return new Point((p.x + q.x)/2,(p.y + q.y)/2,(p.z + q.z)/2);		
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

	// COMPARISON METHODS

	/**
	 * Returns true if this point dominates point q. One point is said to dominate another 
	 * if it has a higher x-coordinate. If two points have identical x-coordinates, the 
	 * y-coordinate is considered and so forth.
	 */
	public boolean dominates(Point q) { 
		if (x > q.x) return true;
		if (x < q.x) return false;
		if (y > q.y) return true;
		if (y < q.y) return false; 
		return z > q.z;
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

	/** Returns true iff o is a point that equals this point. */
	public boolean equals(Object o){
		if(o instanceof Point) return equals((Point)o);
		else return false;
	}
	
	/** Returns true iff this point and point p are overlapping. */
	public boolean equals(Point p) { return x == p.x && y == p.y && z == p.z; }

	/** Return a new object that equals this object. */
	public Point clone(){ return new Point(x,y,z); }
	
	/** Returns the vector from origo to this point. Converts this point to a vector. */
	public Vector toVector() { return new Vector(x, y, z); }

	/** Convert this point to a double-array */
	public double[] toDoubleArray() { 
		return new double[]{x, y, z}; 
	}

	/** Returns a string-representation of this point formatted with two decimals precision. */ 
	public String toString() {	return toString(2); }

	/** Returns a string-representation of this point formatted with <code>dec</code> decimals precision. */ 
	public String toString(int dec) {
		return String.format("Point3d[%."+dec+"f,%."+dec+"f,%."+dec+"f]", x,y,z); 
	}	

	/** Writes this point to <code>System.out</code>. */
	public void toConsole() { System.out.println(toString()); }
	/** Writes this point to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec) { System.out.println(toString(dec)); }


}



