package ProGAL.geom3d;

/** 
 *  A point in (x,y,z)-space represented with double precision. 
 */
public class Point3d {
	protected double x,y,z;

	/** Construct a point with the specified coordinates. */
	public Point3d(double x, double y, double z) { this.x = x; this.y = y; this.z = z; }
	/** Construct a point that is a clone of p. */
	public Point3d(Point3d p) { x = p.x; y = p.y; z = p.z; }
	/** Construct a point at the coordinates of v. */
	public Point3d(Vector3d v) { x = v.x; y = v.y; z = v.z; }

	
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
	public Vector3d vectorTo(Point3d p){
		return new Vector3d(p.x-x, p.y-y, p.z-z);
	}

	/** Returns true if three points are on the same line. */
	public static boolean colinear(Point3d p0, Point3d p1, Point3d p2) {
		Vector3d v1v0 = p1.vectorTo(p0);
		Vector3d v1v2 = p1.vectorTo(p2);
		return v1v0.cross(v1v2).length()==0;
	}

	/** Returns true if four specified points are in the same plane	 */
	public static boolean coplanar(Point3d p0, Point3d p1, Point3d p2, Point3d p3) {
		double ax = p0.x; double ay = p0.y; double az = p0.z;
		double bx = p1.x; double by = p1.y; double bz = p1.z;
		double cx = p2.x; double cy = p2.y; double cz = p2.z;
		double dx = p3.x; double dy = p3.y; double dz = p3.z;
		return   -az*by*cx + ay*bz*cx + az*bx*cy - ax*bz*cy
		-ay*bx*cz + ax*by*cz + az*by*dx - ay*bz*dx
		-az*cy*dx + bz*cy*dx + ay*cz*dx - by*cz*dx
		-az*bx*dy + ax*bz*dy + az*cx*dy - bz*cx*dy
		-ax*cz*dy + bx*cz*dy + ay*bx*dz - ax*by*dz
		-ay*cx*dz + by*cx*dz + ax*cy*dz - bx*cy*dz == 0.0; 
	}

	/** Returns true if the half-plane with normal vector nv and point p0 is beind p. */
	public static boolean isBehind(Point3d p, Point3d p0, Vector3d nv) {
		Vector3d v = p.vectorTo(p0);
		return nv.dot(v) < 0.0;
	}

	/** Returns true if the half-plane through p0, p1, p2 (forming a counterclockwise triangle) is behind p. */
	public static boolean isBehind(Point3d p, Point3d p0, Point3d p1, Point3d p2) {
		return Point3d.getNormal(p0, p1, p2).dot(p0.vectorTo(p)) < 0.0;
	}

	/**
	 * Returns true if the point p is inside the tetrehedron p0,p1,p2,p3.
	 */
	public static boolean isInside(Point3d p, Point3d p0, Point3d p1, Point3d p2, Point3d p3) {
		return isBehind(p,p1,p3,p2) && isBehind(p,p0,p2,p3) && isBehind(p,p0,p3,p1) && isBehind(p,p0,p1,p2);

	}

	/**
	 * Returns true if the specified point p is inside cone with the specified apex
	 */
	public static boolean isInCone(Point3d p0, Point3d[] face, Point3d p) {
		if (Point3d.isBehind(p0, face[0], face[2], face[1])) {
			Point3d temp = face[1];
			face[1]= face[2];
			face[2] = temp;
		}
		//		System.out.println("p0 = " + p0.toString(3));
		//		System.out.println("Face: " + face[0] + " " + face[1] + " " + face[2]);
		for (int i = 0; i < 3 ; i++) {
			//			System.out.println("p = " + p.toString(3));
			if (!Point3d.isBehind(p, face[i], p0, face[(i+1)%3])) return false;
		}
		return true;
	}

	/** Returns vector orthogonal to the triangle formed by p0, p1, p2. */
	public static Vector3d getNormal(Point3d p0, Point3d p1, Point3d p2) {
		return p0.vectorTo(p1).cross(p0.vectorTo(p2));
	}

	/** Translates this point by (x,y,z). */
	public void translate(double x, double y, double z) {
		this.x += x;
		this.y += y;
		this.z += z;
	}
	
	/** Scale this point by a factor s */
	public void scale(double s){
		this.x*=s;
		this.y*=s;
		this.z*=s;
	}
	
	/** Returns p added to this (changing this object). */
	public Point3d addThis(Vector3d p) { translate(p.x,p.y,p.z); return this; }
	
	/** Returns p added to this (without changing this object). */
	public Point3d add(Vector3d p) { return new Point3d(x+p.x, y+p.y, z+p.z); }
	
	/** Returns p subtracted from this (changing this object). */
	public Point3d subtractThis(Vector3d p) { x-=p.x; y-=p.y; z-=p.z; return this;	}
	
	/** Returns p subtracted from this (without changing this object). */
	public Point3d subtract(Vector3d p) {return new Point3d(x-p.x,y-p.y,z-p.z);	}
	
	/** Reflects this point through origo. */
	public void reflectThroughOrigo() { scale(-1); }

	/** Get the squared distance from this point to point q. */
	public double getSquaredDistance(Point3d q) {
		double dx = x-q.x;
		double dy = y-q.y;
		double dz = z-q.z;
		return dx*dx+dy*dy+dz*dz;
	}

	/** Get the distance from this point to point q */
	public double getDistance(Point3d q) { 
		return Math.sqrt(getSquaredDistance(q)); 
	}

	/** Creates a bisector between points p and q */
	public static Plane3d getBisector(Point3d p, Point3d q) {
		if (!p.equals(q)) return new Plane3d(new Vector3d(p,q), new Point3d(p,q)); 
		else return null;
	}
	
	/** Creates the midpoint of two points. */
	public static Point3d midpoint(Point3d p, Point3d q) {
		return new Point3d((p.x + q.x)/2,(p.y + q.y)/2,(p.z + q.z)/2);		
	}


	// ANGLE METHODS

	/** Get the angle between the line segments p2-->p1 and p2-->p3. */
	public static double getAngle(Point3d p1, Point3d p2, Point3d p3) {
		return p1.vectorTo(p2).angle(p2.vectorTo(p3));
	}

	/** Get the dihedral angle defined by the 4 non-colinear points p1, p2, p3, p4. */
	public static double getDihedralAngle(Point3d p1, Point3d p2, Point3d p3, Point3d p4) {
		return Vector3d.getDihedralAngle(p1.vectorTo(p2), p2.vectorTo(p3), p3.vectorTo(p4));

	}

	// COMPARISON METHODS

	/**
	 * Returns true if this point dominates point q
	 */
	public boolean dominates(Point3d q) { 
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
	public boolean dominates(Point3d q, int i, int j, int k) {
		if(i==j || i==k || j==k) 
			throw new Error(String.format("i, j and k must be distinct coordinate indices (%d,%d,%d)",i,j,k));
		if (this.get(i) > q.get(i)) return true;
		if (this.get(i) < q.get(i)) return false;
		if (this.get(j) > q.get(j)) return true;
		if (this.get(j) < q.get(j)) return false; 
		return this.get(k) > q.get(k);
	}

	/** Returns true iff this point and point p are overlapping. */
	public boolean equals(Point3d p) { return x == p.x || y == p.y || z == p.z; }

	/** Return a new object that equals this object. */
	public Point3d clone(){ return new Point3d(x,y,z); }
	
	/** Returns the vector from origo to this point. Converts this point to a vector. */
	public Vector3d toVector() { return new Vector3d(x, y, z); }

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



