package ProGAL.geom3d;

/** 
 *  A vector in (x,y,z)-space represented with double precision.
 *  @todo cache the length so several calls to getLengthSquared and getLength
 *  takes less time. This can be relevant e.g. in Line3d for multiple 
 *  projections, but can also cause serious problems if anyone chooses to extend the 
 *  Vector3d object.
 */
public class Vector3d {
	protected double x,y,z;
		
	/** Construct a vector with the specified coordinates. */
	public Vector3d(double x, double y, double z) { this.x = x; this.y = y; this.z = z; }
	
	/** Construct a vector pointing from origo to p. */
	public Vector3d(Point3d p) { x = p.x; y = p.y; z = p.z; }
	
	/** Construct a vector that is a clone of v. */
	public Vector3d(Vector3d v) { x = v.x; y = v.y; z = v.z; }	
	
	//////// Getters and setters /////////
	
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


	/** Get the squared length of this vector. */
	public double getLengthSquared() { return x*x + y*y + z*z; }

	/** Get the length of this vector. */
	public double getLength() {return Math.sqrt(getLengthSquared()); }

	/** Get the squared length of this vector. */
	public double lengthSquared() {return getLengthSquared(); }
	
	/** Get the length of this vector. */
	public double length() {return getLength(); }

	
	
	/** Return true if the length of this vector is 0. */
	public boolean isZeroVector() { 
		return ((x == 0.0) && (y == 0.0) && (z == 0.0)); 
	}
	
	/** Get the dot-product of this vector and v. */
	public double dot(Vector3d v){ 
		return x*v.x+y*v.y+z*v.z; 
	}

	/** Get the angle between this vector and v. */
	public double angle(Vector3d v) {
		return Math.acos(Math.min(  1, this.dot(v)/(this.getLength()*v.getLength())  ));
	}
	
	/** Add v to this vector and return the result (without changing this object). */
	public Vector3d add(Vector3d v){ return new Vector3d(x+v.x, y+v.y, z+v.z); }
	
	/** Add v to this vector and return the result (changing this object). */ 
	public Vector3d addThis(Vector3d v){ x+=v.x; y+=v.y; z+=v.z; return this; }
	
	/** Multiply this vector by s and return the result (without changing this object). */
	public Vector3d multiply(double s){ return new Vector3d(x*s, y*s, z*s); }
	
	/** Multiply this vector by s and return the result (changing this object). */
	public Vector3d multiplyThis(double s){ x*=s;y*=s;z*=s;return this; }
	
	/** Normalize this vector and return the result (without changing this object). */
	public Vector3d normalize(){ return this.multiply(1/getLength()); }
	
	/** Normalize this vector and return the result (changing this object). */
	public Vector3d normalizeThis(){ return this.multiplyThis(1/getLength()); }

	/** Scale this vector to a certain length (returns new object and does not change this object). */
	public Vector3d scaleToLength(double length) { return multiply(length/getLength()); }
	
	/** Scale this vector to a certain length (changes this object). */
	public Vector3d scaleToLengthThis(double length) { return multiplyThis(length/getLength()); }
	
	/** Get the cross-product of this vector and v (without changing this object). */
	public Vector3d cross(Vector3d v){ return new Vector3d(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x); }
	
	/** Get the cross-product of this vector and v and store the result in this vector (changes this object). */
	public Vector3d crossThis(Vector3d v){ 
		double newX = y*v.z - z*v.y, newY = z*v.x - x*v.z, newZ = x*v.y - y*v.x;
		this.x = newX;this.y = newY;this.z = newZ;
		return this;
	}
	
	/** Convert this vector to a point. */
	public Point3d toPoint() { return new Point3d(x, y, z); }
	
	/** Returns a string-representation of this vector formatted with two decimals precision. */
	public String toString() { return toString(2); }
	
	/** Returns a string-representation of this vector formatted with <code>dec</code> decimals precision. */
	public String toString(int dec) {
		return String.format("Vector3d[%."+dec+"f,%."+dec+"f,%."+dec+"f]",x,y,z);
	}	
	
	/** Writes this point to <code>System.out</code>. */
	public void toConsole() { toConsole(2); }
	
	/** Writes this point to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec) { System.out.println(toString(dec)); }

	/** Return true if this vector equals v. */
	public boolean equals(Vector3d v){ return x==v.x && y==v.y && z==v.z; }
	
	/** Create a clone of this vector. */
	public Vector3d clone(){ return new Vector3d(x, y, z); }

	///////// Static methods and fields ////////

	/** Get the angle between vector u and v. */
	public static double getAngle(Vector3d u, Vector3d v) { return u.angle(v); }
	
	/** Get the dihedral angle between 3 non-colinear vectors v12, v23, v34. */
	public static double getDihedralAngle(Vector3d v12, Vector3d v23, Vector3d v34) {
		return Math.atan2( v12.scaleToLength(v23.getLength()).dot(v23.cross(v34)),
				v12.cross(v23).dot(v23.cross(v34)));
	}
	

	/** An immutable vector pointing in the (1,0,0)-direction. */
	public static Vector3d X = new ImmutableVector3d(1,0,0);
	
	/** An immutable vector pointing in the (0,1,0)-direction. */
		public static Vector3d Y = new ImmutableVector3d(0,1,0);
		
	/** An immutable vector pointing in the (0,0,1)-direction. */
	public static Vector3d Z = new ImmutableVector3d(0,0,1);
	
	///////// Inner classes ////////
	
	/** 
	 * A wrapper class for <code>Vector3d</code> which makes the vector immutable. 
	 * All methods that can change the x, y and z-coordinates are overwritten. If 
	 * e.g. <code>setX(0.1)</code> is called a RuntimeException is thrown. If any of the 
	 * arithmetic methods such as <code>multiplyThis</code> are called, the result of 
	 * their corresponding non-mutating variant (<code>multiply</code>) is returned instead.
	 */
	public static class ImmutableVector3d extends Vector3d{
		public ImmutableVector3d(double x, double y, double z) { super(x,y,z);	}
		public void setCoord(int i, double v){throw new RuntimeException("This vector is immutable");}
		public void set(int i, double v){throw new RuntimeException("This vector is immutable");}
		public void setX(double x) { throw new RuntimeException("This vector is immutable"); }
		public void setY(double y) { throw new RuntimeException("This vector is immutable"); }
		public void setZ(double z) { throw new RuntimeException("This vector is immutable"); }
		public Vector3d addThis(Vector3d v){ return add(v); }
		public Vector3d multiplyThis(double s){ return multiply(s); }
		public Vector3d normalizeThis(){ return multiply(1/getLength()); }
		public Vector3d scaleToLengthThis(double length) { return multiply(length/getLength()); }
		public Vector3d crossThis(Vector3d v){ return cross(v); }
	}
	
}

