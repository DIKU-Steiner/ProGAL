package ProGAL.geom3d;

/** 
 *  A vector in (x,y,z)-space represented with double precision.
 *  @todo cache the length so several calls to getLengthSquared and getLength
 *  takes less time. This can be relevant e.g. in Line3d for multiple 
 *  projections, but can also cause serious problems if anyone chooses to extend the 
 *  Vector3d object.
 */
public class Vector extends ProGAL.geomNd.Vector{
//	protected double x,y,z;
		
	/** Construct a vector with the specified coordinates. */
	public Vector(double x, double y, double z) { 
		super(new double[]{x,y,z});
	}	
	/** Construct a vector pointing from origo to p. */
	public Vector(Point p) { 	super(p);	}
	
	/** Construct a vector that is a clone of v. */
	public Vector(Vector v) {	super(v);	}	
	
	/** Construct a vector using the double-array as coordinates. Note: The array is not cloned */
	public Vector(double[] coords) {	super(coords);	}
	
	
	/** Get the first coordinate. */
	public double getX() { return coords[0]; }
	
	/** Get the second coordinate. */
	public double getY() { return coords[1]; }
	
	/** Get the third coordinate. */
	public double getZ() { return coords[2]; }
	
	/** Set the first coordinate */
	public void setX(double x) { this.coords[0] = x; }
	/** Set the second coordinate */
	public void setY(double y) { this.coords[1] = y; }
	/** Set the third coordinate */
	public void setZ(double z) { this.coords[2] = z; }


	/** Get the dot-product of this vector and v. */
	public double dot(Vector v){ 
		return coords[0]*v.coords[0]+coords[1]*v.coords[1]+coords[2]*v.coords[2]; 
	}

	/** Get the angle between this vector and v. */
	public double angle(Vector v) {
		return Math.acos(Math.min(  1, this.dot(v)/Math.sqrt(this.getLengthSquared()*v.getLengthSquared())  ));
	}
	
	/** Add v to this vector and return the result (without changing this object). */
	public Vector add(Vector v){ return new Vector(coords[0]+v.coords[0],coords[1]+v.coords[1],coords[2]+v.coords[2]); }
	
	/** Add v to this vector and return the result (changing this object). */ 
	public Vector addThis(Vector v){ coords[0]+=v.coords[0]; coords[1]+=v.coords[1]; coords[2]+=v.coords[2]; return this; }
	
	/** Multiply this vector by s and return the result (without changing this object). */
	public Vector multiply(double s){ return new Vector(coords[0]*s, coords[1]*s, coords[2]*s); }
	
	/** Multiply this vector by s and return the result (changing this object). */
	public Vector multiplyThis(double s){ coords[0]*=s;coords[1]*=s;coords[2]*=s;return this; }
	
	/** Divide this vector by s and return the result (without changing this object). */
	public Vector divide(double s){ return new Vector(coords[0]/s, coords[1]/s, coords[2]/s); }
	
	/** Divide this vector by s and return the result (changing this object). */
	public Vector divideThis(double s){ coords[0]/=s;coords[1]/=s;coords[2]/=s;return this; }

	/** Normalize this vector and return the result (without changing this object). */
	public Vector normalize(){ return this.multiply(1/getLength()); }
	
	/** Normalize this vector and return the result (changing this object). */
	public Vector normalizeThis(){ return this.multiplyThis(1/getLength()); }

	/** Scale this vector to a certain length (returns new object and does not change this object). */
	public Vector scaleToLength(double length) { return multiply(length/getLength()); }
	
	/** Scale this vector to a certain length (changes this object). */
	public Vector scaleToLengthThis(double length) { return multiplyThis(length/getLength()); }
	
	/** Get the cross-product of this vector and v (without changing this object). */
	public Vector cross(Vector v){ return new Vector(coords[1]*v.coords[2] - coords[2]*v.coords[1], coords[2]*v.coords[0] - coords[0]*v.coords[2], coords[0]*v.coords[1] - coords[1]*v.coords[0]); }
	
	/** Get the cross-product of this vector and v and store the result in this vector (changes this object). */
	public Vector crossThis(Vector v){ 
		double newX = coords[1]*v.coords[2] - coords[2]*v.coords[1], newY = coords[2]*v.coords[0] - coords[0]*v.coords[2], newZ = coords[0]*v.coords[1] - coords[1]*v.coords[0];
		this.coords[0] = newX;this.coords[1] = newY;this.coords[2] = newZ;
		return this;
	}
	
	/** Convert this vector to a point. */
	public Point toPoint() { return new Point(coords[0], coords[1], coords[2]); }
	
	/** Returns a string-representation of this vector formatted with two decimals precision. */
	public String toString() { return toString(2); }
	
	/** Returns a string-representation of this vector formatted with <code>dec</code> decimals precision. */
	public String toString(int dec) {
		return String.format("Vector3d[%."+dec+"f,%."+dec+"f,%."+dec+"f]",coords[0],coords[1],coords[2]);
	}	
	
	/** Writes this vector to <code>System.out</code>. */
	public void toConsole() { toConsole(2); }
	
	/** Writes this vector to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec) { System.out.println(toString(dec)); }

	/** Create a clone of this vector. */
	public Vector clone(){ return new Vector(coords[0], coords[1], coords[2]); }

	///////// Static methods and fields ////////

	/** Get the angle between vector u and v. */
	public static double getAngle(Vector u, Vector v) { return u.angle(v); }
	
	/** Get the dihedral angle between 3 non-colinear vectors b1, b2, b3. */
	public static double getDihedralAngle(Vector b1, Vector b2, Vector b3) {
		Vector b2xb3 = b2.cross(b3);
		double y = b1.multiply(b2.getLength()).dot(b2xb3);
		double x = b1.cross(b2).dot(b2xb3);
		return Math.atan2(y,x);
	}
	

	/** An immutable vector pointing in the (1,0,0)-direction. */
	public static Vector X = new ImmutableVector3d(1,0,0);
	
	/** An immutable vector pointing in the (0,1,0)-direction. */
		public static Vector Y = new ImmutableVector3d(0,1,0);
		
	/** An immutable vector pointing in the (0,0,1)-direction. */
	public static Vector Z = new ImmutableVector3d(0,0,1);
	
	///////// Inner classes ////////
	
	/** 
	 * A wrapper class for <code>Vector3d</code> which makes the vector immutable. 
	 * All methods that can change the x, y and z-coordinates are overwritten. If 
	 * e.g. <code>setX(0.1)</code> is called a RuntimeException is thrown. If any of the 
	 * arithmetic methods such as <code>multiplyThis</code> are called, the result of 
	 * their corresponding non-mutating variant (<code>multiply</code>) is returned instead.
	 */
	public static class ImmutableVector3d extends Vector{
		public ImmutableVector3d(double x, double y, double z) { super(x,y,z);	}
		public void setCoord(int i, double v){throw new RuntimeException("This vector is immutable");}
		public void set(int i, double v){throw new RuntimeException("This vector is immutable");}
		public void setX(double x) { throw new RuntimeException("This vector is immutable"); }
		public void setY(double y) { throw new RuntimeException("This vector is immutable"); }
		public void setZ(double z) { throw new RuntimeException("This vector is immutable"); }
		public Vector addThis(Vector v){ return add(v); }
		public Vector multiplyThis(double s){ return multiply(s); }
		public Vector normalizeThis(){ return multiply(1/getLength()); }
		public Vector scaleToLengthThis(double length) { return multiply(length/getLength()); }
		public Vector crossThis(Vector v){ return cross(v); }
	}

	
	
}

