package ProGAL.geom3d;

import static java.lang.Math.cos;
import static java.lang.Math.sin;

import java.awt.Color;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.math.Constants;
import ProGAL.math.Matrix;

/** 
 *  A vector in (x,y,z)-space represented with double precision.
 *  @todo cache the length so several calls to getLengthSquared and getLength
 *  takes less time. This can be relevant e.g. in Line3d for multiple 
 *  projections, but can also cause serious problems if anyone chooses to extend the 
 *  Vector3d object.
 */
public class Vector extends ProGAL.geomNd.Vector{
		
	/** Construct a vector with the specified coordinates. */
	public Vector(double x, double y, double z) { 
		super(new double[]{x,y,z});
	}	
	/** Construct a vector pointing from origo to p. */
	public Vector(Point p) { 	super(p);	}
	
	/** Constructs a vector between two points p1 and p2 - added by pawel 12-11-2011 */
	public Vector(Point p1, Point p2) { super(p1, p2); }

	/** Construct a vector that is a clone of v. */
	public Vector(Vector v) {	this(v.get(0), v.get(1), v.get(2)); }	
	
	/** Construct a vector that is a clone of v. */
	public Vector(ProGAL.geomNd.Vector v) {	this(v.get(0), v.get(1), v.get(2)); }	
	
	/** Construct a vector using the double-array as coordinates. Note: The array is not cloned */
	public Vector(double[] coords) {	super(coords);	}
	
	
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

	/** Set all coordinates of this vector equal to those of v */
	public void set(Vector v) {
		this.coords[0] = v.coords[0];
		this.coords[1] = v.coords[1];
		this.coords[2] = v.coords[2];
	}

	/** Get the dot-product of this vector and v. */
	public double dot(Vector v){ 
		return coords[0]*v.coords[0]+coords[1]*v.coords[1]+coords[2]*v.coords[2]; 
	}
	/** Get the dot-product of this vector and position vector to point p. */
	public double dot(Point p){ 
		return coords[0]*p.x()+coords[1]*p.y()+coords[2]*p.z(); 
	}

	/** Get the angle between this vector and v. */
	public double angle(Vector v) {
		return Math.acos(dot(v)/(length()*v.length()));
//		return Math.acos(Math.min(  1, this.dot(v)/Math.sqrt(this.getLengthSquared()*v.getLengthSquared())  ));
	}

	/** Add v to this vector and return the result (without changing this object). */
	public Vector add(Vector v){ return new Vector(coords[0]+v.coords[0],coords[1]+v.coords[1],coords[2]+v.coords[2]); }
	
	/** Add v to this vector and return the result (changing this object). */ 
	public Vector addThis(Vector v){ coords[0]+=v.coords[0]; coords[1]+=v.coords[1]; coords[2]+=v.coords[2]; return this; }

	/** Add p to this vector and return the result (without changing this object). */
	public Vector add(Point p){ return new Vector(coords[0]+p.x(),coords[1]+p.y(),coords[2]+p.z()); }
	
	/** Add p to this vector and return the result (changing this object). */ 
	public Vector addThis(Point p){ coords[0]+=p.x(); coords[1]+=p.y(); coords[2]+=p.z(); return this; }

	/** Subtract v from this vector and return the result (without changing this object). */
	public Vector subtract(Vector v){ return new Vector(coords[0]-v.coords[0],coords[1]-v.coords[1],coords[2]-v.coords[2]); }
	
	/** Subract v from this vector and return the result (changing this object). */ 
	public Vector subtractThis(Vector v){ coords[0]-=v.coords[0]; coords[1]-=v.coords[1]; coords[2]-=v.coords[2]; return this; }
	
	/** Multiply this vector by s and return the result (without changing this object). */
	public Vector multiply(double s){ return new Vector(coords[0]*s, coords[1]*s, coords[2]*s); }
	
	/** Multiply this vector by s and return the result (changing this object). */
	public Vector multiplyThis(double s){ coords[0]*=s;coords[1]*=s;coords[2]*=s;return this; }
	
	//Daisy
	public Vector multiply(Vector v){ 
		double[] ret = new double[dim];
		for(int d=0;d<dim;d++) ret[d] = coords[d]*v.getCoord(d);
		return new Vector(ret); 
	}
	
	/** Divide this vector by s and return the result (without changing this object). */
	public Vector divide(double s){ return new Vector(coords[0]/s, coords[1]/s, coords[2]/s); }
	
	/** Divide this vector by s and return the result (changing this object). */
	public Vector divideThis(double s){ coords[0]/=s;coords[1]/=s;coords[2]/=s;return this; }

	/** Normalize this vector and return the result (without changing this object). */
	public Vector normalize(){ return this.divide(length()); }
	
	/** Normalize this vector and return the result (changing this object). */
	public Vector normalizeThis(){ return this.divideThis(length()); }

	public Vector normalizeFast(){ return this.multiply(invSqrt(getLengthSquared())); }
	public Vector normalizeThisFast(){ return this.multiplyThis(invSqrt(getLengthSquared())); }
	
	/** The fast inverse square root hack from quake 3 */
	private static double invSqrt(double x) {
	    double xhalf = 0.5d*x;
	    long i = Double.doubleToLongBits(x);
	    i = 0x5fe6ec85e7de30daL - (i>>1);
	    x = Double.longBitsToDouble(i);
	    x = x*(1.5d - xhalf*x*x);
	    return x;
	}
	
	
	/** Scale this vector to a certain length (returns new object and does not change this object). */
	public Vector scaleToLength(double length) { return multiply(length/length()); }
	
	/** Scale this vector to a certain length (changes this object). */
	public Vector scaleToLengthThis(double length) { return multiplyThis(length/length()); }
	
	/** Get the cross-product of this vector and v (without changing this object). */
	public Vector cross(Vector v){ 
		return new Vector(coords[1]*v.coords[2] - coords[2]*v.coords[1], 
						  coords[2]*v.coords[0] - coords[0]*v.coords[2], 
						  coords[0]*v.coords[1] - coords[1]*v.coords[0]); }
	
	/** Get the cross-product of this vector and v and store the result in this vector (changes this object). */
	public Vector crossThis(Vector v){ 
		double newX = coords[1]*v.coords[2] - coords[2]*v.coords[1], newY = coords[2]*v.coords[0] - coords[0]*v.coords[2], newZ = coords[0]*v.coords[1] - coords[1]*v.coords[0];
		this.coords[0] = newX;this.coords[1] = newY;this.coords[2] = newZ;
		return this;
	}
	
	/** rotates (clockwise) the vector around the line through the  origo with the direction unit vector v.
	 * For counterclockwise rotation change signs within parentheses in non-diagonal terms. */
	public void rotation(Vector v, double alpha) {
		double c = Math.cos(alpha);
		double d = 1.0-c;
		double s = Math.sin(alpha);
		double vxyd = v.x()*v.y()*d, vxzd = v.x()*v.z()*d, vyzd = v.y()*v.z()*d;
		double vxs = v.x()*s, vys = v.y()*s, vzs = v.z()*s; 
		double xNew = (v.x()*v.x()*d+c)*x() + (vxyd-vzs)*y()    + (vxzd+vys)*z();
		double yNew = (vxyd+vzs)*x()    + (v.y()*v.y()*d+c)*y() + (vyzd-vxs)*z();
		setZ((vxzd-vys)*x()    + (vyzd+vxs)*y()    + (v.z()*v.z()*d+c)*z());
		setX(xNew);
		setY(yNew);
	}

	
	/** 
	 * Perform a right-handed rotation of v around this vector. 
	 * TODO: Test
	 */
	public Vector rotateIn(Vector v, double angle) {
		double l = length();
		if(l==0) throw new Error("Trying to rotate around 0-vector");
		
		double ux = coords[0]/l;
		double uy = coords[1]/l;
		double uz = coords[2]/l;
		double sin = sin(angle);
		double cos = cos(angle);
		
		double a00 = (ux*ux + cos*(1.0-ux*ux));
		double a10 = (ux*uy*(1.0-cos)+uz*sin);
		double a20 = (uz*ux*(1-cos) - uy*sin);
        
		double a01 = (ux*uy*(1-cos) - uz*sin);
		double a11 = (uy*uy + cos*(1.0-uy*uy));
		double a21 = (uy*uz*(1.0-cos) + ux*sin);
        
		double a02 = (uz*ux*(1.0-cos) + uy*sin);
		double a12 = (uy*uz*(1.0-cos) - ux*sin);
		double a22 = (uz*uz + cos*(1.0 - uz*uz));

		double newX = a00*v.coords[0]+a01*v.coords[1]+a02*v.coords[2];
		double newY = a10*v.coords[0]+a11*v.coords[1]+a12*v.coords[2];
		double newZ = a20*v.coords[0]+a21*v.coords[1]+a22*v.coords[2];
		v.setX(newX);
		v.setY(newY);
		v.setZ(newZ);
		return v;
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

	/** Draws a vector starting at the origin */
	public LSS toScene(J3DScene scene, Color clr, double width) {
		LSS lss = new LSS(new Point(0,0,0), toPoint(), width);
		scene.addShape(lss, clr, 3);
		return lss;
	}
	/** Draws a vector starting at point p */
	public LSS toScene(J3DScene scene, Point p, Color clr, double width) {
		LSS lss = new LSS(p, p.add(this), width);
		scene.addShape(lss, clr, 3);
		return lss;
	}

	
	/** Create a clone of this vector. */
	public Vector clone(){ return new Vector(coords[0], coords[1], coords[2]); }

	///////// Static methods and fields ////////

	/** Get the angle between vector u and v. */
	public static double getAngle(Vector u, Vector v) { return u.angle(v); }
	
	/** Get the dihedral angle between 3 non-colinear vectors b1, b2, b3. */
	public static double getDihedralAngle(Vector b1, Vector b2, Vector b3) {
		Vector b2xb3 = b2.cross(b3);
		double y = b1.multiply(b2.length()).dot(b2xb3);
		double x = b1.cross(b2).dot(b2xb3);
		return Math.atan2(y,x);
	}
	
	/** Get a vector (one of many possible) orthonormal to vector v (v is assumed to be unit a vector). */
	public Vector getOrthonormal() {
		if (Math.abs(z()) > Constants.EPSILON) { 
			double ratio = x()/z();
			double length = Math.sqrt(1 + ratio*ratio);
			return new Vector(1/length, 0, -ratio/length);
		}
		else return new Vector(0, 0, 1);
	}

	/** Get a vector (one of many possible) orthogonal (i.e.,  having arbitrary length) to vector v. */
	public Vector getOrthogonal() {
		if (Math.abs(z()) > Constants.EPSILON) return new Vector(1, 0, -x()/z());
		else return new Vector(0, 0, 1); }

	
	/*
	 * added by pawel 12-11-2011
	 */
	public boolean isParallel(Vector v) { return cross(v).isZeroVector(); }
	
	public boolean isSteinerAngle(Vector v) { return dot(v)/(this.length()*v.length()) > -0.5; }

	
  	/*
  	 * returns cosinus of the dihedral angle between three vectors 
  	 * added by pawel on 12-11-2011
  	 */
  	public static double getCosDihedralAngle(Vector u, Vector v, Vector w) {
  		if (u.isParallel(v)) throw new Error("Vectors u and v are colinear");
  		if (v.isParallel(w)) throw new Error("Vectors v and w are colinear");
  		Vector uv = u.cross(u.add(v));
  		Vector vw = v.cross(v.add(w));
  		return (uv.dot(vw)/(uv.length()*vw.length()));
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
		public Vector normalizeThis(){ return multiply(1/length()); }
		public Vector scaleToLengthThis(double length) { return multiply(length/length()); }
		public Vector crossThis(Vector v){ return cross(v); }
	}
	public static void main(String[] args) {
		Vector a = new Vector(0,0,1);
		Vector b = new Vector(1,0,1);
		Point center = new Point(1,1,0);
		Circle C = new Circle(center, 0.2, b);
		Vector rotAxis = b.cross(a);
		double angle = b.angle(a);
//		double angle = b.dot(a.toPoint());
		Matrix rotMatrix = Matrix.createRotationMatrix(angle, rotAxis);
		Vector translate = new Vector(-center.x(),-center.y(),-center.z());
		Vector newVector = rotMatrix.multiply(b);
		Circle newC = new Circle(center.add(translate), 0.2, newVector.normalize());
		System.out.println("newC = "+newC.toString());
		Vector n = new Vector(0, 0, 1);
		Vector u = new Vector(0.3, 0.4, 0.0).normalize(); 
		Vector v = n.cross(u);
		System.out.println(u.x() + " " + v.y());
		System.out.println("cross = "+n.cross(n));
	}

	
}

