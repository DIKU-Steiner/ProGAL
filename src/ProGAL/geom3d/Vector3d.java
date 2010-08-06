package ProGAL.geom3d;

/** 
 *  A vector in (x,y,z)-space represented with double precision. 
 */
public class Vector3d {
	protected double x,y,z;
		
	public Vector3d() { x = 0.0; y = 0.0; z = 0.0; }
	public Vector3d(double x, double y, double z) { this.x = x; this.y = y; this.z = z; }
	public Vector3d(Point3d p, Point3d q) { x = q.x - p.x; y = q.y - p.y; z = q.z - p.z; }
	public Vector3d(Point3d p) { x = p.x; y = p.y; z = p.z; }
	public Vector3d(Segment3d s) { x = s.b.x-s.a.x; y = s.b.y-s.a.y; z = s.b.z - s.a.z; }
	public Vector3d(Vector3d v) { x = v.x; y = v.y; z = v.z; }	
	
	
	public double getX() { return x; }     public void setX(double x) { this.x = x; }
	public double getY() { return y; }     public void setY(double y) { this.y = y; }
	public double getZ() { return z; }     public void setZ(double z) { this.z = z; }
	
	public void set(Point3d p) { x = p.x; y = p.y; z = p.z; }
	
	public double getSquaredLength() { return x*x + y*y + z*z; }
		
	public double getLength() {return Math.sqrt(getSquaredLength()); }
	
	public static double getAngle(Vector3d u, Vector3d v) { return Math.acos(dotProduct(u,v)/(u.getLength()*v.getLength())); }
		
	/*
	 * returns dihedral angle between 3 non-colinear vectors v12, v23, v34.
	 */
	public static double getDihedralAngle(Vector3d v12, Vector3d v23, Vector3d v34) {
		return Math.atan2(dotProduct(v12.createScaledVector3d(v23.getLength()), crossProduct(v23,v34)),
				dotProduct(crossProduct(v12,v23),crossProduct(v23,v34)));
	}

	public boolean isZeroVector() { return ((x == 0.0) && (y == 0.0) && (z == 0.0)); }

	public void negative() { scale(-1); }
	
	/*
	 * scale by factor a
	 */
	//public void scale(double a) { x *= a; y *= a; z *= a; }
	
	/**
	 * scale to specified length
	 */
	public Vector3d scaleToLength(double length) { scale(length/getLength()); return this; }
	
	public void makeUnitVector() { scale(1/getLength()); }
	
	//public void normalize() { makeUnitVector(); }
	
	public void addVector(Vector3d v) { x += v.x; y += v.y; z += v.z; }
	
	public void subtractVector(Vector3d v) { x -= v.x; y -= v.y; z -= v.z; }

	public Vector3d add(Vector3d v){ return new Vector3d(x+v.x, y+v.y, z+v.z); }
	public Vector3d addIn(Vector3d v){ x+=v.x; y+=v.y; z+=v.z; return this; }
	public Vector3d multiply(double s){ return new Vector3d(x*s, y*s, z*s); }
	public Vector3d multiplyIn(double s){ x*=s;y*=s;z*=s;return this; }
	public Vector3d norm(){ return this.multiply(1/length()); }
	public Vector3d normIn(){ return this.multiplyIn(1/length()); }
	
	public Vector3d createNegativeVector3d() { return createScaledVector3d(-1); }
	
	/*
	 * create scaled vector
	 */
	public Vector3d createScaledVector3d(double a) { return new Vector3d(a*x, a*y, a*z); }
	public Vector3d mult(double a) { return new Vector3d(a*x, a*y, a*z); }
	
	/*
	 * creates the vector with the opposite direction
	 */
	public Vector3d createReverseVector3d() { return createScaledVector3d(-1); }
	
	/*
	 * create vector scaled to given length
	 */
	public Vector3d createScaledToLengthVector3d(double length) { return createScaledVector3d(length/getLength()); }
	
	public Vector3d createUnitVector() { return createScaledVector3d(1/getLength()); }
	
	public Point3d toPoint() { return new Point3d(x, y, z); }
	
	public String toString() {
		return toString(2);
		//return "(" + x + " , " + y + " , " + z + ")"; 
	}
	public String toString(int dec) {
		return String.format("Vector3d[%."+dec+"f,%."+dec+"f,%."+dec+"f]",x, y,z);
//		return "(" + Functions.toString(x,dec) + " , " + Functions.toString(y,dec) + " , " + Functions.toString(z,dec) + ")";
	}	
	
	public void toConsole() { System.out.println(toString()); }
	public void toConsole(int dec) { System.out.println(toString(dec)); }

	public static Vector3d createSum(Vector3d u, Vector3d v) { 
		return new Vector3d(u.x+v.x, u.y+v.y, u.z+v.z); }
	
	public static Vector3d createSum(Vector3d u, Vector3d v, Vector3d w) { 
		return new Vector3d(u.x+v.x+w.x, u.y+v.y+w.y, u.z+v.z+w.z); }
	
	public static Vector3d createDiff(Vector3d u, Vector3d v) { 
		return new Vector3d(u.x-v.x, u.y-v.y, u.z-v.z); }
	
	public static Vector3d createUnitVector(Point3d p, Point3d q) {
		double dist = Point3d.getDistance(p, q);
		return new Vector3d((q.getX()-p.getX())/dist, (q.getY()-p.getY())/dist, (q.getZ()-p.getZ())/dist);
	}
	
	public static Vector3d crossProduct(Vector3d u, Vector3d v) {
		return new Vector3d(u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y - u.y*v.x); }
	
	public Vector3d cross(Vector3d v) {
		return crossProduct(this,v);
	}
	
	public static double dotProduct(Vector3d u, Vector3d v) { return u.x*v.x + u.y*v.y + u.z*v.z; }
	public        double dotProduct(Vector3d v)             { return x*v.x   + y*v.y   + z*v.z;   }
	
	public static double dotProduct(Vector3d u, Point3d p)  { return Vector3d.dotProduct(u,new Vector3d(p)); }
	
	public static double dotProduct(Point3d  p, Vector3d u) { return Vector3d.dotProduct(u, p); }
	
	public static double dotProduct(Point3d  p, Point3d q)  { 
		//Edited by Rasmus ... for performance
		//return Vector3d.dotProduct(new Vector3d(p), new Vector3d(q));
		return p.x*q.x + p.y*q.y + p.z*q.z;
	}
	
	public Vector3d clone(){ return new Vector3d(x, y, z); }
	
	public static void main(String[] args) {
	}
	public double angle(Vector3d v2) {
		return Math.acos(Math.min(  1, this.dot(v2)/(this.length()*v2.length())  ));
	}
	public double get(int i) {
		switch(i){
		case 0: return x;
		case 1: return y;
		case 2: return z;
		}
		throw new Error("Illegal coordinate referenced: "+i);
	}
	public Vector3d reverse() {
		return new Vector3d(-x,-y,-z);
	}
	

}

