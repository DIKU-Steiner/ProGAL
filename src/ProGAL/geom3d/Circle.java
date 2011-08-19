package ProGAL.geom3d;

/*
 * If the circle has center c, radius r, and unit-length normal vector n, compute unit-length vectors u and v so that {u,v,n} are mutually
 * orthogonal. The circle is then parameterized by P(t) = c + r*(cos(t)*u + sin(t)*v) for 0 <= t < 2*pi. 
 */

/** 
 * A circle in (x,y,z)-space represented by a center-point, a radius and a normal-vector.
 */
public class Circle implements Shape{
	private Point center;
	private double radius;
	private Vector normal;

	public Circle(Point center, double radius, Vector normalVector) {
		this.center = center;
		this.radius = radius;
		this.normal = normalVector;
	}
	
	/*
	 * Given three points p0, p1, p2 and a vector v, find the circle or just the radius of the circle through the 
	 * projection of the points onto the plane with v as its normal vector. Claim: radius is minimized when of the 
	 * plane goes through p0, p1, p2
	 */
	/*public Circle3d(Point3d p0, Point3d p1, Point3d p2, Vector3d v) {
		// create the plane through the origo with v as its normal vector
		Plane3d plane = new Plane3d(v, new Point3d(0,0,0));
		Point3d q0 = plane.projectPoint(p0);
		Point3d q1 = plane.projectPoint(p1);
		Point3d q2 = plane.projectPoint(p2);
		double a = q0.getSquaredDistance(q1);
		double b = q1.getSquaredDistance(q2);
		double c = q2.getSquaredDistance(q0);
		radius = Math.sqrt(a*b*c)/Math.sqrt(2*a*b + 2*b*c + 2*c*a - a*a - b*b - c*c);
		normalVector = Vector3d.crossProduct(new Vector3d(p0,p1), new Vector3d(p0,p2));
	}*/
	

	/*
	 * returns the radius of the circle through 3 given points (without creating the circle)
	 */
	/*public static double getRadius(Point3d p0, Point3d p1, Point3d p2) {

		// get the plane through q0, q1, q2
		
		Point3d q0 = new Point3d(0,0,0);
		Point3d q1 = new Point3d(p1.x-p0.x, p1.y-p0.y, p1.z-p0.z);
		Point3d q2 = new Point3d(p2.x-p0.x, p2.y-p0.y, p2.z-p0.z);
		
		// get the plane through q0, q1, q2
		
		Plane3d plane = new Plane3d(q0,q1,q2);
		Vector3d verticalVector = new Vector3d(0,0,1);
		double angle = Vector3d.getAngle(plane.n, verticalVector);
		Vector3d rotationVector = Vector3d.crossProduct(plane.n, verticalVector);
		rotationVector.normalizeThis();
		q1.rotation(rotationVector, angle);
		q2.rotation(rotationVector, angle);
		Point2d r0 = new Point2d(0,0);
		Point2d r1 = new Point2d(q1.x, q1.y);
		Point2d r2 = new Point2d(q2.x, q2.y);
		Circle2d circle2 = new Circle2d(r0, r1, r2);
		return circle2.getRadius();
	}*/
	
	/** Get the center of the circle. */
	public Point getCenter() { return center; }
	/** Get the radius of the circle. */
	public double getRadius() { return radius; }
	/** Get the normal of the circle. */
	public Vector getNormalVector() { return normal; }

	/** Create the equilateral circle of two points. */
	public static Circle getEquilateralCircle(Point a, Point b) {
		Point center = Point.getMidpoint(a, b);
		Vector ab = a.vectorTo(b);
		double radius = Math.sqrt(3)*ab.length()/2;
		return new Circle(center, radius, ab);
	}
	

	/** Returns a string-representation of this circle formatted with two decimals precision. */
	public String toString(){
		return toString(2);
	}
	
	/** Returns a string-representation of this circle formatted with <code>dec</code> decimals precision. */
	public String toString(int dec) {
		return String.format("Circle[center=%s,radius=%"+dec+"f,normal=%s]",center.toString(dec),radius,normal.toString(dec));
	}

	/** Writes this circle to <code>System.out</code> with 2 decimals precision. */
	public void toConsole() { toConsole(2); }

	/** Writes this circle to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec) {
		System.out.println(toString(dec)); 
	}
 
}

