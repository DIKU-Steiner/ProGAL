package ProGAL.geom2d;


import java.awt.Color;

import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.geom2d.Line;
import ProGAL.math.Constants;
import ProGAL.math.Functions;

/** 
 *  A point in (x,y)-space represented using double precision. 
 * @author R.Fonseca
 */
public class Point extends ProGAL.geomNd.Point {
	private static final long serialVersionUID = 8095991200265432551L;
	public static Point origo = new Point(0.0, 0.0);

	
	/** Construct a point at (0,0) */
	public Point() { this(0,0); }

	/** Construct a point with the specified coordinates. */
	public Point(double x, double y) { super(new double[]{x,y}); }
	

	/** Construct a point with the specified coordinates. */
	public Point(double[] coords) { super(coords); }
	
	/** Construct a point with the same two first coordinates as p. */
	public Point(ProGAL.geomNd.Point p){
		this(p.getCoords());
	}

	/** Return the x-coordinate */
	public double x(){ return coords[0]; }
	
	/** Return the y-coordinate */
	public double y(){ return coords[1]; }
	
	public Point clone(){
		return new Point(coords[0], coords[1]);
	}
	
	/** Return a vector pointing from this point to <code>p</code> */
	public Vector vectorTo(Point p){
		return new Vector(p.coords[0]-coords[0], p.coords[1]-coords[1]);
	}

	/** Add the vector <code>v</code> to this point and return the result */
	public Point add(Vector v){ return new Point(x()+v.x(), y()+v.y()); }; 
	public Point add(double a, double b) { return new Point(x() + a, y() + b); }
	/** Add the vector <code>v</code> to this point and return the result (changes this object) */
	public Point addThis(Vector v){ coords[0] += v.x(); coords[1] += v.y(); return this; }
	public Point addThis(Point p){ coords[0] += p.x(); coords[1] += p.y(); return this; }
	public Point addThis(double a, double b) { coords[0] += a; coords[1] += b; return this; } 
	/** Subtract the vector <code>v</code> from this point and return the result */
	public Point subtract(Vector v){ return new Point(x()-v.x(), y()-v.y()); }
	/** Subtract the vector <code>v</code> from this point and return the result (changes this object) */
	public Point subtractThis(Vector v){ coords[0] -= v.x(); coords[1] -= v.y(); return this; }
	public Point subtractThis(Point p) { coords[0] -= p.x(); coords[1] -= p.y(); return this; }
	
	/** Returns the signed area of the triangle defined by the three points a, b and c (positive if counterclockwise) */
	public static double area(Point a, Point b, Point c) { 
		return a.x()*(b.y()-c.y()) + a.y()*(c.x()-b.x()) + b.x()*c.y() - c.x()*b.y(); 
	}
	
	/* Returns true if the three points a, b, and c are collinear */
	public static boolean collinear(Point a, Point b, Point c) { return area(a, b, c) < Constants.EPSILON; }
	
	
	/** returns squared distance of this point to the origo. */
	public double getSquaredDistance() { return x()*x() + y()*y(); }
	/** returns squared distance of this point to point q. */
	public double getSquaredDistance(Point q) { return Math.pow(x()-q.x(), 2) + Math.pow(y()-q.y(), 2); }

	/** returns squared distance of this point to segment s */
	public double getSquaredDistance(LineSegment s) {
		double squaredSegmentLength = s.getSquaredLength();
		if (squaredSegmentLength == 0.0) return getSquaredDistance(s.a);
		double sx = s.b.x() - s.a.x(); 
		double sy = s.b.y() - s.a.y();
		double t = ((x() - s.a.x()) * sx + (y() - s.a.y()) * sy)/squaredSegmentLength;
		if (t < 0) return getSquaredDistance(s.a);
		if (t > 1) return getSquaredDistance(s.b);
		return getSquaredDistance(new Point(s.a.x() + t*sx, s.a.y() + t*sy));	
	}

	/** returns distance of this point to segment s */
	public double getDistance(LineSegment s) { return Math.sqrt(getSquaredDistance(s)); }
	

	/** Returns the angle between the points, but if they make a left turn the angle will be negative. */
	public static double getSignedAngle(Point p1, Point p2, Point p3){
		Vector v1 = p2.vectorTo(p1);
		Vector v2 = p2.vectorTo(p3);
		return Math.atan2(v1.x()*v2.y()-v1.y()*v2.x(), v1.dot(v2));
	}
	
	/** Returns a positive double if the point is inside the circle through the 3 specified points (must be in clockwise order). */
	public static double inCircle(Point p, Point q, Point r, Point s) {
		double pp = p.x()*p.x() + p.y()*p.y();
		double qq = q.x()*q.x() + q.y()*q.y();
		double rr = r.x()*r.x() + r.y()*r.y();
		double ss = s.x()*s.x() + s.y()*s.y();
		
		return pp*(q.x()*(r.y()-s.y()) + q.y()*(s.x()-r.x()) + r.x()*s.y() - r.y()*s.x()) -
			   qq*(p.x()*(r.y()-s.y()) + p.y()*(s.x()-r.x()) + r.x()*s.y() - r.y()*s.x()) +
			   rr*(p.x()*(q.y()-s.y()) + p.y()*(s.x()-q.x()) + q.x()*s.y() - q.y()*s.x()) -
			   ss*(p.x()*(q.y()-r.y()) + p.y()*(r.x()-q.x()) + q.x()*r.y() - q.y()*r.x());
	}
	
	/** Returns polar angle of this point */
	public double polarAngle() {
		double angle = Math.acos(polarAngleCos());
		if ((coords[1] < 0) || ((coords[1] == 0.0) && (coords[0] < 0.0))) angle = Constants.TAU - angle;
		return angle;
	}
	
	/** Returns the sinus of the polar angle of this point */
	public double polarAngleSin() { return coords[1]/distance(); }

	/** Returns the cosinus of the polar angle of this point */
	public double polarAngleCos() { return coords[0]/distance(); }

	/** Returns true if points a, b and c make a left turn at b */
	public static boolean leftTurn(Point a, Point b, Point c ) { return area(a,b,c) > 0.0; }

	/** Returns true if points a, b and c make a right turn at b or are colinear */
	public static boolean rightTurn(Point a, Point b, Point c ) { return area(a,b,c) <= 0.0; }


	/** Returns the midpoint of two points. */
	public static Point midPoint(Point p, Point q) { return new Point((p.coords[0] + q.coords[0])/2, (p.coords[1] + q.coords[1])/2); }
	
	/** Creates a bisector line between points p and q */
	public static Line getBisector(Point p, Point q) {
		if (!p.equals(q)) return new Line(midPoint(p,q), p.vectorTo(q)); 
		return null;
	}
	
	/** Returns the line y = ax + b() where the coordinates of this point are (a,b) */
	public Line getDualLine() { return new Line(x(), y()); }


	
	/** Return true iff <code>o</code> is a Point that has the same same coordinates as this. */ 
	public boolean equals(Object o){
		if(o instanceof Point) return equals((Point)o);
		return false;
	}

	/** Return true iff <code>p</code> has the same same coordinates as this. */ 
	public boolean equals(Point p){
		return Math.abs(coords[0]-p.coords[0])<Constants.EPSILON && Math.abs(coords[1]-p.coords[1])<Constants.EPSILON;
	}

	
	/** rotates the point around the origo (counterclockwise) */
	public void rotation(double alpha) {
		double cosA = Math.cos(alpha);
		double sinA = Math.sin(alpha);
		double x = cosA*x() - sinA*y();
		set(1, cosA*y() + sinA*x());
		set(0, x);
	}
	
	public Point rotationClone(double alpha) {
		double cosA = Math.cos(alpha);
		double sinA = Math.sin(alpha);
		return new Point(cosA*x() - sinA*y(), cosA*y() + sinA*x());
	}
	
	public void rotation(double cos, double sin) {
		double xNew =  cos*x() - sin*y();
		set(1, cos*y() + sin*x());
		set(0, xNew);
		
	}
	/** rotates point p around this point by angle */
	public void rotation(Point p, double angle) {
		p.addThis(-x(), -y());
		p.rotation(angle);
		p.addThis(x(), y());
	}
	
	public void toScene(J2DScene scene, double rad, Color col) { scene.addShape(new Circle(this, rad), col); }
	public void toScene(J2DScene scene, double rad, Color col, double width) { scene.addShape(new Circle(this, rad), col, width); }

	public static void main(String[] args) {
		Line l1 = new Line(new LineSegment(new Point(-4,10), new Point(-3,6)));
		Line l2 = new Line(new LineSegment(new Point(-3,6),  new Point(0,4)));
		Line l3 = new Line(new LineSegment(new Point(0,4),   new Point(4,3)));
		Line l4 = new Line(new LineSegment(new Point(4,3),   new Point(7,5)));
		Line l5 = new Line(new LineSegment(new Point(7,5),   new Point(9,9)));
		Line l6 = new Line(new LineSegment(new Point(9,9),   new Point(10,14)));
		Line l7 = new Line(new LineSegment(new Point(-2,-5), new Point(3,-6)));
		Line l8 = new Line(new LineSegment(new Point(6,-7),  new Point(12,4)));

		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		l1.toScene(scene, 50, Color.black,  0.15);
		l3.toScene(scene, 50, Color.black,  0.15);
		l5.toScene(scene, 50, Color.black,  0.15);
		l7.toScene(scene, 50, Color.black,  0.15);
		l2.toScene(scene, 50, Color.black,  0.15);
		l4.toScene(scene, 50, Color.black,  0.15);
		l6.toScene(scene, 50, Color.black,  0.15);
		l8.toScene(scene, 50, Color.black,  0.15);

		
	}

}

