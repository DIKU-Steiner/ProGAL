package ProGAL.geom2d;

import ProGAL.math.Constants;

/** 
 *  A point in (x,y)-space represented using double precision. 
 * @author R.Fonseca
 */
public class Point extends ProGAL.geomNd.Point {
	private static final long serialVersionUID = 8095991200265432551L;
	
	/** Construct a point at (0,0) */
	public Point() { this(0,0); }
	
	/** Construct a point with the specified coordinates. */
	public Point(double x, double y) { super(new double[]{x,y}); }

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
	public Point add(Vector v){ return new Point(x()+v.x(), y()+v.y()); }
	/** Add the vector <code>v</code> to this point and return the result (changes this object) */
	public Point addThis(Vector v){ coords[0] += v.x(); coords[1] += v.y(); return this; }
	/** Subtract the vector <code>v</code> from this point and return the result */
	public Point subtract(Vector v){ return new Point(x()-v.x(), y()-v.y()); }
	/** Subtract the vector <code>v</code> from this point and return the result (changes this object) */
	public Point subtractThis(Vector v){ coords[0] -= v.x(); coords[1] -= v.y(); return this; }

	/** Returns the signed area of the triangle defined by the three points a, b and c (positive if counterclockwise) */
	public static double area(Point a, Point b, Point c) { 
		return a.x()*(b.y()-c.y()) + a.y()*(c.x()-b.x()) + b.x()*c.y() - c.x()*b.y(); 
	}
	
	/** Returns a positive double if the point is inside the circle through the 3 specified points (must be in clockwise order). */
	public static double inCircle(Point p, Point q, Point r, Point s) {
//		double aa = a.x()*a.x() + a.y()*a.y();
//		double bb = b.x()*b.x() + b.y()*b.y();
//		double cc = c.x()*c.x() + c.y()*c.y();
//		double dd = p.x()*p.x() + p.y()*p.y();
//		return a.x()*(b.y()-p.y())*cc + b.x()*(c.y()-a.y())*dd + c.x()*(p.y()-b.y())*aa + p.x()*(a.y()-c.y())*bb;
		
		
		double pos1 = p.x()*q.y()*(r.x()*r.x()+r.y()*r.y());
		double pos2 = p.y()*s.x()*(q.x()*q.x()+q.y()*q.y());
		double pos3 = r.x()*s.y()*(p.x()*p.x()+p.y()*p.y());
		double pos4 = q.x()*r.y()*(s.x()*s.x()+s.y()*s.y());
		
		double neg1 = p.x()*s.y()*(r.x()*r.x()+r.y()*r.y());
		double neg2 = p.y()*q.x()*(s.x()*s.x()+s.y()*s.y());
		double neg3 = q.y()*r.x()*(p.x()*p.x()+p.y()*p.y());
		double neg4 = r.y()*s.x()*(q.x()*q.x()+q.y()*q.y());
		
		return pos1+pos2+pos3+pos4-neg1-neg2-neg3-neg4;
	}
	
	/** Returns true if points a, b and c make a left turn at b */
	public static boolean leftTurn(Point a, Point b, Point c ) { return area(a,b,c) > 0.0; }

	/** Returns true if points a, b and c make a right turn at b or are colinear */
	public static boolean rightTurn(Point a, Point b, Point c ) { return area(a,b,c) <= 0.0; }


	/** Returns the midpoint of two points. */
	public static Point midPoint(Point p, Point q) { return new Point((p.coords[0] + q.coords[0])/2, (p.coords[1] + q.coords[1])/2); }
	
	/** Creates a bisector between points p and q */
	public static Line getBisector(Point p, Point q) {
		if (!p.equals(q)) 
			return new Line(midPoint(p,q), p.vectorTo(q)); 
		return null;
	}
	
	/** Return true iff <code>o</code> is a Point that has the same same coordinates as this. */ 
	public boolean equals(Object o){
		if(o instanceof Point) return equals((Point)o);
		return false;
	}

	/** Return true iff <code>p</code> has the same same coordinates as this. */ 
	public boolean equals(Point p){
		return Math.abs(coords[0]-p.coords[0])<Constants.EPSILON && Math.abs(coords[1]-p.coords[1])<Constants.EPSILON;
	}
	
	
}

