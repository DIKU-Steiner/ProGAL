package ProGAL.geom2d;

import ProGAL.math.Constants;

public class Point {
	protected double x,y;
	
	public Point() { this(0,0); }
	
	public Point(double x, double y) { this.x = x; this.y = y; }
	

	public Point clone(){
		return new Point(x,y);
	}
	
	// GET METHODS
	
	public double getX() { return x; }
	public double getY() { return y; }
	
	public Vector vectorTo(Point p){
		return new Vector(p.x-x, p.y-y);
	}
	
	/** Returns the midpoint of two points. */
	public static Point midPoint(Point p, Point q) { return new Point((p.x + q.x)/2, (p.y + q.y)/2); }
	
//	/*
//	 * creates the equilateral point to the right of pq
//	 */
//	public Point2d getEquilateralPoint(Point2d p, Point2d q) {
//		Point2d e = q.clone();
//		p.rotation(e, Math.PI/3.0);
//		return e;
//	}
//	public static Point2d createEquilateralPoint(Point2d p, Point2d q) {
//		Point2d e = new Point2d(q);
//		p.rotation(e, Math.PI/3.0);
//		return e;
//	}


//	/*
//	 * rotates the point around the origo (counterclockwise)
//	 */
//	public void rotation(double alpha) {
//		double cosA = Math.cos(alpha);
//		double sinA = Math.sin(alpha);
//		double xNew =  cosA*x - sinA*y;
//		y = cosA*y + sinA*x;
//		x = xNew;
//	}
//	
//	/*
//	 * rotates point p around this point by angle
//	 */
//	public void rotation(Point2d p, double angle) {
//		p.translate(-x, -y);
//		p.rotation(angle);
//		p.translate(x,y);
//	}
	
	/** Scales this point by the factor f. */
	public Point scale(double f) { x = x*f; y = y*f; return this; }
	
	/** Returns squared distance of this point to point q. */
	public double squaredDistance(Point q) { 
		double dx = x-q.x;
		double dy = y-q.y;
		return dx*dx+dy*dy; 
	}
	
	/** Returns distance of this point to point q. */
	public double distance(Point q) {
		double dx = x-q.x;
		double dy = y-q.y;
		return Math.sqrt(dx*dx+dy*dy); 
	}
	
//	/*
//	 * returns squared distance of this point to segment s 
//	 */
//	public double getSquaredDistance(Segment2d s) {
//		double squaredSegmentLength = s.getSquaredLength();
//		if (squaredSegmentLength == 0.0) return getSquaredDistance(s.a);
//		double sx = s.b.x - s.a.x, sy = s.b.y - s.a.y;
//		double t = ((s.a.x - x) * sx + (s.a.y - y) * sy)/squaredSegmentLength;
//		if (t < 0) return getSquaredDistance(s.a);
//		if (t > 1) return getSquaredDistance(s.b);
//		return getSquaredDistance(new Point2d(s.a.x + t*sx, s.a.y + t*sy));	
//	}
//	
	/*
	 * creates a bisector between points p and q
	 */
	public static Line getBisector(Point p, Point q) {
		if (!p.equals(q)) return new Line(midPoint(p,q), new Vector(p,q)); 
		return null;
	}
	
	
//	/*
//	 * returns true if points a, b and c make a left turn at b
//	 */
//	public static boolean leftTurn(Point2d a, Point2d b, Point2d c ) { return area2(a,b,c) > 0.0; }
//
//	/*
//	 * returns true if points a, b and c make a right turn at b or are colinear
//	 */
//	public static boolean rightTurn(Point2d a, Point2d b, Point2d c ) { return area2(a,b,c) < 0.0; }
//
//	/*
//	 * returns true if points a, b and c are colinear
//	 */
//	public static boolean colinear(Point2d a, Point2d b, Point2d c ) { return area2(a,b,c) == 0.0; }
//
//	/*
//	 * returns true if this point dominates point q
//	 */
//	public boolean dominates(Point2d q) { return   x > q.x || ( x == q.x && y > q.y); }
//	
//	/*
//	 * returns true if this point dominates point q (i=0,1 is the most important coordinate,
//	 * and j=0,1 is the least important coordinate).
//	 */
//	public boolean dominates(Point2d q, int i, int j) { 
//		double[] pc = new double[2], qc = new double[2];
//		pc[i] = x;   pc[j] = y;  
//		qc[i] = q.x; qc[j] = q.y;
//		if (pc[0] > qc[0]) return true;
//		if (pc[0] < qc[0]) return false; else return pc[1] > qc[1];
//	}
//
//	/*
//	 * returns true if this point yx-dominates point q
//	 */
//	public boolean yxDominates(Point2d q) { return   y > q.y || ( y == q.y && x > q.x); }
//	
//	/*
//	 * returns true if this point and point p are overlapping
//	 */
//	public boolean sameCoordinates(Point2d p) { return Math.abs(x-p.x)<0.000001 && Math.abs(y-p.y)<0.000001; }
//
//	public boolean isBelow(Point2d p) { return y < p.y; }
//
//	/*
//	 * returns true if the point is inside the circle through specified 3 points
//	 */
//	public boolean inCircle(Point2d a, Point2d b, Point2d c) {
//		double aa = a.x*a.x + a.y*a.y, bb = b.x*b.x + b.y*b.y, cc = c.x*c.x + c.y*c.y, dd = x*x + y*y;
//		return a.x*(b.y-y)*cc + b.x*(c.y-a.y)*dd + c.x*(y-b.y)*aa + x*(a.y-c.y)*bb > 0.0;
//	}
//	
//	/*
//	 * returns true if the point is inside (counterclockwise) triangle with the specified corners
//	 */
//	public boolean inTriangle(Point2d a, Point2d b, Point2d c) {
//		return (rightTurn(this,a,b) && rightTurn(this, b, c) && rightTurn(this, c, a)) ||
//		       (leftTurn(this,a,b) && leftTurn(this, b,c) && leftTurn(this, c, a));
//	}
	
	public boolean equals(Object o){
		if(o instanceof Point) return equals((Point)o);
		return false;
	}
	
	public boolean equals(Point p){
		return Math.abs(x-p.x)<Constants.EPSILON && Math.abs(y-p.y)<Constants.EPSILON;
	}
	
	// TRANSFORMATION METHODS

	
//	public Vector2d toVector() { return new Vector2d(x, y); }

	public String toString() { return toString(2); }
	
	public String toString(int dec) { 
		return String.format("Point[%"+dec+"f,%"+dec+"f]",x,y);
	}
	
	public void toConsole() { toConsole(2); }
	public void toConsole(int dec) { System.out.println(toString(dec)); }
	

//	public Line2d bisector(Point2d p) {
//		return Point2d.getBisector(this, p);
//	}

}

//