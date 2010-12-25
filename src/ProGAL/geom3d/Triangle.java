package ProGAL.geom3d;

/**
 * A triangle in (x,y,z)-space represented by the three corner-points.
 */
public class Triangle implements Simplex{
	protected Point p1, p2, p3;

	/** Construct a triangle using the three specified points as corners */
	public Triangle(Point p1, Point p2, Point p3) {
		this.p1 = p1;
		this.p2 = p2;
		this.p3 = p3;
	}

	/** Get the first corner */
	public Point getP1(){ return p1; }
	/** Get the second corner */
	public Point getP2(){ return p2; }
	/** Get the third corner */
	public Point getP3(){ return p3; }
	/** Get the specified corner of this triangle */
	public Point getCorner(int c){ return getPoint(c); }
	/** Get the specified corner-point of this triangle */
	public Point getPoint(int c){
		switch(c){
		case 0: return p1;
		case 1: return p2;
		case 2: return p3;
		}
		throw new Error("Badly specified point number ("+c+"). Should be between 0 and 2");
	}
	
	/** Return the center of the triangle. Here average of the corners is used.*/
	public Point getCenter() { return new Point( (p1.x+p2.x+p3.x)/3, (p1.y+p2.y+p3.y)/3, (p1.z+p2.z+p3.z)/3); }

	/** Return the area of one side of the triangle. */
	public double getArea(){
		return 0.5*p1.vectorTo(p2).crossThis(p1.vectorTo(p3)).length();
	}
	
	/** Return a vector that is normal to this triangle. */
	public Vector getNormal() {
		return p1.vectorTo(p2).crossThis(p1.vectorTo(p3)).normalizeThis();
	}
	
	/** 
	 * Return the circumradius of the triangle. If one side has zero length this method returns 
	 * the length of the two remaining sides.
	 */
	public double getCircumradius(){
		double a = p1.getDistance(p2);
		double b = p1.getDistance(p3);
		double c = p2.getDistance(p3);
		double s = (a+b+c)/2;//Semiperemiter
		return a*b*c/(4*Math.sqrt(s*(a+b-s)*(a+c-s)*(b+c-s)));
	}

	/** Return the circumcenter of the triangle. 
	 * TODO: Test
	 * TODO: Make more efficient (transform to origo with n as z and use 2D formula)	 
	 */
	public Point getCircumcenter() {
		Vector n = getNormal();
		Point m1 = Point.getMidpoint(p1, p2);
		Point m2 = Point.getMidpoint(p1, p3);
		Line l1 = new Line(m1, p1.vectorTo(p2).crossThis(n));
		Line l2 = new Line(m2, p1.vectorTo(p3).crossThis(n));
		return l1.getIntersection(l2);
	}

	/** Returns a string-representation of this triangle formatted with two decimals precision. */
	public String toString(){	return toString(2);	}

	/** Returns a string-representation of this triangle formatted with <code>dec</code> decimals precision. */
	public String toString(int dec) {
		return String.format("Triangle[p1=%s,p2=%s,p3=%s]",p1.toString(dec), p2.toString(dec), p3.toString(dec));
	}

	/** Writes this triangle to <code>System.out</code> with 2 decimals precision. */
	public void toConsole() { toConsole(2); }

	/** Writes this triangle to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec) {
		System.out.println(toString(dec)); 
	}



 
}

