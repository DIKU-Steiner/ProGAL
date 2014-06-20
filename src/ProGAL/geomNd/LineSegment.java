package ProGAL.geomNd;

public class LineSegment {
	protected Point a, b;

	/** Constructs a segment between points a and b. */
	public LineSegment(Point a, Point b) {
		this.a = a;
		this.b = b;
	}
	
	/** Constructs a segment from a to a+v. */
	public LineSegment(Point a, Vector v) {
		this.a = a;
		this.b = a.add(v);
	}

	/** Get the first point spanning the segment. */
	public Point getA() { return a; } 

	/** Get the second point spanning the segment. */
	public Point getB() { return b; }
	
	/** Change the first point spanning the segment. */
	public void setA(Point a) { this.a = a; }
	
	/** Change the second point spanning the segment. */
	public void setB(Point b) { this.b = b; }
	
	/** Get the direction of the segment. This method returns a new object. */
	public Vector getAToB(){ return a.vectorTo(b); }
	
	/** Get the length of the segment. */
	public double getLength() { return Math.sqrt(getLengthSquared()); }
	
	/** Get the squared length of the segment. */
	public double getLengthSquared() { 
		return a.distanceSquared(b);
	}
	
	/** Get the point on the segment closest to a given point q. This method always returns 
	 * a new object.*/
	public Point getClosestPoint(Point q) {
		Vector dir = a.vectorTo(b);
		Vector aq = a.vectorTo(q);
		double t = dir.dot(aq)/dir.getLengthSquared();
		t = Math.min(1, Math.max(0,t));
		return a.add(dir.multiplyThis(t));
	}

	/** Gets the squared distance from q to the nearest point on this segment. */
	public double getSquaredDistance(Point q) { return q.distanceSquared(getClosestPoint(q)); }
	
	/** Gets the distance from q to the nearest point on this segment. */
	public double getDistance(Point q) { return q.distance(getClosestPoint(q)); }
	
	/** Gets the midpoint of the segment. */
	public Point getMidPoint() { 
		return Point.getMidpoint(a, b);
	}

	/**
	 * Returns the mid-point of this line-segment. Since a line-segment can be interpreted as 
	 * a geometric shape (a 1-simplex) the Shape interface requires the getCenter method to be implemented.
	 * TODO: Test
	 */
	public Point getCenter() {
		return getMidPoint();
	}
	
	/** TODO: Comment and test */
	public Point getPoint(int i) {
		if(i<0 || i>1) throw new IllegalArgumentException("Invalid index ("+i+") 1-simplex has two points only");
		if(i==0) return a;
		else return b;
	}
	

	/** Return the 'dimension' of this object. Required by the interface Simplex. */
	public int getDimension() { return 1; }
	
	/** Returns true iff the argument is a line-segment and equals this. */
	public boolean equals(Object o){
		if(o instanceof LineSegment) return equals((LineSegment)o);
		return false;
	}
	
	/** Returns true iff this line-segment and ls are the same. Two line-segments are not 
	 * considered equal if they have the same end-points but in different orders. */
	public boolean equals(LineSegment ls){
		return a.equals(ls.a)&& b.equals(ls.b);
	}
	
	/** Returns a deep clone of this line segment. */
	public LineSegment clone(){
		return new LineSegment(a.clone(),b.clone());
	}
	
	/** Returns a string representation of this segments. */
	public String toString(){ return "Segment["+a+","+b+"]"; } 
	
	/** Returns a string representation of this segments with <code>dec</code> decimals precision. */
	public String toString(int dec){ return "Segment["+a.toString(dec)+","+b.toString(dec)+"]"; } 
	
	/** Writes this segment to <code>System.out</code>. */
	public void toConsole(){ System.out.println(toString()); }
	
	/** Writes this segment to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec){ System.out.println(toString(dec)); }

}
