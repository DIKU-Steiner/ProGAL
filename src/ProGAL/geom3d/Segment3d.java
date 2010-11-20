package ProGAL.geom3d;

/**
 * A line segment spanned by two points, a and b.  
 */
public class Segment3d {
	protected Point3d a, b;

	/** Constructs a segment between points a and b. */
	public Segment3d(Point3d a, Point3d b) {
		this.a = a;
		this.b = b;
	}
	
	/** Constructs a segment from a to a+v. */
	public Segment3d(Point3d a, Vector3d v) {
		this.a = a;
		this.b = a.add(v);
	}

	/** Get the first point spanning the segment. */
	public Point3d getA() { return a; } 

	/** Get the second point spanning the segment. */
	public Point3d getB() { return b; }
	
	/** Change the first point spanning the segment. */
	public void setA(Point3d a) { this.a = a; }
	
	/** Change the second point spanning the segment. */
	public void setB(Point3d b) { this.b = b; }
	
	/** Get the direction of the segment. */
	public Vector3d getAToB(){ return a.vectorTo(b); }
	
	/** Get the length of the segment. */
	public double getLength() { return Math.sqrt(getSquaredLength()); }
	
	/** Get the squared length of the segment. */
	public double getSquaredLength() { 
		double bax = b.x - a.x, bay = b.y - a.y, baz = b.z - a.z;
		return bax*bax + bay*bay + baz*baz;
	}
	
	/** Get the point on the segment closest to a given point q. This method always returns 
	 * a new object.*/
	public Point3d getClosestPoint(Point3d q) {
		Vector3d dir = a.vectorTo(b);;
		Vector3d aq = a.vectorTo(q);;
		double t = dir.dot(aq)/dir.getLengthSquared();
		t = Math.min(1, Math.max(0,t));
		return new Point3d(a.x + t*dir.x, a.y + t*dir.y, a.z + t*dir.z);
	}

	/** Gets the squared distance from q to the nearest point on this segment. */
	public double getSquaredDistance(Point3d q) { return q.getDistanceSquared(getClosestPoint(q)); }
	
	/** Gets the distance from q to the nearest point on this segment. */
	public double getDistance(Point3d q) { return q.getDistance(getClosestPoint(q)); }
	
	/** Gets the midpoint of the segment. */
	public Point3d getMidPoint() { 
		return new Point3d( a.x + (b.x-a.x)/2 , a.y + (b.y-a.y)/2 , a.z + (b.z-a.z)/2 ); 
	}

	/** Returns a string representation of this segments. */
	public String toString(){ return "Segment3d["+a+","+b+"]"; } 
	
	/** Returns a string representation of this segments with <code>dec</code> decimals precision. */
	public String toString(int dec){ return "Segment3d["+a.toString(dec)+","+b.toString(dec)+"]"; } 
	
	/** Writes this segment to <code>System.out</code>. */
	public void toConsole(){ System.out.println(toString()); }
	
	/** Writes this segment to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec){ System.out.println(toString(dec)); }
}
