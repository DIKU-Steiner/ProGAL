package ProGAL.geom3d;

/**
 * A plane represented by a point and a normal. 
 * 
 * Assuming that <i>p</i> is a point on the plane and <i>n</i> is the normal vector, 
 * the half-space <i>{q|pq·n>0}</i> is called the 'upper halfspace' wrt. the plane 
 * and vice versa for the 'lower halfspace'.
 */
public class Plane {
	/** Normal vector of the plane. */
	protected Vector normal; 
	/** Point in the plane. */
	protected Point point;  
	/** Perpendicular distance from the origin. */
	protected double d; 

	/** Constructs a plane with the normal vector n containing point p. */
	public Plane(Point p, Vector n) {
		this.normal = n;
		this.point = p;
		d = -n.x*p.x - n.y*p.y - n.z*p.z;
	}

	/** Constructs a plane through three distinct points. */
	public Plane(Point p, Point q, Point r) {
		normal = p.vectorTo(q).crossThis(p.vectorTo(r)).normalizeThis();
		this.point = p;
		d = -normal.x*p.x - normal.y*p.y - normal.z*p.z;
	}
	
	/** Get the point defining this plane. */
	public Point getPoint(){	return point;	}
	
	/** Return the normal defining this plane. */
	public Vector getNormal(){ 	return normal; 	}

	/** Set the normal to n. */
	public void setNormal(Vector n) { 
		this.normal = n; 
	}

	/** Get the projection of p onto this plane. */
	public Point projectPoint(Point p) {
		double t = normal.x*p.x + normal.y*p.y + normal.z*p.z +d;
		return new Point(p.x - normal.x*t, p.y - normal.y*t, p.z - normal.z*t);
	}

	/** Returns 1/0/-1 if point p is above/on/below this plane. */
	public int above(Point p) {
		double dotP = normal.dot(p.toVector());
		if (dotP > -d) return 1;
		if (dotP < -d) return -1; 
		return 0;
	}

	/** Returns 1/0/-1 if point p is below/on/above this plane */
	public int below(Point p) { return -above(p); }

	/** Gets perpendicular distance of point to this plane */
	public double getDistance(Point p) { 
		return Math.abs(normal.dot(p.toVector()) + d) / normal.getLength(); 
	}

	/** Get the angle between two planes (no sign). */
	public static double dihedralAngle(Plane h1, Plane h2) { 
		return Math.acos(h1.normal.dot(h2.normal)); 
	} 

	/** Get the intersection of a line with the plane. Returns null if line is 
	 * parallel to plane. */
	public Point getIntersection(Line line) {
		double denom = normal.dot(line.getDir());
		if (denom==0) return null;
		else {
			Point a = line.getP();
			Vector pa = point.vectorTo(a);
			double u = normal.dot(pa)/denom;
			return new Point(a.x + u*line.dir.x, a.y + u*line.dir.y, a.z + u*line.dir.z);
		}
	}
	

	/** Get the line-parameter of the intersection between a plane and a line. Returns infinity 
	 * if line is parallel to plane. TODO: Consider moving to Line3d */
	public double getIntersectionParameter(Line line) {
		double denom = normal.dot(line.getDir());
		if (denom == 0) return Double.POSITIVE_INFINITY;
		else {
			Point a = line.getP();
			Vector pa = point.vectorTo(a);
			double u = normal.dot(pa)/denom;
			return u;
		}
	}

	/** Get the intersection of a segment with the plane. Returns null if line is 
	 * parallel to plane.*/
	public Point getIntersection(Segment sgm) {
		Vector dir = sgm.getAToB();
		double denom = normal.dot(dir);
		if (denom == 0) return null;
		else {
			Vector pa = point.vectorTo(sgm.a);
			double u = normal.dot(pa)/denom;
			if ((u < 0) || (u > 1)) return null;
			else return new Point(sgm.a.x + u*dir.x, sgm.a.y + u*dir.y, sgm.a.z + u*dir.z);
		}
	}

}
