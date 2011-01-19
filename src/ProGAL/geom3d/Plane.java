package ProGAL.geom3d;

/**
 * A plane in (x,y,z)-space represented by a point and a normal. 
 * 
 * Assuming that <i>p</i> is a point on the plane and <i>n</i> is the normal vector, 
 * the half-space <i>{q|pq·n>0}</i> is called the 'upper halfspace' wrt. the plane 
 * and vice versa for the 'lower halfspace'.
 */
public class Plane implements Shape{
	/** Normal vector of the plane. */
	protected Vector normal; 
	/** Point in the plane. */
	protected Point point;  

	/** Constructs a plane with the normal vector n containing point p. */
	public Plane(Point p, Vector n) {
		this.normal = n;
		this.point = p;
	}

	/** Constructs a plane with the normal vector n containing the point (0,0,0). */
	public Plane(Vector n) {
		this.normal = n;
		this.point = new Point(0,0,0);
	}

	/** 
	 * Constructs a plane through three points using the first point as defining point. 
	 * The normal of the plane will be decided by the order of p, q and r such that if 
	 * the right hand follows the rotation from p to q to r the thumb points in the 
	 * normals direction. TODO: Test this
	 * 
	 * An error is thrown if the points are collinear.
	 */
	public Plane(Point p, Point q, Point r) {
		if(Point.collinear(p, q, r)) throw new Error("Cant construct plane: Points are collinear");
		normal = p.vectorTo(q).crossThis(p.vectorTo(r)).normalizeThis();
		this.point = p;
	}
	
	private double getD(){
		return -normal.getX()*point.getX() - normal.getY()*point.getY() - normal.getZ()*point.getZ();
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
		double t = normal.getX()*p.getX() + normal.getY()*p.getY() + normal.getZ()*p.getZ() + getD();
		return new Point(p.getX() - normal.getX()*t, p.getY() - normal.getY()*t, p.getZ() - normal.getZ()*t);
	}

	/** Returns 1/0/-1 if point p is above/on/below this plane. */
	public int above(Point p) {
		double dotP = normal.dot(p.toVector());
		double d = getD();
		if (dotP > -d) return 1;
		if (dotP < -d) return -1; 
		return 0;
	}

	/** Returns 1/0/-1 if point p is below/on/above this plane */
	public int below(Point p) { return -above(p); }

	/** Gets perpendicular distance of point to this plane */
	public double getDistance(Point p) {
		double d = getD();
		return Math.abs(normal.dot(p.toVector()) + d) / normal.getLength(); 
	}

	/** Get the unsigned angle between this plane and p. */
	public double getUnsignedDihedralAngle(Plane p){
		return Math.acos(normal.dot(p.normal));
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
			return new Point(a.getX() + u*line.dir.getX(), a.getY() + u*line.dir.getY(), a.getZ() + u*line.dir.getZ());
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
	public Point getIntersection(LineSegment sgm) {
		Vector dir = sgm.getAToB();
		double denom = normal.dot(dir);
		if (denom == 0) return null;
		else {
			Vector pa = point.vectorTo(sgm.a);
			double u = normal.dot(pa)/denom;
			if ((u < 0) || (u > 1)) return null;
			else return new Point(sgm.a.getX() + u*dir.getX(), sgm.a.getY() + u*dir.getY(), sgm.a.getZ() + u*dir.getZ());
		}
	}

	/** Returns the defining point for this plane. The center of a plane is not well-defined, so  
	 * to implement the shape interface the defining point is simply used. */
	public Point getCenter() {
		return point.clone();
	}

}
