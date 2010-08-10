package ProGAL.geom3d;


public class Plane3d {

	protected Vector3d n; // normal vector of the plane - unit vector
	protected Point3d p;  // point in the plane
	protected double d;   // perpendicular distance from the origin
	
	/*
	 * creates a plane with the normal vector n and containing point p
	 */
	public Plane3d(Vector3d n, Point3d p) {
		this.n = n.createUnitVector();
		this.p = p;
		d = -n.x*p.x - n.y*p.y - n.z*p.z;
	}
	
	/*
	 * creates a plane with the normal vector n and trough the tip of vector v 
	 */
	public Plane3d(Vector3d n, Vector3d v) {
		this.n = n.createUnitVector();
		this.p = new Point3d(v.x, v.y, v.z);
		d = -n.x*p.x - n.y*p.y - n.z*p.z;
	}
	
	/*
	 * creates a plane through three distinct points
	 */
	public Plane3d(Point3d p, Point3d q, Point3d r) {
		n = Vector3d.crossProduct(new Vector3d(q,p),new Vector3d(q,r));
//		n.makeUnitVector();
		this.p = q;
		d = -n.x*p.x - n.y*p.y - n.z*p.z;
	}
	
	public void setNormal(Vector3d n) { this.n = n; }
	
	/*
	 * projects point p onto this plane
	 */
	public Point3d projectPoint(Point3d p) {
		double t = n.x*p.x + n.y*p.y + n.z*p.z +d;
		return new Point3d(p.x - n.x*t, p.y - n.y*t, p.z - n.z*t);
	}
	
	/*
	 * returns (1,0,-1) if point p is (above,on,below) this plane
	 */
	public int above(Point3d p) {
		double dotP = Vector3d.dotProduct(n,new Vector3d(p));
		if (dotP > -d) return 1;
		else {
			if (dotP < -d) return -1; else return 0;
		}
	}
	
	/*
	 * returns (1,0,-1) if point p is (below,on,above) this plane
	 */
	public int below(Point3d p) { return -above(p); }
	
	/*
	 * returns perpendicular distance of point p to this plane
	 */
	public double getDistance(Point3d p) { return Math.abs(Vector3d.dotProduct(n,p) + d) / n.getLength(); }
	
	/*
	 * returns the angle between two planes (no sign)
	 */
	public static double dihedralAngle(Plane3d h1, Plane3d h2) { return Math.acos(Vector3d.dotProduct(h1.n, h2.n)); } 
	
	/*
	 * returns intersection of a line with the plane
	 */
	public Point3d getIntersection(Line3d line) {
		double denom = Vector3d.dotProduct(n, line.getDir());
		if (denom < Constants.epsilon) return null;
		else {
			Point3d a = line.getP();
			Vector3d pa = new Vector3d(p, a);
			double u = Vector3d.dotProduct(n, pa)/denom;
			return new Point3d(a.x + u*line.dir.x, a.y + u*line.dir.y, a.z + u*line.dir.z);
		}
	}
	
	/*
	 * returns intersection of a segment with the plane
	 */
	public Point3d getIntersection(Segment3d sgm) {
		Vector3d dir = new Vector3d(sgm.a, sgm.b);
		double denom = Vector3d.dotProduct(n, dir);
		if (denom < Constants.epsilon) return null;
		else {
			Vector3d pa = new Vector3d(p, sgm.a);
			double u = Vector3d.dotProduct(n, pa)/denom;
			if ((u < 0) || (u > 1)) return null;
			else {
				return new Point3d(sgm.a.x + u*dir.x, sgm.a.y + u*dir.y, sgm.a.z + u*dir.z);}
		}
	}
	
	public static void main(String[] args) {
		Plane3d pl = new Plane3d(new Vector3d(1.0,1.0,1.0), new Point3d(0.0,0.0,2.0));
		Point3d p = new Point3d(11,12,13);
		Point3d q = pl.projectPoint(p);
		System.out.println(q.toString(2));
	}

	/**
	 * @author Rasmus 
	 */
	public double intersectionParameter(Line3d line) {
		double denom = Vector3d.dotProduct(n, line.getDir());
		if (denom < Constants.epsilon) throw new RuntimeException("Line is parallel to plane");
		else {
			Point3d a = line.getP();
			Vector3d pa = new Vector3d(p, a);
			double u = Vector3d.dotProduct(n, pa)/denom;
			return u;
		}
	}
}
