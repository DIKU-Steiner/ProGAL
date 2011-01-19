package ProGAL.geom3d;

import ProGAL.math.Constants;

/**
 * A line represented by a point and a direction. Several methods work with a  
 * parameter that can be used to specify a point on the line. For instance 
 * <pre>
 * Line3d l = new Line3d( new Point3d(0,0,0), new Vector3d(2,0,0) );
 * System.out.println( l.getPoint( 0.0 ) );
 * System.out.println( l.getPoint( 1.0 ) );
 * </pre>
 * will print the points (0,0,0) and (2,0,0). Similarly, the lines  
 * <pre>
 * l.orthogonalProjection(new Point3d(1,1,0));
 * l.orthogonalProjectionParameter(new Point3d(1,1,0));
 * </pre>
 * will return the point (1,0,0) and the parameter 0.5 respectively. 
 */
public class Line {
	protected Point p;
	protected Vector dir;
	
	/** Constructs a line through origo with direction d. */
	public Line(Vector d) {
		p = new Point(0, 0, 0);
		dir = d;
	}
	
	/** Constructs a line through p with direction d.*/
	public Line(Point p, Vector d) {
		this.p = p;
		dir = d;
	}
	
	/** 
	 * Constructs a line through the segment s. The point will be s.getA() and the direction 
	 * will be s.getAToB(). Subsequent changes to the line segment should not change this line.
	 */
	public Line(LineSegment s) {
		p = s.getA().clone();
		dir = s.getAToB();
	}
	
	/** Returns the point defining this line. */
	public Point getP()   { return p; }

	/** Returns the direction vector defining this line. */
	public Vector getDir() { return dir; }

	/** Gets the point on the line defined by the specified parameter. If the parameter is 
	 * zero this method will not return a reference to the defining point (but the returned 
	 * point will equal it). */
	public Point getPoint(double t) {
		return new Point(p.getX() + t*dir.getX(), p.getY() + t*dir.getY(), p.getZ() + t*dir.getZ());
	}
	
	/** Returns the othogonal projection of the point q onto this line. */
	public Point orthogonalProjection(Point q) {
		return getPoint(orthogonalProjectionParameter(q));
	}
	
	/** Returns the line-parameter of the othogonal projection of the point q onto this line. */
	public double orthogonalProjectionParameter(Point q) {
		Vector pq = p.vectorTo(q);
		return pq.dot(dir)/dir.getLengthSquared();
	}
	
	/** Returns the smallest segment that contains all orthogonol 
	 * projections of a point set onto this line. */
	public LineSegment orthogonalProjection(PointList points) {
		double minT = Double.POSITIVE_INFINITY, maxT = Double.NEGATIVE_INFINITY;
		for(Point q: points){
			double t = orthogonalProjectionParameter(q);
			minT = Math.min(minT, t);
			maxT = Math.max(maxT, t);
		}
		return new LineSegment(getPoint(minT), getPoint(maxT)); 
	}
	
	/** Returns the line-parameters of the end-points of the smallest segment that 
	 * contains all orthogonol projections of a point set onto this line. 
	 * @return a double-array with two entries containing the smallest and largest 
	 * orthogonal projection line-parameter of the point-set.*/
	public double[] orthogonalProjectionParameters(PointList points) {
		double minT = Double.POSITIVE_INFINITY, maxT = Double.NEGATIVE_INFINITY;
		for(Point q: points){
			double t = orthogonalProjectionParameter(q);
			minT = Math.min(minT, t);
			maxT = Math.max(maxT, t);
		}
		return new double[]{minT, maxT}; 
	}

	/** Gets the squared orthogonal distance to a point. */
	public double getDistanceSquared(Point q) { 
		return dir.cross(q.vectorTo(p)).getLengthSquared()/dir.getLengthSquared();
	}
	
	/** Gets the orthogonal distance to a point. */
	public double getDistance(Point q) { 
		return Math.sqrt(getDistanceSquared(q));
	}

	/** 
	 * Gets the minimum squared distance to another line. 
	 * @see [Ericsson 05, p. 147]
	 * @hop 14-32
	 */
	public double getSquaredDistance(Line l) {
		double a = dir.getLengthSquared();
		double b = dir.dot(l.dir);
		double e = l.dir.getLengthSquared();
		double d = a*e-b*b;
		
		if(Math.abs(d)<Constants.EPSILON) return getDistanceSquared(l.p);
		
		Vector r = l.p.vectorTo(p);
		double c = dir.dot(r);
		double f = l.dir.dot(r);
		
		double s = (b*f-c*e)/d;
		double t = (a*f-b*c)/d;

		double dx = (p.getX()+s*dir.getX()) - (l.p.getX()+t*l.dir.getX());
		double dy = (p.getY()+s*dir.getY()) - (l.p.getY()+t*l.dir.getY());
		double dz = (p.getZ()+s*dir.getZ()) - (l.p.getZ()+t*l.dir.getZ());
		return dx*dx+dy*dy+dz*dz;
		//return getPoint(s).getDistanceSquared(l.getPoint(t));
	}
	
	/** 
	 * Gets the intersection-point of this line with l. If the lines do not intersect 
	 * then null is returned.
	 */
	public Point getIntersection(Line l){
		double a = dir.getLengthSquared();
		double b = dir.dot(l.dir);
		double e = l.dir.getLengthSquared();
		double d = a*e-b*b;
		
		if(Math.abs(d)<Constants.EPSILON) return null;
		
		Vector r = l.p.vectorTo(p);
		double c = dir.dot(r);
		double f = l.dir.dot(r);
		
		double s = (b*f-c*e)/d;
		double t = (a*f-b*c)/d;

		double dx = (p.getX()+s*dir.getX()) - (l.p.getX()+t*l.dir.getX());
		double dy = (p.getY()+s*dir.getY()) - (l.p.getY()+t*l.dir.getY());
		double dz = (p.getZ()+s*dir.getZ()) - (l.p.getZ()+t*l.dir.getZ());
		if(dx*dx+dy*dy+dz*dz>Constants.EPSILON) return null;
		return new Point(p.getX()+s*dir.getX(), p.getY()+s*dir.getY(), p.getZ()+s*dir.getZ());
		
	}

	/** Gets the largest squared distance from the points to the line. */
	public double getMaxDistanceSquared(PointList points) {
		if (points.size() == 0) throw new Error("No point");
		double maxDist = Double.NEGATIVE_INFINITY;
		for (Point q: points) {
			double dist = getDistanceSquared(q);
			if (dist > maxDist) maxDist = dist;
		}
		return maxDist;
	}

	/** Gets the largest distance from the points to the line. */
	public double getMaxDistance(PointList points) {
		return Math.sqrt(getMaxDistanceSquared(points));
	}

	/** Returns a string-representation of this line.*/
	public String toString(){
		return String.format("Line3d[p:%s,dir:%s]", p, dir);
	}
	
	/** Returns a string-representation of this line with <code>dec</code> decimals precision.*/
	public String toString(int dec){
		return String.format("Line3d[p:%s,dir:%s]", p.toString(dec), dir.toString(dec));
	}
}
