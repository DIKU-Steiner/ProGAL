package ProGAL.geom3d;

import java.awt.Color;

import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Cylinder;
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
		dir = d.normalize();
	}
	
	/** Constructs a line through p with direction d.*/
	public Line(Point p, Vector d) {
		this.p = p;
		dir = d.normalize();
	}
		
	/** 
	 * Constructs a line through the segment s. The point will be s.getA() and the direction 
	 * will be s.getAToB(). Subsequent changes to the line segment should not change this line.
	 */
	public Line(LineSegment s) {
		p = s.getA().clone();
		dir = s.getAToB().normalize();
	}
	
	/** Constructs a line through the two specified points. */
	public Line(Point p1, Point p2){
		this(p1, p1.vectorTo(p2));
	}
	
	/** Constructs a line trisecting three points */
	public Line(Point a, Point b, Point c) {
		Circle circle = new Circle(a,b,c);
		p = circle.getCenter();
		dir = new Vector(a,b).cross(new Vector(a,c)).normalizeThis();
	}
	
	/** Construct a line that is a clone of L. */
	public Line clone() { 
		return new Line(new Point(p), new Vector(dir));
	}

	/** Returns the point defining this line. */
	public Point getP()   { return p; }

	/** Returns the direction vector defining this line. */
	public Vector getDir() { return dir; }

	/** Gets the point on the line defined by the specified parameter. If the parameter is 
	 * zero this method will not return a reference to the defining point (but the returned 
	 * point will equal it). */
	public Point getPoint(double t) {
		return new Point(p.x() + t*dir.x(), p.y() + t*dir.y(), p.z() + t*dir.z());
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

		double dx = (p.x()+s*dir.x()) - (l.p.x()+t*l.dir.x());
		double dy = (p.y()+s*dir.y()) - (l.p.y()+t*l.dir.y());
		double dz = (p.z()+s*dir.z()) - (l.p.z()+t*l.dir.z());
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

		double dx = (p.x()+s*dir.x()) - (l.p.x()+t*l.dir.x());
		double dy = (p.y()+s*dir.y()) - (l.p.y()+t*l.dir.y());
		double dz = (p.z()+s*dir.z()) - (l.p.z()+t*l.dir.z());
		if(dx*dx+dy*dy+dz*dz>Constants.EPSILON) return null;
		return new Point(p.x()+s*dir.x(), p.y()+s*dir.y(), p.z()+s*dir.z());
		
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
	
	/** 
	 * Return a rotation of p around this line. The rotation is a right-handed one  
	 * (thumb in the direction of dir)
	 */
	public Point rotate(Point p, double angle){
		Point ret = p.clone();
		return rotateIn(ret, angle);
	}
	
	/** 
	 * Rotate point around this line, store the result in point and return the results. The rotation
	 * is a right-handed one (thumb in the direction of dir)
	 */
	public Point rotateIn(Point point, double angle){
		Vector v = p.vectorTo(point);
		dir.rotateIn(v, angle);
		for(int i=0;i<3;i++)
			point.set(i, v.get(i)+p.get(i));
		return point;
	}
	
	/**
	 * Return the optimal right-hand rotation around the line that brings the m-points as 
	 * close to the f-points as possible. Formally the rotation angle that minimizes 
	 * <code>double S = m[0].distanceSquared(f[0]) + ... + m[m.length-1].distanceSquared(f[m.length-1]);</code>
	 * This method follows the description by Canutescu and Dunbrack 2003.
	 * @param moving Array of effector points
	 * @param target Array of target points 
	 * @return an angle which the effector points should be rotated so they come close to the target points. 
	 */
	public double optimalRotation(Point[] moving, Point[] target){
		Point[] m = moving;
		Point[] f = target;
		if(f.length!=m.length) throw new RuntimeException("Lengths must match");
		
		double tanNum = 0, tanDenom = 0;
		for(int i=0;i<m.length;i++){
			Point O_i = orthogonalProjection(m[i]);
			Vector fVec = O_i.vectorTo(f[i]);
			Vector rVec = O_i.vectorTo(m[i]);
			Vector sVec = dir.cross(rVec).normalizeThis();
			double rLen = rVec.length();
			rVec.multiplyThis(1/rLen);
			tanNum		+= fVec.dot(sVec)*rLen;
			tanDenom 	+= fVec.dot(rVec)*rLen;
		}
		double alpha = Math.atan(tanNum/tanDenom);
		double secDeriv = Math.cos(alpha)*tanDenom + Math.sin(alpha)*tanNum;
		
		if(secDeriv<0) alpha -= Math.signum(alpha)*Math.PI;
		
		return alpha;
	}
	

	/** Returns a string-representation of this line.*/
	public String toString(){
		return String.format("Line3d[p:%s,dir:%s]", p, dir);
	}
	
	/** Returns a string-representation of this line with <code>dec</code> decimals precision.*/
	public String toString(int dec){
		return String.format("Line3d[p:%s,dir:%s]", p.toString(dec), dir.toString(dec));
	}
	
	public Cylinder toScene(J3DScene scene, double rad, Color clr) {
		LineSegment seg = new LineSegment(this.p.add(dir.multiply(1000)), p.add(dir.multiply(-1000)));
		return seg.toScene(scene, rad, clr);
	}
}
