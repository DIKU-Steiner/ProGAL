
package ProGAL.geom2d;

import ProGAL.math.Constants;

public class Line implements Shape {
	protected Point p;     // point on the line
	protected Vector n;    // normal vector of the line, unit length
	
	/** creates a line through a given point and with a given normal vector */
	public Line(Point p, Vector n) { 
		this.p = p;  
		this.n = n.normalize(); 
	}
	
	/** creates a line through 2 given points  */
	public Line(Point p, Point q) {
		this.p = p;
		n = new Vector(p.y()-q.y(), q.x()-p.x());
		n.normalizeThis();
	}
	
	/** creates a line through a given segment */
	public Line(LineSegment seg) {
		p = seg.a;
		n = new Vector(seg.a.y()-seg.b.y(), seg.b.x() - seg.a.x());
		n = n.normalize();
	}
	
	/** creates a line ax + by + c = 0 */
	public Line(double a, double b, double c) {
		n = new Vector(a,b);
		n.normalizeThis();
		if (b != 0.0) p = new Point(0.0,-c/b); else p = new Point(-c/a,0.0);
	}

	public Vector getDirection() { return new Vector(n.y(),-n.x()); }
	
	public double getSlope() { if (!isVertical()) return n.x()/n.y(); else return Double.MAX_VALUE;
	}
	
	public boolean isVertical() { return n.y() == 0.0; }
	
	public boolean isParallelWith(Line l) { return Math.abs(Vector.crossProduct(n, l.n)) < Constants.EPSILON; }
	
	public static boolean areParallel(Line l1, Line l2) {
		return Vector.crossProduct(l1.n, l2.n) == 0.0;
	}
	
	/*
	 * project point p onto this line
	 */
	public Point projectPoint(Point q) {
		double t = n.y()*(q.x()-p.x()) - n.x()*(q.y()-p.y());
		return new Point(n.y()*t+p.x(), p.y()-n.x()*t);
	}
	
	public double projectionParameter(Point q){
		return n.y()*(q.x()-p.x()) - n.x()*(q.y()-p.y());
	}
	
	/*
	 * returns the intersection of two lines using Cramer's rule, see Ericson, Section 3.1.5. 
	 * Runtime exception is thrown if the lines are parallel.
	 * Modified by Pawel on June 23, 2010
	 */
	public static Point getIntersection(Line l1, Line l2) {
		double denom = l1.n.x()*l2.n.y() - l1.n.y()*l2.n.x();
		if (Math.abs(denom) < Constants.EPSILON) throw new RuntimeException("Lines are parallel");
		else {
			double e = l1.n.x()*l1.p.x() + l1.n.y()*l1.p.y();
			double f = l2.n.x()*l2.p.x() + l2.n.y()*l2.p.y();
			return new Point((e*l2.n.y() - f*l1.n.y())/denom, (f*l1.n.x() - e*l2.n.x())/denom);
		}
	}
	
	public String toString(String name) {
		return "Line[" + name + ",point:" + p.toString() + ",normal:" + n.toString()+"]";
	}
	@Override
	public String toString() { return toString(""); }
	
	public void toConsole(String name) { System.out.println(toString(name)); }
	public void toConsole() { System.out.println(toString("")); }

	/** returns the point on the line at distance d from the line-defining point p */
	public Point getPoint(double d) {
		Vector dir = this.getDirection();
		return new Point(p.x()+d*dir.x(), p.y()+d*dir.y());
	}

	/** returns the distance of the point q to the line */
	public double getDistance(Point q) {
		return Math.abs(n.x()*q.x() + n.y()*q.y() - n.x()*p.x() - n.y()*p.y())/Math.sqrt(n.x()*n.x() + n.y()*n.y());
	}

	
	public double intersectionParameter(Line l) {
		Vector dir = getDirection();
		Vector lDir = l.getDirection();
		double denom = lDir.y()*dir.x()-lDir.x()*dir.y();
		Vector c = l.p.vectorTo(p);
		double s = (lDir.x()*c.y()-lDir.y()*c.x())/denom;
		return s;
	}

	@Override
	public Point getCenter() {
		return getPoint(0);
	}
	
}
