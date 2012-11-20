package ProGAL.geom2d;

import java.util.Random;

import ProGAL.dataStructures.Set;

public class PointSet extends Set<Point>{

	public PointSet() { super(); }
	
	/** creates a point set consisting of n uniformly distributed points in the unit square */
	public PointSet(int n) {
		Random random = new Random(3);
		for (int i = 0; i < n; i++) insert(new Point(new Double(random.nextDouble()), new Double(random.nextDouble())));
	}

	public Point getCentroid() {
		double x = 0.0; 
		double y = 0.0;
		for (Point p : this) { x += p.x(); y += p.y(); }
		return new Point(x/getSize(), y/getSize());
	}

	
	/** returns the index of the leftmost point (in case of ties, index of the topmost one is returned) */
	public int leftExtremePointIndx() {
		int k = 0;
		for (int i = 1; i < n; i++)
			if ((get(i).x() < get(k).x()) || ((get(i).x() == get(k).x()) && (get(k).y() > get(i).y()))) k = i;
		return k;
	}
	
	/** returns the leftmost point (in case of ties, the topmost one is returned */
	public Point leftExtremePoint() { return get(leftExtremePointIndx()); }

	/** returns next extreme point of the point set (counterclockwise)  */
	public Point getNextExtremePoint(Point p) {
		Point r = (get(0) != p)? get(0) : get(1);
		for (Point q : this) 
			if ((q != p) && (q != r)) {
				if (Point.rightTurn(p, r, q) ||
				    (Point.collinear(p,  r, q) && (p.distanceSquared(r) < p.distanceSquared(q)))) r = q; 
			}
		return r;
	}
	
	/** translates all points by the specified vector */
	public void translate(Vector v) { for (Point p : this) p.addThis(v); }
	
	/** prints points on the console */
	public void toConsole(int dec) { for (Point p : this) p.toConsole(dec); }
	
}
