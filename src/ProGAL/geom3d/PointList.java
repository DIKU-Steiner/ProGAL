package ProGAL.geom3d;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import ProGAL.math.Matrix3x3;
import ProGAL.math.Randomization;

/**
 * An ArrayList-wrapper for storing 3d points. Adds functionality specific to 
 * point-sets such as finding centroid, diameter or covariance.  
 */
public class PointList extends ArrayList<Point> {
	private static final long serialVersionUID = -4824374877674925546L;

	/** Construct an empty point-list. */
	public PointList() { super(); }
	
	
	/** Construct a point-list from an array of points. */
	public PointList(Point[] elements) { 
		this();
		for(Point p: elements) add(p);
	}
	
	public PointList(Collection<Point> points){
		super(points);
	}
	
	/** Returns the i'th coordinate of k'th point. */
	public double getCoord(int k, int i) { return get(k).getCoord(i); }
	
	/** 
	 * Returns a sub-list of elements from (including) <code>from</code> to
	 * (not including) <code>to</code>. The returned sub-list is not a 'view' 
	 * of the sublist (such as AbstractList.subList), but a shallow copy of 
	 * the specified range.   
	 */
	public PointList getSubList(int from, int to){
		PointList ret = new PointList();
		for(int i=from;i<to;i++) ret.add(get(i));
		return ret;
	}
	
	/** Returns a new point-list containing all the elements in random order. */
	public PointList getRandomPermutation(){
		PointList ret = this.clone();
		Collections.shuffle(ret, Randomization.getGenerator());
		return ret;
	}
		
	/** Get the centroid of the points. */
	public Point getCentroid() {
		double x = 0.0, y = 0.0, z = 0.0;
		for(Point p: this){  
			x += p.x(); y += p.y(); z += p.z();
		}
		int n = size();
		return new Point(x/n, y/n, z/n);
	}
	
	/** Get the variance of the points. */
	public double getVariance() {
		Point c = getCentroid();
		double sum = 0.0;
		for (Point p: this) sum += c.distanceSquared(p);
		return sum/size();
	}
	
	/** Get the covariance of the points. */
	public Matrix3x3 getCovariance() {
		Matrix3x3 cov = new Matrix3x3();
		double cv, ci, cj;
		Point c = getCentroid();
		Point p;
		for (int i = 0; i < 3; i++) {
			ci = c.getCoord(i);
			for (int j = i; j < 3; j++) {
				cj = c.getCoord(j);
				cv = 0.0;
				for (int k = 0; k < size(); k++) {
					p = get(k);
					cv += (p.getCoord(i) - ci) * (p.getCoord(j) - cj);
				}
				cov.set(i,j,cv/size()); 
				cov.set(j,i,cov.get(i,j));
			}
		}
		return cov;
	}
	
	/** Get the standard deviation of the point set. */
	public double getStandardDeviation() { return Math.sqrt(getVariance()); }

	/** Get the most extreme point in the direction specified by the vector. */
	public Point getExtreme(Vector direction) {
		double maxDot = Double.NEGATIVE_INFINITY;
		Point ret = null;
		for (Point p: this) {
			double dot = direction.dot(p.toVector());
			if(dot>maxDot) {
				maxDot = dot;
				ret = p;
			}
		}
		return ret;	
	}
	
	/** Get the most extreme point in the direction specified by the vector. */
	public int getExtremeIndex(Vector direction) {
		double maxDot = Double.NEGATIVE_INFINITY;
		int ret = -1;
		for (Point p: this) {
			double dot = direction.dot(p.toVector());
			if(dot>maxDot) {
				maxDot = dot;
				ret = indexOf(p);
			}
		}
		return ret;	
	}
	
	/**
	 * Get the index of the extreme point. Go along axis ix first, then iy and finally iz. 
	 * If high is true, the index of the point with highest coordinate is returned. Otherwise 
	 * the index of the point with lowest coordinate is returned. For instance
	 * <pre>
	 * PointList3d list = new PointList3d();
	 * list.add( new Point3d(0,0,0) );
	 * list.add( new Point3d(1,0,0) );
	 * list.add( new Point3d(0,1,0) );
	 * list.add( new Point3d(0,0,1) );
	 * System.out.println(list.getExtremeIndex(1,2,0,true));
	 * System.out.println(list.getExtremeIndex(1,2,0,false));
	 * </pre>
	 * will print the indices 2 and 0.
	 */
	public int getExtremeIndex(int ix, int iy, int iz, boolean high) {
		int indx = 0;
		Point q = get(0);
		for (int i = 1; i < size(); i++) {
			Point p = (Point)get(i);
			if (high) { 
				if (p.dominates(q, ix, iy, iz)) { indx = i; q = p; }
			} else {
				if (q.dominates(p, ix, iy, iz)) { indx = i; q = p; }
			}
		}
		return indx;
	}
	
	/** Get the rightmost point (in case of ties the top rightmost point is returned). */
	public Point getExtremeRight() { return get(getExtremeIndex(0,1,2,true)); }

	/** Get the leftmost point (in case of ties the bottom lefttmost point is returned). */
	public Point getExtremeLeft() { return get(getExtremeIndex(0,1,2,false)); }

	/** Get the highest point (in case of ties the left topmost point is returned).	 */
	public Point getExtremeTop() { return (Point)get(getExtremeIndex(1,0,2,true)); }

	/** Get the lowest point (in case of ties the bottom rightmost point is returned). */
	public Point getExtremeBottom() { return (Point)get(getExtremeIndex(1,0,2,false)); }

	/** Get the frontmost point (in case of ties the front rightmost point is returned). */
	public Point getExtremeFront() { return (Point)get(getExtremeIndex(1,2,0,true)); }
	
	/** Get the deepest point (in case of ties the left topmost point is returned). */
	public Point getExtremeBack() { return (Point)get(getExtremeIndex(1,2,0,false)); }
	
	/** Get the diameter of the point set - trivial O(n^2) algorithm.*/
	public LineSegment getDiameter() {
		Point p;
		Point q;
		Point best1=null, best2=null;
		double pq; 
		double best = 0.0;
		for (int i = 0; i < size()-1; i++) {
			p = (Point)get(i);
			for (int j = i+1; j < size(); j++) {
				q = get(j);
				pq = p.distanceSquared(q);
				if (pq > best) { best = pq; best1 = p; best2 = q; }
			}
		}
		return new LineSegment(best1, best2);
	}
	
	/** Get a segment <code>seg</code> between two points in the set such that 
	 * <code>diameter/sqrt(3) <= seg</code>. Requires O(n) time. */
	public LineSegment diameterSqrt3Approx() {
		LineSegment seg0 = new LineSegment(getExtremeLeft(), getExtremeRight());
		LineSegment seg1 = new LineSegment(getExtremeBottom(), getExtremeTop());
		LineSegment seg2 = new LineSegment(getExtremeFront(), getExtremeBack());
		double l0 = seg0.getLengthSquared();
		double l1 = seg1.getLengthSquared();
		double l2 = seg2.getLengthSquared();
		if (l0 < l1) return (l2 < l0)? seg2 : seg0;  else return (l2 < l1)?  seg2 :seg1; 
	}

	/** Get a shallow copy of this list. */
	public PointList clone(){
		PointList ret = new PointList();
		ret.addAll(this);
		return ret;
	}
	
	/** Writes this point-list to <code>System.out</code>. */
	public void toConsole() {
		System.out.println("PointList3d:");
		for (int i = 0; i < size(); i++){
			System.out.print(String.format("%3d> ",i));
			get(i).toConsole();
		}
	}
	
	/** Writes this point-list to <code>System.out</code>. */
	public void toConsole(int dec) {
		System.out.println("PointList3d:");
		for (int i = 0; i < size(); i++){
			System.out.print(String.format("%3d> ",i));
			get(i).toConsole(dec);
		}
	}

	/** Construct a point-list containing n uniformly distributed random points in 
	 * the unit cube. Uses the ProGAL.math.Randomization class.*/
	public static PointList generatePointsInCube(int n) { return generatePointsInCube(n, -1, 1, -1, 1, -1, 1); }

	public static PointList generatePointsInCube(int n, double xL, double xH, double yL, double yH, double zL, double zH) {
		PointList list = new PointList();
		for (int i = 0; i < n; i++) 
			list.add(new Point( 
					Randomization.randBetween(xL,xH),
					Randomization.randBetween(yL,yH),
					Randomization.randBetween(zL,zH) 
					));
		return list;
	}

	/** Construct a point-list of n uniformly distributed random points on the unit sphere.
	 * Uses the ProGAL.math.Randomization class. */
	public static PointList generateRandomPointsOnSphere(int n) {
		PointList list = new PointList();
		for(int i=0;i<n;i++){
				double theta0 = Randomization.randBetween(0, 2*Math.PI);
				double theta1 = Math.acos(Randomization.randBetween(-1.0,1.0));
			list.add(new Point(
					Math.sin(theta0)*Math.cos(theta1),
					Math.sin(theta0)*Math.sin(theta1),
					Math.cos(theta0)
					));
		}
		return list;
	}
	

	/** Construct a point-list of n evenly distributed points (not random) on the unit sphere.
	 * Taken from http://www.cgafaq.info/wiki/Evenly_distributed_points_on_sphere */
	public static PointList generatePointsOnSphere(int n){
		PointList list = new PointList();
		double l = 0;
		double dl = Math.PI*(3-Math.sqrt(5));
		double dz = 2.0/n;
		double z = 1-dz/2;
		for(int k=0;k<n;k++){
			double r = Math.sqrt(1-z*z);
			list.add(new Point(Math.cos(l)*r, Math.sin(l)*r,z));
			z-=dz;
			l+=dl;
		}
		return list;
	}
	
}

