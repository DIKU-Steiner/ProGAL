package ProGAL.geom3d;

import java.util.ArrayList;
import java.util.Collections;

import ProGAL.math.Matrix3x3;
import ProGAL.math.Randomization;

/**
 * An ArrayList-wrapper for storing 3d points. Adds functionality specific to 
 * point-sets such as finding centroid, diameter or covariance.  
 */
public class PointList3d extends ArrayList<Point3d> {
	private static final long serialVersionUID = -4824374877674925546L;

	/** Construct an empty point-list. */
	public PointList3d() { super(); }
	
	
	/** Construc a point-list from an array of points. */
	public PointList3d(Point3d[] elements) { 
		this();
		for(Point3d p: elements) add(p);
	}
	
	/** Returns the i'th coordinate of k'th point. */
	public double getCoord(int k, int i) { return get(k).getCoord(i); }
	
	/** 
	 * Returns a sub-list of elements from (including) <code>from</code> to
	 * (not including) <code>to</code>. The returned sub-list is not a 'view' 
	 * of the sublist (such as AbstractList.subList), but a shallow copy of 
	 * the specified range.   
	 */
	public PointList3d getSubList(int from, int to){
		PointList3d ret = new PointList3d();
		for(int i=from;i<to;i++) ret.add(get(i));
		return ret;
	}
	
	/** Returns a new point-list containing all the elements in random order. */
	public PointList3d getRandomPermutation(){
		PointList3d ret = this.clone();
		Collections.shuffle(ret, Randomization.getGenerator());
		return ret;
	}
	
	/** Get the centroid of the points. */
	public Point3d getCentroid() {
		double x = 0.0, y = 0.0, z = 0.0;
		for(Point3d p: this){  
			x += p.x; y += p.y; z += p.z;
		}
		int n = size();
		return new Point3d(x/n, y/n, z/n);
	}
	
	/** Get the variance of the points. */
	public double getVariance() {
		Point3d c = getCentroid();
		double sum = 0.0;
		for (Point3d p: this) sum += c.getDistanceSquared(p);
		return sum/size();
	}
	
	/** Get the covariance of the points. */
	public Matrix3x3 getCovariance() {
		Matrix3x3 cov = new Matrix3x3();
		double cv, ci, cj;
		Point3d c = getCentroid();
		Point3d p;
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
	public Point3d getExtreme(Vector3d direction) {
		double maxDot = Double.NEGATIVE_INFINITY;
		Point3d ret = null;
		for (Point3d p: this) {
			double dot = direction.dot(p.toVector());
			if(dot>maxDot) {
				maxDot = dot;
				ret = p;
			}
		}
		return ret;	
	}
	
	/** Get the most extreme point in the direction specified by the vector. */
	public int getExtremeIndex(Vector3d direction) {
		double maxDot = Double.NEGATIVE_INFINITY;
		int ret = -1;
		for (Point3d p: this) {
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
		Point3d q = get(0);
		for (int i = 1; i < size(); i++) {
			Point3d p = (Point3d)get(i);
			if (high) { 
				if (p.dominates(q, ix, iy, iz)) { indx = i; q = p; }
			} else {
				if (q.dominates(p, ix, iy, iz)) { indx = i; q = p; }
			}
		}
		return indx;
	}
	
	/** Get the rightmost point (in case of ties the top rightmost point is returned). */
	public Point3d getExtremeRight() { return get(getExtremeIndex(0,1,2,true)); }

	/** Get the leftmost point (in case of ties the bottom lefttmost point is returned). */
	public Point3d getExtremeLeft() { return get(getExtremeIndex(0,1,2,false)); }

	/** Get the highest point (in case of ties the left topmost point is returned).	 */
	public Point3d getExtremeTop() { return (Point3d)get(getExtremeIndex(1,0,2,true)); }

	/** Get the lowest point (in case of ties the bottom rightmost point is returned). */
	public Point3d getExtremeBottom() { return (Point3d)get(getExtremeIndex(1,0,2,false)); }

	/** Get the frontmost point (in case of ties the front rightmost point is returned). */
	public Point3d getExtremeFront() { return (Point3d)get(getExtremeIndex(1,2,0,true)); }
	
	/** Get the deepest point (in case of ties the left topmost point is returned). */
	public Point3d getExtremeBack() { return (Point3d)get(getExtremeIndex(1,2,0,false)); }
	
	/** Get the diameter of the point set - trivial O(n^2) algorithm.*/
	public Segment3d getDiameter() {
		Point3d p;
		Point3d q;
		Point3d best1=null, best2=null;
		double pq; 
		double best = 0.0;
		for (int i = 0; i < size()-1; i++) {
			p = (Point3d)get(i);
			for (int j = i+1; j < size(); j++) {
				q = get(j);
				pq = p.getDistanceSquared(q);
				if (pq > best) { best = pq; best1 = p; best2 = q; }
			}
		}
		return new Segment3d(best1, best2);
	}
	
	/** Get a segment <code>seg</code> between two points in the set such that 
	 * <code>diameter/sqrt(3) <= seg</code>. Requires O(n) time. */
	public Segment3d diameterSqrt3Approx() {
		Segment3d seg0 = new Segment3d(getExtremeLeft(), getExtremeRight());
		Segment3d seg1 = new Segment3d(getExtremeBottom(), getExtremeTop());
		Segment3d seg2 = new Segment3d(getExtremeFront(), getExtremeBack());
		double l0 = seg0.getSquaredLength();
		double l1 = seg1.getSquaredLength();
		double l2 = seg2.getSquaredLength();
		if (l0 < l1) return (l2 < l0)? seg2 : seg0;  else return (l2 < l1)?  seg2 :seg1; 
	}

	/** Get a shallow copy of this list. */
	public PointList3d clone(){
		PointList3d ret = new PointList3d();
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
	 * the unit cube. Uses the math.Randomization class.*/
	public static PointList3d genPointsInCube(int n) {
		PointList3d list = new PointList3d();
		for (int i = 0; i < n; i++) 
			list.add(new Point3d( 
					Randomization.randBetween(-1.0,1.0),
					Randomization.randBetween(-1.0,1.0),
					Randomization.randBetween(-1.0,1.0) 
					));
		return list;
	}
}

