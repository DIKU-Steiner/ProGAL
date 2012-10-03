package ProGAL.geomNd;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import ProGAL.math.Matrix;
import ProGAL.math.Randomization;

/**
 * An ArrayList-wrapper for storing Nd points. Adds functionality specific to 
 * point-sets such as finding centroid, diameter or covariance.  
 */
public class PointList extends ArrayList<Point> {
	private static final long serialVersionUID = -2860146038113339241L;

	private final int dim;
	
	/** Construct an empty point-list. */
	public PointList(int dim) { 
		super(); 
		this.dim = dim;
	}
	
	/** Construc a point-list from an array of points. */
	public PointList(Point[] elements) { 
		this(elements[0].dim);
		for(Point p: elements) add(p);
	}

	public boolean add(Point p){
		if(dim!=p.dim) throw new IllegalArgumentException("Dimension of point must match that of PointList");
		return super.add(p);
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
		PointList ret = new PointList(dim);
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
		double[] cen = new double[dim];
		for(Point p: this)
			for(int d=0;d<dim;d++) cen[d]+=p.get(d);
		
		for(int d=0;d<dim;d++) cen[d]/=size();
		return new Point(cen);
	}
	
	/** Get the variance of the points. */
	public double getVariance() {
		Point c = getCentroid();
		double sum = 0.0;
		for (Point p: this) sum += c.distanceSquared(p);
		return sum/size();
	}
	
	/** Get the covariance of the points. */
	public Matrix getCovariance() {
		Matrix cov = new Matrix(dim, dim);
		double cv, ci, cj;
		Point c = getCentroid();
		Point p;
		for (int i = 0; i < dim; i++) {
			ci = c.getCoord(i);
			for (int j = i; j < dim; j++) {
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


	/** Get a shallow copy of this list. */
	public PointList clone(){
		PointList ret = new PointList(dim);
		ret.addAll(this);
		return ret;
	}
	
	/** Writes this point-list to <code>System.out</code>. */
	public void toConsole() {
		toConsole(2);
	}
	
	/** Writes this point-list to <code>System.out</code>. */
	public void toConsole(int dec) {
		System.out.println("PointList:");
		for (int i = 0; i < size(); i++){
			System.out.print(String.format("%3d> ",i));
			get(i).toConsole(dec);
		}
	}

	public static List<Point> generatePointsInCube(int n, int d) {
		List<Point> ret = new LinkedList<Point>();
		for(int i=0;i<n;i++){
			double[] coords = new double[d];
			for(int j=0;j<d;j++) coords[j] = Randomization.randBetween(0.0, 1.0);
			ret.add(new Point(coords));
		}
		return ret;
	}

}

