package ProGAL.geomNd;

import ProGAL.math.Constants;

/** 
 * A class representing a point in an N-dimensional Euclidean space. The underlying 
 * representation is a double-array. Changes to the elements of this double-array 
 * will also result in a change to the point. The double-array can be specified on  
 * construction and can be retrieved lated. 
 * @author R. Fonseca
 */
public class Point {
	/** The double-array that holds the value of this point */
	protected final double[] coords;
	/** The dimension of this point. */
	protected final int dim;
	
	public Point(double[] coords){
		this.coords = coords;
		this.dim = coords.length;
	}
	
	public Point(int dimensions){
		this.coords = new double[dimensions];
		this.dim = dimensions;
	}
	
	public double get(int d){			return coords[d];		}
	public double getCoord(int d) {		return coords[d];		}	
	public double[] getCoords() {		return coords;			}
	
	public int getDimension() {			return dim;				}

	public void set(int d, double v){		coords[d]=v;			}
	public void setCoord(int d, double v){	coords[d]=v;			}
	
	public double getDistanceSquared(Point p){
		double sum = 0;
		for(int d=0;d<dim;d++) {
			double delta = coords[d]-p.coords[d];
			sum+=delta*delta;
		}
		return sum;
	}
	public double getDistance(Point p){
		return Math.sqrt(getDistanceSquared(p));
	}
	
	/** Creates the midpoint of two points. */
	public static Point getMidpoint(Point p, Point q) {
		if(p.dim!=q.dim) throw new IllegalArgumentException("Dimension of points must match");
		double coords[] = new double[p.dim];
		for(int d=0;d<p.dim;d++) coords[d] = (p.coords[d]+q.coords[d])/2;
		return new Point(coords);
	}
	
	public static boolean collinear(Point p1, Point p2, Point p3) {
		double a = getAngle(p1,p2,p3);
		return ( (Math.abs(a-Math.PI)<Constants.EPSILON) || (Math.abs(a)<Constants.EPSILON) );
	}

	public static double getAngle(Point p1, Point p2, Point p3){
		Vector v1 = p2.vectorTo(p1);
		Vector v2 = p2.vectorTo(p3);
		return Math.acos(v1.dot(v2)/Math.sqrt(v1.getLengthSquared()*v2.getLengthSquared()));
	}
	
	private Vector vectorTo(Point p1) {
		double[] coords = new double[dim];
		for(int d=0;d<dim;d++) coords[d] = p1.coords[d]-coords[d];
		return new Vector(coords);
	}

	public Vector toVector(){
		return new Vector(coords);
	}
	
	public Point clone(){
		double[] newArr = new double[dim];
		for(int d=0;d<dim;d++) newArr[d] = coords[d];
		return new Point(newArr);
	}
	
	public void toConsole() {
		toConsole(2);
	}

	public void toConsole(int dec) {
		System.out.print("Point[");
		for(int d=0;d<dim-1;d++) System.out.printf("%."+dec+"f, ", coords[d]);
		if(dim>0)
			System.out.printf("%."+dec+"f]\n", coords[dim-1]);
	}

}
