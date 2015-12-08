package ProGAL.geomNd;

import java.io.Serializable;
import java.util.Arrays;

import ProGAL.math.Constants;
import ProGAL.math.Randomization;

/** 
 * A class representing a point in an N-dimensional Euclidean space. The underlying 
 * representation is a double-array. Changes to the elements of this double-array 
 * will also result in a change to the point. The double-array can be specified on  
 * construction and can be retrieved lated. 
 * @author R. Fonseca
 */
public class Point implements Serializable{
	private static final long serialVersionUID = -6776647998875885511L;
	
	/** The double-array that holds the value of this point */
	protected final double[] coords;
	
	/** The dimension of this point. */
	protected final int dim;

	public Point(Point p){
		this.dim = p.coords.length;
		this.coords = new double[dim];
		for(int i=0;i<dim;i++) coords[i] = p.coords[i];
	}

	public Point(double[] coords){
		this.coords = coords;
		this.dim = coords.length;
	}
	
	public Point(int dimensions){
		this.coords = new double[dimensions];
		this.dim = dimensions;
	}
	
	/** Get the d'th coordinate. */
	public double get(int d){			return coords[d];		}
	
	/** Get the d'th coordinate. */
	public double getCoord(int d) {		return coords[d];		}
	
	/** Get all coordinates in an array. */
	public double[] getCoords() {		return coords;			}
	
	/** 
	 * Return the dimensionality of this point. If it is a point in xyz-space 
	 * then the dimensionality is 3. The method getDimensions should not be confused 
	 * with getDimension (in Simplex and geom3d.Point) which is always 0 for points. 
	 */
	public int getDimensions() {			return dim;				}

	/** Set the d'th coordinate. */
	public void set(int d, double v){		coords[d]=v;			}
	
	/** Set the d'th coordinate. */
	public void setCoord(int d, double v){	coords[d]=v;			}
	
	/** Set all coordinates to the value v */
	public void fill(double v){				Arrays.fill(coords, v);	}
	
	/** Set the coordinates to be identical to those in the specified point. If the dimensions of 
	 * <code>p</code> and <code>this</code> are not the same, then the largest possible number of 
	 * coordinates are copied */
	public Point set(Point p){		
		for(int i=0;i<Math.min(dim, p.dim);i++) coords[i] = p.coords[i];
		return this;
	}
	
	/** Set the coordinates to be identical to those in the specified point. If the dimensions of 
	 * <code>p</code> and <code>this</code> are not the same, then the largest possible number of 
	 * coordinates are copied. */
	public Point setCoord(Point p){	
		for(int i=0;i<Math.min(dim, p.dim);i++) coords[i] = p.coords[i];
		return this;
	}
	
	/** Add the specified vector and return the result without modifying this point. */
	public Point add(Vector v) {
		return clone().addThis(v);
	}
	
	/** Add the specified vector to this point. */
	public Point addThis(Vector v){
		for(int i=0;i<Math.min(dim, v.dim);i++) coords[i] += v.coords[i];
		return this;
	}
	
	/** Scale this point by the specified value */
	public Point multiplyThis(double s) {
		for(int d=0;d<dim;d++) 
			coords[d]*=s;
		return this;
	}
	
/** Returns squared distance to the point p */
	public double distanceSquared(Point p){
		assert dim==p.dim: "Dimensions don't match";
		double sum = 0;
		for(int d=0;d<dim;d++) {
			double delta = coords[d]-p.coords[d];
			sum+=delta*delta;
		}
		return sum;
	}
	
	/** Returns squared distance to the origo */
	public double distanceSquared(){
		double sum = 0;
		for(int d=0;d<dim;d++) { sum+=coords[d]*coords[d]; }
		return sum;
	}

	/** Returns distance to the point p */
	public double distance(Point p){ return Math.sqrt(distanceSquared(p)); }

	public double dot(Point p){
		double sum = 0;
		for(int d=0;d<dim;d++) {
			sum+=coords[d]*p.coords[d];
		}
		return sum;
	}
	
	/** Returns distance to the origo */
	public double distance() { return Math.sqrt(distanceSquared()); }

	
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
	
	/** Returns the squared distance to the origo of this point projected onto xy plane */
	public double distanceSquaredXY() { return coords[0]*coords[0] + coords[1]*coords[1]; }

	/** Returns the distance to the origo of this point projected onto xy plane */
	public double distanceXY() { return Math.sqrt(distanceSquaredXY()); }

	/** Returns the sinus of the polar angle of this point projected onto xy plane */
	public double polarAngleSinXY() { return coords[1]/(Math.sqrt(coords[0]*coords[0] + coords[1]*coords[1])); }

	/** Returns the cosinus of the polar angle of this point projected onto xy plane */
	public double polarAngleCosXY() { return coords[0]/(Math.sqrt(coords[0]*coords[0] + coords[1]*coords[1])); }

	/** Returns polar angle of this point projected onto xy plane */
	public double polarAngleXY() {
		double angle = Math.acos(polarAngleCosXY());
		if (coords[1] < 0) angle = 2*Math.PI - angle;
		return angle;
	}

	public static Point getRandomPoint(int d, double minCoord, double maxCoord){
		Point ret = new Point(d);
		for(int i=0;i<d;i++){
			ret.coords[i] = Randomization.randBetween(minCoord, maxCoord);
		}
		return ret;
	}
	
	public Vector vectorTo(Point p1) {
		double[] newCoords = new double[dim];
		for(int d=0;d<dim;d++) newCoords[d] = p1.coords[d]-this.coords[d];
		return new Vector(newCoords);
	}

	public Vector toVector(){
		return new Vector(coords);
	}
	
	public Point clone(){
		double[] newArr = new double[dim];
		for(int d=0;d<dim;d++) newArr[d] = coords[d];
		return new Point(newArr);
	}
	
	
	public String toString(){
		return toString(2);
	}
	public String toString(int dec){
		StringBuilder ret = new StringBuilder();
		ret.append("Point[");
		for(int d=0;d<dim-1;d++) ret.append(String.format("%."+dec+"f, ", coords[d]));
		if(dim>0)
			ret.append(String.format("%."+dec+"f]", coords[dim-1]));
		return ret.toString();
	}
	public void toConsole() {
		toConsole(2);
	}

	public void toConsole(int dec) {
		System.out.println(toString(dec));
	}



}
