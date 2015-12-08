package ProGAL.geomNd;

import ProGAL.math.Constants;
import ProGAL.math.Matrix;
import ProGAL.math.Randomization;

/** 
 *  A vector in metric space represented with double precision.
 *  @todo cache the length so several calls to getLengthSquared and getLength
 *  takes less time. This can be relevant e.g. in Line3d for multiple 
 *  projections, but can also cause serious problems if anyone chooses to extend the 
 *  Vector object.
 */
public class Vector {
	protected double[] coords;
	protected final int dim;
	
	/** Construct a vector pointing from origo to p. */
	public Vector(Point p) { 
		dim = p.dim;
		coords = new double[dim];
		for(int d=0;d<dim;d++) coords[d] = p.coords[d];
	}
	
	/** Constructs a vector between two points p1 and p2 - added by pawel 12-11-2011 */
	public Vector(Point p1, Point p2) {
		dim = p1.dim;
		coords = new double[dim];
		for(int d=0;d<dim;d++) coords[d] = p2.coords[d] - p1.coords[d];
	}
	
	/** Construct a vector that is a clone of v. */
	public Vector(Vector v) { 
		dim = v.dim;
		coords = new double[dim];
		for(int d=0;d<dim;d++) coords[d] = v.coords[d];
	}	
	
	/** Construct a vector with the specified coordinates */
	public Vector(double[] coords) {
		dim = coords.length;
		this.coords = coords;
	}

	public Vector(int dim) {
		this.dim = dim;
		this.coords = new double[dim];
	}

	/** Get the i'th coordinate. */
	public double getCoord(int i) { return coords[i]; 	}
	
	/** Get the i'th coordinate. */
	public double get(int i) { return getCoord(i); }

	/** Get the dimensionality of this vector */
	public int getDimensions() {	return coords.length;	}

	/** Set the i'th coordinate to v */
	public void setCoord(int i, double v){ coords[i] = v; }

	/** Set the coordinates */
	public void setCoords(double[] coords){
//		if(this.coords.length!=coords.length) throw new Error("Wrong dimension of coordinates");
		for(int i=0;i<Math.min(coords.length, this.coords.length);i++) this.coords[i] = coords[i]; 
	}
	
	/** Set the i'th coordinate to v */
	public void set(int i, double v){ setCoord(i,v); }

	/** Get the squared length of this vector. */
	public double getLengthSquared() { 
		double sum = 0;
		for(int d=0;d<dim;d++) sum+=coords[d]*coords[d];
		return sum;
	}

	/** Get the length of this vector. */
	public double length() {return Math.sqrt(getLengthSquared()); }

	
	/** Return true if the length of this vector is 0. */  // rettet af Pawel 12-11-2011
	public boolean isZeroVector() {
		for(int d=0;d<dim;d++) if(Math.abs(coords[d])>Constants.EPSILON) return false;
		return true;
	}
	
	/** Get the dot-product of this vector and v. */
	public double dot(Vector v){
		double sum = 0;
		for(int d=0;d<dim;d++) sum+=coords[d]*v.coords[d];
		return sum; 
	}

	/** Get the angle between this vector and v. */
	public double angle(Vector v) {
		return Math.acos(Math.min(  1, this.dot(v)/Math.sqrt(this.getLengthSquared()*v.getLengthSquared())  ));
	}

	/** Add v to this vector and return the result (without changing this object). */
	public Vector add(Vector v){ 
		double[] sum = new double[dim];
		for(int d=0;d<dim;d++) sum[d] = coords[d]+v.coords[d];
		return new Vector(sum); 
	}
	
	/** Add v to this vector and return the result (changing this object). */ 
	public Vector addThis(Vector v){ 
		for(int d=0;d<dim;d++) coords[d] += v.coords[d];
		return this; 
	}

	/** Subtract v from this vector and return the result (without changing this object). */
	public Vector subtract(Vector v){ 
		double[] sum = new double[dim];
		for(int d=0;d<dim;d++) sum[d] = coords[d]-v.coords[d];
		return new Vector(sum); 
	}
	
	/** Subtract v from this vector and return the result (changing this object). */ 
	public Vector subtractThis(Vector v){ 
		for(int d=0;d<dim;d++) coords[d] -= v.coords[d];
		return this; 
	}
	
	/** Multiply this vector by s and return the result (without changing this object). */
	public Vector multiply(double s){ 
		double[] ret = new double[dim];
		for(int d=0;d<dim;d++) ret[d] = coords[d]*s;
		return new Vector(ret); 
	}
	
	/** Multiply this vector by s and return the result (changing this object). */
	public Vector multiplyThis(double s){ 
		for(int d=0;d<dim;d++) coords[d]*=s;
		return this; 
	}
	
	/** Divide this vector by s and return the result (without changing this object). */
	public Vector divide(double s){ 
		double[] ret = new double[dim];
		for(int d=0;d<dim;d++) ret[d] = coords[d]/s;
		return new Vector(ret); 
	}
	
	/** Divide this vector by s and return the result (changing this object). */
	public Vector divideThis(double s){
		for(int d=0;d<dim;d++) coords[d]/=s;
		return this; 
	}
	
	public Matrix outerProduct(Vector v) {
		Matrix M = new Matrix(dim, dim);
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				M.set(i, j, coords[i]*coords[j]);
		return M;
	}

	/** Normalize this vector and return the result (without changing this object). */
	public Vector normalize(){ return this.multiply(1/length()); }
	
	/** Normalize this vector and return the result (changing this object). */
	public Vector normalizeThis(){ return this.multiplyThis(1/length()); }

	/** Scale this vector to a certain length (returns new object and does not change this object). */
	public Vector scaleToLength(double length) { return multiply(length/length()); }
	
	/** Scale this vector to a certain length (changes this object). */
	public Vector scaleToLengthThis(double length) { return multiplyThis(length/length()); }
	
	/** Convert this vector to a point. */
	public Point toPoint() { return new Point(clone().coords); }
	
	/** Returns a string-representation of this vector formatted with two decimals precision. */
	public String toString() { return toString(2); }
	
	/** Returns a string-representation of this vector formatted with <code>dec</code> decimals precision. */
	public String toString(int dec) {
		StringBuilder sb = new StringBuilder();
		sb.append("Vector[");
		for(int d=0;d<dim-1;d++) sb.append(String.format("%."+dec+"f, ",coords[d]));
		sb.append(String.format("%."+dec+"f]",coords[dim-1]));
		return sb.toString();
	}	
	
	/** Writes this vector to <code>System.out</code>. */
	public void toConsole() { toConsole(2); }
	
	/** Writes this vector to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec) { System.out.println(toString(dec)); }

	/** Return true if this vector equals v. */
	public boolean equals(Vector v){
		for(int d=0;d<dim;d++) 
			if(Math.abs(coords[d]-v.coords[d])>Constants.EPSILON) return false;
		return true;
	}
	public boolean equals(Object v){
		if(v instanceof Vector) return equals((Vector)v);
		return false;
	}

	/** Create a clone of this vector. */
	public Vector clone(){ 
		double[] newCoords = new double[dim];
		for(int d=0;d<dim;d++) newCoords[d] = coords[d];
		return new Vector(newCoords); 
	}

	///////// Static methods and fields ////////

	/** Get the angle between vector u and v. */
	public static double getAngle(Vector u, Vector v) { return u.angle(v); }

	/** Create a vector of length <code>len</code> pointing a random direction */
	public static Vector randomVector(int dim, double len) {
		Vector ret = new Vector(dim);
		do{
		for(int d=0;d<dim;d++)
			ret.coords[d] = Randomization.randBetween(-1.0, 1.0);
		}while(ret.length()>1.0);
		
		ret.multiplyThis(len/ret.length());
		return ret;
	}

	

}

