package ProGAL.geom2d;

public class Vector extends ProGAL.geomNd.Vector{
//	private double x,y;
		
	/** Creates a vector with the given coordinates */
	public Vector(double x, double y) { super(new double[]{x,y}); }
	/** Creates a vector with the given coordinates given in an array */
	public Vector(double[] coords) { super(coords); }
	/** Creates a vector from point p to point q */
	public Vector(Point p, Point q) { this(q.x()-p.x(), q.y() - p.y()); }
	/** Creates a position vector (from origo to p) */
	public Vector(Point p) { this(p.x(), p.y()); }
	
	public double x() { return coords[0]; }
	public double y() { return coords[1]; }
	public double[] getCoords(){ return coords; }
		
	public double getSquaredLength() { return coords[0]*coords[0] + coords[1]*coords[1]; }
	public double length() { return Math.sqrt(getSquaredLength()); }
	
	public double getSlope() {
		if (coords[0] == 0) return 9999.0; else return coords[1]/coords[0];
	}

	public void negative() { scale(-1); }
	public Vector scale(double a) { coords[0]*=a; coords[1]*=a; return this; }
	
	/** Perform a counter-clock-wise rotation by a radians. */
	public Vector rotateThis(double a) {
		double cosA = Math.cos(a);
		double sinA = Math.sin(a);
		double xOld = get(0);
		coords[0] = cosA*coords[0] - sinA*coords[1];
		coords[1] = sinA*xOld + cosA*coords[1];
		return this;
	}
	public Vector add(Vector v) { return new Vector(coords[0]+v.coords[0],coords[1]+v.coords[1]); }
	public Vector addThis(Vector v) { coords[0]+=v.coords[0]; coords[1]+=v.coords[1]; return this; }
	public Vector addThis(Point p) { coords[0]+=p.x(); coords[1]+=p.y(); return this; }
//	public void subtractThis(Vector v) { x()=x()-v.x(); y=y-v.y; }
//	
//	public Vector reverse() { return multiply(-1); }
	public Vector multiply(double a) { return new Vector(a*coords[0], a*coords[1]); }
	public Vector multiplyThis(double a) { coords[0]*=a;coords[1]*=a;return this; }
	public Vector normalize() { return multiply(1/length()); }
	public Vector normalizeThis() { return multiply(1/length()); }
	public Vector createRotatedVector(double a) {
		return new Vector(Math.cos(a)*coords[0] - Math.sin(a)*coords[1], Math.sin(a)*coords[0] + Math.cos(a)*coords[1]);
	}
	public Vector rotate90() { return new Vector(-coords[1],coords[0]); }
	public Vector rotate90This() { double tmp = coords[0]; coords[0] = -coords[1]; coords[1]=tmp; return this; }
	public Point toPoint() { return new Point(coords[0],coords[1]); }
	
	
	public static Vector createSum(Vector  u, Vector v) { return new Vector(u.coords[0]+v.coords[0], u.coords[1]+v.coords[1]); }
	public static Vector createDiff(Vector u, Vector v) { return new Vector(u.coords[0]-v.coords[0], u.coords[1]-v.coords[1]); }
	
	public static double crossProduct(Vector u, Vector v) { return u.coords[0]*v.coords[1] - u.coords[1]*v.coords[1]; }
	public static double dotProduct(Vector u, Vector v) { return u.coords[0]*v.coords[0] + u.coords[1]*v.coords[1]; }
	
	public static boolean leftTurn(Vector u, Vector v)  { return Vector.crossProduct(u, v) >  0.0; }
	public static boolean rightTurn(Vector u, Vector v) { return Vector.crossProduct(u, v) <= 0.0; }

	
	@Override
	public String toString() { return toString(2); }
	public String toString(int dec) { return String.format("Vector[%"+dec+"f,%f"+dec+"]",coords[0],coords[1]); }
	public void toConsole() { toConsole(2); }
	public void toConsole(int dec) { System.out.println(toString(dec));
	
	}
	public Vector scaleToLength(double l) {
		double l2 = this.length();
		return new Vector(l*coords[0]/l2, l*coords[1]/l2);
	}
	public double dot(Vector v) {
		return coords[0]*v.coords[0] + coords[1]*v.coords[1];
	}

	public Vector clone(){
		return new Vector(coords[0],coords[1]);
	}
}

