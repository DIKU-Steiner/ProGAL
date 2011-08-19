package ProGAL.geom2d;

public class Vector {
	protected double x,y;
		
	public Vector(double x, double y) { this.x = x; this.y = y; }
	public Vector(Point p, Point q) { x = q.x - p.x; y = q.y - p.y; }
	public Vector(Point p) { x = p.x; y = p.y; }
		
	public double getX() { return x; }
	public double getY() { return y; }
		
	public double getSquaredLength() { return x*x + y*y; }
		
	public double getLength() { return Math.sqrt(getSquaredLength()); }
	
	public double getSlope() {
		if (x == 0) return 9999.0; else return y/x;
	}

	public void negative() { scale(-1); }
	public Vector scale(double a) { x=a*x; y=a*y; return this; }
	public void rotate(double a) {
		double cosA = Math.cos(a);
		double sinA = Math.sin(a);
		double xOld = x;
		x = cosA*x - sinA*y;
		y = sinA*xOld + cosA*y;
	}
	public void rotate90() { 
		double xOld = x;
		x = -y; 
		y = xOld;
	}
	public void makeUnitVector() { scale(1/getLength()); }
	public Vector addThis(Vector v) { x=x+v.x; y=y+v.y; return this; }
	public Vector addThis(Point p) { x=x+p.x; y=y+p.y; return this; }
	public void subtractVector(Vector v) { x=x-v.x; y=y-v.y; }
	
	public Vector createNegativeVector2d() { return createScaledVector2d(-1); }
	public Vector createScaledVector2d(double a) { return new Vector(a*x, a*y); }
	public Vector createUnitVector() { return createScaledVector2d(1/getLength()); }
	public Vector createRotatedVector(double a) {
		return new Vector(Math.cos(a)*x - Math.sin(a)*y, Math.sin(a)*x + Math.cos(a)*y);
	}
	public Vector createRotatedVector90() { return new Vector(-y,x); }
	public Point toPoint() { return new Point(x,y); }
	
	public static Vector createSum(Vector  u, Vector v) { return new Vector(u.x+v.x, u.y+v.y); }
	public static Vector createDiff(Vector u, Vector v) { return new Vector(u.x-v.x, u.y-v.y); }
	
	public static double crossProduct(Vector u, Vector v) { return u.x*v.y - u.y*v.x; }
	public static double crossProduct(Point u, Point v) { return u.x*v.y - u.y*v.x; }
	public static double dotProduct(Vector u, Vector v) { return u.x*v.x + u.y*v.y; }
	public static double dotProduct(Vector u, Point p)  { return Vector.dotProduct(u,new Vector(p)); }
	public static double dotProduct(Point  p, Vector u) { return Vector.dotProduct(u, p); }
	public static double dotProduct(Point  p, Point q)  { return Vector.dotProduct(new Vector(p), new Vector(q)); }
	
	public static boolean leftTurn(Vector u, Vector v)  { return Vector.crossProduct(u, v) >  0.0; }
	public static boolean rightTurn(Vector u, Vector v) { return Vector.crossProduct(u, v) <= 0.0; }

	
	@Override
	public String toString() { return toString(2); }
	public String toString(int dec) { return String.format("Vector[%"+dec+"f,%f"+dec+"]",x,y); }
	public void toConsole() { toConsole(2); }
	public void toConsole(int dec) { System.out.println(toString(dec));
	
	}
	public Vector scaleToLength(double l) {
		double l2 = this.getLength();
		return new Vector(l*x/l2, l*y/l2);
	}
	public double dot(Vector v) {
		return x*v.x + y*v.y;
	}

}

