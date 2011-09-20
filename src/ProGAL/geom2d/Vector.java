package ProGAL.geom2d;

public class Vector {
	protected double x,y;
		
	public Vector(double x, double y) { this.x = x; this.y = y; }
		
	public double x() { return x; }
	public double y() { return y; }
		
	public double getSquaredLength() { return x*x + y*y; }
		
	public double length() { return Math.sqrt(getSquaredLength()); }
	
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
	public Vector add(Vector v) { return new Vector(x+v.x,y+v.y); }
	public Vector addThis(Vector v) { x=x+v.x; y=y+v.y; return this; }
	public Vector addThis(Point p) { x=x+p.x(); y=y+p.y(); return this; }
	public void subtractThis(Vector v) { x=x-v.x; y=y-v.y; }
	
	public Vector reverse() { return multiply(-1); }
	public Vector multiply(double a) { return new Vector(a*x, a*y); }
	public Vector multiplyThis(double a) { x=a*x;y=a*y;return this; }
	public Vector normalize() { return multiply(1/length()); }
	public Vector normalizeThis() { return multiply(1/length()); }
	public Vector createRotatedVector(double a) {
		return new Vector(Math.cos(a)*x - Math.sin(a)*y, Math.sin(a)*x + Math.cos(a)*y);
	}
	public Vector rotate90() { return new Vector(-y,x); }
	public Vector rotate90This() { double tmp = x; x = -y; y=tmp; return this; }
	public Point toPoint() { return new Point(x,y); }
	
	public static Vector createSum(Vector  u, Vector v) { return new Vector(u.x+v.x, u.y+v.y); }
	public static Vector createDiff(Vector u, Vector v) { return new Vector(u.x-v.x, u.y-v.y); }
	
	public static double crossProduct(Vector u, Vector v) { return u.x*v.y - u.y*v.x; }
	public static double dotProduct(Vector u, Vector v) { return u.x*v.x + u.y*v.y; }
	
	public static boolean leftTurn(Vector u, Vector v)  { return Vector.crossProduct(u, v) >  0.0; }
	public static boolean rightTurn(Vector u, Vector v) { return Vector.crossProduct(u, v) <= 0.0; }

	
	@Override
	public String toString() { return toString(2); }
	public String toString(int dec) { return String.format("Vector[%"+dec+"f,%f"+dec+"]",x,y); }
	public void toConsole() { toConsole(2); }
	public void toConsole(int dec) { System.out.println(toString(dec));
	
	}
	public Vector scaleToLength(double l) {
		double l2 = this.length();
		return new Vector(l*x/l2, l*y/l2);
	}
	public double dot(Vector v) {
		return x*v.x + y*v.y;
	}

}

