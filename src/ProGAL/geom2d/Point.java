package ProGAL.geom2d;

import ProGAL.math.Constants;

public class Point extends ProGAL.geomNd.Point {
	private static final long serialVersionUID = 8095991200265432551L;
	
	public Point() { this(0,0); }
	
	public Point(double x, double y) { super(new double[]{x,y}); }

	public double x(){ return coords[0]; }
	public double y(){ return coords[1]; }
	
	public Point clone(){
		return new Point(coords[0], coords[1]);
	}
	
	public Vector vectorTo(Point p){
		return new Vector(p.coords[0]-coords[0], p.coords[1]-coords[1]);
	}

	public Point add(Vector v){ return new Point(x()+v.x(), y()+v.y()); }
	public Point addThis(Vector v){ coords[0] = x()+v.x(); coords[1] = y()+v.y(); return this; }
	public Point subtract(Vector v){ return new Point(x()-v.x(), y()-v.y()); }
	public Point subtractThis(Vector v){ coords[0] = x()-v.x(); coords[1] = y()-v.y(); return this; }
	
	/** Returns the midpoint of two points. */
	public static Point midPoint(Point p, Point q) { return new Point((p.coords[0] + q.coords[0])/2, (p.coords[1] + q.coords[1])/2); }
	
	/*
	 * creates a bisector between points p and q
	 */
	public static Line getBisector(Point p, Point q) {
		if (!p.equals(q)) return new Line(midPoint(p,q), p.vectorTo(q)); 
		return null;
	}
	
	
	public boolean equals(Object o){
		if(o instanceof Point) return equals((Point)o);
		return false;
	}
	
	public boolean equals(Point p){
		return Math.abs(coords[0]-p.coords[0])<Constants.EPSILON && Math.abs(coords[0]-p.coords[0])<Constants.EPSILON;
	}
	
	
}

