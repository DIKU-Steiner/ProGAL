package ProGAL.geom2d;

import java.awt.Color;

import ProGAL.geom2d.viewer.J2DScene;


public class LineSegment implements Shape{
	protected Point a,b;
	
	public LineSegment(Point a, Point b){
		this.a = a;
		this.b = b;
	}

	public Point getA(){ return a; }
	public Point getB(){ return b; }
	
	public void setA(Point p) { a = p; }
	public void setB(Point p) { b = p; }
	
	/** returns the length of the line segment */
	public double getLength() { return a.distance(b); }
	
	/** returns the squared length of the segment */
	public double getSquaredLength() { 
		double bax = b.x() - a.x();
		double bay = b.y() - a.y();
		return bax*bax + bay*bay;
	}
	
	public LineSegment clone() { return new LineSegment(a.clone(), b.clone()); }

	public Point getCenter() { return Point.midPoint(a, b); }
	
	/** swaps the end-points of this segment */
	public LineSegment reverse() {
		Point c = a;
		a = b;
		b = c;
		return this;
	}
	
	/** returns TRUE if this line segment properly intersects line segment s. */
	public boolean intersects(LineSegment s) {
		return (((Point.leftTurn(a,   b,   s.a) && Point.leftTurn(b,   a,   s.b)) || 
				 (Point.leftTurn(a,   b,   s.b) && Point.leftTurn(b,   a,   s.a))) &&
				((Point.leftTurn(s.a, s.b, a)   && Point.leftTurn(s.b, s.a, b)) || 
				 (Point.leftTurn(s.a, s.b, b)   && Point.leftTurn(s.b, s.a, a))));
	}
	public String toString() { return "[" + a.toString() + b.toString() + "]"; }
	public void toConsole() { System.out.println("[" + a.toString() + b.toString() + "]"); }
	
	public void toScene(J2DScene scene) { scene.addShape(this, Color.black); }
	public void toScene(J2DScene scene, Color clr) { scene.addShape(this, clr); }
	public void toScene(J2DScene scene, Color clr, double width) { scene.addShape(this, clr, width); }
	

	
	public double distance(Point p){
		Vector v = a.vectorTo(b);
		Vector vP = a.vectorTo(p);
		double t = v.dot(vP)/v.getSquaredLength();
		return a.add(v.multiplyThis(Math.min(Math.max(0,t), 1))).distance(p);
	}

	@Override
	public boolean contains(Point p) {
		return false;
	}
}
