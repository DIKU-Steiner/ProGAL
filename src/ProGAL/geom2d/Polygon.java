package ProGAL.geom2d;

import java.awt.Color;
import java.util.ArrayList;

import ProGAL.geom2d.viewer.J2DScene;

public class Polygon extends ArrayList<Point> {
	private static final long serialVersionUID = 1L;

	public Polygon() { super(); }

	public Polygon(Point p0, Point p1, Point p2) { 
		add(p0);
		if (Point.leftTurn(p0,  p1, p2)) { add(p1); add(p2); } else { add(p2); add(p1); }
	}
		
	public Polygon(PointSet points) { for (Point p : points) add(p); }
	
	public Point getCorner(int i) { return get(i); }
	
	public void setCorner(Point p, int i) { set(i, p); }
		
	/** inserts point p into the polygon after the corner with specified index */
	public void insertAfter(Point p, int index) { this.add(index, p); }
	
	/** deletes last corner of the polygon */
	public void deleteLast() { remove(size()-1); }
	
	/** returns the index of the leftmost point (in case of ties, index of the bottommost one is returned) */
	public int leftExtremePointIndx() {
		int k = 0;
		for (int i = 1; i < size(); i++)
			if ((get(i).x() < get(k).x()) || ((get(i).x() == get(k).x()) && (get(k).y() > get(i).y()))) k = i;
		return k;
	}

	/** returns the index of the rightmost point (in case of ties, index of the topmost one is returned) */
	public int rightExtremePointIndx() {
		int k = 0;
		for (int i = 1; i < size(); i++)
			if ((get(i).x() > get(k).x()) || ((get(i).x() == get(k).x()) && (get(k).y() < get(i).y()))) k = i;
		return k;
	}

	/** first corner of the polygon is moved by shiftStep positions */
	public void shift(int shiftStep) {
		Point[] front = new Point[shiftStep];
		for (int i = 0; i < shiftStep; i++) front[i] = get(i);
		for (int i = shiftStep; i < size(); i++) set(i-shiftStep, get(i));
		for (int i = 0; i < shiftStep; i++) set(size()-shiftStep+i, front[i]);
	}
	
	/** draws the polygon */
	public void draw(J2DScene scene, Color clr) {
		for (int i = 1; i < size(); i++) scene.addShape(new LineSegment(get(i-1), get(i)), clr);
		scene.addShape(new LineSegment(get(size()-1), get(0)), clr);
	}
	public void draw(J2DScene scene) { draw(scene, Color.black); }
	
	/** returns convex hull of a simple polygon in O(n) time */
	public ConvexPolygon getConvexPolygon() { return new ConvexPolygon(this); }
	
	public static void main(String[] args) {
		Polygon pol = new Polygon();
		pol.add(new Point(0,0));
		pol.add(new Point(1,-1));
		pol.add(new Point(2,-1));		
		pol.add(new Point(3,2));
		pol.add(new Point(2,2));
		pol.add(new Point(1, 1.5));
		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		pol.draw(scene);
		ConvexPolygon cPol = pol.getConvexPolygon();
		cPol.draw(scene);
		System.out.println(cPol.farthestVertex(0, 3));
		cPol.getDiameter();
	}
	

}
