package ProGAL.geom2d;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import ProGAL.geom2d.viewer.J2DScene;

public class Polygon extends ArrayList<Point> implements Shape {
	private static final long serialVersionUID = 1L;
	

	public Polygon(){	super();	}
	
	public Polygon(List<Point> corners){
		super(corners);
	}

	public Polygon(Point p0, Point p1, Point p2) {
		add(p0);
		if (Point.leftTurn(p0,  p1, p2)) { add(p1); add(p2); } 
		else { add(p2); add(p1); }
	}
		
	public Polygon(PointSet points) { for (Point p : points) add(p); }
	
	public Polygon(Point[] points) {
		for(Point p: points) add(p);
	}

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
	
		public boolean isConvex(){
		if(size()<4) return true;
		
		Point p0=get(0);
		Point p1=get(1);
		Point p2=get(2);
		boolean ccw = Point.leftTurn(p0,p1,p2);
		
		for(int i=1;i<size();i++){
			p0=p1;
			p1=p2;
			p2=get((i+2)%size());
			if(ccw!=Point.leftTurn(p0, p1, p2)) return false;
		}
		return true;
	}
	
	@Override
	public Point getCenter() {
		Vector v = new Vector(0,0);
		for(Point p: this) v.addThis(p);
		return new Point(v.x()/size(), v.y()/size());
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
//		pol.add(new Point(0,0));
//		pol.add(new Point(1,-1));
//		pol.add(new Point(2,-1));		
//		pol.add(new Point(3,2));
//		pol.add(new Point(2,2));
//		pol.add(new Point(1, 1.5));
//		J2DScene scene = J2DScene.createJ2DSceneInFrame();
//		pol.draw(scene);
//		ConvexPolygon cPol = pol.getConvexPolygon();
//		cPol.draw(scene);
//		System.out.println(cPol.farthestVertex(0, 3));
//		cPol.getDiameter();
		
		pol.add(new Point(0,0));
		pol.add(new Point(0,1));
		pol.add(new Point(2,0));
		pol.add(new Point(2,1));
		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		scene.addShape(pol, Color.BLUE, 0, true);
	}

	@Override
	public boolean contains(Point p) {
		//From http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
		boolean result = false;
		for(int i=0;i<size();i++){
			Point p0 = get(i);
			Point p1 = get( (i+1)%size() );
			if( (p0.y() > p.y()) != (p1.y() > p.y()) && (p.x() < (p1.x() - p0.x()) * (p.y() - p0.y()) / (p1.y()-p0.y()) + p0.x()) ) 
				result = !result;
		}
		
		return result;
	}
	
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append("Polygon[");
		for(Point p: this) {sb.append(p.toString());sb.append(","); }
		sb.deleteCharAt(sb.length()-1);
		sb.append("]");
		return sb.toString();
	}

}
