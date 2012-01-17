package ProGAL.geom2d.delaunay;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import ProGAL.geom2d.Circle;
import ProGAL.geom2d.LineSegment;
import ProGAL.geom2d.Point;
import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.math.Randomization;

/** 
 * Not working yet
 * @author ras
 */
public class DelaunayTriangulation {
	final List<Point> points;
	final List<Triangle> triangles;
	final Point[] bigPoints;
	
	public DelaunayTriangulation(List<Point> points){
		scene = J2DScene.createJ2DSceneInFrame();
		this.points = new ArrayList<Point>(points.size());
		this.triangles = new ArrayList<Triangle>(points.size());
		bigPoints = new Point[3];
		bigPoints[0] = new Point(-15,-15); 
		bigPoints[1] = new Point( 15,-15);
		bigPoints[2] = new Point(    0, 15);
		triangles.add(new Triangle(bigPoints[0], bigPoints[1], bigPoints[2]));
		for(Point p: points) 
			addPoint(p);
	}
	
	public void addPoint(Point p){
		System.out.println(p);
		points.add(p);
		Triangle old = locate(p);
		split(old, p);
		draw();
		fixDelaunay();
		draw();
		System.out.println();
	}
	
	public Triangle locate(Point point){
		for(Triangle t: triangles){
			if(t.circumsphere.contains(point) && t.contains(point)){
				return t;
			}
		}
		return null;
	}
	private void split(Triangle t, Point p){
		triangles.remove(t);
		Triangle t0 = new Triangle(p, t.getCorner(1), t.getCorner(2));
		Triangle t1 = new Triangle(p, t.getCorner(2), t.getCorner(0));
		Triangle t2 = new Triangle(p, t.getCorner(0), t.getCorner(1));
		t0.setNeighbors(t.getNeighbor(0), t2, t1);
		t1.setNeighbors(t.getNeighbor(1), t2, t0);
		t2.setNeighbors(t.getNeighbor(2), t0, t1);
		triangles.add(t0);
		triangles.add(t1);
		triangles.add(t2);
	}
	
	/** 
	 * Ensure that the delaunay criterion is maintained for the last three triangles 
	 * and fix it if not.
	 */
	private void fixDelaunay(){
		for(int i=triangles.size()-3;i<triangles.size();i++){
			Triangle t = triangles.get(i);
			for(int n=0;n<3;n++) fixDelaunay(t,n);
		}
	}
	private void fixDelaunay(Triangle t, int n){
		if(t.neighbors[n]!=null && t.circumsphere.contains(t.getOpposite(n))){
			flip(t,n);
			draw();
			fixDelaunay(t,1);
			fixDelaunay(t,2);
			fixDelaunay(t.neighbors[0],1);
			fixDelaunay(t.neighbors[0],2);
		}
	}
	private void flip(Triangle t, int n){
		System.out.println("Flipping "+n);
		Triangle t0 = t;
		Triangle t1 = t.getNeighbor(n);
		System.out.println("t0: "+t0);
		System.out.println("t1: "+t1);
		
		Triangle n0 = t0.getNeighbor( (n+1)%3 );
		Triangle n1 = t0.getNeighbor( (n+2)%3 );
		
		int m = 0; 
		while(t1.getNeighbor(m)!=t0) m++;
		Triangle n2 = t1.getNeighbor( (m+1)%3 );
		Triangle n3 = t1.getNeighbor( (m+2)%3 );
		System.out.println("n0: "+n0);
		System.out.println("n1: "+n1);
		System.out.println("n2: "+n2);
		System.out.println("n3: "+n3);
		
		Point p0 = t0.getCorner(n);
		Point p1 = t0.getCorner( (n+1)%3 );
		Point p2 = t1.getCorner(m);
		Point p3 = t1.getCorner( (m+1)%3 );
		
		
		t0.setCorner(p3, 0);t0.neighbors[0] = t1;
		t0.setCorner(p0, 1);t0.neighbors[1] = n3;
		t0.setCorner(p2, 2);t0.neighbors[2] = n0;
		t1.setCorner(p1, 0);t1.neighbors[0] = t0;
		t1.setCorner(p2, 1);t1.neighbors[1] = n1;
		t1.setCorner(p0, 2);t1.neighbors[2] = n2;
		
		if(n1!=null){
			m = 0;while(n1.neighbors[m]!=t0) m++;
			n1.neighbors[m] = t1;
		}
		if(n3!=null){
			m = 0;while(n3.neighbors[m]!=t1) m++;
			n3.neighbors[m] = t0;
		}
		System.out.println("After flipping: ");
		System.out.println("t0: "+t0);
		System.out.println("t1: "+t1);
		System.out.println("n0: "+n0);
		System.out.println("n1: "+n1);
		System.out.println("n2: "+n2);
		System.out.println("n3: "+n3);
	}
	
	public boolean inCH(Triangle t){
		for(Point bp: bigPoints){
			for(int c=0;c<3;c++) if(t.getCorner(c)==bp) return false;
		}
		return true;
	}
	
	public void draw() {			
		scene.removeAllShapes();
		for(Triangle t: triangles){
			scene.addShape(new LineSegment(t.getCorner(0), t.getCorner(1)), Color.GRAY);
			scene.addShape(new LineSegment(t.getCorner(1), t.getCorner(2)), Color.GRAY);
			scene.addShape(new LineSegment(t.getCorner(2), t.getCorner(0)), Color.GRAY);
		}

	}
	
	J2DScene scene;

	public static void main(String[] args){
		Randomization.seed(1);
		List<Point> points = new ArrayList<Point>();
		points.add(new Point(Randomization.randBetween(-4.0, 4.0), Randomization.randBetween(-3.0, 3.0)));
		points.add(new Point(Randomization.randBetween(-4.0, 4.0), Randomization.randBetween(-3.0, 3.0)));
		points.add(new Point(Randomization.randBetween(-4.0, 4.0), Randomization.randBetween(-3.0, 3.0)));
		points.add(new Point(Randomization.randBetween(-4.0, 4.0), Randomization.randBetween(-3.0, 3.0)));
		DelaunayTriangulation dt = new DelaunayTriangulation(points);
		dt.scene.addShape(new LineSegment(new Point(0,0), new Point(0,1)), Color.black);
		dt.scene.addShape(new LineSegment(new Point(0,0), new Point(1,0)), Color.black);
		for(Point p: points) dt.scene.addShape(new Circle(p,0.03), Color.BLACK);
		dt.scene.addShape(new Circle(new Point(Randomization.randBetween(-4.0, 4.0), Randomization.randBetween(-3.0, 3.0)),0.03), Color.BLACK, 15, true);
		for(Triangle t: dt.triangles){
//			if(t.inCH()){
				dt.scene.addShape(new LineSegment(t.getCorner(0), t.getCorner(1)), Color.GRAY);
				dt.scene.addShape(new LineSegment(t.getCorner(1), t.getCorner(2)), Color.GRAY);
				dt.scene.addShape(new LineSegment(t.getCorner(2), t.getCorner(0)), Color.GRAY);
//			}
		}
//		scene.centerCamera();
		dt.scene.autoZoom();
	}
}
