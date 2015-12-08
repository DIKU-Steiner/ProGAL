package ProGAL.geom2d.delaunay;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;

import ProGAL.geom2d.Circle;
import ProGAL.geom2d.Point;
import ProGAL.geom2d.Vector;
import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.geom3d.predicates.ExactJavaPredicates;
import ProGAL.io.IOToolbox;
import ProGAL.math.Randomization;

public class DTNoBigPoints {
	
	
	public static void main(String[] args){
		
		//Boundary problem at internal point
		List<Point> points = new ArrayList<Point>();
		Randomization.seed(2);
//		for(int i=0;i<300;i++) points.add(new Point(Randomization.randBetween(0.0, 1.0), Randomization.randBetween(0.0, 1.0)));
		for(String line: IOToolbox.readFromFile("/Users/rfonseca/Downloads/punkter.txt").split("\n")){
			String[] coords = line.split(" ");
			points.add(new Point(Double.parseDouble(coords[0]), Double.parseDouble(coords[1])));
		}
		Collections.shuffle(points, Randomization.getGenerator());
		while(points.size()>1000) { points.remove(points.size()-1); }
		
		long start = System.nanoTime();
		DTNoBigPoints dt = new DTNoBigPoints(points);
		long end = System.nanoTime();
		System.out.printf("Took %.2fms.\n", (end-start)/1000000.0);
		
		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		float delta = 1.0f/(dt.triangles.size()+1);
		float sum = 0;
		for(Triangle t: dt.triangles){
			Color c = new Color(0.6f+sum*0.1f, 0.8f+sum*0.19f, 0.7f);
			Color c1 = new Color(100, 180, 180);
			Color c2 = new Color(150, 250, 180);
//			Color c = new Color(0.5f-sum*0.2f, sum*0.9f, 0.7f-sum*0.4f);
//			Color c = new Color(150, 250, 180);
			sum+=delta;
//			boolean fill = t.neighbors[0]==null||t.neighbors[1]==null||t.neighbors[2]==null;
//			boolean fill = t.getCircumCircle().getRadius()<2000;
//			boolean fill = false;
			if(t.getCircumCircle().getRadius()<2000)
			{
				scene.addShape(t, c, 0, true);
				scene.addShape(t, Color.GRAY.darker(), 40, false);
			}
		}
		for(Point p: points){
			scene.addShape(new Circle(p, 140.0), Color.RED, 0, true);
			scene.addShape(new Circle(p, 140.0), Color.BLACK, 30, false);
		}
		
		scene.centerCamera();
		scene.autoZoom();
	}

//		J2DScene scene = J2DScene.createJ2DSceneInFrame();

	public final List<Vertex> vertices;
	public final List<Triangle> triangles = new ArrayList<Triangle>();
	
	public DTNoBigPoints(List<Point> points){
//		scene.frame.setSize(500, 800);
		//Sort from left to right
		vertices = new ArrayList<Vertex>();
		for(Point p: points) vertices.add(new Vertex(p)); 
		Collections.sort(vertices, new Comparator<Point>(){
			public int compare(Point p0, Point p1) {
				int ret = Double.compare(p0.x(), p1.x());
				if(ret==0) return Double.compare(p0.y(), p1.y());
				else return ret;
			}});

		int c = 0;
		for(Vertex v: vertices){
			v.id = c++; 
//			scene.addShape(new Circle(v,0.005));
//			scene.addShape(new TextShape(sortedPoints.indexOf(v)+"", v.add(new Vector(0.005, 0)), 0.034));
		}
//		scene.centerCamera();
//		scene.autoZoom();

		//Initializes
		Vertex p0 = vertices.get(0);
		Vertex p1 = vertices.get(1);
		Vertex p2 = vertices.get(2);
//		scene.addShape(new TextShape(""+p0.id, p0, 0.02));
//		scene.addShape(new TextShape(""+p1.id, p1, 0.02));
//		scene.addShape(new TextShape(""+p2.id, p2, 0.02));
		Triangle t = new Triangle(p0,p1,p2);
		p0.first = p0.last = t;
		p1.first = p1.last = t;
		p2.first = p2.last = t;
		triangles.add(t);
//		scene.addShape(t,Color.GRAY, 0.002);
		
		//Iterate
		for(int i=3;i<vertices.size();i++){
			if(i%10000==0) System.out.println(i);
			Vertex p = vertices.get(i);
//			scene.addShape(new Circle(p, 0.001));
//			scene.addShape(new TextShape(""+p.id, p, 0.02));

			Vertex q = vertices.get(i-1);
			t = q.first;
			Vertex r = t.corners[(t.indexOf(q)+1)%3];
			while(Point.rightTurn(q, r, p)){
				Triangle tNew = new Triangle(p,r,q);//TODO optimize the rightturn in triangle constructor
				tNew.neighbors[0] = t;
				t.neighbors[(t.indexOf(q)+2)%3] = tNew;
				if(p.last!=null) {
					tNew.neighbors[1] = p.first;
					p.first.neighbors[ (p.first.indexOf(p)+2)%3 ] = tNew;
					q.last = tNew;
				}else{
					p.last = tNew;
				}
				q.first = r.last = p.first = tNew;
				triangles.add(tNew);
//				scene.addShape(tNew, Color.GRAY, 0.002);
				legalizeEdge(tNew, 0);
				
				q = r;
				t = q.first;
				r = t.corners[(t.indexOf(q)+1)%3];
			}
			
			q = vertices.get(i-1);
			t = q.last;
			r = t.corners[(t.indexOf(q)+2)%3];
			while(Point.leftTurn(q, r, p)){ 
				Triangle tNew = new Triangle(p,q,r);//TODO optimize the rightturn in triangle constructor
				tNew.neighbors[0] = t;
				t.neighbors[(t.indexOf(q)+1)%3] = tNew;
				if(p.last!=null) {
					tNew.neighbors[2] = p.last;
					p.last.neighbors[ (p.last.indexOf(p)+1)%3 ] = tNew;
					q.first = tNew;
				}else{
					p.first = tNew;
				}
				q.last = r.first = p.last = tNew;
				triangles.add(tNew);
//				scene.addShape(tNew, Color.GRAY, 0.002);
				legalizeEdge(tNew, 0);
				
				q = r;
				t = q.last;
				r = t.corners[(t.indexOf(q)+2)%3];
			}
		}
	}
	
	private final ExactJavaPredicates pred = new ExactJavaPredicates();
	
	void legalizeEdge(Triangle t, int e){
		if(t.neighbors[e]==null) return;
		
		Triangle u = t.neighbors[e];
		int f = (u.indexOf(t.corners[(e+1)%3])+1)%3;
//		Circle c = new Circle(t.corners[0], t.corners[1], t.corners[2]);
//		scene.addShape(c, Color.BLUE, 0.001);
		double inc = pred.incircle(t.corners[0].getCoords(),t.corners[1].getCoords(),t.corners[2].getCoords(), u.corners[f].getCoords());
		boolean illegal = inc>0;
		if(illegal){
//			scene.removeShape(c);
			flip(t,e,u,f);
//			scene.repaint();
			legalizeEdge(t, e);
			legalizeEdge(u, (f+2)%3);
		}else{
//			scene.removeShape(c);
		}
		
	}

	void flip(Triangle t, int e, Triangle u, int f){
		Vertex p0 = t.corners[e];
		Vertex p1 = t.corners[(e+1)%3];
		Vertex p2 = u.corners[f];
		Vertex p3 = t.corners[(e+2)%3];
		Triangle n01 = t.neighbors[(e+2)%3];
		Triangle n12 = u.neighbors[(f+1)%3];
		Triangle n23 = u.neighbors[(f+2)%3];
		Triangle n34 = t.neighbors[(e+1)%3];

		t.corners[(e+2)%3] = p2; t.setCorner(p2, (e+2)%3);
		u.corners[(f+2)%3] = p0; u.setCorner(p0, (f+2)%3);
		
		t.neighbors[e] = n12;
		t.neighbors[(e+1)%3] = u;
		t.neighbors[(e+2)%3] = n01;
		u.neighbors[f] = n34;
		u.neighbors[(f+1)%3] = t;
		u.neighbors[(f+2)%3] = n23;

		if(n12!=null) n12.neighbors[(n12.indexOf(p1)+1)%3] = t;
		if(n34!=null) n34.neighbors[(n34.indexOf(p3)+1)%3] = u;

		if(p0.onBoundary() && p0.last==t) 	p0.last = u;
		if(p1.onBoundary() && p1.first==u) 	p1.first = t;
		if(p2.onBoundary() && p2.last==u) 	p2.last = t;
		if(p3.onBoundary() && p3.first==t) 	p3.first = u;
	}
	
	
	void rotateVertex(int pointId, Point axis, double angle){
		Vector translation = new Vector(axis.x(), axis.y());
		for(Vertex v: vertices) v.subtractThis(translation);
		
		PriorityQueue<FlipEvent> heap = new PriorityQueue<FlipEvent>();
		
		
		for(Vertex v: vertices) v.addThis(translation);
	}
	
	
	static class FlipEvent implements Comparable<FlipEvent>{
		double angle;
		Triangle tri;
		int e;
		FlipEvent(double angle, Triangle tri, int e){
			this.angle = angle;
			this.tri = tri;
		}
		
		public int compareTo(FlipEvent e) {
			return Double.compare(angle, e.angle);
		}
	}
	
}
