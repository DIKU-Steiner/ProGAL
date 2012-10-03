package ProGAL.geom2d.delaunay;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;

import ProGAL.geom2d.Circle;
import ProGAL.geom2d.Point;
import ProGAL.geom2d.Vector;
import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.geom2d.viewer.TextShape;
import ProGAL.math.Randomization;

public class KineticTessellation extends DTWithBigPoints{

	public KineticTessellation(List<Point> points){
		super(points);
	}
	volatile boolean acceptEvent = false;

	void rotateVertex(int pointId, Point axis, double angle){
		//		Vector translation = new Vector(axis.x(), axis.y());
		//		for(Vertex v: vertices) v.subtractThis(translation);

		Vertex p = super.vertices.get(pointId);
		double startAngle = Point.getSignedAngle(axis.add(new Vector(1,0)), axis, p);
		double curAngle = startAngle;
		List<Triangle> star = star(p);

		PriorityQueue<FlipEvent> heap = initializeEventQueue(p, star, axis, curAngle, angle);
		while(curAngle<Math.signum(angle)*(startAngle+angle))
		{
			while(!acceptEvent) Thread.yield();
			FlipEvent evt = heap.poll();
			curAngle += (evt.angle-curAngle);
			angle -= (evt.angle-curAngle);
			p.set(rotateAroundAxis(p, axis, evt.angle-curAngle));
			if(evt.tri.corners[evt.e]==p)
				increaseStar(p,evt,star,heap,axis,angle);
			else
				decreaseStar(p,evt,star,axis,angle);
		}
		//		for(Vertex v: vertices) v.addThis(translation);
	}

	//	J2DScene scene = J2DScene.createJ2DSceneInFrame();

	private static void increaseStar(Vertex p, FlipEvent evt, List<Triangle> star, PriorityQueue<FlipEvent> events, Point axis, double angle){
		Triangle t0 = evt.tri;
		int pIdx = evt.e;
		if(t0.corners[pIdx]!=p) return;//The event has become obsolete
		
		Triangle t1 = evt.tri.neighbors[evt.e];
		int oppositeP = (t1.indexOf(t0.corners[(pIdx+1)%3])+1)%3;
		
		Triangle[] ns = {t1.neighbors[(oppositeP+1)%3], t1.neighbors[(oppositeP+2)%3]};
		
		flip(t0,pIdx, t1,oppositeP);
		
		star.add(star.indexOf(t0)+1, t1);
		
		double rSq = axis.distanceSquared(p);
		double r = Math.sqrt(rSq);
		for(int n=0;n<2;n++){
			Circle c = new Circle(ns[n]);
			double dSq = axis.distanceSquared(c.center());
			double d = Math.sqrt(dSq);
			if(d<r+c.getRadius() && d>Math.abs(r-c.getRadius())){ 
				double x = (dSq-c.getRadius()*c.getRadius()+rSq)/(2*d);
				double theta = Math.acos(x/r);
				double psi = Point.getSignedAngle(p, axis, c.center());
				double alpha = psi-Math.signum(angle)*theta;
				evt = new FlipEvent(p,alpha, t0, pIdx, c);
				events.add(evt);
			}
		}
		
		Triangle tPrev = t0.neighbors[(pIdx+2)%3];
		int pIdxPrev = tPrev.indexOf(p); 

		Point p0 = tPrev.getCorner((pIdxPrev+2)%3);
		Point p1 = t0.getCorner((pIdx+1)%3);
		Point p2 = t0.getCorner((pIdx+2)%3);
		Circle c = new Circle(p0,p1,p2);
		if(c.contains(p)){
			double dSq = axis.distanceSquared(c.center());
			double d = Math.sqrt(dSq);
			if(d<r+c.getRadius() && d>Math.abs(r-c.getRadius())){ 
				double x = (dSq-c.getRadius()*c.getRadius()+rSq)/(2*d);
				double theta = Math.acos(x/r);
				double psi = Point.getSignedAngle(p, axis, c.center());
				double alpha = psi+Math.signum(angle)*theta;
				evt = new FlipEvent(p, alpha, tPrev, (pIdxPrev+2)%3, c);
				events.add(evt);
			}
		}		
		Triangle tNext = t1.neighbors[oppositeP];
		int pIdxNext = tNext.indexOf(p); 

		p0 = t1.getCorner(oppositeP);
		p1 = t1.getCorner((oppositeP+1)%3);
		p2 = tNext.getCorner((pIdxNext+2)%3);
		c = new Circle(p0,p1,p2);
		if(c.contains(p)){
			double dSq = axis.distanceSquared(c.center());
			double d = Math.sqrt(dSq);
			if(d<r+c.getRadius() && d>Math.abs(r-c.getRadius())){ 
				double x = (dSq-c.getRadius()*c.getRadius()+rSq)/(2*d);
				double theta = Math.acos(x/r);
				double psi = Point.getSignedAngle(p, axis, c.center());
				double alpha = psi+Math.signum(angle)*theta;
				evt = new FlipEvent(p, alpha, t1, (oppositeP+2)%3, c);
				events.add(evt);
			}
		}
		
	}
	private static void decreaseStar(Vertex p, FlipEvent evt, List<Triangle> star, Point axis, double angle){

	}

	private static PriorityQueue<FlipEvent> initializeEventQueue(Vertex v, List<Triangle> star, Point axis, double curAngle, double angle){
		PriorityQueue<FlipEvent> ret = new PriorityQueue<FlipEvent>();
		//		for(Triangle t: triangles){
		//			scene.addShape(t, Color.BLACK, 0.001);
		//		}
		double rSq = axis.distanceSquared(v);
		double r = Math.sqrt(rSq);
		for(Triangle t: star){
			//Check for increase-star event
			int vIdx = t.indexOf(v);
			if(t.neighbors[vIdx]!=null){
				Circle c = new Circle(t.neighbors[vIdx]);
				double dSq = axis.distanceSquared(c.center());
				double d = Math.sqrt(dSq);
				if(d<r+c.getRadius() && d>Math.abs(r-c.getRadius())){ 
					double x = (dSq-c.getRadius()*c.getRadius()+rSq)/(2*d);
					double theta = Math.acos(x/r);
					double psi = Point.getSignedAngle(v, axis, c.center());
					double alpha = psi-Math.signum(angle)*theta;
					FlipEvent evt = new FlipEvent(v,alpha, t, vIdx, c);
					ret.add(evt);
				}
			}

			//Check for decrease-star event
			Triangle tNext = t.neighbors[(vIdx+1)%3];
			int vIdxNext = tNext.indexOf(v); 

			Point p0 = t.getCorner((vIdx+1)%3);
			Point p1 = t.getCorner((vIdx+2)%3);
			Point p2 = tNext.getCorner((vIdxNext+2)%3);
			Circle c = new Circle(p0,p1,p2);
			if(c.contains(v)){
				double dSq = axis.distanceSquared(c.center());
				double d = Math.sqrt(dSq);
				if(d<r+c.getRadius() && d>Math.abs(r-c.getRadius())){ 
					double x = (dSq-c.getRadius()*c.getRadius()+rSq)/(2*d);
					double theta = Math.acos(x/r);
					double psi = Point.getSignedAngle(v, axis, c.center());
					double alpha = psi+Math.signum(angle)*theta;
					FlipEvent evt = new FlipEvent(v, alpha, t, (vIdx+1)%3, c);
					ret.add(evt);
				}
			}
		}

		return ret;
	}

	private static List<Triangle> star(Vertex v){
		List<Triangle> ret = new ArrayList<Triangle>();

		Triangle t = v.first;
		do{
			ret.add(t);
			int vIdx = t.indexOf(v);
			t = t.neighbors[(vIdx+1)%3];
		}while(t!=v.first);

		return ret;
	}

	static class FlipEvent implements Comparable<FlipEvent>{
		Circle circle;
		double angle;
		Triangle tri;
		int e;
		Vertex p;
		FlipEvent(Vertex p, double angle, Triangle tri, int e, Circle c){
			this.p = p;
			this.angle = angle;
			this.tri = tri;
			this.e = e;
			this.circle = c;
		}

		public int compareTo(FlipEvent e) {
			return Double.compare(Math.abs(angle), Math.abs(e.angle));
		}
	}

	static Point rotateAroundAxis(Point point, Point axis, double angle){
		Vector p = axis.vectorTo(point);
		return axis.add(p.rotateThis(angle));
	}

	public static void main(String[] args) {
		Randomization.seed(2);
		List<Point> points = new ArrayList<Point>();
		for(int i=0;i<10;i++) points.add(new Point(Randomization.randBetween(-1.0, 1.0), Randomization.randBetween(-1.0, 1.0)));
		KineticTessellation dt = new KineticTessellation(points);

		Vertex v = dt.vertices.get(7);
		List<Triangle> star = star(v);
		Point axis = new Point(0.5,0.5);
		PriorityQueue<FlipEvent> events = dt.initializeEventQueue(v, star, axis, 0, Math.PI);

		J2DScene scene = J2DScene.createJ2DSceneInFrame();

		for(Triangle t: dt.triangles) {
			scene.addShape(t, Color.GRAY, 0.001, false);
		}

		for(Triangle t: star) scene.addShape(t, new Color(100,100,100,40), 0.001, true); 
		scene.addShape(new Circle(v,0.02), Color.BLACK, 0, true);

		Circle path =  new Circle(axis,axis.distance(v));
		scene.addShape(path, Color.RED.darker());
		//		for(FlipEvent e: events)
		FlipEvent e = events.poll();
		{
			if(e.tri.corners[e.e]==v){//Increase star event
				Triangle link = e.tri.neighbors[e.e];
				scene.addShape(link, new Color(100,100,100,80),0,true);
				//				scene.addShape(link.getCircumCircle(), new Color(0,100,0,100));
			}
			scene.addShape(e.circle, new Color(0,100,0,100));

			Circle intersection = new Circle(rotateAroundAxis(v, axis, e.angle),0.04);
			scene.addShape(intersection, Color.GREEN.darker(), 0,true);
		}


		for(Vertex p: dt.vertices) scene.addShape(new TextShape(p.toString(),p, 0.08));
	}

}
