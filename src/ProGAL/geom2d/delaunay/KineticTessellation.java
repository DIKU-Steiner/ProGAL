package ProGAL.geom2d.delaunay;

import java.util.List;
import java.util.PriorityQueue;

import ProGAL.geom2d.Point;
import ProGAL.geom2d.Vector;

public class KineticTessellation extends DTWithBigPoints{
	
	public KineticTessellation(List<Point> points){
		super(points);
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
	
	public static void main(String[] args) {
		
	}

}
