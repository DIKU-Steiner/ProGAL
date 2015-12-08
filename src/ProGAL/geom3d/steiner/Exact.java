package ProGAL.geom3d.steiner;

import java.awt.Color;
import java.util.List;

import ProGAL.geom3d.Circle;
import ProGAL.geom3d.LineSegment;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.kineticDelaunay.Tet;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.math.Randomization;

public class Exact {

	public Exact(List<Point> points) {
		
	}
	
	public static void main(String[] args){

		int n = 4;
		Randomization.seed(2);
		List<Point> points = PointList.generatePointsInCube(n);
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		for(Point p: points) p.toScene(scene, 0.03, Color.blue);
	
		// find optimal solution for the first three points 
		
		Circle c = new Circle(points.get(0), points.get(1), points.get(2));
		
//		c.toScene(scene, 0.005, 32);
//		c.getCenter().toScene(scene, 0.03, Color.red);
//		LineSegment seg = new LineSegment(c.getCenter(), c.getCenter().add(c.getNormal()));
//		seg.toScene(scene, 0.02, Color.green);
		
		Point s = Point.getSteinerPoint(points.get(0), points.get(1), points.get(2));
		s.toScene(scene, 0.03, Color.yellow);
		LineSegment sp0 = new LineSegment(s, points.get(0));
		LineSegment sp1 = new LineSegment(s, points.get(1));
		LineSegment sp2 = new LineSegment(s, points.get(2));
		sp0.toScene(scene, 0.02, Color.pink);
		sp1.toScene(scene, 0.02, Color.pink);
		sp2.toScene(scene, 0.02, Color.pink);
		
		
		for (int i = 4; i < n; i++) {
			points.get(i).toScene(scene, 0.03, Color.blue);
			int j = 2*i - 5;
			for (int k = 0; k < j; k++) {
				
			}
		}
	}
}
