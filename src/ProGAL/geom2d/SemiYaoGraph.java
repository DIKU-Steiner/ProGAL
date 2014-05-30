package ProGAL.geom2d;

import java.awt.Color;

import ProGAL.dataStructures.Graph;
import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.math.Constants;

public class SemiYaoGraph {
	private Graph G;
	private boolean testing = true;
	private J2DScene scene;
	
	public SemiYaoGraph(PointSet points) {
		if (testing) scene = J2DScene.createJ2DSceneInFrame();
		G = new Graph(points.getSize(), true);
		Line[] bisector = new Line[6];
		Line[] wedgeLine = new Line[6];
		for (int i = 0; i < 6; i++) {
			wedgeLine[i] = new Line(Point.origo, new Point(Math.cos(-Math.PI/6 + i*Math.PI/3), Math.sin(-Math.PI/6 + i*Math.PI/3)));
			bisector[i]  = new Line(Point.origo, new Point(Math.cos(i*Math.PI/3), Math.sin(i*Math.PI/3)));
		}
		for (int i = 0; i < G.V(); i++) {
			Point p = points.get(i);
			for (int dir = 0; dir < 6; dir++) {
				wedgeLine[dir].translateTo(p);
				wedgeLine[(dir+1)%6].translateTo(p);
				bisector[dir].translateTo(p);
				int k = -1;
				double dist = Constants.bigDouble;
				for (int j = 0; j < G.V(); j++) {
					Point r = points.get(j);
					if ((r != p) && wedgeLine[dir].isAbove(r) && wedgeLine[(dir+1)%6].isBelow(r)) {
						double cDist = bisector[dir].projectionParameter(r);
						if (cDist < dist) {
							dist = cDist;
							k = j;
						}
					}
				}
				if (k != -1) {
					G.addEdge(k, i);
					if (testing) new LineSegment(points.get(k), p).toScene(scene, Color.red);
				}
			}
		}			
		if (testing) for (Point p : points) p.toScene(scene, 0.01, Color.blue);		
	}
	public Graph getSemiYaoGraph() { return G; }
	
	public static void main(String[] args) {
		PointSet points = new PointSet(100);
		SemiYaoGraph G = new SemiYaoGraph(points);
		System.out.println(G.toString());
	}

}
