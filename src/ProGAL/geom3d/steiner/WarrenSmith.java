package ProGAL.geom3d.steiner;

import java.awt.Color;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.math.Randomization;

public class WarrenSmith {

	
	public static void main(String[] args){

		int n = 5;
		int d = 3;
		
		Randomization.seed(2);
		List<Point> points = PointList.generatePointsInCube(n);
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		for(Point p: points) p.toScene(scene, 0.03, Color.blue);

		if (n == 3) {
			Point.getSteinerPoint(points.get(0), points.get(1), points.get(2));
		}
		else {
			int[] topology = new int[n-3];
		}
		
	}
	
}
