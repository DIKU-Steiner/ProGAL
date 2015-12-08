package ProGAL.geom3d.complex.delaunayComplex.tests;

import java.awt.Color;
import java.util.LinkedList;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointWeighted;
import ProGAL.geom3d.complex.CEdge;
import ProGAL.geom3d.complex.delaunayComplex.DelaunayComplex;
import ProGAL.geom3d.complex.delaunayComplex.RegularComplex;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.math.Matrix;
import ProGAL.math.Matrix3x3;

public class RegularComplexVisualTest {

	static void display1(){
		
		Matrix m = new Matrix3x3(new double[][]{
				{1,5,1},
				{2,5,0},
				{2,7,1}
		});
		m.invert().toConsole();
		
		List<PointWeighted> points = new LinkedList<PointWeighted>();
		points.add(new PointWeighted(1.2,0,0, 1));
		points.add(new PointWeighted(0,1,0, 1));
		points.add(new PointWeighted(0,0,1.3, 1));
		points.add(new PointWeighted(1.2,1,1, 1));
		points.add(new PointWeighted(5,0,0, 1));
		points.add(new PointWeighted(0,0,0, 1));
//		RegularComplex rc = new RegularComplex(points);
//
//		J3DScene scene = J3DScene.createJ3DSceneInFrame();
////		scene.setAxisEnabled(true);
//		for(CEdge e: rc.getEdges()){
//			scene.addShape(new LSS(e,0.05), new Color(100,100,100));
//		}
//		for(PointWeighted pw: points){
////			scene.addShape(new Sphere(pw), new Color(100,200,50,100));
//		}
		
//		List<Point> dpoints = new LinkedList<Point>();
//		dpoints.add(new Point(0,0,0));
//		dpoints.add(new Point(1,0,0));
//		dpoints.add(new Point(0,1,0));
//		dpoints.add(new Point(0,0,1.1));
//		dpoints.add(new Point(1,1,1));
//		dpoints.add(new Point(5,0,0));
//		DelaunayComplex dc = new DelaunayComplex(dpoints);
//		
//		J3DScene scene = J3DScene.createJ3DSceneInFrame();
////		scene.setAxisEnabled(true);
//		for(CEdge e: dc.getEdges()){
//			scene.addShape(new LSS(e,0.05), new Color(100,100,100));
//		}
//		for(Point pw: dpoints){
//			scene.addShape(new Sphere(pw, 0.1), new Color(100,200,50,100));
//		}
		
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		display1();
	}

}
