package ProGAL.geom3d.complex.alphaComplex;

import java.util.ArrayList;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.complex.alphaComplex.VoidTree;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;

public class EiffelPointList extends PointList{
	private static final long serialVersionUID = 1L;
	

	public static void main(String[] args){
		List<Point> points = new EiffelPointList();
//		J3DScene scene = J3DScene.createJ3DSceneInFrame();
//		for(Point p: points){
//			scene.addShape(new Sphere(p,0.6), new java.awt.Color(255,0,0,200));
//		}
		VoidTree vt = new VoidTree(points,0);
		new ProGAL.dataStructures.viewer.BinaryTreePainter(vt.root);
	}
	
	public EiffelPointList(){
		super();
		for(Point p: points) add(p);
		for(int[] edge: edges){
			Point p1 = points[(edge[0])];
			Point p2 = points[(edge[1])];
			double delta = 5/p1.distance(p2);
			for(double t=delta;t<1;t+=delta){
				Point p = p1.add(p1.vectorTo(p2).multiplyThis(t));
				add(p);
			}
		}
	}
	
	
	private final Point[] points = {
			new Point(-2.11, -19.54, 24.38),
			new Point(-22.92, 0.86, 24.38),
			new Point(-23.17, 0.85, 24.38),
			new Point(-2.11, -19.29, 24.38),
			new Point(-2.78, 21.92, 24.38),
			new Point(-23.17, 0.85, 21.29),
			new Point(18.03, 1.52, 24.38),
			new Point(-2.11, -19.29, 24.02),
			new Point(-2.78, 21.67, 24.38),
			new Point(-2.11, -19.54, 21.29),
			new Point(18.29, 1.53, 24.38),
			new Point(-22.92, 0.86, 24.02),
			new Point(-2.78, 21.67, 24.02),
			new Point(-2.78, 21.92, 21.29),
			new Point(-23.48, 0.85, 19.01),
			new Point(18.29, 1.53, 21.29),
			new Point(18.03, 1.52, 24.02),
			new Point(-20.12, 0.90, 24.02),
			new Point(-2.78, 22.23, 19.01),
			new Point(-2.10, -19.85, 19.01),
			new Point(18.60, 1.53, 19.01),
			new Point(15.23, 1.48, 24.02),
			new Point(-2.16, -16.49, 24.02),
			new Point(-33.03, 0.69, 6.71),
			new Point(-2.94, 31.78, 6.71),
			new Point(-2.73, 18.87, 24.02),
			new Point(-2.27, -9.59, 41.72),
			new Point(-16.01, 0.97, 24.02),
			new Point(28.14, 1.69, 6.71),
			new Point(-1.95, -29.39, 6.71),
			new Point(11.13, 1.41, 24.02),
			new Point(-2.62, 11.97, 41.72),
			new Point(-2.66, 14.76, 24.02),
			new Point(-13.22, 1.02, 41.72),
			new Point(-2.22, -12.38, 24.02),
			new Point(-37.37, 0.62, 0.28),
			new Point(-3.01, 36.12, 0.28),
			new Point(32.49, 1.76, 0.28),
			new Point(8.34, 1.37, 41.72),
			new Point(-3.14, -10.49, 43.89),
			new Point(-1.88, -33.74, 0.28),
			new Point(9.00, 0.25, 43.89),
			new Point(-1.51, 12.63, 43.89),
			new Point(-14.12, 1.88, 43.89),
			new Point(-14.09, 0.12, 43.89),
			new Point(-1.37, -10.46, 43.89),
			new Point(8.96, 2.50, 43.89),
			new Point(-3.75, 12.60, 43.89),
			new Point(-3.14, -10.49, 44.62),
			new Point(-1.37, -10.46, 44.62),
			new Point(9.00, 0.25, 44.62),
			new Point(8.96, 2.50, 44.62),
			new Point(-3.75, 12.60, 44.62),
			new Point(-14.12, 1.88, 44.62),
			new Point(-14.09, 0.12, 44.62),
			new Point(-1.51, 12.63, 44.62),
			new Point(-12.95, 1.02, 44.62),
			new Point(-2.27, -9.31, 44.62),
			new Point(8.06, 1.36, 44.62),
			new Point(-2.61, 11.70, 44.62),
			new Point(-2.27, -9.31, 45.27),
			new Point(-2.61, 11.70, 45.27),
			new Point(-12.95, 1.02, 45.27),
			new Point(8.06, 1.36, 45.27),
			new Point(-12.00, 1.04, 45.27),
			new Point(7.12, 1.35, 45.27),
			new Point(-2.60, 10.75, 45.27),
			new Point(-2.29, -8.37, 45.27),
			new Point(-2.34, -4.99, 63.96),
			new Point(-7.47, 1.11, 45.27),
			new Point(2.59, 1.27, 45.27),
			new Point(-2.54, 7.37, 63.96),
			new Point(-2.52, 6.22, 45.27),
			new Point(-8.63, 1.09, 63.96),
			new Point(-2.36, -3.84, 45.27),
			new Point(3.74, 1.29, 63.96),
			new Point(-2.36, -3.74, 82.79),
			new Point(-2.52, 6.13, 82.79),
			new Point(-7.38, 1.11, 82.79),
			new Point(2.49, 1.27, 82.79),
			new Point(-2.38, -2.62, 96.18),
			new Point(-2.51, 5.00, 96.18),
			new Point(-6.25, 1.13, 96.18),
			new Point(1.37, 1.25, 96.18),
			new Point(-2.38, -2.62, 96.90),
			new Point(-2.51, 5.00, 96.90),
			new Point(-6.25, 1.13, 96.90),
			new Point(1.37, 1.25, 96.90),
			new Point(-2.35, -4.39, 99.56),
			new Point(3.13, 1.28, 99.56),
			new Point(-2.53, 6.77, 99.56),
			new Point(-8.02, 1.10, 99.56),
			new Point(-2.35, -4.39, 100.75),
			new Point(-2.53, 6.77, 100.75),
			new Point(-8.02, 1.10, 100.75),
			new Point(3.13, 1.28, 100.75),
			new Point(-2.37, -3.40, 103.74),
			new Point(-2.52, 5.79, 103.74),
			new Point(-7.04, 1.12, 103.74),
			new Point(2.15, 1.27, 103.74),
			new Point(-3.25, 1.18, 106.42),
			new Point(-1.64, 1.20, 106.42),
			new Point(-2.46, 2.00, 106.42),
			new Point(-2.43, 0.39, 106.42),
			new Point(-2.43, 0.39, 114.02),
			new Point(-2.46, 2.00, 114.02),
			new Point(-3.25, 1.18, 114.02),
			new Point(-1.64, 1.20, 114.02)
	};
	static int[][] edges = {
		{73,68},
		{68,75},
		{75,71},
		{71,73},
		{64,78},
		{78,77},
		{77,79},
		{79,76},
		{76,78},
		{76,68},
		{79,75},
		{77,71},
		{78,73},
		{82,78},
		{80,76},
		{83,79},
		{77,81},
		{83,81},
		{83,80},
		{82,80},
		{81,82},
		{67,68},
		{65,75},
		{71,66},
		{74,70},
		{72,69},
		{74,69},
		{72,70},
		{38,31},
		{33,26},
		{26,38},
		{33,31},
		{65,70},
		{67,74},
		{69,64},
		{72,66},
		{38,21},
		{31,25},
		{33,17},
		{22,26},
		{22,21},
		{17,25},
		{17,22},
		{21,25},
		{20,18},
		{ 4,10},
		{ 0,10},
		{20,19},
		{ 2,0 },
		{19,14},
		{12,2 },
		{14,18},
		{28,37},
		{28,20},
		{36,24},
		{24,18},
		{35,23},
		{14,23},
		{40,29},
		{19,29},
	};		
}