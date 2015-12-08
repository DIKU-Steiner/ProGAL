package ProGAL.geom3d.kineticDelaunay;

import java.awt.Color;

import ProGAL.geom3d.Circle;
import ProGAL.geom3d.Line;
import ProGAL.geom3d.Plane;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.geom3d.volumes.Torus;

public class Tester {

	public static Color blueTransp = new Color(0,0,255,100);
	public static Color redTransp = new Color(255,0,0,100);
	
	private static void edgeCase1() {
		double alpha = 0.5;
		J3DScene scene1 = J3DScene.createJ3DSceneInFrame();
		Point origo = new Point(0,0,0);
		Vector zUnitVector = new Vector(0,1,0);
		Line zAxis = new Line(origo, zUnitVector);
		zAxis.toScene(scene1, 0.01, Color.black);
		Point rotPoint = new Point(Math.cos(Math.PI/36), 0, Math.sin(Math.PI/36));
		Sphere rotSphere = new Sphere(rotPoint, alpha);
		Point stPoint = new Point(0.4, 0.5, 0.7);
		Sphere stSphere = new Sphere(stPoint, alpha);
		rotPoint.toScene(scene1, 0.03, Color.blue, 16);
		stPoint.toScene(scene1, 0.03, Color.red, 16);
		scene1.addShape(stSphere, redTransp, 64);
		Circle rotOrbit = new Circle(origo, 1, zUnitVector);
		rotOrbit.toScene(scene1, 0.01, Color.blue);
		scene1.addShape(rotSphere, blueTransp, 64);
	}	

	
	private static void edgeCase2() {
		double alpha = 0.5;
		J3DScene scene2 = J3DScene.createJ3DSceneInFrame();
		Point origo = new Point(0,0,0);
		Vector zUnitVector = new Vector(0,1,0);
		Line zAxis = new Line(origo, zUnitVector);
		zAxis.toScene(scene2, 0.01, Color.black);
		Point rotPoint = new Point(Math.cos(Math.PI/36), 0, Math.sin(Math.PI/36));
		Point stPoint = new Point(0.4, 0.5, 0.7);
		Sphere stSphere = new Sphere(stPoint, 2*alpha);
		rotPoint.toScene(scene2, 0.03, Color.blue, 16);
		stPoint.toScene(scene2, 0.03, Color.red, 16);
		scene2.addShape(stSphere, redTransp, 64);
		Circle rotOrbit = new Circle(origo, 1, zUnitVector);
		rotOrbit.toScene(scene2, 0.01, Color.blue);
	}	

	private static void triangleCase() {
		double alpha = 0.5;
		J3DScene scene3 = J3DScene.createJ3DSceneInFrame();
		Point origo = new Point(0,0,0);
		Vector zUnitVector = new Vector(0,0, 1);
		Line zAxis = new Line(origo, zUnitVector);
		zAxis.toScene(scene3, 0.01, Color.black);
		
		// rotating point, its a-sphere and torus
		Point rotPoint = new Point(Math.cos(Math.PI/36), 0, Math.sin(Math.PI/36));
		rotPoint.toScene(scene3, 0.03, Color.blue, 16);
		Sphere S = new Sphere(rotPoint, alpha);
		S.toScene(scene3, blueTransp, 64);
		Torus torus = new Torus(origo, zUnitVector, 1, alpha);
		torus.toScene(scene3, blueTransp, 36, 72);
		
		// stationary points, their bisecting plane and circle
		Point stPoint1 = new Point(0.4, 0.5, 0.7);
		Point stPoint2 = new Point(0.2, 0.3, 0.2);
		stPoint1.toScene(scene3, 0.03, Color.red, 16);
		stPoint2.toScene(scene3, 0.03, Color.red, 16);
		Point midPoint = Point.getMidpoint(stPoint1, stPoint2);
		Vector normal = new Vector(midPoint, stPoint2).normalize();
		Plane P = new Plane(midPoint, normal);
		P.toScene(scene3, Color.red, 3);
		double radius = Math.sqrt(alpha*alpha - midPoint.distanceSquared(stPoint2));
		Circle redCircle = new Circle(midPoint, radius, normal);
		redCircle.toScene(scene3, 0.01, Color.red);
		Circle rotOrbit = new Circle(origo, 1, zUnitVector);
		rotOrbit.toScene(scene3, 0.01, Color.blue);
	}	

	
	public static void main(String[] args){
//		edgeCase1();
//		edgeCase2();
		triangleCase();
	}
}
