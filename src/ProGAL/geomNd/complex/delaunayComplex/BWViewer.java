package ProGAL.geomNd.complex.delaunayComplex;

import ProGAL.geom2d.Circle;
import ProGAL.geom2d.Triangle;
import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.geom3d.volumes.Tetrahedron;
import ProGAL.geomNd.Point;
import ProGAL.geomNd.complex.Tessel;

public class BWViewer {

	public static void displayTessel(BowyerWatson bw) {
		int d = bw.getDim();
		if (d==2) {
			display2D(bw);
		}
		if (d==3) {
			display3D(bw);
		}
	}
	
	public static void display2D(BowyerWatson bw) {
		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		for (Tessel t: bw.getTessels()) {
			Triangle tri = new Triangle(new ProGAL.geom2d.Point(t.getPoint(0).getCoords()),
					new ProGAL.geom2d.Point(t.getPoint(1).getCoords()),
					new ProGAL.geom2d.Point(t.getPoint(2).getCoords()));
			scene.addShape(tri, java.awt.Color.GREEN.darker());
		}
		for (Point p: bw.getVertices()) {
			scene.addShape(new Circle(new ProGAL.geom2d.Point(p.getCoords()), 0.03), java.awt.Color.BLUE);
		}
	}
	
	public static void display3D(BowyerWatson bw) {
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		for (Tessel t: bw.getTessels()) {
			Tetrahedron tet = new Tetrahedron(new ProGAL.geom3d.Point(t.getPoint(0).getCoords()),
					new ProGAL.geom3d.Point(t.getPoint(1).getCoords()),
					new ProGAL.geom3d.Point(t.getPoint(2).getCoords()),
					new ProGAL.geom3d.Point(t.getPoint(3).getCoords()));
			scene.addShape(tet, new java.awt.Color(0,255,0,80));
			LSS line1 = new LSS(new ProGAL.geom3d.Point(t.getPoint(0).getCoords()), new ProGAL.geom3d.Point(t.getPoint(1).getCoords()),0.003);
			scene.addShape(line1, java.awt.Color.BLACK, 3);
			LSS line2 = new LSS(new ProGAL.geom3d.Point(t.getPoint(0).getCoords()), new ProGAL.geom3d.Point(t.getPoint(2).getCoords()),0.003);
			scene.addShape(line2, java.awt.Color.BLACK, 3);
			LSS line3 = new LSS(new ProGAL.geom3d.Point(t.getPoint(0).getCoords()), new ProGAL.geom3d.Point(t.getPoint(3).getCoords()),0.003);
			scene.addShape(line3, java.awt.Color.BLACK, 3);
			LSS line4 = new LSS(new ProGAL.geom3d.Point(t.getPoint(1).getCoords()), new ProGAL.geom3d.Point(t.getPoint(2).getCoords()),0.003);
			scene.addShape(line4, java.awt.Color.BLACK, 3);
			LSS line5 = new LSS(new ProGAL.geom3d.Point(t.getPoint(1).getCoords()), new ProGAL.geom3d.Point(t.getPoint(3).getCoords()),0.003);
			scene.addShape(line5, java.awt.Color.BLACK, 3);
			LSS line6 = new LSS(new ProGAL.geom3d.Point(t.getPoint(2).getCoords()), new ProGAL.geom3d.Point(t.getPoint(3).getCoords()),0.003);
			scene.addShape(line6, java.awt.Color.BLACK, 3);
		}
		for (Point p: bw.getVertices()) {
			scene.addShape(new Sphere(new ProGAL.geom3d.Point(p.getCoords()), 0.03), java.awt.Color.BLUE);
		}
	}

}
