package ProGAL.geom3d.complex.delaunayComplex;

import java.awt.Color;
import java.util.ArrayList;
//import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
//import java.util.PriorityQueue;

import ProGAL.geom3d.Circle;
import ProGAL.geom3d.Line;
import ProGAL.geom3d.LineSegment;
import ProGAL.geom3d.Plane;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.complex.CTriangle;
import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.predicates.ExactJavaPredicates;
import ProGAL.geom3d.predicates.Predicates;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;

import ProGAL.dataStructures.*;

public class DelaunayComplexRotation {
	J3DScene scene;
	List<Point> points;
	DelaunayComplex dt;
	
	Point p;
	Point q;
	Line rotationAxis;
	HeapItemComp heapItemComp = new HeapItemComp();
	Heap heap;
//	PriorityQueue<HeapItem> heap;
	HashMap<CTetrahedron, HeapItem> map;
	private final Predicates predicates;
//	private final Walk walk;
//	private final Flip14 f14;
	private final Flips flips;

	public DelaunayComplexRotation(List<Point> points) {
		this.predicates = new ExactJavaPredicates();
		this.flips = new Flips(predicates);

		this.points = points;
		dt = new DelaunayComplex(points);
		heap = new Heap(10, heapItemComp);
		map = new HashMap<CTetrahedron, HeapItem>();
	}
	
	private class HeapItem {
		
		CTetrahedron tetr;
		CTetrahedron oppTetr;
		CVertex v;
		CVertex oppV;
		double angle;
		boolean alive;
		
		public HeapItem(Double angle, CTetrahedron tetr, CTetrahedron oppTetr, CVertex v, CVertex oppV) {
			this.angle = angle;
			this.tetr = tetr;
			this.oppTetr = oppTetr;
			this.v = v;
			this.oppV = oppV;
			alive = true;
		}
		
		public CTetrahedron getTetrahedron() { return tetr; }
		public CTetrahedron getOppTetrahedron() { return oppTetr; }
		public CVertex getVertex() { return v; }
		public CVertex getOppVertex() { return oppV; }
		public double getAngle() { return angle; }
		public void setAlive(boolean flag) { alive = flag; }
		public boolean isAlive() { return alive; }
	}

	public class HeapItemComp implements SortTool { 
		public int compare(Object item1, Object item2) {
			double dif = ((HeapItem)item1).getAngle() - ((HeapItem)item2).getAngle();
			if (dif < 0.0) return -1; 
			else {
				if (dif == 0.0) return 0; else return 1;
			}
		}
	}

	
	private class Flipper {

		List<CTetrahedron> bornTetrahedra = new ArrayList<CTetrahedron>();
		List<CTetrahedron> deadTetrahedra = new ArrayList<CTetrahedron>();
		
		public Flipper() {
		}
		
		
		private void flip(HeapItem item) {
			CTetrahedron tetr = item.getTetrahedron();
			CTetrahedron oppTetr = item.getOppTetrahedron();
			CVertex v = item.getVertex();
			int vID = tetr.getID(v);
			CVertex oppV = item.getOppVertex();
			int oppVID = oppTetr.getID(oppV);
			
			scene.removeAllShapes();
			for (CTetrahedron t : dt.getAllTetrahedra()) t.toScene(scene, 0.01, Color.black);        

			tetr.toScene(scene, 0.05, Color.red);
			oppTetr.toScene(scene, 0.05, Color.blue);
			LineSegment seg = new LineSegment(v, oppV);
			seg.toScene(scene, 0.02, Color.green);
			v.toScene(scene, 0.1, Color.yellow);
			
			Flips.ApexConfig flipcase = flips.apexConfig(tetr, v, vID, oppTetr, oppV);
			if(flipcase == Flips.ApexConfig.CONVEX) flips.getFlip23().flip23(tetr, vID, oppVID);
			else if(flipcase== Flips.ApexConfig.CONCAVE){
				if(flips.getFlip23().getT3()!=null){
					flips.getFlip32().flip32(tetr, oppTetr,flips.getFlip23().getT3(), vID, oppVID);
					flips.getFlip23().setT3(null);
				}					
			}
		}
	}
	
	public void rotate(int rotatingBondNumber, double rotationAngle) {
		p = points.get(rotatingBondNumber);
		q = points.get(rotatingBondNumber+1);
		rotationAxis = new Line(p, new Vector(p, q));
		rotateOneAtTime(rotatingBondNumber,  rotationAngle);
//		rotateAllAtOnce(rotatingBondNumber, rotationAngle);
//		rotateAllAtOnceBySmallAngle(rotatingBondNumber, rotationAngle);
	}
	
	public void rotateOneAtTime(int rotatingBondNumber, double rotationAngle) {
		Vector dir = new Vector(points.get(rotatingBondNumber), points.get(rotatingBondNumber+1)).normalize();
		rotationAxis = new Line(points.get(rotatingBondNumber), dir);
		
		for (int i = rotatingBondNumber+2; i < points.size(); i++) {
			CVertex v = dt.getVertex(i);                                   			
			Point o = rotationAxis.orthogonalProjection(v);
			double rad = v.distance(o);
			Circle circle = new Circle(o, rad, dir);     
			double angle;
			double accAngle = 0.0;
			Flipper flipper = new Flipper();
			J3DScene scene1 = J3DScene.createJ3DSceneInFrame();
			J3DScene scene2 = J3DScene.createJ3DSceneInFrame();
			J3DScene scene3 = J3DScene.createJ3DSceneInFrame();
			
			drawScene(scene1, v, circle);
			for (CTetrahedron tetr : dt.getAllTetrahedra())  addToHeap(tetr, v, rotationAngle, circle);        

			while (!heap.isEmpty()) {
				HeapItem item = (HeapItem)heap.extract();
				if (item.isAlive()) {
					scene2.removeAllShapes();
					scene3.removeAllShapes();
					v.toScene(scene2, 0.1, Color.pink);
					circle.toScene(scene2, 0.01, 72);
					CTetrahedron tetr = item.getTetrahedron();
					tetr.toScene(scene2, 0.05, Color.yellow);
					CTetrahedron oppTetr = item.getOppTetrahedron();
					oppTetr.toScene(scene2, 0.05, Color.blue);
					angle = item.getAngle();

					flipper.flip(item);

					rotationAxis.rotateIn(v, angle - accAngle);					
					v.toScene(scene3, 0.1, Color.pink);
					circle.toScene(scene3, 0.01, 72);
					tetr.toScene(scene3, 0.05, Color.yellow);
					oppTetr.toScene(scene3, 0.05, Color.blue);

					accAngle = accAngle + angle;
					for (CTetrahedron t : flipper.deadTetrahedra) {
						HeapItem delItem = map.get(t);
						delItem.setAlive(false);
					}
					for (CTetrahedron t : flipper.bornTetrahedra) {
						addToHeap(t, v, rotationAngle, circle);
					}
				}
			}
			J3DScene scene4 = J3DScene.createJ3DSceneInFrame();
			drawScene(scene4, v, circle);
		}
	}
	
	public void rotateAllAtOnce(int rotatingBondNumber, double rotationAngle) {
		
	}
	
	public void rotateAllAtOnceBySmallAngle(int rotatingBondNumber, double rotationAngle) {
		int numberSteps = 10;
		double stepAngle = rotationAngle/numberSteps;
		for (int i = 0; i < numberSteps; i++) {
			rotateOneAtTime(rotatingBondNumber, stepAngle);
		}
	}
	
	
	public void addToHeap(CTetrahedron tetr, CVertex v, double rotationAngle, Circle circle) {
		int id = tetr.getID(v);
		if (id != -1) {
			scene.removeAllShapes();                                
			                                                        drawScene(scene,v,circle); tetr.toScene(scene, 0.07, Color.yellow);
			CTetrahedron oppTetr = tetr.getNeighbour(id); 			
																	oppTetr.toScene(scene, 0.07, Color.blue); 																	
			CVertex oppV = oppTetr.findVertex(tetr);                
			                                                        LineSegment seg = new LineSegment(v, oppV); 
			                                                        seg.clone().toScene(scene, 0.07, Color.green);
			                                                        CVertex[] vert = tetr.getCommonVertices(oppTetr);
			                                                        CTriangle tr = new CTriangle(vert[0],vert[1],vert[2], tetr, oppTetr);
			                                                        tr.toScene(scene, new Color(0,0,255,100));
			Vector dir = circle.getNormal();
			Double angle = null;
			if (rotationAngle < 0) dir = dir.multiplyThis(-1);
			if (oppTetr.containsBigPoint()) {
				Plane plane = tetr.getPlane(oppTetr);
				                                                    plane.toScene(scene, Color.red, 5);
				angle = plane.getIntersectionAngle(circle, v , dir, scene);
			}
			else {
				Sphere sphere = new Sphere(oppTetr);                         scene.addShape(sphere, new Color(255,0,0,100), 72);
				angle = sphere.getIntersectionAngle(circle, v, dir);
			}
			if ((angle != null) && (angle < rotationAngle)) {
																				  	Point intPoint = rotationAxis.rotate(v, angle); intPoint.toScene(scene, 0.05, Color.red);
				HeapItem item = new HeapItem(angle, tetr, oppTetr, v, oppV);
				heap.insert(item);
				System.out.println("heap size = " + heap.size());
			}
		}
	}
	
	public void drawScene(J3DScene scene, CVertex v, Circle circle) {
		for (CTetrahedron t : dt.getAllTetrahedra()) t.toScene(scene, 0.01, Color.black);        
		for (CTetrahedron t : dt.getTetrahedra()) t.toScene(scene, 0.05, Color.black);   
		Vector transl = new Vector(0.1,0.1,0.1);
		for (int i = 0; i < 5; i++) {
			Integer j = new Integer(i);
			scene.addText(j.toString(), points.get(i).add(transl));
		}
		v.toScene(scene, 0.1, Color.yellow);
		circle.toScene(scene, 0.01, 72);
	}
	
	public void setupRotationHeap(CVertex v, double rotationAngle, Circle circle) {
		// go through all tetrahedra and determine the angle they are hit by v that rotates on the circle
		drawScene(scene, v, circle);
		
		for (CTetrahedron tetr : dt.getAllTetrahedra()) {
			addToHeap(tetr, v, rotationAngle, circle);        
		}
	}
	
	public static void main(String[] args) {
		List<Point> points = new ArrayList<Point>();
	//	points.add(new Point(0,0,0));
//		points.add(new Point(4,0,0));
//		points.add(new Point(0,3,0));
//		points.add(new Point(1,1,3));
//		points.add(new Point(0,0,3));
//		points.add(new Point(4,2,5));

		points.add(new Point(0,0,0));
		points.add(new Point(4,0,0));
		points.add(new Point(0,4,0));
		points.add(new Point(-3.1,0.1,8));
		points.add(new Point(1,2,-8));
		DelaunayComplexRotation dcr = new DelaunayComplexRotation(points);
		dcr.scene = J3DScene.createJ3DSceneInFrame();

		dcr.rotateOneAtTime(2, Math.PI);

	}

}

