package ProGAL.geom3d.complex.delaunayComplex;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.PriorityQueue;

import ProGAL.geom3d.Circle;
import ProGAL.geom3d.Line;
import ProGAL.geom3d.Plane;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;

public class DelaunayComplexRotation_Stavefejl {
//	List<Point> points;
//	DelaunayComplex dt;
//	
//	Point p;
//	Point q;
//	Line rotationAxis;
//	PriorityQueue<HeapItem> heap;
//	HashMap<CTetrahedron, HeapItem> map;
//	enum ApexConfig { CONVEX, CONCAVE, COPLANAR };
//
//	public DelaunayComplexRotation_Stavefejl(List<Point> points) {
//		this.points = points;
//		dt = new DelaunayComplex(points);
//		heap = new PriorityQueue<HeapItem>();
//		map = new HashMap<CTetrahedron, HeapItem>();
//	}
//	
//	private class HeapItem {
//		CTetrahedron tetr;
//		CTetrahedron oppTetr;
//		CVertex v;
//		CVertex oppV;
//		double angle;
//		boolean alive;
//		
//		public HeapItem(Double angle, CTetrahedron tetr, CTetrahedron oppTetr, CVertex v, CVertex oppV) {
//			this.angle = angle;
//			this.tetr = tetr;
//			this.oppTetr = oppTetr;
//			this.v = v;
//			this.oppV = oppV;
//			alive = true;
//		}
//		
//		public CTetrahedron getTetrahedron() { return tetr; }
//		public CTetrahedron getOppTetrahedron() { return oppTetr; }
//		public void setAlive(boolean flag) { alive = flag; }
//		public boolean isAlive() { return alive; }
//	}
//	
//	private class Flipper {
//
//		List<CTetrahedron> bornTedrahedra = new ArrayList<CTetrahedron>();
//		List<CTetrahedron> deadTedrahedra = new ArrayList<CTetrahedron>();
//		
//		public Flipper () {
//			
//		}
//		
//		private void flip(HeapItem item) {
//			CTetrahedron tetr = item.getTetrahedron();
//			CTetrahedron oppTetr = item.getOppTetrahedron();
//			
//			ApexConfig flipcase = apexConfig(tetr, p, pid, oppTetr, d);
//
//			if(flipcase==ApexConfig.CONVEX){
//
//				next_t = f23.flip23(t, pid, did);
//			}
//			else if(flipcase==ApexConfig.CONCAVE){
//
//				if(f23.getT3()!=null){
//
//					next_t = f32.flip32(t,t2,f23.getT3(), pid, did);
//
//					f23.setT3(null);
//
//				}						
//
//			
//			
//			
//			
//		}
//	}
//	
//	public void rotate(int rotatingBondNumber, double rotationAngle) {
//		p = points.get(rotatingBondNumber);
//		q = points.get(rotatingBondNumber+1);
//		rotationAxis = new Line(p, new Vector(p, q));
//		rotateOneAtTime(rotatingBondNumber,  rotationAngle);
////		rotateAllAtOnce(rotatingBondNumber, rotationAngle);
////		rotateAllAtOnceBySmallAngle(rotatingBondNumber, rotationAngle);
//	}
//	
//	public void rotateOneAtTime(int rotatingBondNumber, double rotationAngle) {
//		Vector dir = rotationAxis.getDir();
//		for (int i = rotatingBondNumber+2; i < points.size(); i++) {
//			CVertex v = dt.getVertex(i);
//			Point o = rotationAxis.orthogonalProjection(v);
//			double rad = v.distance(o);
//			Circle circle = new Circle(o, rad, dir);
//
//			setupRotationHeap(v, rotationAngle, circle);
//			while (!heap.isEmpty()) {
//				HeapItem item = heap.poll();
//				if (item.isAlive()) {
//					flipper.flip(item);
//					for (CTetrahedron tetr : flipper.deadTetrahedra) {
//						HeapItem delItem = map.get(tetr);
//						delItem.setAlive(false);
//					}
//					for (CTetrahedron tetr : flipper.bornTetrahedra) {
//						addToHeap(tetr, v, rotationAngle, circle);
//					}
//				}
//			}
//		}
//	}
//	
//	public void rotateAllAtOnce(int rotatingBondNumber, double rotationAngle) {
//		
//	}
//	
//	public void rotateAllAtOnceBySmallAngle(int rotatingBondNumber, double rotationAngle) {
//		int numberSteps = 10;
//		double stepAngle = rotationAngle/numberSteps;
//		for (int i = 0; i < numberSteps; i++) {
//			rotateOneAtTime(rotatingBondNumber, stepAngle);
//		}
//	}
//	
//	
//	public void addToHeap(CTetrahedron tetr, CVertex v, double rotationAngle, Circle circle) {
//		int id = tetr.getID(v);
//		CTetrahedron oppTetr = tetr.getNeighbour(id);
//		Vector dir = circle.getNormal();
//		Double angle = null;
//		if (rotationAngle < 0) dir = dir.multiplyThis(-1);
//		if (oppTetr.containsBigPoint()) {
//			Plane plane = tetr.getPlane(oppTetr);
//			angle = plane.getIntersectionAngle(circle, v , dir);
//		}
//		else {
//			Sphere sphere = new Sphere(oppTetr);
//			angle = sphere.getIntersectionAngle(circle, v, dir);
//		}
//		if ((angle != null) && (angle < rotationAngle)) {
//			CVertex oppV = oppTetr.findVertex(tetr);
//			HeapItem item = new HeapItem(angle, tetr, oppTetr, v, oppV);
//			heap.add(item);
//		}
//	}
//	
//	public void setupRotationHeap(CVertex v, double rotationAngle, Circle circle) {
//		// go through all tetrahedra and determine the angle they are hit by v[i]
//		for (CTetrahedron tetr : dt.getTetrahedra()) addToHeap(tetr, v, rotationAngle, circle);
//	}
//	
//	public static void main(String[] args) {
//		List<Point> points = new ArrayList<Point>();
//		points.add(new Point(0,0,0));
//		points.add(new Point(4,0,0));
//		points.add(new Point(0,3,0));
//		points.add(new Point(1,1,3));
//		points.add(new Point(0,0,3));
//		points.add(new Point(4,2,5));
//		DelaunayComplexRotation dcr = new DelaunayComplexRotation(points);
//		dcr.rotateOneAtTime(3, Math.PI);
//
//	}

}
