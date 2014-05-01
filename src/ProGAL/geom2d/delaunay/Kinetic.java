package ProGAL.geom2d.delaunay;


//import java.awt.Color; 
//import java.util.ArrayList;
//import java.util.List;
//import ProGAL.geom2d.Vector;
//
//import ProGAL.dataStructures.Heap;
//import ProGAL.dataStructures.SortTool;
//import ProGAL.geom2d.Circle;
//import ProGAL.geom2d.Line;
//import ProGAL.geom2d.Point;
//import ProGAL.geom2d.PointSet;
//import ProGAL.geom2d.Shape;
//import ProGAL.geom2d.Triangulation;
//import ProGAL.geom2d.Triangulation.TriangulationAlgorithm;
//import ProGAL.geom2d.TriangulationFace;
//import ProGAL.geom2d.TriangulationVertex;
//import ProGAL.geom2d.viewer.J2DScene;
//import ProGAL.math.Randomization;




public class Kinetic {
// Commented out because there were compilation problems
//	private class HeapItem {
//		private double key;
//		private TriangulationVertex v;
//		private TriangulationFace face;
//		private Object oObject;
//		Circle c;
//		int type;
//
//		public HeapItem(double key, TriangulationVertex v, TriangulationFace vFace, Object oObject, int type) {
//			this.key = key;
//			this.v = v;
//			this.face = face;
//			this.oObject = oObject;
//			this.type = type;
//		}
//		
//		private double getKey() { return key; }
//		private TriangulationVertex getVertex() { return v; }
//		private TriangulationFace getFace() { return face; }
//		private Object getOObject() { return oObject; }
//		private int getType() { return type; }
//	}
//	
//	private class SortToolHeapItem implements ProGAL.dataStructures.SortTool {
//		public int compare(Object p1, Object p2) {
//			if ((p1 instanceof HeapItem) && (p2 instanceof HeapItem)) {
//				double x1 = ((HeapItem)p1).getKey();
//				double x2 = ((HeapItem)p2).getKey();
//				if (x1 < x2) return COMP_LESS;
//				else { 
//					if (x1 > x2) return COMP_GRTR; else return COMP_EQUAL; 
//				}
//			}
//			else throw SortTool.err1;
//		}
//	}
//
//	private Triangulation DT;
//	SortToolHeapItem sorter = new SortToolHeapItem();
//	Heap heap = new Heap(10, sorter);
//	J2DScene scene = J2DScene.createJ2DSceneInFrame();
//	int rotationPointIndx;
//	TriangulationVertex rotationPoint;
//	boolean ccw;
//	double rotationAngle;
//	double rotatedSoFar;
//
//	private void support0(TriangulationFace face, TriangulationFace oppFace) {
//		TriangulationVertex v = face.getCorner(0);
//		Circle c = oppFace.getCircumCircle();
//		scene.addShape(c, Color.blue, 0.001);
//
//		Double angle = v.getOrbit(rotationPoint).enteringAngle(v, c, ccw);
//		if (angle != null) {  // intersection found
//			angle =+ rotatedSoFar;
//			if ((ccw && (angle < rotationAngle)) || (!ccw && (angle > rotationAngle))) 
//				heap.insert(new HeapItem(angle, v, face, oppFace, 0));
//		}
//		scene.removeShape(c);
//	}
//
//	private void support1(TriangulationFace face, Line oppLine) {
//		TriangulationVertex v = face.getCorner(0);
//		Double angle = v.getOrbit(rotationPoint).enteringAngle(face.getCorner(0), oppLine, ccw);
//		if (angle != null) {
//			angle =+ rotatedSoFar;
//			if ((ccw && (angle < rotationAngle)) || (!ccw && (angle > rotationAngle))) 
//				heap.insert(new HeapItem(angle, v, face, oppLine, 1));
//		}
//	}					
//
//	private void support4(TriangulationVertex v, TriangulationFace face, TriangulationFace nextFace) {
//		TriangulationVertex u = face.getCorner((face.getIndex(v)+1)%3);
//		TriangulationVertex w = face.getCorner((face.getIndex(v)+2)%3);
//		TriangulationVertex z = nextFace.getCorner((nextFace.getIndex(v)+2)%3);
//		if (Point.leftTurn(u, w, z))  {
//			Circle C = new Circle(u,w,z);
//			scene.addShape(C, Color.blue, 0.001);
//			Double angle = v.getOrbit(rotationPoint).exitingAngle(v, C, ccw);
//			if (angle != null) {
//				angle =+ rotatedSoFar;
//				if ((ccw && (angle < rotationAngle)) || (!ccw && (angle > rotationAngle))) {
//					HeapItem heapItem = new HeapItem(angle, v, face, C, 4);
//					heap.insert(heapItem);
//				}
//			}
//			scene.removeShape(C);
//		}					
//	}
//
//	
//	public void rotate(Point rotationPoint, List<TriangulationVertex> rotVertices, double rotAngle, boolean ccw) {
//
//		Circle rotationPointCircle = new Circle(rotationPoint, 0.01);
//		scene.addShape(rotationPointCircle, Color.red, 0.01);	
//		rotatedSoFar = 0.0;
//		for (TriangulationVertex v : rotVertices) {
//			Point p = (Point)v;
//			scene.addShape(new Circle(rotationPoint, p.distance(rotationPoint)), Color.green, 0.004);
//			scene.addShape(new Circle(p, 0.01), Color.blue, 0.01);
//		}
//		// building the heap
//		for (TriangulationVertex v : rotVertices) {
//			Circle rotatingPointCircle = new Circle(v, 0.01);
//			Circle orbit = new Circle(rotationPoint, ((Point)v).distance(rotationPoint));
//			scene.addShape(rotatingPointCircle, Color.green, 0.01);
//			scene.addShape(orbit, Color.green, 0.005);
//			TriangulationFace firstFace = ((TriangulationVertex)v).getFace();
//			TriangulationFace face = firstFace;   
//			TriangulationFace prevOppFace = null;
//			do {
//				ProGAL.geom2d.Triangle faceTr = face.getTriangle();
//				scene.addShape(faceTr, Color.red, 0.01);
//				TriangulationFace oppFace = face.getOppFace(v);
//				if (oppFace != null) { // oppFace exists, circle-circle in-case
//					ProGAL.geom2d.Triangle oppFaceTr = oppFace.getTriangle();
//					scene.addShape(oppFaceTr, Color.pink, 0.01);
//					if (oppFace != prevOppFace) support0(face, oppFace);
//					scene.removeShape(oppFaceTr);
//					prevOppFace = oppFace;
//				}
//				else support1(face, face.getOppLine(v));
//				// out-case
//				TriangulationFace nextFace = v.getNextFace(face);
//				if (nextFace != null) support4(v, face, nextFace);
//				scene.removeShape(faceTr);
//				face = nextFace;
//			} while ((face != null) && (face != firstFace));
//			if (face == null) { 
//				// first face boundary in-case, 
//				TriangulationFace vFace = v.getFace();
//				TriangulationVertex vNext = vFace.getCorner((vFace.getIndex(v)+1)%3);
//				TriangulationFace fNext = vNext.getFace();
//				Line oppLine = fNext.getOppLine(fNext.getCorner((fNext.getIndex(vNext)+2)%3));			
//				Double angle = v.getOrbit(rotationPoint).enteringAngle(v, oppLine, ccw);
//				if (angle != null) {
//					angle =+ rotatedSoFar;
//					if ((ccw && (angle < rotAngle)) || (!ccw && (angle > rotAngle))) {
//						HeapItem heapItem = new HeapItem(angle, v, vFace, oppLine, 2);
//						heap.insert(heapItem);
//					}
//				}
//				// last face boundary in-case
//				vFace = v.getLastFace();
//				vNext = vFace.getCorner((vFace.getIndex(v)+2)%3);
//				fNext = vNext.getLastFace();
//				oppLine = fNext.getOppLine(fNext.getCorner((fNext.getIndex(vNext)+2)%3));			
//				angle = v.getOrbit(rotationPoint).enteringAngle(v, oppLine, ccw);
//				if (angle != null) {
//					angle =+ rotatedSoFar;
//					if ((ccw && (angle < rotAngle)) || (!ccw && (angle > rotAngle))) {
//						HeapItem heapItem = new HeapItem(angle, v, vFace, oppLine, 3);
//						heap.insert(heapItem);
//					}
//				}
//				scene.removeShape(orbit);
//				scene.removeShape(rotatingPointCircle);
//			}
//		}
//		// moving vertices
//		while (!heap.isEmpty()) {
//			HeapItem heapItem = (HeapItem)heap.extract();
//			TriangulationVertex v0 = heapItem.getVertex();
//			TriangulationFace face = heapItem.getFace();
//			if (face.isAlive()) {
//				Double angle = heapItem.getKey() - rotatedSoFar;
//				switch (heapItem.getType()) {
//				case 0:
//					TriangulationFace oppFace = (TriangulationFace)heapItem.getOObject();
//					if (oppFace.isAlive()) {
//						for (TriangulationVertex u : rotVertices) rotationPoint.rotation(u, angle);
//						DT.flip(face, oppFace, false);
//						TriangulationFace t12 = face.getNeighbor(0);
//						if (t12 != null) support0(face, t12); else support1(face, face.getOppLine(v0));
//						TriangulationFace t23 = oppFace.getNeighbor(0);
//						if (t23 != null) support0(face, t23); else support1(face, face.getOppLine(v0));
//						TriangulationVertex v1 = face.getCorner(1);
//						support0(face, oppFace);
//						TriangulationVertex v2 = face.getCorner(2);
//						TriangulationFace t01 = face.getNeighbor(2);
//						if (t01 != null) support0(face, t01); else support1(face, face.getOppLine(v2));
//						TriangulationFace t30 = oppFace.getNeighbor(1);
//						if (t30 != null) support0(oppFace, t30); else support1(oppFace, oppFace.getOppLine(v2));
//						TriangulationVertex v3 = oppFace.getCorner(1);
//						support0(oppFace, face);
//						if (t01 != null) support4(v0, t01, face);
//						support4(v0, face, oppFace);
//						if (t30 != null) support4(v0, oppFace, t30);
//						if (rotVertices.contains(v1)) {
//							if (t12 != null) support4(v1, t12, face);
//							if (t01 != null) support4(v1, face, t01);
//						}
//						if (rotVertices.contains(v2)) {
//							if (t23 != null) support4(v2, t23, oppFace);
//							support4(v2, oppFace, face);
//							if (t12 != null) support4(v2, face, t12);
//						}
//						if (rotVertices.contains(v3)) {
//							if (t30 != null) support4(v3, t30, oppFace);
//							if (t23 != null) support4(v3, oppFace, t23);
//						}
//					}
//					break;
//				case 1:
////					DT.boundaryFlip(v0, face);//Error
//					break;
//				case 2:
//					int indx = face.getIndex(v0);
//					TriangulationVertex a = face.getCorner((indx+2)%3);
//					if (a.getFace() == face) {
//						face.setAlive(false);
//						DT.triangulationFaces.remove(face);
//						TriangulationFace newFace = new TriangulationFace(v0, a, face.getCorner((indx+2)%3));
//						DT.triangulationFaces.add(newFace);
//						TriangulationFace nextFace = face.getNeighbor((indx+1)%3);
//						TriangulationFace firstNextFace = nextFace;
//						while (nextFace != null) {
//							indx = nextFace.getIndex(v0);
//							DT.triangulationFaces.remove(nextFace);
//							newFace = new TriangulationFace(a, nextFace.getCorner((indx+1)%3), nextFace.getCorner((indx+2)%3));
//							DT.triangulationFaces.add(newFace);
//							nextFace = nextFace.getNeighbor((indx+1)%3);
//						}
//						if (firstNextFace != null) {
//							newFace = new TriangulationFace(a, nextFace.getCorner(1), v0);
//							DT.triangulationFaces.add(newFace);
//						}
//					}
//					break;
//				case 3: 
//					break;
//				case 4:
//					break;
//				case 5:
//					break;
//				}
//			}
//		}
//	}
//	
//	public Kinetic(PointSet points) {
//		DT = new Triangulation(points, TriangulationAlgorithm.Delaunay);
//		
//	}
//	
//	public static void main(String[] args) {
//		Randomization.seed(1);
//		int n = 15;
//		PointSet points = new PointSet(n);
//		int k = n-4;
//		Point rotationPoint = points.get(k);
//		points.translate(new Vector(-rotationPoint.x(), -rotationPoint.y()));
//		points.toConsole(3);
//		Kinetic kineticDT = new Kinetic(points);
//		kineticDT.DT.draw(kineticDT.scene);
//		
//		kineticDT.rotationPoint = kineticDT.DT.vertices.get(k);                                        // points will rotate around this point
//		List<TriangulationVertex> rotatingVertices = new ArrayList<TriangulationVertex>();             // this subset of vertices will rotate
//		for (int i = k+1; i < n; i++) rotatingVertices.add(kineticDT.DT.vertices.get(i));
//		kineticDT.rotationAngle = Math.PI;                                                             // rotation stops when this angle is reached
//		kineticDT.rotatedSoFar = 0.0;
//		kineticDT.ccw = true;
//		kineticDT.rotate(kineticDT.rotationPoint,  rotatingVertices, kineticDT.rotationAngle, kineticDT.ccw);
//	}
}
