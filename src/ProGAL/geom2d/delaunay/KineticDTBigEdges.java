package ProGAL.geom2d.delaunay;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import ProGAL.dataStructures.Heap;
import ProGAL.dataStructures.SortTool;
import ProGAL.dataStructures.SortToolPoint2dDistance;
import ProGAL.dataStructures.Sorter;
import ProGAL.dataStructures.SorterQuick;
import ProGAL.geom2d.Circle;
import ProGAL.geom2d.Point;
import ProGAL.geom2d.PointSet;
import ProGAL.geom2d.Shape;
import ProGAL.geom2d.Triangulation;
import ProGAL.geom2d.TriangulationFace;
import ProGAL.geom2d.TriangulationVertex;
import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.geom2d.viewer.TextShape;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.math.Constants;
import ProGAL.math.Functions;
import ProGAL.math.Trigonometry;

public class KineticDTBigEdges extends Triangulation {
	
	private static enum Direction { CW, CCW }
	
	private Point rotationPoint;                                  // rotation center
	private List<Integer> rotIndx = new ArrayList<Integer>();     // indicies of vertices that move
	private Direction rotDir;
	private Double[] angles = new Double[2];
	private List<Shape> markShapes = new ArrayList<Shape>();

	private double angleTotal = 0.0;
	private int flipNr = 0;
	private J2DScene scene = J2DScene.createJ2DSceneInFrame();
	private Heap heap = new Heap(this.vertices.size(), new SortToolHeapItems());
	
	private boolean testing = false;
	private boolean circleAnimation = false;
	private boolean printLabels = false;
	
	private class HeapItem {
		private Double[] angles;
		private TriangulationFace t;
		private TriangulationFace oppT;
		private List<Shape> stops;
		private List<Shape> labels;
		private HeapItem createdBy;
		
		private HeapItem(Double[] angles, TriangulationFace t, TriangulationFace oppT, List<Shape> stops, List<Shape> labels, HeapItem createdBy) {
			this.angles = angles;
			this.t = t;
			this.oppT = oppT;
			this.stops = stops;
			this.labels = labels;
			this.createdBy = createdBy;
		}
				
		private Double getAngle() { return angles[0];} 
		private TriangulationFace getT() { return t; }
		private TriangulationFace getOppT() { return oppT; }
		private List<Shape> getStops() { return stops; }
		private List<Shape> getLabels() { return labels; }
		private Double getAngle1() { return angles[1]; }
	}
	
	private class SortToolHeapItems implements SortTool {
		public int compare(Object x1, Object x2) {
			if ((x1 instanceof HeapItem) && (x2 instanceof HeapItem)) {
				double d1 = ((HeapItem)x1).getAngle();
				double d2 = ((HeapItem)x2).getAngle();
				if (d1 < d2) return COMP_LESS;
				else { if (d1 > d2) return COMP_GRTR; else return COMP_EQUAL; }
			}
			else throw SortTool.err1;
		}

	}
	public KineticDTBigEdges(PointSet points) {
		super(points, TriangulationAlgorithm.Delaunay);	
		Circle c1 = new Circle(new Point(15,15),0.1);
		Circle c2 = new Circle(new Point(-5,-5),0.1);
		scene.addShape(c1);
		scene.addShape(c2);
		scene.centerCamera();
		scene.autoZoom();
		scene.removeAllShapes();
		
		draw(scene); 
		scene.frame.setLocation(700, 0);
		scene.frame.setSize(750,900);		
	}

	public Point getRotationPoint() { return rotationPoint; }
	public void setRotationPoint(Point rotationPoint) { 
		if (this.rotationPoint != null) for (TriangulationVertex v : vertices) v.addThis(this.rotationPoint);
		this.rotationPoint = rotationPoint; 
		for (TriangulationVertex v : vertices) v.subtractThis(rotationPoint);	
	}
	
	public Direction getDirection() { return rotDir; }
	public void setDirection(Direction rotDir) { this.rotDir = rotDir; }
	
	/* specify what is rotating - 3 methods to do it */
	public void setRotVertices(int a) { setRotVertices(a, a); }
	public void setRotVertices(int a, int b) { 
		for (int i = a; i <= b; i++) rotIndx.add(i); 
		for (TriangulationVertex v : vertices) v.setType(TriangulationVertex.VertexType.S); 
		for (int i : rotIndx) vertices.get(i).setType(TriangulationVertex.VertexType.R);
	}
	public void setRotVertices(List<Integer> rotList) { 
		rotIndx = rotList; 
		for (TriangulationVertex v : vertices) v.setType(TriangulationVertex.VertexType.S); 
		for (int i : rotIndx) vertices.get(i).setType(TriangulationVertex.VertexType.R);
	}
	
	public void animateCircles() {
		for (TriangulationFace face : triangulationFaces) {
			if (!face.isFlat()) {
				Circle c = face.getCircumCircle();
				if (c == null) {
					face.setCircumCircle(new Circle(face.getCorner(0), face.getCorner(1), face.getCorner(2)));
					scene.addShape(face.getCircumCircle(), Color.red);
				}
				else {
					Circle d = new Circle(face.getCorner(0), face.getCorner(1), face.getCorner(2));
					c.setCenter(d.center());
					c.setRadius(d.getRadius()); 
				}
			}
		}
	}


	public void animate(double alpha) {
		int steps = 1;
		double angleStep = alpha/steps;
		double sin = Math.sin(angleStep);
		double cos = Math.cos(angleStep);

		for (int k = 0; k < steps; k++) {
			for (int i : rotIndx) vertices.get(i).rotation(cos, sin);
			if (true) {
				try {
					if (vertices.size() < 10000) Thread.sleep(2);
				} catch (InterruptedException e) {}
			}
			if (circleAnimation) animateCircles();
			scene.repaint();
		}
	}

	
	public Double[] getRoot(TriangulationFace face, TriangulationVertex oppVertex) {
		int count = face.getCount();
		if (oppVertex.getType() == TriangulationVertex.VertexType.R) count = count + 8;	
		
		if ((count == 0) || (count == 15)) return null;
		if (count <= 7) {
			if (count <= 3) {
				if (count == 1) return getRootSSSR(face.getCorner(1), face.getCorner(2), oppVertex, face.getCorner(0), 0);
				if (count == 2) return getRootSSSR(face.getCorner(0), face.getCorner(2), oppVertex, face.getCorner(1), 0);
				return getRootSSRR(face.getCorner(2), oppVertex, face.getCorner(0), face.getCorner(1));
			}
			if (count <= 5) {
				if (count == 4) return getRootSSSR(face.getCorner(0), face.getCorner(1), oppVertex, face.getCorner(2), 0);
				return getRootSSRR(face.getCorner(1), oppVertex, face.getCorner(0), face.getCorner(2));
			}
			if (count == 6) return getRootSSRR(face.getCorner(0), oppVertex, face.getCorner(1), face.getCorner(2));					
			return getRootSSSR(face.getCorner(0), face.getCorner(1), face.getCorner(2), oppVertex, 1);
		}
		if (count <= 11) {
			if (count <= 9) {
				if (count == 8) return getRootSSSR(face.getCorner(0), face.getCorner(1), face.getCorner(2), oppVertex, 0);
				return getRootSSRR(face.getCorner(1), face.getCorner(2), face.getCorner(0), oppVertex);
			}
			if (count == 10) return getRootSSRR(face.getCorner(0), face.getCorner(2), face.getCorner(1), oppVertex);
			return getRootSSSR(face.getCorner(0), face.getCorner(1), oppVertex, face.getCorner(2), 1);
		}
		if (count <= 13) {
			if (count == 12) return getRootSSRR(face.getCorner(0), face.getCorner(1), face.getCorner(2), oppVertex);
			return getRootSSSR(face.getCorner(0), face.getCorner(2), oppVertex, face.getCorner(1), 1);
		}
		return getRootSSSR(face.getCorner(1), face.getCorner(2), oppVertex, face.getCorner(0), 1);
	}
		
		
	private Double[] getRootSSRR(TriangulationVertex A, TriangulationVertex B, TriangulationVertex C, TriangulationVertex D) {	
		double aa  = A.getSquaredPolarRadius(); 
		double bb  = B.getSquaredPolarRadius(); 
		double cc  = C.getSquaredPolarRadius(); 
		double c = C.getPolarRadius();
		double dd  = D.getSquaredPolarRadius(); 
		double d = D.getPolarRadius();
	
		double m11 = aa - bb;
		double m12 = A.y() - B.y();
		double m13 = A.y()*bb - B.y()*aa;
		double m22 = A.x() - B.x();
		double m23 = A.x()*bb - B.x()*aa;
		double m33 = A.x()*B.y() - A.y()*B.x();

		double cosAlpha = C.getCosAngle();
		double sinAlpha = C.getSinAngle();
		double cosBeta  = D.getCosAngle();
		double sinBeta  = D.getSinAngle();
		double m13ccm12 = m13 - cc*m12;
		double m23ccm22 = m23 - cc*m22;
		double m13ddm12 = m13 - dd*m12;
		double m23ddm22 = m23 - dd*m22;
		double coefSin = d*(cosBeta*m23ccm22 + sinBeta*m13ccm12) - c*(cosAlpha*m23ddm22 + sinAlpha*m13ddm12);
		double coefCos = d*(sinBeta*m23ccm22 - cosBeta*m13ccm12) + c*(cosAlpha*m13ddm12 - sinAlpha*m23ddm22);
		double coef    = d*c*m11*Math.sin(D.getPolarAngle()-C.getPolarAngle()) + (cc - dd)*m33;
	
		angles = Trigonometry.solveAsinXPlusBcosXplusC(coefSin, coefCos, coef);
		return getRotAngle(angles, 0);
	}
		
	
	public Double[] getRotAngle(Double[] angles, int dir) {
		if ((angles == null) || (angles[1] == null)) return null;
		if (dir == 1) { 
			angles[0] = Constants.TAU - angles[0] + angleTotal;
			if (angles[0] > Constants.TAU) angles[0] = angles[0] - Constants.TAU;
			angles[1] = Constants.TAU - angles[1] + angleTotal;
			if (angles[1] > Constants.TAU) angles[1] = angles[1] - Constants.TAU;
		}
		if (angles[1] < angles[0]) { double tmp = angles[1]; angles[1] = angles[0]; angles[0] = tmp; }
		if (rotDir == Direction.CCW) {
			if (angles[0] < angleTotal + Constants.EPSILON) {
				if (angles[1] < angleTotal + Constants.EPSILON) angles[0] = angles[1] = Constants.TAU;
				else {
					angles[0] = angles[1];
					angles[1] = Constants.TAU;
				}
			}
		}
		else {
			if (Constants.TAU - angles[1] < angleTotal + Constants.EPSILON) {
				if (Constants.TAU - angles[0] < angleTotal + Constants.EPSILON) angles[0] = angles[1] = Constants.TAU; 
				else { 
					angles[0] = Constants.TAU - angles[0];
					angles[1] = Constants.TAU;
				}
			}
			else {
				double tmp = angles[0];
				angles[0] = Constants.TAU - angles[1];
				angles[1] = Constants.TAU - tmp;
			}
		}
		return angles;
	}	
		
	private Double[] getRootSSSR(TriangulationVertex A, 
							   TriangulationVertex B, 
							   TriangulationVertex C, 
							   TriangulationVertex D, int dir) {
		double aa  = A.getSquaredPolarRadius(); 
		double bb  = B.getSquaredPolarRadius(); 
		double cc  = C.getSquaredPolarRadius();
		double dd  = D.getSquaredPolarRadius(); 
		double d   = D.getPolarRadius();
	
		double m41 = A.y()*(bb-cc) + B.y()*(cc-aa) + C.y()*(aa-bb);
		double m42 = A.x()*(bb-cc) + B.x()*(cc-aa) + C.x()*(aa-bb);
		double m43 = A.x()*(B.y()-C.y()) + B.x()*(C.y()-A.y()) + C.x()*(A.y()-B.y());
		double m44 = A.x()*(B.y()*cc -C.y()*bb) + B.x()*(C.y()*aa - A.y()*cc) + C.x()*(A.y()*bb - B.y()*aa);

		double cosBeta = D.getCosAngle();
		double sinBeta = D.getSinAngle();

		double coefSin = d*(sinBeta*m41 + cosBeta*m42);
		double coefCos = d*(sinBeta*m42 - cosBeta*m41);
		double coef    = m44 - dd*m43;
		angles = Trigonometry.solveAsinXPlusBcosXplusC(coefSin, coefCos, coef);
		return getRotAngle(angles, dir);
	}
	
	private void printHeap() {
		for (Object obj : heap.getObjects()) {
			System.out.println(((HeapItem)obj).getT().toString() + ((HeapItem)obj).getOppT().toString() + " " + ((HeapItem)obj).getAngle());
		}
	}

	private List<Shape> markAllStops(double angle) {
		List<Shape> shapes = new ArrayList<Shape>();
		double rotAngle = angle;
		if (rotDir == Direction.CW) rotAngle = - angle;
		for (int k : rotIndx) {
			Shape mark0 = new Circle(vertices.get(k).rotationClone(rotAngle), 0.1);
			scene.addShape(mark0, Color.red);
			markShapes.add(mark0);
			shapes.add(mark0);
		}
		return shapes;
	}

	private List<Shape> addLabels(TriangulationFace t, TriangulationFace oppT, double angle) {
		List<Shape> labels = new ArrayList<Shape>();
		double rotAngle = angle;
		if (rotDir == Direction.CW) rotAngle = -angle;
		for (int j : rotIndx) {
			Shape text = new TextShape(t.toString() + oppT.toString(), vertices.get(j).rotationClone(rotAngle).add(0.1,0.0), 0.3);
			scene.addShape(text);
			labels.add(text);
		}
		return labels;
	}

	private void addToHeap(Double[] angles, TriangulationFace t, TriangulationFace oppT, HeapItem createdBy) {
		if (angles[0] < Constants.TAU) {
			List<Shape> stops = null;
			List<Shape> labels = null;
			if (testing) {
				System.out.println(t.toString() + oppT.toString() + ", rotate to = " + Functions.toDeg(angles[0])); 
				stops = markAllStops(angles[0] - angleTotal);
					if (printLabels) labels = addLabels(t, oppT, angles[0] - angleTotal);
			}
			heap.insert(new HeapItem(angles, t, oppT, stops, labels, createdBy));
		}
		else if (testing) System.out.println("Faces " + t.toString() + " and " + oppT.toString() + ", skipped ");					
		if (testing) showCircles(t, oppT, angles);
	}
	
	private void showCircles(TriangulationFace t, TriangulationFace oppT, Double[] angles) {
		
		if (angles[0] != Constants.TAU) {
			Shape circle0 = t.drawRotatedCircumCircle(scene, Color.blue, angles[0]-angleTotal, rotIndx);
			Shape oppCircle0 = oppT.drawRotatedCircumCircle(scene, Color.magenta, angles[0]-angleTotal, rotIndx);
			if ((angles[1] != null) && (angles[1] != Constants.TAU)) {
				Shape circle1 = t.drawRotatedCircumCircle(scene, Color.cyan, angles[1]-angleTotal, rotIndx);
				Shape oppCircle1 = oppT.drawRotatedCircumCircle(scene, Color.pink, angles[1]-angleTotal, rotIndx);
				scene.removeShape(circle1);
				scene.removeShape(oppCircle1);
			}
			scene.removeShape(circle0);
			scene.removeShape(oppCircle0);
		}
	}
	
	private void initializeRotation() {
		heap.clear();
		for (TriangulationVertex v : vertices) {
			if (testing) System.out.print(v.getId() + ": " + v.toString(2));
			if (v.getSquaredDistance() < Constants.EPSILON) {                //takes care of the special case when the rotation center overlaps with one of the given points
				v.setSquaredPolarRadius(0.0);                            
				v.setPolarRadius(0.0);                        
				v.setPolarAngle(0.0);			
				if (testing) System.out.println(", polar angle: 0.0"); 
				v.setCosAngle(1.0);
				v.setSinAngle(0.0);
				v.setType(TriangulationVertex.VertexType.S);
				rotIndx.remove((Integer)v.getId());
			}
			else {
				v.setSquaredPolarRadius(v.getSquaredDistance());             // remains unchanged when rotating around the same point
				v.setPolarRadius(Math.sqrt(v.getSquaredDistance()));         // remains unchanged when rotating around the same point
				v.setPolarAngle(v.polarAngle());			
				if (testing) System.out.println(", polar angle: " + Functions.toDeg(v.getPolarAngle())); 
				v.setCosAngle(v.polarAngleCos());
				v.setSinAngle(v.polarAngleSin());
			}
		}
		
		if (testing) for (int i = 2; i < vertices.size()-1; i++) 
			scene.addShape(new TextShape(Integer.toString(i), vertices.get(i), 0.3));
		
		// draw the orbits of rotating points
		if (testing) for (int i : rotIndx) scene.addShape(new Circle(new Point(0,0), vertices.get(i).distance()));

		for (int k : rotIndx) {
			for (TriangulationFace t : vertices.get(k).getFaces()) {
				if (t.isAlive()) {
					t.setAlive(false);				
					for (int i = 0; i < 3; i++) {
						TriangulationFace oppT = t.getNeighbor(i);
						if (oppT != null) {
							TriangulationVertex oppV = t.getThirdVertex(oppT);
							if ((oppT.getCount() == 0) || (t.getCorner(i).getId() < oppV.getId())) {
								if (testing) System.out.print(t.toString() + oppT.toString() + " initialization ");
								angles = getRoot(t, oppV);
								if (angles != null) addToHeap(angles, t, oppT, null);
								else if (testing) System.out.println(t.toString() + oppT.toString() + " no solution");
							}
							else if (testing) System.out.println(t.toString() + oppT.toString() + " " + t.getCorner(i).getId() + " > " + oppV.getId() + " and second triangle is not static");
						}
						else if (testing) System.out.println(t.toString()  + " no triangle opposite to vertex " + t.getCorner(i).getId());
					}
				}
				else if (testing) System.out.println(t.toString() + " has been processed earlier");
			}
		}
		for (TriangulationFace t : triangulationFaces) t.setAlive(true);
	}

	private void traceHeapItems(TriangulationFace t, TriangulationFace oppT, HeapItem heapItem) {
		System.out.println("Tracing " + t.toString() + oppT.toString());
		HeapItem item = heapItem;
		while (item != null) {
			System.out.println("created by " + item.t + item.oppT);
			item = item.createdBy;
		}
	}
	
	public void flip(TriangulationFace t013, TriangulationFace t123, double rotAngle, Direction rotDir, HeapItem heapItem){

		flipNr++;

		int p0Indx = t013.getIndex(t123);
		int p2Indx = t123.getIndex(t013);
		TriangulationVertex p0 = t013.getCorner(p0Indx);
		TriangulationVertex p1 = t013.getCorner((p0Indx+1)%3);
		TriangulationVertex p2 = t123.getCorner(p2Indx);
		TriangulationVertex p3 = t123.getCorner((p2Indx+1)%3);

		if (testing) System.out.println(" Flip "+ flipNr + ", Animation angle: " + Functions.toDeg(rotAngle-angleTotal));
		animate(rotAngle-angleTotal);
		angleTotal = rotAngle;

		t013.setAlive(false);
		t013.hide(scene, testing);
		if (circleAnimation) t013.hideCircumCircle(scene);
		triangulationFaces.remove(t013);
		t123.setAlive(false);
		t123.hide(scene, testing);
		if (circleAnimation) t123.hideCircumCircle(scene);
		triangulationFaces.remove(t123);

		TriangulationFace t012 = new TriangulationFace(p0, p1, p2); 
		triangulationFaces.add(t012);
		t012.draw(scene, testing);
		TriangulationFace t023 = new TriangulationFace(p0, p2, p3); 
		triangulationFaces.add(t023);
		t023.draw(scene, testing);

		if (circleAnimation) animateCircles();
		
		TriangulationFace t01 = t013.getNeighbor((p0Indx+2)%3);
		TriangulationFace t12 = t123.getNeighbor((p2Indx+1)%3);
		TriangulationFace t23 = t123.getNeighbor((p2Indx+2)%3);
		TriangulationFace t30 = t013.getNeighbor((p0Indx+1)%3);
	
		t012.setNeighbor(0, t12);
		t012.setNeighbor(1, t023);
		t012.setNeighbor(2, t01);
		t023.setNeighbor(0, t23);
		t023.setNeighbor(1, t30);
		t023.setNeighbor(2, t012);                               
		if (t01 != null) t01.setNeighbor(t01.getIndex(t013), t012);
		if (t12 != null) t12.setNeighbor(t12.getIndex(t123), t012);
		if (t23 != null) t23.setNeighbor(t23.getIndex(t123), t023);
		if (t30 != null) t30.setNeighbor(t30.getIndex(t013), t023);
	                                                         
		if (p0.getFace() == t013) p0.setFace(t012);
		if ((p1.getFace() == t123) || (p1.getFace() == t013)) p1.setFace(t012);
		if (p2.getFace() == t123) p2.setFace(t023);
		if ((p3.getFace() == t013) || (p3.getFace() == t123)) p3.setFace(t023);
	
		t012.setId(t013.getId());
		t023.setId(t123.getId());
			
		if (t01 != null) {
			angles = getRoot(t012, t012.getThirdVertex(t01));
			if (angles != null) addToHeap(angles, t012, t01, heapItem);
		}
		
		if (t12 != null) {
			angles = getRoot(t012, t012.getThirdVertex(t12));
			if (angles != null) addToHeap(angles, t012, t12, heapItem);
		}
		
		angles = getRoot(t012, p3);
		if ((heapItem.getAngle1() != null) && (heapItem.getAngle1() < Constants.TAU)) {
			Double[] angles1 = new Double[2];
			angles1[0] = heapItem.getAngle1();
			angles1[1] = null;
			addToHeap(angles1, t012, t023, heapItem);
		}

		if (t23 != null) {
			angles = getRoot(t023, t023.getThirdVertex(t23));
			if (angles != null) addToHeap(angles, t023, t23, heapItem);
		}
		
		if (t30 != null) {
			angles = getRoot(t023, t023.getThirdVertex(t30));
			if (angles != null) addToHeap(angles, t023, t30, heapItem);
		}
	}
	
	private void rotate() {	
		HeapItem heapItem;
		TriangulationFace t, oppT;
		initializeRotation();
		while (!heap.isEmpty()) {
			heapItem = (HeapItem)heap.extract();
			t = heapItem.getT();
			oppT = heapItem.getOppT();
			if (testing) System.out.print(t.toString() + oppT.toString());
			if (t.isAlive() && oppT.isAlive()) flip(t, oppT, heapItem.getAngle(), rotDir, heapItem);
			if (testing) {
				for (Shape stop : heapItem.getStops()) scene.removeShape(stop);
				if (printLabels) for (Shape label : heapItem.getLabels()) scene.removeShape(label);	
			}
		}
	}
	
	public static void main(String[] args) {
		PointSet points = new PointSet(3000);
		Sorter sort = new SorterQuick();
		sort.Sort(points, new SortToolPoint2dDistance());
		for (int i = 1; i < points.getSize(); i++) {
			System.out.println(points.get(i-1).distance());
			if (points.get(i).distance() - points.get(i-1).distance() < Constants.EPSILON) System.out.println(i);
		}
/*		PointSet points = new PointSet();
		points.insert(new Point(0, 0));
		points.insert(new Point(1, 0));
		points.insert(new Point(0, 1));
		points.insert(new Point(1, 1));
*/
		// add big points
		double t = 5000.0;
		Point big0 = new Point(-t/Math.sqrt(2), -t/Math.sqrt(2));
		Point big1 = new Point(-Constants.EPSILON, t);
		Point bigN = new Point(t, 0.0);
		
		points.insert(big0);
		points.insert(big1);
		points.insert(bigN);
		// magnify 
		for (Point p : points) for (int i = 0; i < 2; i++) p.setCoord(i, 10*p.getCoord(i));

		// compute initial DT
		KineticDTBigEdges kDT = new KineticDTBigEdges(points);
		// initialize
		kDT.setRotationPoint(new Point(0,0));
		kDT.setDirection(KineticDTBigEdges.Direction.CCW);

		double fraction = 0.05;
		List<Integer> rotList = new ArrayList<Integer>();
		Random random = new Random(2);
		for (int i = 2; i < points.getSize()-1; i++) 
			if (random.nextDouble() < fraction) {
				rotList.add(i);
			}
		kDT.setRotVertices(rotList);
		
//		kDT.setRotVertices(150, 350);
		long start = System.nanoTime();
		kDT.rotate();
		long end = System.nanoTime();
		System.out.printf(" time in miliseconds %.2f\n", (end - start)/1000000.0);

	}
}
