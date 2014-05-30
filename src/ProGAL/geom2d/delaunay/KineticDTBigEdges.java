package ProGAL.geom2d.delaunay;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import ProGAL.dataStructures.DLCyclicList;
import ProGAL.dataStructures.Heap;
import ProGAL.dataStructures.SortTool;
import ProGAL.dataStructures.SortToolPoint2dDistance;
import ProGAL.dataStructures.Sorter;
import ProGAL.dataStructures.SorterQuick;
import ProGAL.geom2d.Circle;
import ProGAL.geom2d.Line;
import ProGAL.geom2d.LineSegment;
import ProGAL.geom2d.Point;
import ProGAL.geom2d.PointSet;
import ProGAL.geom2d.Shape;
import ProGAL.geom2d.Triangulation;
import ProGAL.geom2d.TriangulationFace;
import ProGAL.geom2d.TriangulationVertex;
import ProGAL.geom2d.Vector;
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
	private int itNr = 0;
	private J2DScene scene = J2DScene.createJ2DSceneInFrame();
	
	private Heap heap = new Heap(this.vertices.size(), new SortToolHeapItems());
	
//	private boolean testing = true;
	private boolean testingPrint = false;
	private boolean testingScene = true;
	private boolean showOrbits = true;
	private boolean circleAnimation = false;
	private boolean printLabels = true;
	
	private double alpha = 2000000;
//	private double alpha = 0.35;
	
	private class HeapItem {
		private Double[] angles;
		private TriangulationFace t;
		private Object oppT;
		private List<Shape> stops;
		private List<Shape> labels;
		private HeapItem createdBy;
		
		private HeapItem(Double[] angles, TriangulationFace t, Object oppT, List<Shape> stops, List<Shape> labels, HeapItem createdBy) {
			this.angles = angles;
			this.t = t;
			this.oppT = oppT;
			this.stops = stops;
			this.labels = labels;
			this.createdBy = createdBy;
		}
				
		private Double getAngle() { return angles[0];} 
		private TriangulationFace getT() { return t; }
		private Object getOppT() { return oppT; }
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
	}

	public Point getRotationPoint() { return rotationPoint; }
	public void setRotationPoint(Point rotationPoint) { 
		if (this.rotationPoint != null) for (TriangulationVertex v : vertices) v.addThis(this.rotationPoint);
		this.rotationPoint = rotationPoint; 
		for (TriangulationVertex v : vertices) v.subtractThis(rotationPoint);	
	}
	
	public Direction getDirection() { return rotDir; }
	public void setDirection(Direction rotDir) { this.rotDir = rotDir; }
	
	/* specify what is rotating - 4 methods to do it */
	public void setRotVertices(int a) { setRotVertices(a, a); }
	public void setRotVertices(int a, int b) { 
		for (int i = a; i <= b; i++) rotIndx.add(i); 
		setRotVertices();
	}
	public void setRotVertices(double fraction) {
		Random random = new Random(2);
		for (int i = 0; i < vertices.size(); i++) 
			if (!vertices.get(i).isBigPoint() && (random.nextDouble() < fraction)) {
				rotIndx.add(i);
			}
		setRotVertices();	
	}
	public void setRotVertices() { 
		for (TriangulationVertex v : vertices) v.setType(TriangulationVertex.VertexType.S); 
		for (int i : rotIndx) vertices.get(i).setType(TriangulationVertex.VertexType.R);
	}
	
	public double getAlpha() { return alpha; }
	public void setAlpha(double alpha) { this.alpha = alpha; }
	
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


	public void animateEvent(double angle) {
		int steps = 1;
		double angleStep = angle/steps;
		double sin = Math.sin(angleStep);
		double cos = Math.cos(angleStep);

		for (int k = 0; k < steps; k++) {
//			if ((steps > 1) && (k == steps/2) && (!isDelaunay())) { System.out.println("not Delaunay"); System.exit(0); }
			for (int i : rotIndx) vertices.get(i).rotation(cos, sin);
			if (true) { try { if (vertices.size() < 1000) Thread.sleep(1); } catch (InterruptedException e) {} }
			if (circleAnimation) animateCircles();
		}
		int count;
		for (TriangulationFace t : triangulationFaces) {
			count = t.nrRotatingCorners();
			if ((count > 0) && (count < 3)) t.setCircumCircleRadius();
		}
	}
	
	/** returns up two a pair of angles when two edge-sharing faces have their 4 vertices on a common circle */
	private Double[] getRoot(TriangulationFace face, TriangulationVertex oppVertex) {
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
		
	/** returns up two a pair of angles when the circumcircle of a face has radius alpha */
	private Double[] getRoot(TriangulationFace face) {
		int count = face.getCount();
		if ((count == 0) || (count == 7)) return null;
		if (count <= 3) {
			if (count == 1) return getRootSSR(face.getCorner(1), face.getCorner(2), face.getCorner(0), 0);
			if (count == 2) return getRootSSR(face.getCorner(0), face.getCorner(2), face.getCorner(1), 0);
			return getRootSSR(face.getCorner(0), face.getCorner(1), face.getCorner(2), 1);
		}
		else {
			if (count == 4) return getRootSSR(face.getCorner(0), face.getCorner(1), face.getCorner(2), 0);
			if (count == 5) return getRootSSR(face.getCorner(0), face.getCorner(2), face.getCorner(1), 1);
			return getRootSSR(face.getCorner(1), face.getCorner(2), face.getCorner(0), 1);
		}
	}
	
	
	private double getRotAngle(double polarAngleU, double polarAngleC) {
		if (polarAngleC > angleTotal) {
			if (polarAngleU > polarAngleC+Constants.EPSILON) return polarAngleU - polarAngleC + angleTotal;
			if (polarAngleU+Constants.EPSILON < polarAngleC - angleTotal) return Constants.TAU -polarAngleC + polarAngleU + angleTotal;
			return Constants.TAU;
		}
		if (polarAngleU > polarAngleC+Constants.EPSILON) {
			if (polarAngleU > Constants.TAU+ Constants.EPSILON - angleTotal + polarAngleC) return Constants.TAU;
			return polarAngleU - polarAngleC + angleTotal;
		}
		return Constants.TAU;
	}
	// returns up to 2 rotation angles needed to bring vertex C to be 2*alpha away from the vertex u 
	private Double[] getRoot(TriangulationVertex u, TriangulationVertex C) {
		Circle uCircle = new Circle(u, 2*getAlpha());
		Point[] intersections = uCircle.intersections(C.getOrbit());
		if (intersections == null) return null;
		Double[] angles = new Double[2];
		angles[0] = getRotAngle(intersections[0].polarAngle(), C.polarAngle());
		if (intersections[1] != null) {
			angles[1] = getRotAngle(intersections[1].polarAngle(), C.polarAngle()); 
			if (angles[1] < angles[0]) angles[0] = angles[1];
		}
		return angles;
	}
	
	private Double[] getRootSSRR(TriangulationVertex A, TriangulationVertex B, TriangulationVertex C, TriangulationVertex D) {	
		double aa  = A.getSquaredPolarRadius(); 
		double bb  = B.getSquaredPolarRadius(); 
		double cc  = C.getSquaredPolarRadius(); 
		double c   = C.getPolarRadius();
		double dd  = D.getSquaredPolarRadius(); 
		double d   = D.getPolarRadius();
	
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
	
		double bb_cc = bb - cc;
		double cc_aa = cc - aa;
		double aa_bb = aa - bb;
		double m41 = A.y()*bb_cc + B.y()*cc_aa + C.y()*aa_bb;
		double m42 = A.x()*bb_cc + B.x()*cc_aa + C.x()*aa_bb;
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
	
	private Double[] getRootSSR(TriangulationVertex A, TriangulationVertex B, TriangulationVertex C, int dir) {
		double AB2_4 = A.getSquaredDistance(B)/4;
		if (AB2_4 > alpha*alpha) return null;
		Line bisector = Line.getBisectorLine(A, B);
		double a = Math.sqrt(alpha*alpha - AB2_4);
		Vector dVector = bisector.getDirection().scale(a);
		Point center1 = bisector.getPoint().add(dVector);
		Point center2 = bisector.getPoint().subtract(dVector);
		Circle circum1 = new Circle(center1, alpha);
		Circle circum2 = new Circle(center2, alpha);
		
		Point[] intersect1 = circum1.intersections(bisector);
		Point[] intersect2 = circum2.intersections(bisector);
		
		if (circum2.contains(intersect1[0])) intersect1[0] = intersect1[1]; 
		if (circum1.contains(intersect2[0])) intersect2[0] = intersect2[1];
 
		TriangulationVertex C1 = new TriangulationVertex(intersect1[0]);
		C1.setType(TriangulationVertex.VertexType.S);
		C1.setSquaredPolarRadius(C1.getSquaredDistance());

		TriangulationVertex C2 = new TriangulationVertex(intersect2[0]);
		C2.setType(TriangulationVertex.VertexType.S);
		C2.setSquaredPolarRadius(C2.getSquaredDistance());

		Double[] angles1 = getRootSSSR(A, B, C1, C, dir);
		Double[] angles2 = getRootSSSR(A, B, C2, C, dir);
		if (angles1 == null) {
			if (angles2 != null) return angles2;else return null;
		}
		else {
			if (angles2 == null) return angles1;
			else {
				if (angles1[0] < angles2[0]) {
					if (angles2[0] < angles1[1]) angles1[1] = angles2[0];
					return angles1;
				}
				else {
					if (angles1[0] < angles2[1]) angles2[1] = angles1[0]; 
					return angles2;
				}
			}
		}
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

	private void addToHeap(Double[] angles, TriangulationFace t, Object oppT, HeapItem createdBy) {
		if (angles[0] < Constants.TAU) {
			List<Shape> stops = null;
			List<Shape> labels = null;
			if (testingPrint) {
				if (oppT != null) {
					if (oppT instanceof TriangulationFace) System.out.println("   " + t.toString() + oppT.toString() + ", new flip event added, rotate to = " + Functions.toDeg(angles[0])); 
					else System.out.println("   " + t.toString() + " " + ((TriangulationVertex)oppT).getId() + ", new edge length event added, rotate to = " + Functions.toDeg(angles[0]));
				}
				else System.out.println("   " + t.toString() + ", new circum radius event added, rotate to = " + Functions.toDeg(angles[0])); 
			}
			if (showOrbits) stops = markAllStops(angles[0] - angleTotal);
//			if (printLabels) labels = addLabels(t, oppT, angles[0] - angleTotal);	
			heap.insert(new HeapItem(angles, t, oppT, stops, labels, createdBy));
		}
		else if (testingPrint) 
			if (oppT != null) System.out.println("Faces " + t.toString() + " and " + oppT.toString() + ", skipped ");					
//		if (testing) showCircles(t, oppT, angles);
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
			if (testingPrint) System.out.print(v.getId() + ": " + v.toString(2));
			if (v.getSquaredDistance() < Constants.EPSILON) {                //takes care of the special case when the rotation center overlaps with one of the given points
				v.setSquaredPolarRadius(0.0);                            
				v.setPolarRadius(0.0);                        
				v.setPolarAngle(0.0);			
				if (testingPrint) System.out.println(", polar angle: 0.0"); 
				v.setCosAngle(1.0);
				v.setSinAngle(0.0);
				v.setType(TriangulationVertex.VertexType.S);
				rotIndx.remove((Integer)v.getId());
			}
			else {
				v.setSquaredPolarRadius(v.getSquaredDistance());             // remains unchanged when rotating around the same point
				v.setPolarRadius(Math.sqrt(v.getSquaredDistance()));         // remains unchanged when rotating around the same point
				v.setPolarAngle(v.polarAngle());			
				if (testingPrint) System.out.println(", polar angle: " + Functions.toDeg(v.getPolarAngle())); 
				v.setCosAngle(v.polarAngleCos());
				v.setSinAngle(v.polarAngleSin());
			}
		}
		
		// identifies short triangles
		for (TriangulationFace t : triangulationFaces) t.setShort(t.getCircumRadius() <= alpha);

		// draws the initial alha complex
		if (testingScene) {
			if (vertices.size() > 100) printLabels = false;
			draw(scene, getAlpha(), printLabels);
		}
		
		// draw the orbits of rotating points
		if (testingScene) {
			if (vertices.size() > 50) showOrbits = false; 
			if (showOrbits) for (int i : rotIndx) scene.addShape(new Circle(new Point(0,0), vertices.get(i).distance()));
		}

		// add flip events to the heap
		for (int k : rotIndx) {
			for (TriangulationFace t : vertices.get(k).getFaces()) {
				if (t.isAlive()) {
					t.setAlive(false);				
					for (int i = 0; i < 3; i++) {
						TriangulationFace oppT = t.getNeighbor(i);
						if (oppT != null) {
							TriangulationVertex oppV = t.getThirdVertex(oppT);
							if ((oppT.getCount() == 0) || (t.getCorner(i).getId() < oppV.getId())) {
								if (testingPrint) System.out.print(t.toString() + oppT.toString() + " initialization, ");
								angles = getRoot(t, oppV);
								if (angles != null) addToHeap(angles, t, oppT, null);
								else if (testingPrint) System.out.println(t.toString() + oppT.toString() + " no solution");
							}
							else if (testingPrint) System.out.println(t.toString() + oppT.toString() + " " + t.getCorner(i).getId() + " > " + oppV.getId() + " and second triangle is not static");
						}
						else if (testingPrint) System.out.println(t.toString()  + " no triangle opposite to vertex " + t.getCorner(i).getId());
					}
				}
				else if (testingPrint) System.out.println(t.toString() + " has been processed earlier");
			}
		}
		// clean up
		for (int k: rotIndx) {
			for (TriangulationFace t : vertices.get(k).getFaces()) t.setAlive(true);
		}
		
		// add circumradius length events
		for (int k : rotIndx) {
			TriangulationVertex v = vertices.get(k);
			for (TriangulationFace t : v.getFaces()) {
				if (t.isAlive()) {
					t.setAlive(false);
					angles = getRoot(t);
					if (angles != null) {
						if (testingPrint) {
							System.out.println(t.toString() + ": circumradius length event added to heap. Current radius of circumcircle is " +
											   t.getCircumRadius() + ". Circumradius will be " +
											   getAlpha() + " by rotating by angle " + Functions.toDeg(angles[0]));
						}
						addToHeap(angles, t, null, null);
					}
				}
			}
		}
		// clean up
		for (int k: rotIndx) 
			for (TriangulationFace t : vertices.get(k).getFaces()) t.setAlive(true);

		// add edge length events
		for (int k : rotIndx) {
			TriangulationVertex v = vertices.get(k);
			for (TriangulationFace t : v.getFaces()) {
				int indx = t.getIndex(v);
				TriangulationVertex vk = t.getCorner((indx+1)%3);
				addLengthEvent(t, v, vk);
			}
		}
	}

	private void addLengthEvent(TriangulationFace t, TriangulationVertex v, TriangulationVertex vk) {
		Double [] angles;
		if (!vk.isBigPoint() && (vk.getType() == TriangulationVertex.VertexType.S)) {
			angles = getRoot(vk, v);
			if (angles != null) {
				if (testingPrint) {
					System.out.println("[" + vk.getId() + "," + v.getId() + "]: length event added to heap. Current edge length is " +
							v.distance(vk) + " Edge length will be " +
							(2*getAlpha()) + " by rotating " + v.getId() + " by angle " + Functions.toDeg(angles[0]));
				}
				addToHeap(angles, t, v, null);
			}
		}

	}
	private void traceHeapItems(TriangulationFace t, TriangulationFace oppT, HeapItem heapItem) {
		System.out.println("Tracing " + t.toString() + oppT.toString());
		HeapItem item = heapItem;
		while (item != null) {
			System.out.println("created by " + item.t + item.oppT);
			item = item.createdBy;
		}
	}
	
	/** rotate to next circumradius length event and update the alpha complex if necessary */
	private void lengthEvent(TriangulationFace t, double rotAngle) {
		if (testingPrint) System.out.println(", Animation angle: " + Functions.toDeg(rotAngle-angleTotal));
		animateEvent(rotAngle-angleTotal);
		angleTotal = rotAngle;
	
		t.setShort(!t.isShort());	
		if (testingScene)
			if (t.isShort()) t.draw(scene, Color.black); else t.draw(scene, Color.red);
		
		angles = getRoot(t);
		if (angles != null) addToHeap(angles, t, null, null);

	}

	/** rotate to next edge-length event and update the alpha complex if necessary */
	private void lengthEvent(TriangulationFace t, TriangulationVertex v, TriangulationVertex vOpp, double rotAngle) {
		if (testingPrint) System.out.println(", Animation angle: " + Functions.toDeg(rotAngle-angleTotal));
		animateEvent(rotAngle-angleTotal);
		angleTotal = rotAngle;
	
		int indxT = t.getIndex(v);
		TriangulationFace oppT = t.getNeighbor((indxT+2)%3);
		int indxOppT = oppT.getIndex(vOpp);
		if (testingScene) {
			Shape shape = t.getEdgeShape(indxT);
			if (shape != null) {
				scene.removeShape(shape);
				scene.removeShape(oppT.getEdgeShape(indxOppT));
			}
			else {
				t.setEdgeShape(indxT, new LineSegment(t.getCorner(indxT), t.getCorner((indxT+1)%3)));
				scene.addShape(t.getEdgeShape(indxT), Color.red, 0.01);
				oppT.setEdgeShape(indxOppT, new LineSegment(oppT.getCorner(indxOppT), oppT.getCorner((indxOppT+1)%3)));
				scene.addShape(oppT.getEdgeShape(indxOppT), Color.red, 0.01);
			}
		}
		
		angles = getRoot(vOpp, v);
		if (angles != null) addToHeap(angles, t, v, null);

	}

	
	
	/** rotate to next flip event, flip, and identify new events */
	private void flipEvent(TriangulationFace t013, TriangulationFace t123, double rotAngle, Direction rotDir, HeapItem heapItem){
		int p0Indx = t013.getIndex(t123);
		int p2Indx = t123.getIndex(t013);
		TriangulationVertex p0 = t013.getCorner(p0Indx);
		TriangulationVertex p1 = t013.getCorner((p0Indx+1)%3);
		TriangulationVertex p2 = t123.getCorner(p2Indx);
		TriangulationVertex p3 = t123.getCorner((p2Indx+1)%3);

		// rotate points to the flip event
		if (testingPrint) System.out.print(", Animation angle: " + Functions.toDeg(rotAngle-angleTotal));
		animateEvent(rotAngle-angleTotal);
		
		angleTotal = rotAngle;

		// remove flip-out triangles
		t013.setAlive(false);
		t123.setAlive(false);
		if (testingScene) { 
			t013.hide(scene);
			if (circleAnimation) t013.hideCircumCircle(scene);
			t123.hide(scene);
			if (circleAnimation) t123.hideCircumCircle(scene);
		}
		triangulationFaces.remove(t013);
		triangulationFaces.remove(t123);

		// add flip-in triangles
		TriangulationFace t012 = new TriangulationFace(p0, p1, p2); 
		triangulationFaces.add(t012);
		TriangulationFace t023 = new TriangulationFace(p0, p2, p3); 
		triangulationFaces.add(t023);
		
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
			
		if (testingScene) {
			t012.draw(scene, getAlpha(), testingScene);
			t023.draw(scene, getAlpha(), testingScene);
			if (circleAnimation) animateCircles();
		}

		
		// identify new flip events
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
		

		
		// identify new circum radius length events
		angles = getRoot(t012);
		if (angles != null) addToHeap(angles, t012, null, null);
		angles = getRoot(t023);
		if (angles != null) addToHeap(angles, t023, null, null);
		
		// identify new edge length events
		if (p0.getType() == TriangulationVertex.VertexType.R) addLengthEvent(t012, p0, p1);
		if (p1.getType() == TriangulationVertex.VertexType.R) addLengthEvent(t012, p1, p2);
		if (p2.getType() == TriangulationVertex.VertexType.R) addLengthEvent(t012, p2, p0);
		if (p0.getType() == TriangulationVertex.VertexType.R) addLengthEvent(t023, p0, p2);
		if (p2.getType() == TriangulationVertex.VertexType.R) addLengthEvent(t023, p2, p3);
		if (p3.getType() == TriangulationVertex.VertexType.R) addLengthEvent(t023, p3, p0);
	}
	
	
	/** keep rotating the points until the heap is empty */
	private void rotate() {	
		HeapItem heapItem;
		Object oppO;
		TriangulationFace t, oppT;
		TriangulationVertex v, oppV;
		initializeRotation();
		itNr = 0;
		while (!heap.isEmpty()) {
			heapItem = (HeapItem)heap.extract();
			t = heapItem.getT();
			if (t.isAlive()) {
				oppO = heapItem.getOppT();
				if (oppO != null) {
					if (oppO instanceof TriangulationFace) {
						oppT = (TriangulationFace)oppO;
						if (oppT.isAlive()) {
							if (testingPrint) System.out.print("iteration nr.: " + ++itNr + ": flip event: " + t.toString() + oppT.toString());
							flipEvent(t, oppT, heapItem.getAngle(), rotDir, heapItem);
						}
					}
					else {
						v = (TriangulationVertex)oppO;
						oppV = t.getCorner((t.getIndex(v)+1)%3);
						if (testingPrint) System.out.print("iteration nr.: " + ++itNr + ": edge length event: [" + oppV.getId() + "," + v.getId() + "]");
						lengthEvent(t, v, oppV, heapItem.getAngle());
					}
				}
				else {
					if (testingPrint) System.out.print("iteration nr.: " + ++itNr + ": circum radius length event: " + t.toString());
					lengthEvent(t, heapItem.getAngle());
				}
				if (testingScene) {
					if (showOrbits) for (Shape stop : heapItem.getStops()) scene.removeShape(stop);
				}
			}
		}
	}
	
	public static void main(String[] args) {
		PointSet points = new PointSet(1000);
		
		// points are sorted by their distance from the origo
		Sorter sort = new SorterQuick();
		sort.Sort(points, new SortToolPoint2dDistance());
		
		// check for overlapping points
		int size = points.getSize();
		for (int i = 0; i < size; i++) {
//			System.out.println((i) + ": " + points.get(i).distance());
			if ((i < size-1) && (points.get(i+1).distance() - points.get(i).distance() < Constants.EPSILON)) System.out.println(i);
		}
		
		// add three big points
		double t = 5000.0;		
		points.insert(new Point(-t/Math.sqrt(2), -t/Math.sqrt(2)));
		points.insert(new Point(-Constants.EPSILON, t));
		points.insert(new Point(t, 0.0));
		
		// magnify 
		for (Point p : points) for (int i = 0; i < 2; i++) p.setCoord(i, 10*p.getCoord(i));

		// compute initial DT
//		J2DScene scene =  J2DScene.createJ2DSceneInFrame();
		KineticDTBigEdges kDT = new KineticDTBigEdges(points);

		// initialize
		kDT.setRotationPoint(new Point(0,0));
		kDT.setDirection(KineticDTBigEdges.Direction.CCW);		
		
		// points that will rotate are picked randomly. Big points cannot rotate.
		kDT.setRotVertices(0.25);
		
		// rotate
		long start = System.nanoTime();
		kDT.rotate();
		long end = System.nanoTime();
		System.out.printf(" time in miliseconds %.2f\n", (end - start)/1000000.0);

	}
}
