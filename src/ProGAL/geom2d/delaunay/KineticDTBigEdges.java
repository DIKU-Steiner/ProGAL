package ProGAL.geom2d.delaunay;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import ProGAL.dataStructures.Heap;
import ProGAL.dataStructures.SortTool;
import ProGAL.geom2d.Circle;
import ProGAL.geom2d.Point;
import ProGAL.geom2d.PointSet;
import ProGAL.geom2d.Polygon;
import ProGAL.geom2d.Shape;
import ProGAL.geom2d.Triangulation;
import ProGAL.geom2d.TriangulationFace;
import ProGAL.geom2d.TriangulationVertex;
import ProGAL.geom2d.TriangulationVertex.VertexType;
import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.geom2d.viewer.TextShape;
import ProGAL.math.Constants;
import ProGAL.math.Functions;
import ProGAL.math.Trigonometry;




public class KineticDTBigEdges extends Triangulation{
	
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
		private double angle;
		private TriangulationFace t;
		private TriangulationFace oppT;
		private List<Shape> stops;
		private List<Shape> labels;
		
		private HeapItem(double angle, TriangulationFace t, TriangulationFace oppT, List<Shape> stops, List<Shape> labels) {
			this.angle = angle;
			this.t = t;
			this.oppT = oppT;
			this.stops = stops;
			this.labels = labels;
		}
		
		
		
		private Double getAngle() { return angle;} 
		private TriangulationFace getT() { return t; }
		private TriangulationFace getOppT() { return oppT; }
		private List<Shape> getStops() { return stops; }
		private void addStop(Shape stop) { stops.add(stop); }
		private List<Shape> getLabels() { return labels; }
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
		draw(scene); 
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


	public void animate(double alpha, Direction dir) {
		int steps = 10;
		double angleStep = alpha/steps;
		double sin = Math.sin(angleStep);
		if (dir == Direction.CW) sin = -sin;
		double cos = Math.cos(angleStep);

		for (int k = 0; k < steps; k++) {
			for (int i : rotIndx) vertices.get(i).rotation(cos, sin);
			if (true) {
				try {
					if (vertices.size() < 10000) Thread.sleep(1);
				} catch (InterruptedException e) {}
			}
			if (circleAnimation) animateCircles();
			scene.repaint();
		}
	}

	
	public double getRoot(TriangulationFace face, TriangulationFace oppFace, TriangulationVertex oppVertex) {
		int dirCCW, dirCW;
		int count;
		if (rotDir == Direction.CCW) { dirCCW = 0; dirCW = 1; } else { dirCCW = 1; dirCW = 0; }

		count = face.getCount();
		if (oppVertex.getType() == TriangulationVertex.VertexType.R) count = count + 8;	
		
		if ((count == 0) || (count == 15)) return Constants.TAU;
		if (count <= 7) {
			if (count <= 3) {
				if (count == 1) return getRootSSSR(face.getCorner(1), face.getCorner(2), oppVertex, face.getCorner(0), dirCCW);
				if (count == 2) return getRootSSSR(face.getCorner(0), face.getCorner(2), oppVertex, face.getCorner(1), dirCCW);
				return getRootSSRR(face.getCorner(2), oppVertex, face.getCorner(0), face.getCorner(1));
			}
			else {
				if (count <= 5) {
					if (count == 4) return getRootSSSR(face.getCorner(0), face.getCorner(1), oppVertex, face.getCorner(2), dirCCW);
					return getRootSSRR(face.getCorner(1), oppVertex, face.getCorner(0), face.getCorner(2));
				}
				else {
					if (count == 6) return getRootSSRR(face.getCorner(0), oppVertex, face.getCorner(1), face.getCorner(2));
					return getRootSSSR(face.getCorner(0), face.getCorner(1), face.getCorner(2), oppVertex, dirCW);
				}
			}
		}
		else {
			if (count <= 11) {
				if (count <= 9) {
					if (count == 8) return getRootSSSR(face.getCorner(0), face.getCorner(1), face.getCorner(2), oppVertex, dirCCW);
					return getRootSSRR(face.getCorner(1), face.getCorner(2), face.getCorner(0), oppVertex);
				}
				else {
					if (count == 10) return getRootSSRR(face.getCorner(0), face.getCorner(2), face.getCorner(1), oppVertex);
					return getRootSSSR(face.getCorner(0), face.getCorner(1), oppVertex, face.getCorner(2), dirCW);
				}
			}
			else {
				if (count <= 13) {
					if (count == 12) return getRootSSRR(face.getCorner(0), face.getCorner(1), face.getCorner(2), oppVertex);
					return getRootSSSR(face.getCorner(0), face.getCorner(2), oppVertex, face.getCorner(1), dirCW);
				}
				else return getRootSSSR(face.getCorner(1), face.getCorner(2), oppVertex, face.getCorner(0), dirCW);
			}
		}
	}
		
	private double getRootSSRR(TriangulationVertex A, TriangulationVertex B, TriangulationVertex C, TriangulationVertex D) {	
		double aa  = A.getSquaredPolarRadius(); 
		double bb  = B.getSquaredPolarRadius(); 
		double cc  = C.getSquaredPolarRadius(); double c = C.getPolarRadius();
		double dd  = D.getSquaredPolarRadius(); double d = D.getPolarRadius();
	
		double m11 = aa - bb;
		double m12 = A.y() - B.y();
		double m13 = A.y()*bb - B.y()*aa;
		double m22 = A.x() - B.x();
		double m23 = A.x()*bb - B.x()*aa;
		double m33 = A.x()*B.y() - A.y()*B.x();

		double cosAlpha = C.getCosAngle();
		double sinAlpha = C.getSinAngle();
		double cosBeta = D.getCosAngle();
		double sinBeta = D.getSinAngle();
		double coefSin = d*(cosBeta*(m23 - cc*m22)  + sinBeta*(m13 - cc*m12)) + 
				         c*(cosAlpha*(dd*m22 - m23) + sinAlpha*(dd*m12 - m13));
		double coefCos = d*(cosBeta*(cc*m12 - m13)  + sinBeta*(m23 - cc*m22)) + 
				         c*(cosAlpha*(m13 - dd*m12) + sinAlpha*(dd*m22 - m23));
		double coef    = d*c*m11*Math.sin(D.getPolarAngle()-C.getPolarAngle()) + (cc - dd)*m33;
	
		angles = Trigonometry.solveAsinXPlusBcosXplusC(coefSin, coefCos, coef);
		return getRotAngle(angles, 0);
	}
		
	
	public double getRotAngle(Double[] angles, int dir) {
		if (angles == null) return Constants.TAU;
		if (angles[1] == null) {
			if (dir == 1) {
				angles[0] = Constants.TAU - angles[0] + angleTotal;
				if (angles[0] > Constants.TAU) angles[0] = angles[0] - Constants.TAU;
			}
			if (rotDir == Direction.CCW) {
				if (angles[0] < angleTotal + Constants.EPSILON) return Constants.TAU; else return angles[0];
			}
			else {
				if (Constants.TAU - angles[0] < angleTotal + Constants.EPSILON) return Constants.TAU; else return Constants.TAU - angles[0];
			}
		}
		else {
			if (dir == 1) { 
				angles[0] = Constants.TAU - angles[0] + angleTotal;
				if (angles[0] > Constants.TAU) angles[0] = angles[0] - Constants.TAU;
				angles[1] = Constants.TAU - angles[1] + angleTotal;
				if (angles[1] > Constants.TAU) angles[1] = angles[1] - Constants.TAU;
			}
			if (angles[1] < angles[0]) { double tmp = angles[1]; angles[1] = angles[0]; angles[0] = tmp; }
			if (rotDir == Direction.CCW) {
				if (angles[0] < angleTotal + Constants.EPSILON) {
					if (angles[1] < angleTotal + Constants.EPSILON) return Constants.TAU; else return angles[1];
				}
				else return angles[0];
			}
			else {
				if (Constants.TAU - angles[1] < angleTotal + Constants.EPSILON) {
					if (Constants.TAU - angles[0] < angleTotal + Constants.EPSILON) return Constants.TAU; else return Constants.TAU - angles[0];
				}
				else return Constants.TAU - angles[1];
			}
		}	
	}	
		
	private double getRootSSSR(TriangulationVertex A, 
							   TriangulationVertex B, 
							   TriangulationVertex C, 
							   TriangulationVertex D, int dir) {
		double aa  = A.getSquaredDistance(); 
		double bb  = B.getSquaredDistance(); 
		double cc  = C.getSquaredDistance();
		double dd  = D.getSquaredDistance(); double d = Math.sqrt(dd);
	
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
		if (testing && (dir == 1) && (angles != null)) {
			Point A0 = A.rotationClone(Constants.TAU-angles[0]);
			Point B0 = B.rotationClone(Constants.TAU-angles[0]);
			Point C0 = C.rotationClone(Constants.TAU-angles[0]);
			Shape markA0 = new Circle(A0, 0.1);
			Shape markB0 = new Circle(B0, 0.1);
			Shape markC0 = new Circle(C0, 0.1);
			scene.addShape(markA0, Color.red);
			scene.addShape(markB0, Color.red);
			scene.addShape(markC0, Color.red);
			Shape circle0 = new Circle(A0, B0, C0);
			scene.addShape(circle0, Color.green);
			Shape markA1 = null;
			Shape markB1 = null;
			Shape markC1 = null;
			Shape circle1 = null;
			if (angles[1] != null) {
				Point A1 = A.rotationClone(Constants.TAU-angles[1]);
				Point B1 = B.rotationClone(Constants.TAU-angles[1]);
				Point C1 = C.rotationClone(Constants.TAU-angles[1]);
				markA1 = new Circle(A1, 0.1);
				markB1 = new Circle(B1, 0.1);
				markC1 = new Circle(C1, 0.1);
				scene.addShape(markA1, Color.red);
				scene.addShape(markB1, Color.red);
				scene.addShape(markC1, Color.red);			
				circle1 = new Circle(A1, B1, C1);
				scene.addShape(circle1, Color.green);
			}
			scene.removeShape(markA0);
			scene.removeShape(markB0);
			scene.removeShape(markC0);
			scene.removeShape(circle0);
			scene.removeShape(markA1);
			scene.removeShape(markB1);
			scene.removeShape(markC1);
			scene.removeShape(circle1);
		}
		return getRotAngle(angles, dir);
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

	private void addToHeap(double angle, TriangulationFace t, TriangulationFace oppT) {
		if (angle < Constants.TAU) {
			List<Shape> stops = null;
			List<Shape> labels = null;
			if (testing) {
				System.out.println(t.toString() + oppT.toString() + ", rotate to = " + Functions.toDeg(angle)); 
				stops = markAllStops(angle - angleTotal);
				if (printLabels) labels = addLabels(t, oppT, angle - angleTotal);
			}
			heap.insert(new HeapItem(angle, t, oppT, stops, labels));
		}
		else if (testing) System.out.println("Faces " + t.toString() + " and " + oppT.toString() + ", skipped ");					

	}
	
	private void initializeRotation() {
		heap.clear();
		for (TriangulationVertex v : vertices) {
			if (testing) System.out.print(v.getId() + ": " + v.toString(2));
			v.setSquaredPolarRadius(v.distanceSquared());          // remains unchanged when rotating around the same point
			v.setPolarRadius(Math.sqrt(v.getSquaredDistance()));   // remains unchanged when rotating around the same point
			v.setPolarAngle(v.polarAngle());			
			if (testing) System.out.println(", polar angle: " + Functions.toDeg(v.getPolarAngle())); 
			v.setCosAngle(v.polarAngleCos());
			v.setSinAngle(v.polarAngleSin());
		}
		
		if (testing) for (int i = 2; i < vertices.size()-1; i++) 
			scene.addShape(new TextShape(Integer.toString(i), vertices.get(i), 0.3));
		
		// draw the orbits of rotating points
		if (testing) for (int i : rotIndx) scene.addShape(new Circle(new Point(0,0), vertices.get(i).distance()));

		for (int k : rotIndx) {
			TriangulationVertex v = vertices.get(k);
			for (TriangulationFace t : v.getFaces()) {
				if (t.isAlive()) {
					t.setAlive(false);				
					for (int i = 0; i < 3; i++) {
						TriangulationFace oppT = t.getNeighbor(i);
						if (oppT != null) {
							TriangulationVertex oppV = t.getThirdVertex(oppT);
							if ((oppT.getCount() == 0) || (t.getCorner(i).getId() < oppV.getId())) {
								if (testing) System.out.print("Initialization: " + t.toString() + oppT.toString());
								double angle = getRoot(t, oppT, oppV);
								addToHeap(angle, t, oppT);
								if (testing && (angle !=Constants.TAU)) {
									Shape circle = t.drawRotatedCircumCircle(scene, Color.blue, angle, rotIndx);
									Shape oppCircle = oppT.drawRotatedCircumCircle(scene, Color.magenta, angle, rotIndx);
									scene.removeShape(circle);
									scene.removeShape(oppCircle);
								}							
							}
						}
					}
				}
			}
		}
		for (TriangulationFace t : triangulationFaces) t.setAlive(true);
	}
	
	private boolean convex(TriangulationFace t013, TriangulationFace t123){
		int p0Indx = t013.getIndex(t123);
		int p2Indx = t123.getIndex(t013);
		TriangulationVertex p0 = t013.getCorner(p0Indx);
		TriangulationVertex p1 = t013.getCorner((p0Indx+1)%3);
		TriangulationVertex p2 = t123.getCorner(p2Indx);
		TriangulationVertex p3 = t123.getCorner((p2Indx+1)%3);

		if (Point.rightTurn(p0, p1, p2) || Point.leftTurn(p0, p3, p2) || Point.rightTurn(p2, p3, p0) || Point.leftTurn(p2, p1, p0)) {
			System.out.println(flipNr + " Back animation angle: " + Functions.toDeg(angleTotal));
			return false;
//			System.exit(0);
		}
		return true;
	}
	
	public void flip(TriangulationFace t013, TriangulationFace t123, double rotAngle, Direction rotDir){

		flipNr++;
		if (testing && (flipNr == 11)) {
			System.out.println("Stopping here");
		}
		int p0Indx = t013.getIndex(t123);
		int p2Indx = t123.getIndex(t013);
		TriangulationVertex p0 = t013.getCorner(p0Indx);
		TriangulationVertex p1 = t013.getCorner((p0Indx+1)%3);
		TriangulationVertex p2 = t123.getCorner(p2Indx);
		TriangulationVertex p3 = t123.getCorner((p2Indx+1)%3);

		if (testing) System.out.println(" Flip "+ flipNr + ", Animation angle: " + Functions.toDeg(rotAngle-angleTotal));
		animate(rotAngle-angleTotal, rotDir);

		
		if (Point.rightTurn(p0, p1, p2) || Point.leftTurn(p0, p3, p2) || Point.rightTurn(p2, p3, p0) || Point.leftTurn(p2, p1, p0)) {
			animate(angleTotal-rotAngle, rotDir);
			System.out.println(flipNr + " Back animation angle: " + Functions.toDeg(angleTotal-rotAngle));
			
			J2DScene scene = J2DScene.createJ2DSceneInFrame();
			scene.addShape(new Polygon(t013.getCorners()),Color.BLACK, 0.0004);
			scene.addShape(new Polygon(t123.getCorners()),Color.BLACK, 0.0004);
			scene.addShape(new Polygon(t123.getNeighbor(0).getCorners()),Color.GRAY, 0.0004);
			scene.addShape(new Polygon(t123.getNeighbor(1).getCorners()),Color.GRAY, 0.0004);
			scene.addShape(new Polygon(t123.getNeighbor(2).getCorners()),Color.GRAY, 0.0004);
			scene.addShape(new Polygon(t013.getNeighbor(0).getCorners()),Color.GRAY, 0.0004);
			scene.addShape(new Polygon(t013.getNeighbor(1).getCorners()),Color.GRAY, 0.0004);
			scene.addShape(new Polygon(t013.getNeighbor(2).getCorners()),Color.GRAY, 0.0004);
			scene.addShape(new Circle(t013.getCorner(0),0.0006), t013.getCorner(0).getType()==VertexType.R?Color.RED:Color.GRAY, 0, true);
			scene.addShape(new Circle(t013.getCorner(1),0.0006), t013.getCorner(1).getType()==VertexType.R?Color.RED:Color.GRAY, 0, true);
			scene.addShape(new Circle(t013.getCorner(2),0.0006), t013.getCorner(2).getType()==VertexType.R?Color.RED:Color.GRAY, 0, true);
			scene.addShape(new Circle(t123.getCorner(0),0.0006), t123.getCorner(0).getType()==VertexType.R?Color.RED:Color.GRAY, 0, true);
			scene.addShape(new Circle(t123.getCorner(1),0.0006), t123.getCorner(1).getType()==VertexType.R?Color.RED:Color.GRAY, 0, true);
			scene.addShape(new Circle(t123.getCorner(2),0.0006), t123.getCorner(2).getType()==VertexType.R?Color.RED:Color.GRAY, 0, true);
			for(TriangulationVertex v: t013.getCorners()){
				if(v.getType()==VertexType.R){
					scene.addShape(v.getOrbit(new Point(0,0)), Color.BLUE, 0.0002);
				}
			}
			for(TriangulationVertex v: t123.getCorners()){
				if(v.getType()==VertexType.R){
					scene.addShape(v.getOrbit(new Point(0,0)), Color.BLUE, 0.0002);
				}
			}
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			
			
			
			throw new RuntimeException("No msg now");
//			System.exit(0);
		} else {
			
			
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
			t013 = t012;
			t123 = t023;
				
			TriangulationVertex v;
			double angle;

			if (t01 != null) {
				v = t012.getThirdVertex(t01);
				angle = getRoot(t012, t01, v);
				addToHeap(angle, t012, t01);
				if (testing && (angle != Constants.TAU)) {
					Shape circle = t012.drawRotatedCircumCircle(scene, Color.blue, angle-angleTotal, rotIndx);
					Shape oppCircle = t01.drawRotatedCircumCircle(scene, Color.magenta, angle-angleTotal, rotIndx);
					scene.removeShape(circle);
					scene.removeShape(oppCircle);
				}
			}
			
			if (t12 != null) {
				v = t012.getThirdVertex(t12);
				angle = getRoot(t012, t12, v);
				addToHeap(angle, t012, t12);
				if (testing && (angle != Constants.TAU)) {
					Shape circle = t012.drawRotatedCircumCircle(scene, Color.blue, angle-angleTotal, rotIndx);
					Shape oppCircle = t12.drawRotatedCircumCircle(scene, Color.magenta, angle-angleTotal, rotIndx);
					scene.removeShape(circle);
					scene.removeShape(oppCircle);
				}

			}
			
			angle = getRoot(t012, t023, p3);
			addToHeap(angle, t012, t023);
			System.out.println(angle);
			if (testing && (angle != Constants.TAU)) {
				Shape circle = t012.drawRotatedCircumCircle(scene, Color.blue, angle-angleTotal, rotIndx);
				Shape oppCircle = t023.drawRotatedCircumCircle(scene, Color.magenta, angle-angleTotal, rotIndx);
				scene.removeShape(circle);
				scene.removeShape(oppCircle);
			}

			
			if (t23 != null) {
				v = t023.getThirdVertex(t23);
				angle = getRoot(t023, t23,  v);
				addToHeap(angle, t023, t23);
				if (testing && (angle != Constants.TAU)) {
					Shape circle = t023.drawRotatedCircumCircle(scene, Color.blue, angle-angleTotal, rotIndx);
					Shape oppCircle = t23.drawRotatedCircumCircle(scene, Color.magenta, angle-angleTotal, rotIndx);
					scene.removeShape(circle);
					scene.removeShape(oppCircle);
				}

			}
			
			if (t30 != null) {
				v = t023.getThirdVertex(t30);
				angle = getRoot(t023,t30, v);
				addToHeap(angle, t023, t30);
				if (testing && (angle != Constants.TAU)) {
					Shape circle = t023.drawRotatedCircumCircle(scene, Color.blue, angle-angleTotal, rotIndx);
					Shape oppCircle = t30.drawRotatedCircumCircle(scene, Color.magenta, angle-angleTotal, rotIndx);
					scene.removeShape(circle);
					scene.removeShape(oppCircle);
				}
		
			}
		}
	}


	
	private void rotate() {		
		initializeRotation();
		if (testing && circleAnimation) animateCircles();
		while (!heap.isEmpty()) {
			HeapItem heapItem = (HeapItem)heap.extract();
			double angle = heapItem.getAngle();
			TriangulationFace t = heapItem.getT();
			TriangulationFace oppT = heapItem.getOppT();
			if (testing) System.out.print(t.toString() + oppT.toString());
			if (t.isAlive() && oppT.isAlive()) {
				
				flip(t, oppT, angle, rotDir);
			}
			if (testing) {
				List<Shape> stops = heapItem.getStops();
				for (Shape stop : stops) scene.removeShape(stop);
				if (printLabels) {
					List<Shape> labels = heapItem.getLabels();
					for (Shape label : labels) scene.removeShape(label);	
				}
			}
		}
	}
	
	public static void main(String[] args) {
		PointSet points = new PointSet(7000);
/*		PointSet points = new PointSet();
		points.insert(new Point(0.411, 0.745));
		points.insert(new Point(0.586, 0.3));
		points.insert(new Point(0.62, 0.28));
		points.insert(new Point(0.723, 0.206));
*/
		// add big points
		double t = 1000.0;
		Point big0 = new Point(-t/Math.sqrt(2), -t/Math.sqrt(2));
		Point big1 = new Point(-1/t, t);
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

/*		List<Integer> rotList = new ArrayList<Integer>();
		rotList.add(3);
		rotList.add(5);
		kDT.setRotVertices(rotList);
*/
		kDT.setRotVertices(110,126);
		kDT.rotate();

	}
}
