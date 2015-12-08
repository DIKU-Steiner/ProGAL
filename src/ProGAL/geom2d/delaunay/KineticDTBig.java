package ProGAL.geom2d.delaunay;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ProGAL.geom2d.Circle;
import ProGAL.geom2d.Point;
import ProGAL.geom2d.PointSet;
import ProGAL.geom2d.Shape;
import ProGAL.geom2d.Triangulation;
import ProGAL.geom2d.TriangulationFace;
import ProGAL.geom2d.TriangulationVertex;
import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.geom2d.viewer.TextShape;
import ProGAL.math.Constants;
import ProGAL.math.Functions;
import ProGAL.math.Trigonometry;



public class KineticDTBig extends Triangulation {

	List<Integer> rotIndx = new ArrayList<Integer>();
	J2DScene scene = J2DScene.createJ2DSceneInFrame();
	final double angleStep = 0.0001;
	int nrRounds;
	Double[] angles = new Double[2];
	final double cos = Math.cos(angleStep);
	final double sin = Math.sin(angleStep);
	int k; // index of the first rotating vertex
	public static enum Direction  { CW, CCW };
	private boolean testing = true;
	List<Shape> markShapes = new ArrayList<Shape>();
	Direction rotDir;

	
	public KineticDTBig(PointSet points) {
		super(points, TriangulationAlgorithm.Delaunay);	
		scene.removeAllShapes();
	}
		
	public void rotate(Point rotPoint, Direction dir) {
		initializeBeforeFirstRotations(rotPoint);
		double angle = 0.0;
		double alpha;
		nrRounds = 0;
		while (angle < 2*Math.PI) {
			if (testing) testPrint();
			alpha = getNextFlip(dir);
			System.out.println("Round " + ++nrRounds + " alpha = " + alpha);
			angle = angle + alpha;
			updateAfterRotation();
		}
	}
	
	public void testPrint() {
		System.out.print("Faces: ");
		for (TriangulationFace face : triangulationFaces) System.out.print(face.toString());
		System.out.println();
		System.out.print("Opposite faces: ");
		for (TriangulationFace face : triangulationFaces) {
			System.out.print("[");
			for (int i = 0; i < 3; i++) {
				TriangulationFace neighbor = face.getNeighbor(i); 
				if (neighbor == null) System.out.print(" null "); 
				else System.out.print(face.getNeighbor(i).toString());
			}
			System.out.print("] ");
		}
		System.out.println();
		System.out.print("First faces: ");
		for (TriangulationVertex v : vertices) System.out.print(v.getFace().toString());
		System.out.println();
	}
	
	private void initializeBeforeFirstRotations(Point rotPoint) {
		// translate so rotPoint is in the origo
		for (TriangulationVertex v : vertices) v.subtractThis(rotPoint);
		// vertices are either stationary or rotating
		for (TriangulationVertex v : vertices) v.setType(TriangulationVertex.VertexType.S); 
		for (int i : rotIndx) vertices.get(i).setType(TriangulationVertex.VertexType.R);
		// compute polar coordinates
		for (TriangulationVertex v : vertices) {
			v.setSquaredPolarRadius(v.distanceSquared());          // remains unchanged when rotating around the same point
			v.setPolarRadius(Math.sqrt(v.getSquaredDistance()));   // remains unchanged when rotating around the same point
			v.setPolarAngle(v.polarAngle());
			v.setCosAngle(v.polarAngleCos());
			v.setSinAngle(v.polarAngleSin());
		}

		for (TriangulationVertex v : vertices) {
			for (TriangulationVertex u : vertices) {
				if ((v != u) && Math.abs(v.getPolarRadius() - u.getPolarRadius()) < Constants.EPSILON) System.out.println(v.getId() + " " + u.getId());
			}
		}
		// draw DT (with circumcircles and vertex numbers if TESTING == TRUE)
		if (scene != null) {
			draw(scene);
			scene.autoZoom();
			scene.centerCamera();
			if (testing) {
				drawNewCircles(); 
				if (vertices.size() < 200) {
					for (int i = 0; i < vertices.size(); i++) scene.addShape(new TextShape(Integer.toString(i), vertices.get(i), 0.3));
				}
			}
			// draw the orbits of rotating points
			if (testing) for (int i : rotIndx) scene.addShape(new Circle(new Point(0,0), vertices.get(i).distance()));
		}
	}

	private void updateAfterRotation() {
		for (int i : rotIndx)  {
			TriangulationVertex v = vertices.get(i);
			v.setPolarAngle(v.polarAngle());
			v.setCosAngle(v.polarAngleCos());
			v.setSinAngle(v.polarAngleSin());
		
		}
		if (testing) for (Shape shape : markShapes) scene.removeShape(shape);
	}
	
	public void drawNewCircles() {
		for (TriangulationFace face : triangulationFaces) {
			scene.removeShape(face.getCircumCircle());
			if (face.isFlat()) face.setCircumCircle(null);
			else {
				face.setCircumCircle(new Circle(face.getCorner(0), face.getCorner(1), face.getCorner(2)));
				if (face.getCircumCircle().getRadius() < 100.0) scene.addShape(face.getCircumCircle(), Color.red);
			}
		}
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
		double angleStep = alpha/100;
		double sin = Math.sin(angleStep);
		double cos = Math.cos(angleStep);

		for (int k = 0; k < 100; k++) {
			if (dir == Direction.CCW) for (int i : rotIndx) vertices.get(i).rotation(cos, sin);
			else for (int i : rotIndx) vertices.get(i).rotation(cos, -sin);
			if (testing) animateCircles();
			if (true) {
				try {
					if (vertices.size() < 10000) Thread.sleep(5);
				} catch (InterruptedException e) {}
			}
			scene.repaint();
		}
	}
	
	public double getNextFlip(Direction dir) {
		double delta = Constants.TAU;
		TriangulationFace selectedFace = null;
		TriangulationFace selectedOppFace = null;
		// visit faces one by one
		for (TriangulationFace face : triangulationFaces) {
			int fcount = face.getCount();
			for (int i = 0; i < 3; i++) {
				TriangulationVertex v  = face.getCorner(i);
				TriangulationFace oppFace = face.getNeighbor(i);
				if (oppFace != null) {
					TriangulationVertex oppVertex = face.getThirdVertex(oppFace);
					if (oppVertex.getId() > v.getId()) {
						int vcount = fcount;
						if (oppVertex.getType() == TriangulationVertex.VertexType.R) vcount = fcount + 8;						
						double angle = getRoot(face, oppVertex, vcount, dir);
						if (testing && (angle < Constants.TAU)) {
							if ((vcount == 7) || (vcount == 11) || (vcount == 13) || (vcount == 14)) markAllStops(false); else markAllStops(true);
							System.out.println("Face: " + face.toString() + ", vertex: " + face.getCorner(i).getId() + ", entering vertex: " + oppVertex.getId() + 
										", rotation angle = " + Functions.toDeg(angle) + ", type 0");
						}
						if (angle < delta) {
							delta = angle;
							selectedFace = face;
							selectedOppFace = oppFace;
						}
					}
				}
			}
		}
		if (testing) System.out.println("Selected face: " + selectedFace.toString() + 
				", Selected opposite face: " + selectedOppFace.toString() + ", delta = " + Functions.toDeg(delta));
		if (delta <= Constants.EPSILON) System.out.println("ZERO");
		animate(delta, dir);
		if (testing) {
			Shape sh = new Circle(selectedFace.getCorner(0),selectedFace.getCorner(1),selectedFace.getCorner(2));
			scene.addShape(sh, Color.green);
			scene.removeShape(sh);
		}
//		TriangulationFace oppFace = selectedFace.getOppFace(selectedVertex);
		flip(selectedFace, selectedOppFace, false, scene, testing);
		return delta;
	}
	
	
	public double getRoot(TriangulationFace face, TriangulationVertex oppVertex, int count, Direction dir) {
		int dirCCW, dirCW;
		if (dir == Direction.CCW) { dirCCW = 0; dirCW = 1; } else { dirCCW = 1; dirCW = 0; }

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
		return getRotAngle(angles);
	}
	
	public double getRotAngle(Double[] angles) {
		if (angles == null) return Constants.TAU;
		if ((angles[0] != null)  && (angles[0] > Constants.TAU - Constants.EPSILON)) angles[0] = 0.0;
		if ((angles[1] != null)  && (angles[1] > Constants.TAU - Constants.EPSILON)) angles[1] = 0.0;
		
		if (angles[1] == null) {
			if (angles[0] < Constants.EPSILON) return Constants.TAU; else {
				if (rotDir == Direction.CCW) return angles[0]; else return Constants.TAU - angles[0];
			}
		}
		if (angles[0] < Constants.EPSILON) {
			if (angles[1] < Constants.EPSILON) return Constants.TAU; else {
				if (rotDir == Direction.CCW) return angles[1]; else return Constants.TAU - angles[1];
			}
		}
		if (angles[1] < Constants.EPSILON) return angles[0];
		if (angles[0] < angles[1]) {
			if (rotDir == Direction.CCW) return angles[0]; else return Constants.TAU - angles[1]; 
		}
		else {	
			if (rotDir == Direction.CCW) return angles[1]; else return Constants.TAU - angles[0];
		}
	}
	
	public double getRotAngle(Double[] angles, int dir) {
		if (angles == null) return Constants.TAU;
		if ((angles[0] != null)  && (angles[0] > Constants.TAU - Constants.EPSILON)) angles[0] = 0.0;
		if ((angles[1] != null)  && (angles[1] > Constants.TAU - Constants.EPSILON)) angles[1] = 0.0;

		if (angles[1] == null) { 
			if (angles[0] < Constants.EPSILON) return Constants.TAU; else return angles[0]; 
		}
		if (dir == 0) {
			if (angles[0] < Constants.EPSILON) {
				if (angles[1] < Constants.EPSILON) return Constants.TAU; else return angles[1];
			}
			if (angles[0] < angles[1]) return angles[0]; 
			else {
				if (angles[1] < Constants.EPSILON) return Constants.TAU; else return angles[1];
			}
		}
		else {
			if (angles[0] < Constants.EPSILON) {
				if (angles[1] < Constants.EPSILON) return Constants.TAU; else return Constants.TAU - angles[1];
			}
			if (angles[0] < angles[1]) return Constants.TAU - angles[1]; else return Constants.TAU - angles[0];			
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
		return getRotAngle(angles, dir);
	}
	

	public void setRotVertices(int a) { setRotVertices(a, a); }
	public void setRotVertices(int a, int b) {
		for (int i = a; i <= b; i++) rotIndx.add(i);
	}
	public void setRotVertices(List<Integer> rotList) {
		rotIndx = rotList;
	}
	
	private void markAllStops(boolean direction) {
		double angle;
		if (angles[0] != null) {
			if (direction) angle = angles[0]; else angle = Constants.TAU - angles[0];
			for (int k : rotIndx) {
				Shape mark0 = new Circle(vertices.get(k).rotationClone(angle), 0.1);
				scene.addShape(mark0, Color.blue);
				markShapes.add(mark0);
			}
		}
		if (angles[1] != null) {
			if (direction) angle = angles[1]; else angle = Constants.TAU - angles[1];
			for (int k : rotIndx) {
				Shape mark1 = new Circle(vertices.get(k).rotationClone(angle), 0.1);
				scene.addShape(mark1, Color.magenta);
				markShapes.add(mark1);
			}
		}
	}
	
	public static void main(String[] args) {
		PointSet points = new PointSet();
		points.insert(new Point(0.02, 0.005));
//		points.insert(new Point(0, 0.5));
		points.insert(new Point(-0.01, 1.02));
//		points.insert(new Point(0.5, 0));
//		points.insert(new Point(0.5, 0.5));
//		points.insert(new Point(0.5, 1));
		points.insert(new Point(1.03, -0.03));
//		points.insert(new Point(1, 0.5));
		points.insert(new Point(1.04, 1.05));
		
		points.insert(new Point(100.0, 0.0));
		points.insert(new Point(-100.0,-100));
		points.insert(new Point(-1.0,100));
		
		for (Point p : points) for (int i = 0; i < 2; i++) p.setCoord(i, 10*p.getCoord(i));
		KineticDTBig kDT = new KineticDTBig(points);
		if (kDT.testing) {  //checks if there are points at the same distance from the origo
			points.sort();
			Point q = null;
			int i = -1;
			for (Point p : points) {
				i++;
				System.out.print(i + " ");
				if ((q != null) && (p.distance() - q.distance() < Constants.EPSILON)) System.out.println("cocircircular");
				q = p;
			}
		}
		List<Integer> rotList = Arrays.asList(2,6);
		kDT.setRotVertices(3, 3);
		Point rotPoint = new Point(4,2);  // rotation point
		kDT.rotDir = Direction.CCW;
		kDT.rotate(rotPoint, kDT.rotDir);

	}
}
