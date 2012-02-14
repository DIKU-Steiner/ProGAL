package ProGAL.geom3d.complex.delaunayComplex;

import java.awt.Color;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import com.sun.org.apache.xml.internal.serialize.OutputFormat.DTD;



import ProGAL.geom3d.Circle;
import ProGAL.geom3d.Line;
import ProGAL.geom3d.LineSegment;
import ProGAL.geom3d.Plane;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Shape;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.complex.CEdge;
import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.complex.CTriangle;
import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.complex.delaunayComplex.Flips.ApexConfig;
import ProGAL.geom3d.predicates.ExactJavaPredicates;
import ProGAL.geom3d.predicates.Predicates;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Cylinder;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.geom3d.volumes.Tetrahedron;
import ProGAL.math.Constants;
import ProGAL.math.Matrix;
import ProGAL.math.Matrix3x3;



public class DelaunayComplexKinetic {
//	private Vector textVector = new Vector(0.05, 0.05, 0.05);
//	private Color tetraColor = new Color(0,255,0,100);
//	private final Flips flips;
//	private final Predicates predicates;
//
//
//	private DelaunayComplex dt;
//
//	public DelaunayComplexKinetic() {
//		this.predicates = new ExactJavaPredicates();
//		this.flips = new Flips(predicates);
//	} 
//	
//	public DelaunayComplexKinetic(DelaunayComplex dt) {
//		this.dt = dt;
//		this.predicates = new ExactJavaPredicates();
//		this.flips = new Flips(predicates);
//
//	}
//	
//	public static double toDegrees(double alpha) { return 180*alpha/Math.PI; }
//	public static double toRadians(double alpha) { return Math.PI*alpha/180; }
//
//
//	public double flipRotation(Line line, CTetrahedron tetr) {
//		double alpha = 99999.0;
//		CVertex a = tetr.getPoint(0);
//		CVertex b = tetr.getPoint(1);
//		CVertex c = tetr.getPoint(2);
//		CVertex d = tetr.getPoint(3);
//		for (int i = 0; i < 4; i++) {
//			CTetrahedron tetrNeighbour = tetr.getNeighbour(i);
//			if (!tetrNeighbour.containsBigPoint()) {
//				CTriangle tr = tetr.getTriangle(i);
//				CVertex e = tetrNeighbour.oppositeVertex(tr);
//				double c1[][] = {{b.y()-e.y(),b.z()-e.z(),b.dot(b)-e.dot(e)},
//						         {c.y()-e.y(),c.z()-e.z(),c.dot(c)-e.dot(e)},
//						         {d.y()-e.y(),d.z()-e.z(),d.dot(d)-e.dot(e)}}; 
//				double c2[][] = {{b.x()-e.x(),b.z()-e.z(),b.dot(b)-e.dot(e)},
//				                 {c.x()-e.x(),c.z()-e.z(),c.dot(c)-e.dot(e)},
//				                 {d.x()-e.x(),d.z()-e.z(),d.dot(d)-e.dot(e)}}; 
//				double c3[][] = {{b.x()-e.x(),b.z()-e.z(),b.dot(b)-e.dot(e)},
//				                 {c.x()-e.x(),c.z()-e.z(),c.dot(c)-e.dot(e)},
//				                 {d.x()-e.x(),d.z()-e.z(),d.dot(d)-e.dot(e)}}; 
//				double c4[][] = {{b.x()-e.x(),b.y()-e.y(),b.z()-e.z()},
//				                 {c.x()-e.x(),c.y()-e.y(),c.z()-e.z()},
//				                 {d.x()-e.x(),d.y()-e.y(),d.z()-e.z()}}; 
//				Matrix3x3 matr1 = new Matrix3x3(c1);
//				Matrix3x3 matr2 = new Matrix3x3(c2);
//				Matrix3x3 matr3 = new Matrix3x3(c3);
//				Matrix3x3 matr4 = new Matrix3x3(c4);
//				double minor1 = matr1.determinant();
//				double minor2 = matr2.determinant();
//				double minor3 = matr3.determinant();
//				double minor4 = matr4.determinant();
//				double rhs = (e.x()*minor1 - e.y()*minor2 + e.z()*minor3 + (circle.getRadius()*circle.getRadius() - e.dot(e))*minor4)/circle.getRadius();
//				double dd = 4*minor2*minor2 - 4*(-minor1-rhs)*(minor1-rhs);
//				if (dd > 0) {
//					double t1 = (2*minor2+Math.sqrt(dd))/2*(-minor1-rhs);
//					double t2 = (2*minor2-Math.sqrt(dd))/2*(-minor2-rhs);
//					double a1 = Math.atan2((1-t1*t1)/(1+t1*t1), 2*t1/(1+t1*t1));
//					double a2 = Math.atan2((1-t2*t2)/(1+t2*t2), 2*t2/(1+t2*t2));
//					if (a1 < a2) {
//						if (a1 < alpha) alpha = a1;
//					}
//					else {
//						if (a2 < alpha) alpha = a2;
//					}
//				}
//				else {
//					if (Math.abs(dd) < Constants.EPSILON) {
//						double t1 = minor2/(-minor1-rhs);
//						double a1 = Math.atan2((1-t1*t1)/(1+t1*t1), 2*t1/(1+t1*t1));
//						if (a1 < alpha) alpha = a1;
//					}
//					else alpha = 9999.0;
//				}
//			}
//		}
//		return alpha;
//	}
//	
//	class NextListener1 implements KeyListener{
//		public void keyPressed(KeyEvent e) {		}
//		public void keyReleased(KeyEvent e) {	
//			if (e.getKeyChar()=='n'){
//				rotatingVector = normalToCircle.rotateIn(rotatingVector, alpha);
//				Point movedPoint = circle.getCenter().add(rotatingVector);
//				movingVertex.setCoord(movedPoint);
//				
//				List<CEdge> edges = tetr.getPoint(0).getAdjacentEdges();
//				for (CEdge edge : edges) {
//					edge.setA(circle.getCenter().add(rotatingVector));
//				}
////				tetr.setPoint(0, circle.getCenter().add(rot));
////				tetr.getPoint(0).set(circle.getCenter().add(rot));
////				sphere. = new Sphere(tetr.getCorners());
//				sphere.setCenter(tetr.getCenter());
//				sphere.setRadius(sphere.getCenter().distance(movedPoint));
//				sphere.getCenter().toScene(scene, 0.05, Color.yellow);
//
//
//				sphere.toConsole(4);
//				scene.autoZoom();
//			}	
//		}
//		public void keyTyped(KeyEvent e) {		}
//	}
//
//	public void drawKineticEdges() {
//		movingEdges.clear();
//		for (CEdge edge : dt.getEdges()) {
//			CVertex u = edge.getPoint(0);
//			CVertex w = edge.getPoint(1);
//			Cylinder cyl = new Cylinder(edge, 0.05);
//			if (u.equals(v[0])) {
//				cyl = new Cylinder(edge, 0.05);
//				movingEdges.add(cyl);
//			}
//			else {
//				if (w.equals(v[0])) {
//					LineSegment seg = new LineSegment(v[0], edge.getPoint(0));
//					cyl = new Cylinder(seg, 0.05);
//					movingEdges.add(cyl);
//				}
//				else cyl = new Cylinder(edge, 0.05);
//			}
//			scene.addShape(cyl, tetraColor);
//		}
//	}
//	
//	class NextListener implements KeyListener{
//		HashMap<Object, Point> map = new HashMap<Object, Point>();
//		Plane circlePlane; 
//		Sphere spheres[] = new Sphere[6];
//
//		
//		public void keyPressed(KeyEvent e) {		}
//		public void keyReleased(KeyEvent e) {	
//			if (e.getKeyChar() =='n'){
//				scene.removeAllShapes();
//
//				
//				
//				// construct Delaunay complex
//				dt = new DelaunayComplex(p); 
//				System.out.println("dt size before flip = " + dt.getTetrahedra().size() );
//				for (int i = 0; i < p.size(); i++) v[i] = dt.getVertex(i);
//				
//				drawKineticEdges();
//
//				for (int i = 1; i < p.size(); i++) {
//					spheres[i] = new Sphere(v[i], 0.05);
//					scene.addShape(spheres[i], Color.blue);
//					scene.addText(Integer.toString(i), v[i].add(textVector));
//				}
//				spheres[0] = new Sphere(v[0], 0.1); // this point rotates
// 				scene.addShape(spheres[0], Color.yellow);
//				scene.addText("0", v[0].add(textVector));
//				
//				
//				// construct the orbit of the rotating point
//				rotatingVector = new Vector(v[3], v[0]);
//				rr = rotatingVector.dot(rotatingVector);
//				normalToCircle = rotatingVector.getOrthonormal(); // rotation around normalToCircle rooted at v[3]
//				circle = new Circle(v[3], v[0], normalToCircle);
//				circlePlane = new Plane(circle.getCenter(), circle.getNormal());
//				circle.toScene(scene,0.02, 72);
//				
//				// identify triangles bounding the rotating point
//				List<CTriangle> boundingTriangles = v[0].getOppositeTriangles();
//				List<CTetrahedron> interiorBoundingTetrahedra = new ArrayList<CTetrahedron>();
//				List<CTriangle> exteriorBoundingTriangles = new ArrayList<CTriangle>();
//				for (CTriangle triangle : boundingTriangles) {
//					CTetrahedron tetr0 = triangle.getAdjacentTetrahedron(0);   
//					CTetrahedron tetr1 = triangle.getAdjacentTetrahedron(1);
//					if (tetr0.containsBigPoint() || tetr1.containsBigPoint()) {
//						exteriorBoundingTriangles.add(triangle);
//						triangle.toScene(scene, new Color(0,0,255));
//					}
//					else {
//						if (tetr0.containsPoint(v[0])) interiorBoundingTetrahedra.add(tetr1);
//						else interiorBoundingTetrahedra.add(tetr0);
//						triangle.toScene(scene, new Color(0,255,0));
//					}
//				}
//				
//				// identify intersections with spheres circumscribing surrounding tetrahedra
//				int k = 0;
//				for (CTetrahedron tetr : interiorBoundingTetrahedra) {
//					sphere = new Sphere(tetr);
//					scene.addShape(sphere, new Color(255,0,0,100));
//					Point[] intPoints = sphere.getIntersections(circle);
//					if (intPoints == null) System.out.println("no intersections");
//					else {
//						if (intPoints.length == 1) System.out.println("touch-point");
//						else {
//							k++;
//							for (int i = 0; i < 2; i++) {
//								if (i == 0) intPoints[i].toScene(scene, 0.1, Color.red); else intPoints[i].toScene(scene, 0.1, Color.pink);
//								scene.addText(Integer.toString(k).toString(), intPoints[i].add(textVector));
//								map.put(tetr, intPoints[i]);
//							}
//						}
//					}
//				}
//				
//				// identify intersections with planes through exterior surrounding triangles
//				for (CTriangle triangle : exteriorBoundingTriangles) {
//					Plane plane = new Plane(triangle.getP1(), triangle.getP2(), triangle.getP3());
//					Point[] intPoints = plane.getIntersection(circle);
//					if (intPoints == null) System.out.println("no intersections");
//					else {
//						if (intPoints.length == 1) System.out.println("touch-point");
//						else {
//							k++;
//							for (int i = 0; i < 2; i++) {
//								if (i == 0) intPoints[i].toScene(scene, 0.1, Color.magenta); else intPoints[i].toScene(scene, 0.1, Color.cyan);
//								scene.addText(Integer.toString(k).toString(), intPoints[i].add(textVector));
//								map.put(triangle, intPoints[i]);
//							}
//						}
//					}
//				}
////				scene.repaint();
//				scene.autoZoom();
//				
//				
//
//			}		
//
//			if (e.getKeyChar() == 'l') {
//				Object object = null;
//				Point point = null;
//				double bestAlpha = 7.0;
//				for (Object obj : map.keySet()) {
//					Point p = map.get(obj);
//					Vector v = new Vector(circle.getCenter(), map.get(obj));
//					double alpha = Math.atan2(rotatingVector.cross(v).length()/rr, rotatingVector.dot(v)/rr);
//					if (alpha < 0) alpha = 2*Math.PI + alpha;
//					System.out.println(alpha);
//					if (alpha < bestAlpha) {
//						bestAlpha = alpha;
//						object = obj;
//						point = p;
//					}
//				}
//
//				int did = -1;
//				int pid = -1;
//				CTetrahedron tetr = null;
//				CTetrahedron oppTetr = null;
//				ApexConfig flipcase;
//				if (object instanceof CTetrahedron) {
//					tetr = (CTetrahedron)object;
//					oppTetr = tetr.findNeighbour(v[0]);  // returns the neighbour containing v[0]
//					CVertex p = tetr.findVertex(oppTetr);
//					pid = tetr.getID(p);
//					did = oppTetr.getID(v[0]); // this does not change - should be moved out
//					flipcase = flips.apexConfig(tetr, p, pid, oppTetr, v[0]);
//				}
//
//				
//				Sphere targetSphere = new Sphere(point, 0.3);
//				scene.addShape(targetSphere, new Color(100,100,100,100));
//				
//				for (int i = 0; i <= 1002; i++) {
//					scene.removeShape(spheres[0]);
//					rotatingVector = normalToCircle.rotateIn(rotatingVector, bestAlpha/1000);
//					v[0].setCoord(v[3].add(rotatingVector));
//					scene.addShape(spheres[0], Color.yellow);
//					for (Cylinder cyl : movingEdges) {
//						scene.removeShape(cyl);
//						LineSegment seg = cyl.getSegment();
//						seg.getA().setCoord(spheres[0].getCenter());
//						scene.addShape(cyl, tetraColor);
//					}
//				}
//				flips.getFlip23().flip23(tetr, pid, did);
//				System.out.println(flips.getFlipstack().size());
//				CTetrahedron tmp = flips.fixDelaunay();
//				dt.completeComplex();
//				System.out.println("dt size after flip = " + dt.getTetrahedra().size() );
//
//				p.get(0).setCoord(v[0]);
//			}
//		}
//		public void keyTyped(KeyEvent e) {		}
//	}
//
//	
//	Vector normalToCircle;
//	Vector rotatingVector;
//	double rr;
//	Vector rotateToVector;
//	Circle circle;
//	J3DScene scene;
//	Line rotationLine;
//	Sphere sphere;
//	double alpha = Math.PI/36;
//	CVertex movingVertex;
//	List<Point> p = new PointList(); 
//	CVertex[] v = new CVertex[6];
//	List<Cylinder> movingEdges = new ArrayList<Cylinder>();
//
//
//	public static void main(String[] args) {
//		System.out.println(Math.atan2(Math.sin(2*Math.PI-Math.PI/4),Math.cos(2*Math.PI-Math.PI/4)));
//		DelaunayComplexKinetic dck = new DelaunayComplexKinetic();
//		dck.p.add(new Point(0,0,0));
//		dck.p.add(new Point(4,0,0));
//		dck.p.add(new Point(0,3,0));
//		dck.p.add(new Point(1,1,3));
//		dck.p.add(new Point(0,0,3));
//		dck.p.add(new Point(4,2,5));
//
//		dck.scene = J3DScene.createJ3DSceneInFrame();
//		dck.scene.getCanvas().addKeyListener(dck.new NextListener());
//
//	}
//	
//	
//	
//	
//	
//	
//	public static void main1(String[] args) {
//		// get points and compute DT
////		List<Point> points = PointList.generatePointsInCube(6);
//		
//		Point p1 = new Point(1,0,0);
//		Point p2 = new Point(1,1,0);
//		Point p3 = new Point(2,2,3);
//		Point p4 = new Point(3,4,2);
//		Point p5 = new Point(0,0,3);
//		Point p6 = new Point(2,2,2);
//		
//		List<Point> points = new PointList();
//		points.add(p1); points.add(p2); points.add(p3); points.add(p4); points.add(p5); points.add(p6);
//		
//		DelaunayComplexKinetic dck = new DelaunayComplexKinetic(new DelaunayComplex(points));
//		dck.scene = J3DScene.createJ3DSceneInFrame();
//		dck.scene.getCanvas().addKeyListener(dck.new NextListener());
//		
//		// draw all points
//		for (int i = 1; i < dck.dt.getVertices().size(); i++) dck.dt.getVertices().get(i).toScene(dck.scene, 0.02, Color.blue);	
//		
//		// define rotation line
//		dck.normalToCircle = new Vector(0,0,1);
//		dck.rotationLine = new Line(dck.normalToCircle);
//		
//		// pick tetrahedron		
//		dck.tetr = dck.dt.getTetrahedra().get(1);
//		dck.scene.addShape(dck.tetr, Color.yellow);
////		dck.tetr.toScene(dck.scene, 0.02, Color.gray);
//		dck.tetr.toConsole();
//
//		// get sphere
//		dck.sphere = new Sphere(dck.tetr.getCorners());
//		dck.sphere.getCenter().toScene(dck.scene, 0.05, Color.yellow);
//		dck.sphere.toConsole();
//		dck.scene.addShape(dck.sphere, new Color(0,255,255,200));
//
//		// get point
//		dck.movingVertex = dck.tetr.getPoint(0);
//		dck.movingVertex.toScene(dck.scene, 0.07, Color.red);
//		
//		// get circle through the selected point
//		Point center = dck.rotationLine.orthogonalProjection(dck.movingVertex);
//		dck.circle = new Circle(center, center.distance(dck.movingVertex), dck.normalToCircle);
//		dck.circle.toScene(dck.scene, 0.02, 64);
//		dck.rotatingVector = new Vector(dck.circle.getCenter(), dck.movingVertex);
////		dck.alpha = dck.flipRotation(dck.rotationLine, dck.tetr);
//	}
}
