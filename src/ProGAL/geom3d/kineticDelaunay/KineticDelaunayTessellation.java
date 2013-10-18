package ProGAL.geom3d.kineticDelaunay; 

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;

import ProGAL.dataStructures.Heap;
import ProGAL.dataStructures.SortTool;
import ProGAL.geom3d.Circle;
import ProGAL.geom3d.LineSegment;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.kineticDelaunay.Hole;
import ProGAL.geom3d.kineticDelaunay.Hole.Face;
import ProGAL.geom3d.kineticDelaunay.Vertex.VertexType;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.geom3d.volumes.Tetrahedron;
import ProGAL.geom3d.Vector;
import ProGAL.math.Constants;
import ProGAL.math.Functions;
import ProGAL.math.Matrix;
import ProGAL.math.Randomization;
import ProGAL.math.Trigonometry;
import ProGAL.proteins.PDBFile;



public class KineticDelaunayTessellation {
	private static enum Direction { CW, CCW }
	public static enum ProblemInstanceType { random, pdb, pdbNCC, toy }

	private ProblemInstanceType instanceType;
	private List<Vertex> vertices = new ArrayList<Vertex>();
	private List<Tet> tets = new LinkedList<Tet>();
	
	private Tet lastTet;
	private Double[] angles = new Double[2];
	private double angleTotal = 0.0;
	private double angleLimit = Constants.TAU;        
	
	private Point rotationPoint;   // rotation center
	private Vector rotationAxis; // rotation vector
	private Direction rotDir;
	private List<Integer> rotIndx = new ArrayList<Integer>();     // indicies of vertices that move

	private Heap heap = new Heap(this.vertices.size(), new SortToolHeapItems());
//	private boolean testing = true;
	private boolean testingPrint = false;
	private boolean testingScreen = true;
	private boolean sphereAnimation = false;
	J3DScene scene;

	Triangle commonTriangle = null;
	Triangle commonTriangle02 = null;
	Triangle commonTriangle12 = null;
	Tetrahedron commonTriangleShape = null;
	Tetrahedron commonTriangleShape02 = null;
	Tetrahedron commonTriangleShape12 = null;
	LineSegment apexSegment = null;
	LSS apexLSS = null;

	private double alpha;
	
	private class HeapItem {
		private Double[] angles;
		private Tet t;
		private Tet nt;
		
		private HeapItem(Double[] angles, Tet t, Tet nt) {
			this.angles = angles;
			this.t = t;
			this.nt = nt;
		}
				
		private Double[] getAngles() { return angles;} 
		private Tet getT()  { return t; }
		private Tet getNT() { return nt; }
	}
	
	private class SortToolHeapItems implements SortTool {
		public int compare(Object x1, Object x2) {
			if ((x1 instanceof HeapItem) && (x2 instanceof HeapItem)) {
				double d1 = ((HeapItem)x1).getAngles()[0];
				double d2 = ((HeapItem)x2).getAngles()[0];
				if (d1 < d2) return COMP_LESS;
				else { if (d1 > d2) return COMP_GRTR; else return COMP_EQUAL; }
			}
			else throw SortTool.err1;
		}
	}

	private class HeapItemDelete<T> {
		private double power;
		private Face f;
		private Face oppFace;
		
		private HeapItemDelete(Face f, Face oppFace, double power) {
			this.power = power;
			this.f = f;
			this.oppFace = oppFace;
		}
				
		private double getPower() { return power;} 
		private Face getF() { return f; }
		private Face getOppFace() { return oppFace; }
	}
	private class SortToolHeapItemsDelete implements SortTool {
		public int compare(Object x1, Object x2) {
			if ((x1 instanceof HeapItemDelete) && (x2 instanceof HeapItemDelete)) {
				double d1 = ((HeapItemDelete)x1).getPower();
				double d2 = ((HeapItemDelete)x2).getPower();
				if (d1 < d2) return COMP_LESS;
				else { if (d1 > d2) return COMP_GRTR; else return COMP_EQUAL; }
			}
			else throw SortTool.err1;
		}
	}

	
	public KineticDelaunayTessellation(double alpha, ProblemInstanceType type) {
		setInstanceType(type);
		List<Point> points = getInstance(getInstanceType());
		if (testingScreen) scene  = J3DScene.createJ3DSceneInFrame();

		Tetrahedron bigT = Tetrahedron.regularTetrahedron();
		bigT.blowUp(1000);
		lastTet = new BigTet(bigT, vertices);
		int j = 0;
		for (Point p: points) {
//			if (j++ > 12) break;
			insertPoint(p);	
		}
		this.setAlpha(alpha);
		
		/** Each vertex gets a pointer to one of its faces */
		for (Tet tet : tets) {
			for (int i = 0; i < 4; i++) tet.getCorner(i).setTet(tet);
		}
		
		if (getInstanceType() == ProblemInstanceType.random) {
			setRotationPoint(new Point(0.5,0.5,0.5));
			setRotationAxis(new Vector(0, 0, 1));
			setDirection(KineticDelaunayTessellation.Direction.CCW);
			setRotVertices(vertices.size()/2,vertices.size()-1);
		}
		else {
			Vector dir = new Vector(vertices.get(10), vertices.get(11));
			double scaleFactor = 1/dir.length();
			this.setAlpha(alpha*scaleFactor);
			dir.scaleToLengthThis(1);
			Vector newDir = new Vector(0,0,1);
			Vector cross = dir.cross(newDir).scaleToLength(1);
			double angle = Vector.getAngle(dir, newDir);
			Point vOld = vertices.get(10).clone();
			Vertex v;
			for (int i = 0; i < vertices.size(); i++) {
				v = vertices.get(i);
				v.translateThis(vOld);
				v.scaleThis(scaleFactor);
				v.rotation(cross, angle);
			}
			setRotationPoint(vertices.get(10));
			setRotationAxis(newDir);
			setDirection(KineticDelaunayTessellation.Direction.CCW);
			setRotVertices(12, 13);
		}
	}
	
	/** Adds new event to the heap */
	private void addToHeap(Double[] angles, Tet t, Tet nt) {
		heap.insert(new HeapItem(angles, t, nt));
	}
	
	/** Rotates and draws kinetic DT */
	public void animate(J3DScene scene, double alpha) {
		int steps = 1;
		double angleStep = alpha/steps;
		for (int k = 0; k < steps; k++) {
			for (int i = 0; i < rotIndx.size(); i++) {
				Vertex v = vertices.get(rotIndx.get(i));
				v.rotation(rotationAxis, angleStep);
			}
			try { Thread.sleep(1000); } catch (InterruptedException e) {}
			if (testingScreen && sphereAnimation) animateSpheres();			
		}
	}

	
	/** Deletes vertex v and all incident tetrahedra. Determines DT of the remaining vertices. */
	public void delete(Vertex u, J3DScene scene) {
		boolean testing = true;
		Face f, oppFace, thirdFace;
		Tet newTet, tet, faceOppTet, oppFaceOppTet;
		Vertex a, b, c, d;
		double gamma;
		Color clr;
		Sphere sphere = null;
		
		Hole hole = new Hole(this, u, scene, testing);
		if (testing)
			for (Face face : hole.faces) face.tet.fromSceneFace(scene, face.tet.indexOf(u));
		
		Heap heap = new Heap(vertices.size(), new SortToolHeapItemsDelete());

		// heap containing pairs of faces surrounding vertex u is set up
		for (Face face: hole.faces) {
			if (testing) face.shape = face.tet.toSceneFace(scene, face.tet.indexOf(u), Color.pink);
			for (int i = 0; i < 3; i++) {
				oppFace = face.neighbors[i];
				if (testing) {
					oppFace.shape = oppFace.tet.toSceneFace(scene, oppFace.tet.indexOf(u), Color.yellow);
				}
				if (!oppFace.processed) {
					gamma = Point.orientation(face.vertices[0], face.vertices[1], face.vertices[2], oppFace.getFreeVertex(face));
					if (gamma > 0.0) {
						double beta = Point.inSphere(face.vertices[0], face.vertices[1], face.vertices[2], oppFace.getFreeVertex(face), u);
						heap.insert(new HeapItemDelete(face, oppFace, -beta/gamma));
						if (testing) { 
							System.out.println(face.toString() + oppFace.toString() + " " + (beta/gamma));
							sphere = new Sphere(face.vertices[0], face.vertices[1], face.vertices[2], oppFace.getFreeVertex(face));
							sphere.toScene(scene, new Color(255,0,0,50));
						}
					}
				}
				if (testing) {
					scene.removeShape(oppFace.shape);
					scene.removeShape(sphere);
				}
			}
			face.processed = true;
			if (testing) scene.removeShape(face.shape);
		}

		if (testing)
			for (Face face : hole.faces) {
				face.oppTet.toSceneFace(scene, face.oppTet.indexOf(face.vertices[0], face.vertices[1], face.vertices[2]), Color.red);
			}
		
		for (Face face: hole.faces) face.processed = false;
		
		Face face;
		while (hole.faces.size() != 4) {
			
			HeapItemDelete item = null;
			
			item = (HeapItemDelete)heap.extract();
			f = item.getF();
			oppFace = item.getOppFace();
			faceOppTet = f.oppTet;
			if (!f.processed && !oppFace.processed) {
				oppFaceOppTet = oppFace.oppTet;
				thirdFace = f.commonNeighbor(oppFace);
				if ((thirdFace != null) && !thirdFace.processed) {
					a = f.getFreeVertex(oppFace);
					b = oppFace.getFreeVertex(thirdFace);
					c = thirdFace.getFreeVertex(f);
					d = f.getCommonVertex(oppFace, thirdFace);
					hole.vertices.remove(d);
					hole.faces.remove(f); 
					f.processed = true;
					hole.faces.remove(oppFace);
					oppFace.processed = true;
					hole.faces.remove(thirdFace);
					thirdFace.processed = true;
					if (testing) {
						f.oppTet.fromSceneFace(scene, f.oppTet.indexOf(f.vertices[0], f.vertices[1], f.vertices[2]));
						oppFace.oppTet.fromSceneFace(scene, oppFace.oppTet.indexOf(oppFace.vertices[0], oppFace.vertices[1], oppFace.vertices[2]));
						thirdFace.oppTet.fromSceneFace(scene, thirdFace.oppTet.indexOf(thirdFace.vertices[0], thirdFace.vertices[1], thirdFace.vertices[2]));
					}
					newTet = new Tet(a, b, c, d);
					tets.add(newTet);
					if (testing) newTet.toSceneEdges(scene, Color.blue, 0.001);
					if (testing) if (!newTet.getCircumSphere().containsNoneButAtMostOne(u, vertices)) System.out.println("Not Delaunay");

					// neighbor pointers are updated
					newTet.neighbors[newTet.indexOf(a)] = oppFace.oppTet;
					newTet.neighbors[newTet.indexOf(b)] = thirdFace.oppTet;
					newTet.neighbors[newTet.indexOf(c)] = f.oppTet;
					newTet.neighbors[newTet.indexOf(d)] = null;
					oppFace.oppTet.neighbors[oppFace.oppTet.indexOf(b, c, d)] = newTet;
					thirdFace.oppTet.neighbors[thirdFace.oppTet.indexOf(c, a, d)] = newTet;
					f.oppTet.neighbors[f.oppTet.indexOf(a, b, d)] = newTet;
					
					// hole is updated
					Face newFace = hole.new Face(a, b, c, d, null, newTet);
					hole.faces.add(newFace);
					if (testing) newTet.toSceneFace(scene, newTet.indexOf(a, b, c), Color.red);
					face = oppFace.neighbors[oppFace.indexOf(d)];
					newFace.neighbors[newFace.indexOf(a)] = face;
					face.neighbors[face.IndexOfFreeVertex(oppFace)] = newFace;
					
					// new face pairs are identified
					gamma = Point.orientation(newFace.vertices[0], newFace.vertices[1], newFace.vertices[2], face.getFreeVertex(newFace));
					if (gamma > 0.0) {
						double beta = Point.inSphere(newFace.vertices[0], newFace.vertices[1], newFace.vertices[2], face.getFreeVertex(newFace), u);
						heap.insert(new HeapItemDelete(face, oppFace, -beta/gamma));
					}
					face = thirdFace.neighbors[thirdFace.indexOf(d)];
					newFace.neighbors[newFace.indexOf(b)] = face;
					face.neighbors[face.IndexOfFreeVertex(thirdFace)] = newFace;
					gamma = Point.orientation(newFace.vertices[0], newFace.vertices[1], newFace.vertices[2], face.getFreeVertex(newFace));
					if (gamma > 0.0) {
						double beta = Point.inSphere(newFace.vertices[0], newFace.vertices[1], newFace.vertices[2], face.getFreeVertex(newFace), u);
						heap.insert(new HeapItemDelete(face, thirdFace, -beta/gamma));
					}
					face = f.neighbors[f.indexOf(d)];
					newFace.neighbors[newFace.indexOf(c)] = face;
					face.neighbors[face.IndexOfFreeVertex(f)] = newFace;
					gamma = Point.orientation(newFace.vertices[0], newFace.vertices[1], newFace.vertices[2], face.getFreeVertex(newFace));
					if (gamma > 0.0) {
						double beta = Point.inSphere(newFace.vertices[0], newFace.vertices[1], newFace.vertices[2], face.getFreeVertex(newFace), u);
						heap.insert(new HeapItemDelete(face, f, -beta/gamma));
					}
				}
				else {
					a = f.getFreeVertex(oppFace);
					d = oppFace.getFreeVertex(f);
					b = f.vertices[0];
					if (b == a) {
						b = f.vertices[1];
						c = f.vertices[2];
					}
					else {
						c = f.vertices[1];
						if (c == a) c = f.vertices[2];
					}
					hole.faces.remove(f);
					f.processed = true;
					hole.faces.remove(oppFace);
					oppFace.processed = true;
					if (testing) {
						f.oppTet.fromSceneFace(scene, f.oppTet.indexOf(f.vertices[0], f.vertices[1], f.vertices[2]));
						oppFace.oppTet.fromSceneFace(scene, oppFace.oppTet.indexOf(oppFace.vertices[0], oppFace.vertices[1], oppFace.vertices[2]));
					}
	
					newTet = new Tet(f.vertices[0], f.vertices[1], f.vertices[2], d);
					tets.add(newTet);
					if (testing) newTet.toSceneEdges(scene, Color.blue, 0.001);
					if (testing) if (!newTet.getCircumSphere().containsNoneButAtMostOne(u, vertices)) System.out.println("Not Delaunay");

					// neighbor pointers are updated
					newTet.neighbors[newTet.indexOf(a)] = oppFace.oppTet;
					newTet.neighbors[newTet.indexOf(d)] = f.oppTet;
					oppFace.oppTet.neighbors[oppFace.oppTet.indexOf(b, c, d)] = newTet;
					f.oppTet.neighbors[f.oppTet.indexOf(a, b, c)] = newTet;
					
					// hole is updated, first new face is added
					Face newFace1 = hole.new Face(a, b, d, c, null, newTet);
					hole.faces.add(newFace1);
					if (testing) newTet.toSceneFace(scene, newTet.indexOf(a, b, d), Color.red);
	
					// new face pairs are identified
					face = oppFace.neighbors[oppFace.indexOf(c)];
					newFace1.neighbors[newFace1.indexOf(a)] = face;
					face.neighbors[face.IndexOfFreeVertex(oppFace)] = newFace1;
					gamma = Point.orientation(newFace1.vertices[0], newFace1.vertices[1], newFace1.vertices[2], face.getFreeVertex(newFace1));
					if (gamma > 0.0) {
						double beta = Point.inSphere(newFace1.vertices[0], newFace1.vertices[1], newFace1.vertices[2], face.getFreeVertex(newFace1), u);
						heap.insert(new HeapItemDelete(face, newFace1, -beta/gamma));
					}
	
					face = f.neighbors[f.indexOf(c)];
					newFace1.neighbors[newFace1.indexOf(d)] = face;
					face.neighbors[face.IndexOfFreeVertex(f)] = newFace1;
					gamma = Point.orientation(newFace1.vertices[0], newFace1.vertices[1], newFace1.vertices[2], face.getFreeVertex(newFace1));
					if (gamma > 0.0) {
						double beta = Point.inSphere(newFace1.vertices[0], newFace1.vertices[1], newFace1.vertices[2], face.getFreeVertex(newFace1), u);
						heap.insert(new HeapItemDelete(face, newFace1, -beta/gamma));
					}
					
					// hole is updated, second new face is added
					Face newFace2 = hole.new Face(a, c, d, b, null, newTet);
					hole.faces.add(newFace2);
					if (testing) newTet.toSceneFace(scene, newTet.indexOf(a, c, d), Color.red);
					
					face = oppFace.neighbors[oppFace.indexOf(b)];
					newFace2.neighbors[newFace2.indexOf(a)] = face;
					face.neighbors[face.IndexOfFreeVertex(oppFace)] = newFace2;
					gamma = Point.orientation(newFace2.vertices[0], newFace2.vertices[1], newFace2.vertices[2], face.getFreeVertex(newFace2));
					if (gamma > 0.0) {
						double beta = Point.inSphere(newFace2.vertices[0], newFace2.vertices[1], newFace2.vertices[2], face.getFreeVertex(newFace2), u);
						heap.insert(new HeapItemDelete(face, newFace2, -beta/gamma));
					}
					
					face = f.neighbors[f.indexOf(b)];
					newFace2.neighbors[newFace2.indexOf(d)] = face;
					face.neighbors[face.IndexOfFreeVertex(f)] = newFace2;
					gamma = Point.orientation(newFace2.vertices[0], newFace2.vertices[1], newFace2.vertices[2], face.getFreeVertex(newFace2));
					if (gamma > 0.0) {
						double beta = Point.inSphere(newFace2.vertices[0], newFace2.vertices[1], newFace2.vertices[2], face.getFreeVertex(newFace2), u);
						heap.insert(new HeapItemDelete(face, newFace2, -beta/gamma));
					}
					
					newFace1.neighbors[newFace1.indexOf(b)] = newFace2;
					newFace2.neighbors[newFace2.indexOf(c)] = newFace1;
				}
			}
		}
		
		heap.clear();
		
		Face face1 = hole.faces.remove(0);
		Face face2 = hole.faces.remove(0);
		Face face3 = hole.faces.remove(0);
		Face face4 = hole.faces.remove(0);
		if (testing) {
			face1.oppTet.fromSceneFace(scene, face1.oppTet.indexOf(face1.vertices[0], face1.vertices[1], face1.vertices[2]));
			face2.oppTet.fromSceneFace(scene, face2.oppTet.indexOf(face2.vertices[0], face2.vertices[1], face2.vertices[2]));
			face3.oppTet.fromSceneFace(scene, face3.oppTet.indexOf(face3.vertices[0], face3.vertices[1], face3.vertices[2]));
			face4.oppTet.fromSceneFace(scene, face4.oppTet.indexOf(face4.vertices[0], face4.vertices[1], face4.vertices[2]));
		}		
		a = face1.vertices[0];
		b = face1.vertices[1];
		c = face1.vertices[2];
		d = face2.getFreeVertex(face1);
		newTet = new Tet(a, b, c, d);
		tets.add(newTet);
		if (testing) newTet.toSceneEdges(scene, Color.blue, 0.001);
		
		// looking for face that does not contain a
		face = face1;
		if (face.containsVertex(a)) {
			face = face2;
			if (face.containsVertex(a)) {
				face = face3;
				if (face.containsVertex(a)) face = face4;
			}
		}
		newTet.neighbors[newTet.indexOf(a)] = face.oppTet;
		face.oppTet.neighbors[face.oppTet.indexOf(b, c, d)] = newTet;

		// looking for face that does not contain b
		face = face1;
		if (face.containsVertex(b)) {
			face = face2;
			if (face.containsVertex(b)) {
				face = face3;
				if (face.containsVertex(b)) face = face4;
			}
		}
		newTet.neighbors[newTet.indexOf(b)] = face.oppTet;
		face.oppTet.neighbors[face.oppTet.indexOf(a, c, d)] = newTet;
		if (testing) {
			System.out.println("Neighbor of " + newTet.toString() + " facing " + b + " is " + face.oppTet.toString() );
			System.out.println("Neighbor of " + face.oppTet.toString() + " facing " + face.oppTet.getCorner(face.oppTet.indexOf(a,c,d)) + " is " + newTet.toString() );
		}
		
		// looking for face that does not contain c
		face = face1;
		if (face.containsVertex(c)) {
			face = face2;
			if (face.containsVertex(c)) {
				face = face3;
				if (face.containsVertex(c)) face = face4;
			}
		}
		newTet.neighbors[newTet.indexOf(c)] = face.oppTet;
		face.oppTet.neighbors[face.oppTet.indexOf(a, b, d)] = newTet;
		
		// looking for face that does not contain d
		face = face1;
		if (face.containsVertex(d)) {
			face = face2;
			if (face.containsVertex(d)) {
				face = face3;
				if (face.containsVertex(d)) face = face4;
			}
		}
		newTet.neighbors[newTet.indexOf(d)] = face.oppTet;
		face.oppTet.neighbors[face.oppTet.indexOf(a, b, c)] = newTet;
	}

	
	/** Returns the rotation direction (clockwise or counterclockwise */
	public Direction getDirection() { return rotDir; }

	/** Reads or creates a problem instance */
	private List<Point> getInstance(ProblemInstanceType instanceType) {
		if (instanceType == ProblemInstanceType.pdbNCC) {
			// reads a pdb-file - at the moment reads the backbone atoms only
			PDBFile f = new PDBFile("/Users/pawel/Downloads/1X0O.pdb");
			return f.getAtomCoords("N,Ca,C");	
		}
		if (instanceType == ProblemInstanceType.random) {
		// Generates a set of uniformly distributed points inside a unit cube
			Randomization.seed(3);
			return PointList.generatePointsInCube(10, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5);
		}
		if (instanceType == ProblemInstanceType.toy) {
		// Generates a toy example
			List<Point> points = new java.util.LinkedList<Point>();
			points.add(new Point(0,0,0));
			points.add(new Point(1,0,0));
			points.add(new Point(0,1,0));
			points.add(new Point(0,0,1));
			points.add(new Point(1.1,1.1,1.1));
			return points;
		}
		return null;
	}
	
	/** Returns the number of vertices in the DT (without 4 big points*/
	public int getNrVertices() { return vertices.size() - 4; }
	
	/** Used when deleting vertices from static DT */
	private double getPower(Point a, Point b, Point c, Point d, Point p) {
		double delta;
		double gamma = Point.orientation(a, b, c, p);
		if (gamma > 0.0) delta = Point.orientation(a, b, c, d); 
		else delta = Point.orientation(a, c, b, d);		
		if (delta <= 0.0) return Constants.bigDouble;
		return Point.inSphere(a, b, c, d, p)/delta;
	}

	/** Returns the direction vector of the rotation axis */
	public Vector getRotationAxis() { return rotationAxis; }

	/** Returns the rotation point */
	public Point getRotationPoint() { return rotationPoint; }
	
	
	/** Returns tetrahedra sharing a face and a vertex */
	public Tet[] getTetrahedra(Vertex v, Tet tet) {
		int indxV = tet.indexOf(v);
		Tet[] nTet = new Tet[3]; 
		for (int i = 0; i < 3; i++) {
			nTet[i] = tet.neighbors[(indxV+i+1)/3];
		}
		return nTet;
	}
	
	/** Returns a tetrahedron containing vertex v */
	public Tet getTetrahedron(Vertex v) {
		for (Tet tet : tets) if (tet.hasVertex(v)) return tet;
		return null;
	}
	
	/** Returns the list of tetrahedra */
	public List<Tet> getTets() { return tets; }

	/** Returns i-th vertex */
	public Vertex getVertex(int i) { return vertices.get(i); }
	
	/** Returns the list of vertices */
	public List<Vertex> getVertices() { return vertices; }
	
	/** Inserts new point into DT */
	public void insertPoint(Point p){
		Vertex v = new Vertex(p);
		vertices.add(v);
		Tet c = walk(p);
		
		//Corresponds to findNTes
		List<Tet> newTets = new LinkedList<Tet>();
		HashSet<Tet> processed = new HashSet<Tet>();
		processed.add(null);
		Stack<Tet> fringe = new Stack<Tet>();
		fringe.add(c);

		while(!fringe.isEmpty()){
			c = fringe.pop();
			if(processed.contains(c)) continue;
			for(int f=0;f<4;f++){
				Tet neigh = c.neighbors[f];
				if(neigh==null || !neigh.insideCircumsphere(p)){
					//-3 -2  0  1 .. c
					//-2 -1  0  1 .. neigh .. neigh.apex(c) should be 0
					//Create new cell
					Vertex[] corners = new Vertex[4];
					corners[3] = v;
					for(int i=1;i<4;i++) corners[i-1] = c.getCorner((f+i)%4);
					Tet newTet = new Tet(corners);
					newTet.neighbors[3] = neigh;
					if(neigh!=null) neigh.neighbors[neigh.apex(c)] = newTet;
					newTets.add(newTet);
				}else if(!processed.contains(neigh)){
					fringe.add(neigh);
				}
			}
			processed.add(c);
			tets.remove(c);
		}
		
		lastTet = newTets.get(0); 
		restoreNeighborhood(newTets);
		tets.addAll(newTets);
	}

	
	/** Removes a tetrahedron from the list of tetrahedra */
	public void removeTetrahedron(Tet tet) { tets.remove(tet); }

	/** Restore neighborhood information after insertion and flips */
	private void restoreNeighborhood(List<Tet> newCells){
		for(Tet c1: newCells){
			cellLoop: for(Tet c2: newCells){
				if(c1==c2) break;
				
				//Check if c1 and c2 share a face based on corner vertices
				//Both will contain the last inserted vertex (highest index)
				//i and j run over the arrays, I and J record the location of differences
				int i=0,j=0, I=-1, J=-1;
				while(i<=3&&j<=3){
					if(c1.getCorner(i)==c2.getCorner(j)) {
						i++;j++;
					}
					else 
						if(i<3 && c1.getCorner(i+1)==c2.getCorner(j)) {
							if(I>=0) continue cellLoop;
							I=i;i++;
						}
						else {
							if(J>=0) continue cellLoop;
							J=j;j++;
						}
				}
				c1.neighbors[I] = c2;
				c2.neighbors[J] = c1;
			}
		}
	}

	/** Sets the direction of the rotation (clockwise or counterclockwise) */
	public void setDirection(Direction rotDir) { this.rotDir = rotDir; }

	/** Sets the direction vector of the rotation axis */
	public void setRotationAxis(Vector v)  { rotationAxis = v; }

	/** Sets the rotation point and translates all vertices so that the rotation point is in the origo.  */
	public void setRotationPoint(Point rotationPoint) { 
		if (this.rotationPoint != null) for (Vertex v : vertices) v.addThis(this.rotationPoint);
		this.rotationPoint = rotationPoint; 
		for (Vertex v : vertices) v.subtractThis(rotationPoint);	
	}

	// Identifies tetrahedron into which new point has to be inserted */
	private Tet walk(Point p){
		Tet t = lastTet;
		mainWalk: while(true){
			for(int f=0;f<4;f++){
				if(!t.insideFace(f, p)){
					t = t.neighbors[f];
					lastTet = t;
					continue mainWalk;
				}
			}
			return t;
		}
	}
	

	public void animateSpheres() {
		Sphere d;
		for (Tet t : tets) {
			if (!t.isFlat()) {
				Sphere c = t.getCircumSphere();
				if (c != null) {
					d = new Sphere(t.getCorner(0), t.getCorner(1), t.getCorner(2), t.getCorner(3));
					c.setCenter(d.getCenter());
					c.setRadius(d.getRadius()); 
				}
			}
		}
		scene.repaint();
	}

	
	private Double[] getRootSSSSR(Vertex A, Vertex B, Vertex C, Vertex D, Vertex E, int dir) {
		double aa = A.getSquaredPolarRadius(); 
		double bb = B.getSquaredPolarRadius(); 
		double cc = C.getSquaredPolarRadius();
		double dd = D.getSquaredPolarRadius(); 
		double ee = E.getSquaredPolarRadius();
	    double e  = Math.sqrt(E.x()*E.x() + E.y()*E.y());

	    double m[][] = { {A.x(), A.y(), A.z(), aa, 1.0}, 
	    				 {B.x(), B.y(), B.z(), bb, 1.0},
	    			     {C.x(), C.y(), C.z(), cc, 1.0},
	    				 {D.x(), D.y(), D.z(), dd, 1.0},
	    			     {0, 0, 0, 0, 0}};
	    Matrix M = new Matrix(m);
	    double det40 = M.minor(4, 0).determinant();
	    double det41 = M.minor(4, 1).determinant();
	    double det42 = M.minor(4, 2).determinant();
	    double det43 = M.minor(4, 3).determinant();
	    double det44 = M.minor(4, 4).determinant();

	    double coefSin = -e*(E.getSinAngle()*det40 + E.getCosAngle()*det41);
	    double coefCos = e*(E.getCosAngle()*det40 - E.getSinAngle()*det41);
	    double coef    = E.z()*det42 - ee*det43 + det44; 
	    angles = Trigonometry.solveAsinXPlusBcosXplusC(coefSin, coefCos, coef);
	    return getRotAngle(angles, dir);
	}

	private Double[] getRootSSSRR(Vertex A, Vertex B, Vertex C, Vertex D, Vertex E, int dir) {
		double aa = A.getSquaredPolarRadius(); 
		double bb = B.getSquaredPolarRadius(); 
		double cc = C.getSquaredPolarRadius();
		double dd = D.getSquaredPolarRadius(); 
		double d  = Math.sqrt(D.x()*D.x() + D.y()*D.y());
		double ee = E.getSquaredPolarRadius();
	    double e  = Math.sqrt(E.x()*E.x() + E.y()*E.y());
	    
	    double dSin = D.getSinAngle();
	    double dCos = D.getCosAngle();
	    double eSin = E.getSinAngle();
	    double eCos = E.getCosAngle();
	    
	    double aa_bb = aa - bb;
	    double bb_cc = bb - cc;
	    double cc_aa = cc - aa;
	    double M12 = A.z()*bb_cc                     + C.z()*aa_bb                     + B.z()*cc_aa;
	    double M13 = A.y()*bb_cc                     + C.y()*aa_bb                     + B.y()*cc_aa;
	    double M14 = A.y()*(B.z()-C.z())             + C.y()*(A.z()-B.z())             + B.y()*(C.z()-A.z());
	    double M15 = A.y()*(B.z()*cc-C.z()*bb)       + C.y()*(A.z()*bb-B.z()*aa)       + B.y()*(C.z()*aa-A.z()*cc);
	    double M23 = A.x()*bb_cc                     + C.x()*aa_bb                     + B.x()*cc_aa;
	    double M24 = A.x()*(B.z()-C.z())             + C.x()*(A.z()-B.z())             + B.x()*(C.z()-A.z());
	    double M25 = A.x()*(B.z()*cc-C.z()*bb)       + C.x()*(A.z()*bb-B.z()*aa)       + B.x()*(C.z()*aa-A.z()*cc);																																		;
	    double M34 = A.x()*(B.y()-C.y())             + C.x()*(A.y()-B.y())             + B.x()*(C.y()-A.y());
	    double M35 = A.x()*(B.y()*cc-C.y()*bb)       + C.x()*(A.y()*bb-B.y()*aa)       + B.x()*(C.y()*aa-A.y()*cc);
	    double M45 = A.x()*(B.y()*C.z()-C.y()*B.z()) + C.x()*(A.y()*B.z()-B.y()*A.z()) + B.x()*(C.y()*A.z()-A.y()*C.z());

	    double coefSin = M13*(E.z()*d*dSin - e*D.z()*eSin) + M14*(e*dd*eSin - ee*d*dSin) + M15*(d*dSin - e*eSin) + M23*(E.z()*d*dCos - e*D.z()*eCos) + M24*(e*dd*eCos - ee*d*dCos) + M25*(d*dCos - e*eCos);
	    double coefCos = M13*(e*D.z()*eCos - E.z()*d*dCos) + M14*(ee*d*dCos - e*dd*eCos) + M15*(e*eCos - d*dCos) + M23*(E.z()*d*dSin - e*D.z()*eSin) + M24*(e*dd*eSin - ee*d*dSin) + M25*(d*dSin - e*eSin);
	    double coef    = d*e*M12*(eSin*dCos - eCos*dSin) + M34*(ee*D.z() -E.z()*dd) + M35*(E.z() - D.z()) + M45*(dd - ee);
	    angles = Trigonometry.solveAsinXPlusBcosXplusC(coefSin, coefCos, coef);
	    return getRotAngle(angles, dir);
	}

	
	public Double[] getRotAngle(Double[] angles, int dir) {
		if (angles == null || angles[1] == null) return null;
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

	public Double[] getRoot(Tet t, Tet oppT) {
		int count = t.getCount();
		Vertex oppV = oppT.getCorner(oppT.apex(t));
		if (oppV.getType() == Vertex.VertexType.R) count = count + 16;
		return getRoot(t, oppV, count);
	}
	
	public Double[] getRoot(Tet t, Vertex oppV, int count) {		
		if (count == 0 || count == 31) return null;
		if (count <= 15) {	
			if (count <= 7) {
				if (count <= 3) {
					if (count == 1) return getRootSSSSR(t.getCorner(1), t.getCorner(2), t.getCorner(3), oppV, t.getCorner(0), 0);
					if (count == 2) return getRootSSSSR(t.getCorner(0), t.getCorner(2), t.getCorner(3), oppV, t.getCorner(1), 0);
					return getRootSSSRR(t.getCorner(2), t.getCorner(3), oppV, t.getCorner(0), t.getCorner(1), 0);
				}
				if (count <= 5) {
					if (count == 4) return getRootSSSSR(t.getCorner(0), t.getCorner(1), t.getCorner(3), oppV, t.getCorner(2), 0);
					return getRootSSSRR(t.getCorner(1), t.getCorner(3), oppV, t.getCorner(0), t.getCorner(2), 0);
				}
				if (count == 6) return getRootSSSRR(t.getCorner(0), t.getCorner(3), oppV, t.getCorner(1), t.getCorner(2), 0);					
				return getRootSSSRR(t.getCorner(0), t.getCorner(1), t.getCorner(2), t.getCorner(3), oppV, 1);
			}
			if (count <= 11) {
				if (count <= 9) {
					if (count == 8) return getRootSSSSR(t.getCorner(0), t.getCorner(1), t.getCorner(2), oppV, t.getCorner(3), 0);
					return getRootSSSRR(t.getCorner(1), t.getCorner(2), oppV, t.getCorner(0), t.getCorner(3), 0);
				}
				if (count == 10) return getRootSSSRR(t.getCorner(0), t.getCorner(2), oppV, t.getCorner(1), t.getCorner(3), 0);
				return getRootSSSRR(t.getCorner(0), t.getCorner(1), t.getCorner(3), t.getCorner(2), oppV, 1);
			}
			if (count <= 13) {
				if (count == 12) return getRootSSSRR(t.getCorner(0), t.getCorner(1), oppV, t.getCorner(2), t.getCorner(3), 0);
				return getRootSSSRR(t.getCorner(0), t.getCorner(2), t.getCorner(3), t.getCorner(1), oppV, 1);
			}
			if (count == 14) return getRootSSSRR(t.getCorner(1), t.getCorner(2), t.getCorner(3), t.getCorner(0), oppV, 1);
			return getRootSSSSR(t.getCorner(0), t.getCorner(1), t.getCorner(2), t.getCorner(3), oppV, 1);
		}
		if (count <= 23) {
			if (count <= 19) {
				if (count <= 17) {
					if (count == 16) return getRootSSSSR(t.getCorner(0), t.getCorner(1), t.getCorner(2), t.getCorner(3), oppV, 0);
					return getRootSSSRR(t.getCorner(1), t.getCorner(2), t.getCorner(3), t.getCorner(0), oppV, 0);
				}
				if (count == 18) return getRootSSSRR(t.getCorner(0), t.getCorner(2), t.getCorner(3), t.getCorner(1), oppV, 0);
				return getRootSSSRR(t.getCorner(0), t.getCorner(1), oppV, t.getCorner(2), t.getCorner(3), 1);
			}
			if (count <= 21) {
				if (count == 20) return getRootSSSRR(t.getCorner(0), t.getCorner(1), t.getCorner(3), t.getCorner(2), oppV, 0);
				return getRootSSSRR(t.getCorner(1), t.getCorner(3), oppV, t.getCorner(0), t.getCorner(2), 0);
			}
			if (count == 22) return getRootSSSRR(t.getCorner(1), t.getCorner(2), oppV, t.getCorner(0), t.getCorner(3), 1);					
			return getRootSSSSR(t.getCorner(0), t.getCorner(1), t.getCorner(2), oppV, t.getCorner(3), 1);
		}
		if (count <= 27) {
			if (count <= 25) {
				if (count == 24) return getRootSSSRR(t.getCorner(0), t.getCorner(1), t.getCorner(2), t.getCorner(3), oppV, 0);
				return getRootSSSRR(t.getCorner(0), t.getCorner(3), oppV, t.getCorner(1), t.getCorner(2), 1);
			}
			if (count == 26) return getRootSSSRR(t.getCorner(1), t.getCorner(3), oppV, t.getCorner(0), t.getCorner(2), 1);
			return getRootSSSSR(t.getCorner(0), t.getCorner(1), t.getCorner(3), oppV, t.getCorner(2), 1);
		}
		if (count <= 29) {
			if (count == 28) return getRootSSSRR(t.getCorner(2), t.getCorner(3), oppV, t.getCorner(0), t.getCorner(1), 1);
			return getRootSSSSR(t.getCorner(0), t.getCorner(2), t.getCorner(3), oppV, t.getCorner(1), 1);
		}
		return getRootSSSSR(t.getCorner(1), t.getCorner(2), t.getCorner(3), oppV, t.getCorner(0), 1);
	}
	
	private void testInfo(Tet t, Color clr, boolean second) {
		if (second) System.out.print(" ");
		System.out.println(t.toString()); 
		t.circumSphere = new Sphere(t.getCorner(0), t.getCorner(1), t.getCorner(2), t.getCorner(3));
		scene.addShape(t.circumSphere, new Color(255, 255, 0, 100), 32);
	}

	
	private void initializeRotation() {
		heap.clear();
		// set constants associated with vertices - does not include vertices of big points
		for (Vertex v : vertices) {
			if (testingPrint) System.out.print(v.getId() + ": " + v.toString(2));
			if (v.distanceSquared() < Constants.EPSILON) { //takes care of the special case when the rotation center overlaps with one of the given points
				v.setSquaredPolarRadius(0.0);                            
				v.setPolarRadius(0.0);                        
				v.setPolarAngle(0.0);			
				if (testingPrint) System.out.println(", polar angle: 0.0"); 
				v.setCosAngle(1.0);
				v.setSinAngle(0.0);
				v.setType(Vertex.VertexType.S);
				rotIndx.remove((Integer)v.getId());
			}
			else {
				v.setSquaredPolarRadius(v.distanceSquared());       // remains unchanged when rotating around the same point
				v.setPolarRadius(v.distance());                     // remains unchanged when rotating around the same point
				v.setPolarAngle(v.polarAngleXY());     			
				if (testingPrint) System.out.println(", initial polar angle: " + Functions.toDeg(v.getPolarAngle())); 
				v.setCosAngle(v.polarAngleCosXY());
				v.setSinAngle(v.polarAngleSinXY());
			}
		}
		int count = 0;
		int count4 = 0;
		Tet nt = null;
		Vertex oppV;
		for (Tet t : tets) {
			count4 = t.getCount();
			for (int i = 0; i < 4; i++) {
				nt = t.neighbors[i];
				if (nt != null) { 
					oppV = nt.getCorner(nt.apex(t));  
					if (t.getCorner(t.apex(nt)).getId() < oppV.getId()) {
						count = count4;
						if (oppV.getType() == Vertex.VertexType.R) count = count + 16;	
						angles = getRoot(t, oppV, count);
						if ((angles != null) && (angles[0] < angleLimit)) {
							addToHeap(angles, t, nt);
							if (testingPrint) System.out.println(t + " " + nt + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
						}
					}
				}
			}
		}
	}

	
	/**
	 * Performs a 2-3-flip of the two specified tetrahedra. The two existing tetrahedra are preserved and 
	 * a new is created and returned. 
	 * @requires 	Arrays.equals(t0.corners, Arrays.sort(t0.corners) ) &&
	 * 				Arrays.equals(t1.corners, Arrays.sort(t1.corners) ) &&
	 * 				convex(t0,t1)
	 * @return The newly created tetrahedron.
	 */
	private Tet[] flip23(Tet t0, Tet t1){

		Vertex vns0 = t0.getCorner(t0.apex(t1));
		Vertex vns1 = t1.getCorner(t1.apex(t0));
		
		// identify three shared vertices
		Vertex[] vs = new Vertex[3];
		int k = 0;
		for (int i = 0; i < 4; i++) if (t0.getCorner(i) != vns0) vs[k++] = t0.getCorner(i);
		
		Tet nt0 = new Tet(vns0, vns1, vs[0], vs[1]);
		Tet nt1 = new Tet(vns0, vns1, vs[1], vs[2]);
		Tet nt2 = new Tet(vns0, vns1, vs[2], vs[0]);
		
		// vertices get new tetrahedron pointers
		vns0.setTet(nt0);
		vns1.setTet(nt0);
		vs[0].setTet(nt0);
		vs[1].setTet(nt0);
		vs[2].setTet(nt1);
		
		// nt0 contain vs[0]Êand vs[1]
		
		if (!nt1.hasVertex(vs[0])) nt0.neighbors[nt0.indexOf(vs[0])] = nt1; else nt0.neighbors[nt0.indexOf(vs[0])] = nt2;
		if (!nt1.hasVertex(vs[1])) nt0.neighbors[nt0.indexOf(vs[1])] = nt1; else nt0.neighbors[nt0.indexOf(vs[1])] = nt2;
		Tet t = t1.neighbors[t1.indexOf(vs[2])];
		nt0.neighbors[nt0.indexOf(vns0)] = t;
		if (t != null) t.neighbors[t.apex(t1)] = nt0;
		t = t0.neighbors[t0.indexOf(vs[2])];
		nt0.neighbors[nt0.indexOf(vns1)] = t;
		if (t != null) t.neighbors[t.apex(t0)] = nt0;
		
		// nt1 contains vs[1] and vs[2]
		
		if (!nt0.hasVertex(vs[1])) nt1.neighbors[nt1.indexOf(vs[1])] = nt0; else nt1.neighbors[nt1.indexOf(vs[1])] = nt2;
		if (!nt0.hasVertex(vs[2])) nt1.neighbors[nt1.indexOf(vs[2])] = nt0; else nt1.neighbors[nt1.indexOf(vs[2])] = nt2;
		t = t1.neighbors[t1.indexOf(vs[0])];
		nt1.neighbors[nt1.indexOf(vns0)] = t;
		if (t != null) t.neighbors[t.apex(t1)] = nt1;
		t = t0.neighbors[t0.indexOf(vs[0])];
		nt1.neighbors[nt1.indexOf(vns1)] = t;
		if (t != null) t.neighbors[t.apex(t0)] = nt1;
		
		// nt2 contains vs[2] and vs[0]
		
		if (!nt0.hasVertex(vs[2])) nt2.neighbors[nt2.indexOf(vs[2])] = nt0; else nt2.neighbors[nt2.indexOf(vs[2])] = nt1;
		if (!nt0.hasVertex(vs[0])) nt2.neighbors[nt2.indexOf(vs[0])] = nt0; else nt2.neighbors[nt2.indexOf(vs[0])] = nt1;
		t = t1.neighbors[t1.indexOf(vs[1])];
		nt2.neighbors[nt2.indexOf(vns0)] = t;
		if (t != null) t.neighbors[t.apex(t1)] = nt2;
		t = t0.neighbors[t0.indexOf(vs[1])];
		nt2.neighbors[nt2.indexOf(vns1)] = t;
		if (t != null) t.neighbors[t.apex(t0)] = nt2;
		
		tets.remove(t0);
		t0.setAlive(false);
		if (testingScreen) t0.fromSceneEdges(scene);
		tets.remove(t1);
		t1.setAlive(false);		
		if (testingScreen) t1.fromSceneEdges(scene);
		Tet[] newTets = new Tet[3];
		newTets[0] = nt0;
		newTets[1] = nt1;
		newTets[2] = nt2;
		if (testingScreen) 
			for (int i = 0; i < 3; i++) newTets[i].toSceneEdges(scene, Color.black, 0.001, 0.0001);
		tets.add(nt0);
		tets.add(nt1);
		tets.add(nt2);
		angles = getRoot(nt0, nt1);
		if (angles[0] < angleLimit) {
			addToHeap(angles, nt0, nt1);
			if (testingPrint) System.out.println(nt0 + " " + nt1 + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
		}
		angles = getRoot(nt0, nt2);
		if (angles[0] < angleLimit) {
			addToHeap(angles, nt0, nt2);
			if (testingPrint) System.out.println(nt0 + " " + nt2 + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
		}
		angles = getRoot(nt1, nt2);
		if (angles[0] < angleLimit) {
			addToHeap(angles, nt1, nt2);
			if (testingPrint) System.out.println(nt1 + " " + nt2 + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
		}
		
		// check for new events
		Tet nti;
		for (int i = 0; i < 4; i++) {
			nti = nt0.neighbors[i];
			if ((nti != null) && (nti != nt1) && (nti != nt2)) {
				angles = getRoot(nt0, nti);
				if ((angles != null) && (angles[0] < angleLimit)) {
					addToHeap(angles, nt0, nti);
					if (testingPrint) System.out.println(nt0 + " " + nti + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
				}
			}
		}
		
		for (int i = 0; i < 4; i++) {
			nti = nt1.neighbors[i];
			if ((nti != null) && (nti != nt0) && (nti != nt2)) {
				angles = getRoot(nt1, nti);
				if ((angles != null) && (angles[0] < angleLimit)) {
					addToHeap(angles, nt1, nti);
					if (testingPrint) System.out.println(nt1 + " " + nti + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
				}
			}
		}

		for (int i = 0; i < 4; i++) {
			nti = nt2.neighbors[i];
			if (nti != null && nti != nt0 && nti != nt1) {
				angles = getRoot(nt2, nti);
				if (angles != null && angles[0] < angleLimit) {
					addToHeap(angles, nt2, nti);
					if (testingPrint) System.out.println(nt2 + " " + nti + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
				}
			}
		}
		return newTets;
	}

	
	/**
	 * Performs a 3-2-flip of the three specified tetrahedra. Two of these tetrahedra are preserved 
	 * and the third (deleted) is returned. 
	 * @requires 	Arrays.equals(t0.corners, Arrays.sort(t0.corners) ) &&
	 * 				Arrays.equals(t1.corners, Arrays.sort(t1.corners) ) &&
	 * 				Arrays.equals(t1.corners, Arrays.sort(t1.corners) ) &&
	 * 				_two vertices shared by all three tetrahedra_ 
	 * @return The deleted tetrahedron.
	 */
	private Tet[] flip32(Tet t0, Tet t1, Tet t2){

		//Locate three non-shared vertices
		Vertex[] vns = new Vertex[3];
		vns[0] = t0.getCorner(t0.apex(t1));
		vns[1] = t1.getCorner(t1.apex(t2));
		vns[2] = t2.getCorner(t2.apex(t0));

		//Locate two shared vertices
		Vertex[] vs = new Vertex[2];
		int nrSharedFound = 0;
		for (int i = 0; (i < 4) && (nrSharedFound < 2); i++) {
			Vertex v = t0.getCorner(i);
			if (t1.hasVertex(v) && t2.hasVertex(v)) vs[nrSharedFound++] = v;
		}
		
		Tet nt0 = new Tet(vns[0], vns[1], vns[2], vs[0]);
		Tet nt1 = new Tet(vns[0], vns[1], vns[2], vs[1]);

		// vertices get new tetrahedron pointers
		vs[0].setTet(nt0);
		vs[1].setTet(nt1);
		vns[0].setTet(nt0);
		vns[1].setTet(nt0);
		vns[2].setTet(nt0);
		

		
		// update neighbor information
		
		nt0.neighbors[nt0.indexOf_slow(vs[0])] = nt1;      
		nt1.neighbors[nt1.indexOf_slow(vs[1])] = nt0;

		Tet tn0 = t2;
		if (!t0.hasVertex(vns[0])) tn0 = t0; else { if (!t1.hasVertex(vns[0])) tn0 = t1; }
		Tet tn1 = t2;
		if (!t0.hasVertex(vns[1])) tn1 = t0; else { if (!t1.hasVertex(vns[1])) tn1 = t1; }
		Tet tn2 = t2;
		if (!t0.hasVertex(vns[2])) tn2 = t0; else { if (!t1.hasVertex(vns[2])) tn2 = t1; }
		
		for (int i = 0; i < 4; i++) {
			Tet n0 = tn0.neighbors[i];
			if ((n0 != null) && (n0 != tn1) && (n0 != tn2)) {
				if (n0.hasVertex(vs[0])) {
					n0.neighbors[n0.apex(tn0)] = nt0;
					nt0.neighbors[nt0.indexOf_slow(vns[0])] = n0;
				}
				else {
					n0.neighbors[n0.apex(tn0)] = nt1;
					nt1.neighbors[nt1.indexOf_slow(vns[0])] = n0;
				}
			}
		}

		for (int i = 0; i < 4; i++) {
			Tet n1 = tn1.neighbors[i];
			if ((n1 != null) && (n1 != tn0) && (n1 != tn2)) {
				if (n1.hasVertex(vs[0])) {
					n1.neighbors[n1.apex(tn1)] = nt0;
					nt0.neighbors[nt0.indexOf_slow(vns[1])] = n1;
				}
				else {
					n1.neighbors[n1.apex(tn1)] = nt1;
					nt1.neighbors[nt1.indexOf_slow(vns[1])] = n1;
				}
			}
		}

		for (int i = 0; i < 4; i++) {
			Tet n2 = tn2.neighbors[i];
			if ((n2 != null) && (n2 != tn0) && (n2 != tn1)) {
				if (n2.hasVertex(vs[0])) {
					n2.neighbors[n2.apex(tn2)] = nt0;
					nt0.neighbors[nt0.indexOf_slow(vns[2])] = n2;
				}
				else {
					n2.neighbors[n2.apex(tn2)] = nt1;
					nt1.neighbors[nt1.indexOf_slow(vns[2])] = n2;
				}
			}
		}
		
		tets.remove(t0);
		t0.setAlive(false);
		tets.remove(t1);
		t1.setAlive(false);		
		tets.remove(t2);
		t2.setAlive(false);
		if (testingScreen) {
			t0.fromSceneEdges(scene);
			t1.fromSceneEdges(scene);
			t2.fromSceneEdges(scene);
		}
		Tet[] newTets = new Tet[2];
		newTets[0] = nt0;
		newTets[1] = nt1;
		if (testingScreen) {
			newTets[0].toSceneEdges(scene, Color.black, 0.001, 0.0001);
			newTets[1].toSceneEdges(scene, Color.black, 0.001, 0.0001);
		}
		tets.add(nt0);
		tets.add(nt1);
		angles = getRoot(nt0, nt1);
		if (angles[0] < angleLimit) {
			addToHeap(angles, nt0, nt1);
			if (testingPrint) System.out.println(nt0 + " " + nt1 + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
		}

		// check for new events
		
		Tet nti;
		for (int i = 0; i < 4; i++) {
			nti = nt0.neighbors[i];
			if ((nti != null) && (nti != nt1)) {
				angles = getRoot(nt0, nti);
				if ((angles != null) && (angles[0] < angleLimit)) {
					addToHeap(angles, nt0, nti);
					if (testingPrint) System.out.println(nt0 + " " + nti + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
				}
			}
		}

		for (int i = 0; i < 4; i++) {
			nti = nt1.neighbors[i];
			if ((nti != null) && (nti != nt0)) {
				angles = getRoot(nt1, nti);
				if ((angles != null) && (angles[0] < angleLimit)) {
					addToHeap(angles, nt1, nti);
					if (testingPrint) System.out.println(nt1 + " " + nti + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
				}
			}
		}	
		return newTets;
	}
	
	/** Returns TRUE id the tetrahedron is Delaunay (circumscribing sphere contains no point) */
	public boolean isDelaunay() {
		boolean cont = true;
		for (Tet t : tets) {
			t.setCircumSphere();
			if (!t.circumSphere.isEmpty(vertices, Constants.EPSILON)) {
				System.out.print(t + " is not empty: ");
				t.circumSphere.contains(vertices, Constants.EPSILON);
				cont = false;
			}
		}
		return cont;
	}
	
	public void showPair(Tet t, Tet nt, Tet tt) {
		commonTriangle = t.commonFace(nt);
		commonTriangleShape = commonTriangle.toScene(scene, Color.blue);
		if (tt != null) {
			commonTriangle02 = t.commonFace(tt);
			commonTriangleShape02 =commonTriangle02.toScene(scene, Color.red);
			commonTriangle12 = nt.commonFace(tt);
			commonTriangleShape12 =commonTriangle12.toScene(scene, Color.green);
		}
		apexSegment = new LineSegment(t.getCorner(t.apex(nt)), nt.getCorner(nt.apex(t)));
		apexLSS = new LSS(apexSegment, 0.02);
		scene.addShape(apexLSS, Color.magenta);
		scene.repaint();
	}
	
	public void hidePair(Tet tt) {
		scene.removeShape(apexLSS);
		scene.removeShape(commonTriangleShape);
		if (tt != null) {
			scene.removeShape(commonTriangleShape02);
			scene.removeShape(commonTriangleShape12);
		}
		scene.repaint();
	}
	
	public void rotate() {
		HeapItem heapItem;
		Tet t, nt, tt;
		Tet tet1 = null; 
		Tet tet2 = null;
		Tet[] newTets = new Tet[3];
		int nrFlips = 0;
		double rotAngle;
//		Color red_tr = new Color(255,0,0,50);
//		Color blue_tr = new Color(0,0,255,50);
//		Color green_tr = new Color(0,128,0,50);
		initializeRotation();
		while (!heap.isEmpty()) {
			heapItem = (HeapItem) heap.extract();
			t = heapItem.getT();   
			nt = heapItem.getNT(); 
			if (testingPrint) {
				System.out.print(t.toString());  if (!t.isAlive())  System.out.print(" ");
				System.out.print(nt.toString()); if (!nt.isAlive()) System.out.print(" ");
			}
			angles = heapItem.getAngles();
			rotAngle = angles[0];
			if (t.isAlive() && nt.isAlive() && (rotAngle < Constants.TAU)) {
				animate(scene, rotAngle-angleTotal);
//				isDelaunay();
//				animate(scene, (rotAngle-angleTotal)/2.0);	
				angleTotal = rotAngle;
				
				if (!t.isConvex(nt)) {
					tt = KineticToolbox.getThirdTet(t, nt);
					if (tt != null) {
						if (testingPrint) {
							System.out.print(tt.toString()); 
							if (!tt.isAlive())  System.out.print(" ");
							System.out.print(" rotated to angle = " + Functions.toDeg(heapItem.getAngles()[0]));						
							System.out.println();
						}
						newTets = flip32(t, nt, tt);
						if (testingPrint) System.out.println(++nrFlips + ". flip.");
					}
					else { System.out.println(" not convex but tt = null"); break; }
				}
				else {
					if (testingPrint) {
						System.out.print(" rotated to angle = " + Functions.toDeg(heapItem.getAngles()[0]));
						System.out.println();
					}
					if (testingScreen) {
						t.fromSceneFaces(scene);
						nt.fromSceneFaces(scene);
					}
					newTets = flip23(t, nt);
					if (testingPrint) System.out.println(++nrFlips + ". flip.");
				}
			}
			else System.out.println();
		}
		animate(scene, Constants.TAU-angleTotal);

	}
	/* specify what is rotating - 3 methods to do it */
	public void setRotVertices() { for (Vertex v : vertices) v.setType(Vertex.VertexType.S); }
	public void setRotVertices(int a) { setRotVertices(a, a); }
	public void setRotVertices(int a, int b) { 
		for (int i = a; i <= b; i++) rotIndx.add(i); 
		for (Vertex v : vertices) v.setType(Vertex.VertexType.S); 
		for (int i : rotIndx) vertices.get(i).setType(Vertex.VertexType.R);
	}
	public void setRotVertices(List<Integer> rotList) { 
		rotIndx = rotList; 
		for (Vertex v : vertices) v.setType(Vertex.VertexType.S); 
		for (int i : rotIndx) vertices.get(i).setType(Vertex.VertexType.R);
	}

	/** returns the scene of this kDT */
	public J3DScene getScene() { return scene; }
	
	/** returns the list of tetrahedra in this kinetic Dalaunay tessellation */
	public List<Tet> getTetrahedra(){ return tets; }
	
	/** Displays Delaunay tesselation */
	public void toScene(J3DScene scene) { toScene(scene, Constants.bigDouble); }
	
	// Displays the backbone of the protein
	public void toSceneBackbone(J3DScene scene, double width, Color clr) {
		int sz = vertices.size();
		for (int i = 4; i < sz-2; i++)
			new LineSegment(vertices.get(i), vertices.get(i+1)).toScene(scene, width, clr);
	}
	
	// Displays triangles with exactly one tetrahedron in the alpha complex 
	public void toSceneSurface(J3DScene scene) {
		Tet nTet;
		for (Tet tet : tets) {
			tet.getCircumSphere();
			if (tet.circumRadius() <= getAlpha()) {
				for (int i = 0; i < 4; i++) {
					nTet = tet.neighbors[i];
					nTet.getCircumSphere();
					if (nTet.circumRadius() > getAlpha()) tet.toSceneFace(scene, i, Color.red);
				}	
			}
		}
	}
	
	/** Displays alpha complex for given alpha */
	public void toScene(J3DScene scene, double alpha) {
		scene.removeAllShapes();
		Color blue_tr = new Color(0,0,255,50);
		Color green_tr = Color.green;
//		Color green_tr = new Color(0,128,0,50);
		Color red_tr = new Color(255,0,0,50);
		Vector tr = new Vector (0.02, 0.02, 0.02);
		int nrTetrahedra = 0;
		for (Tet tet : tets) { 
			nrTetrahedra++;
			System.out.println(nrTetrahedra);
			tet.getCircumSphere();
			if (tet.circumRadius() <= alpha) tet.toSceneFaces(scene, Color.red);
			else {
/*				Circle c;
				c = new Circle(tet.corners[0], tet.corners[1], tet.corners[2]);
				if ((c.getRadius() <= alpha) && (tet.neighbors[3].getCircumSphereRadius() > alpha)) tet.toSceneFace(scene, 3, green_tr);
				c = new Circle(tet.corners[0], tet.corners[1], tet.corners[3]);
				if ((c.getRadius() <= alpha) && (tet.neighbors[2].getCircumSphereRadius() > alpha)) tet.toSceneFace(scene, 2, green_tr);
				c = new Circle(tet.corners[0], tet.corners[2], tet.corners[3]);
				if ((c.getRadius() <= alpha) && (tet.neighbors[1].getCircumSphereRadius() > alpha)) tet.toSceneFace(scene, 1, green_tr);
				c = new Circle(tet.corners[1], tet.corners[2], tet.corners[3]);
				if ((c.getRadius() <= alpha) && (tet.neighbors[3].getCircumSphereRadius() > alpha)) tet.toSceneFace(scene, 0, green_tr);
				double alpha2 = 2*alpha;
				if (tet.corners[0].distance(tet.corners[1]) <= alpha2) new LineSegment(tet.corners[0], tet.corners[1]).toScene(scene, 0.003, Color.black, 1);
				if (tet.corners[0].distance(tet.corners[2]) <= alpha2) new LineSegment(tet.corners[0], tet.corners[2]).toScene(scene, 0.003, Color.black, 1);
				if (tet.corners[0].distance(tet.corners[3]) <= alpha2) new LineSegment(tet.corners[0], tet.corners[3]).toScene(scene, 0.003, Color.black, 1);
				if (tet.corners[1].distance(tet.corners[2]) <= alpha2) new LineSegment(tet.corners[1], tet.corners[2]).toScene(scene, 0.003, Color.black, 1);
				if (tet.corners[1].distance(tet.corners[3]) <= alpha2) new LineSegment(tet.corners[1], tet.corners[3]).toScene(scene, 0.003, Color.black, 1);
				if (tet.corners[2].distance(tet.corners[3]) <= alpha2) new LineSegment(tet.corners[2], tet.corners[3]).toScene(scene, 0.003, Color.black, 1);
*/			}
		}
		for (Vertex v : vertices) {
			if (v.getType() == VertexType.S) v.toScene(scene, 0.01, Color.red);
			else {
				v.toScene(scene, 0.01, Color.blue, 1);
				Circle c = new Circle(new Point(0, 0, v.z()), v, getRotationAxis());
			}
			if (this.vertices.size() < 20) scene.addText(String.valueOf(v.getId()), v);
		}
	}
	
	public static void main(String[] args){
		System.out.println(java.lang.Runtime.getRuntime().maxMemory());
		KineticDelaunayTessellation kDT = new KineticDelaunayTessellation(1000000, ProblemInstanceType.random);
		
/*		if (kDT.getInstanceType() == ProblemInstanceType.pdbNCC) {
			if (kDT.testingScreen) {
				kDT.toSceneBackbone(kDT.scene, 0.1, Color.blue);
//				kDT.toSceneSurface(kDT.scene);
			}
		}
		if (kDT.getInstanceType() == ProblemInstanceType.random) {
			if (kDT.testingScreen) {
				for (Vertex v : kDT.vertices) {
					if (v.getType() == VertexType.S) v.toScene(kDT.scene, 0.01, Color.red);
					else {
						v.toScene(kDT.scene, 0.01, Color.blue);
//						Circle c = new Circle(new Point(0, 0, v.z()), v, kDT.getRotationAxis());
					}
				}
			}
		}	
		HashSet<Tet> alphaTets = new HashSet<Tet>();
		HashSet<Tet> processedTets = new HashSet<Tet>();
		HashSet<Vertex> processedVers = new HashSet<Vertex>();
		long start = System.nanoTime();
		Vertex a = kDT.getVertex(40);
		a.getTet().getVertexSharingTetrahedra(a, kDT.alpha, alphaTets, processedTets, processedVers);
		long totalTime = System.nanoTime() - start;
		a.toScene(kDT.scene, 0.05, Color.yellow, 32);
		for (Tet nTet : processedTets) {
			if (nTet.isAlpha(kDT.alpha)) {
				nTet.toSceneEdges(kDT.scene, Color.black, 0.003);
				nTet.toSceneFaces(kDT.scene, Color.red); //new Color(255,0,0,50));
			}
		}
		
		start = System.nanoTime();
		HashSet<Tet> processed = (HashSet<Tet>)processedTets.clone();
		HashSet<Vertex> processedVersB = new HashSet<Vertex>();
		for (Vertex b : processedVers) {
			b.getTet().getVertexSharingTetrahedra(b, kDT.alpha, alphaTets, processedTets, processedVersB);
		}
		totalTime += System.nanoTime() - start;
		System.out.printf(" time in miliseconds %.2f\n", totalTime/1000000.0);
		for (Tet nTet : processedTets) {
			if (!processed.contains(nTet) && alphaTets.contains(nTet)) {
				nTet.toSceneEdges(kDT.scene, Color.black, 0.003);
				nTet.toSceneFaces(kDT.scene, Color.green); //new Color(255,0,0,50));
			}
		}
		
		
*/		
		long start = System.nanoTime();

		kDT.rotate();
		System.out.printf(" time in miliseconds %.2f\n", (System.nanoTime() - start)/1000000.0);
	}

	public ProblemInstanceType getInstanceType() {
		return instanceType;
	}

	public void setInstanceType(ProblemInstanceType instanceType) {
		this.instanceType = instanceType;
	}

	public double getAlpha() {
		return alpha;
	}

	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}
	
}
