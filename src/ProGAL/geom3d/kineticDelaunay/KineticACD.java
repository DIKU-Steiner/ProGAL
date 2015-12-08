package ProGAL.geom3d.kineticDelaunay; 

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import ProGAL.dataStructures.Heap;
import ProGAL.dataStructures.SortTool;
import ProGAL.geom3d.Circle;
import ProGAL.geom3d.LineSegment;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.kineticDelaunay.Vertex.VertexType;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.geom3d.volumes.Tetrahedron;
import ProGAL.geom3d.volumes.Torus; //
import ProGAL.geom3d.Vector;
import ProGAL.math.Constants; //
import ProGAL.math.Functions;
import ProGAL.math.Matrix;
import ProGAL.math.Randomization;
import ProGAL.math.Trigonometry;
import ProGAL.proteins.PDBFile; //



public class KineticACD {
	private static enum Direction { CW, CCW }
	public static enum ProblemInstanceType { random, pdb1, pdb1_AA, pdb2, chain, chain2, pdbNCC, toy, cro }

	private ProblemInstanceType instanceType;
	private List<Vertex> vertices = new ArrayList<Vertex>();
	private Set<Tet> tets = new HashSet<Tet>();
//	private Map<TetPoints,Tet> mapTets = new HashMap<TetPoints,Tet>();
//	private List<Tri> tris = new LinkedList<Tri>();
	private Map<TrianglePoints,Tri> mapTris = new HashMap<TrianglePoints,Tri>();
//	private List<Edge> edges = new LinkedList<Edge>();
	private Map<EdgePoints,Edge> mapEdges = new HashMap<EdgePoints,Edge>();
	private Set<Edge> alphaEdges = new HashSet<Edge>();
	private Set<Tri> alphaTris = new HashSet<Tri>();
	private Set<Tet> alphaTets = new HashSet<Tet>();
	private Set<Tri> initDTtriangles = new HashSet<Tri>();
	private Set<Edge> initDTedges = new HashSet<Edge>();
	private double alphaVal = Double.POSITIVE_INFINITY;
	public int nrFlips = 0;
	private PDBFile f;
	boolean clashFirst = false;
	boolean clash = false;
	List<Double> clashes = new ArrayList<Double>();
	double shortestEdge = Double.POSITIVE_INFINITY;
	int tet4Rot = 0;
	int tri3Rot = 0;
	int edg2Rot = 0;
	int tet2Rot = 0;
	public double error = 0;
	public double maxError = 0;
	public int nrErrors = 0;
	public double numericalTime = 0;
	public double analyticalTime = 0;
	public double flipTime = 0;
	
	private Tet lastTet;
	private Double[] angles = new Double[4]; // Daisy: changed from 2 to 4
	public double angleTotal = 0.0;
	private double angleLimit = Constants.TAU;        
	
	private Point rotationPoint;   // rotation center
	private Vector rotationAxis; // rotation vector
	private Direction rotDir;
	public List<Integer> rotIndx = new ArrayList<Integer>();     // indicies of vertices that move

	private Heap heap = new Heap(this.vertices.size(), new SortToolHeapItems());
//	private boolean testing = true;
	private boolean testingPrint = false;
	private boolean testingScreen = false;
	private boolean screenAlpha = false;
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

//	private double alpha;
	
	private class HeapItem {
		private Double[] angles;
		private Tet t = null;
		private Tet nt = null;
		private Tri tri = null;
		private Edge edg = null;
		
		private HeapItem(Double[] angles, Tet t, Tet nt, Tri tri, Edge e) {
			this.angles = angles;
			this.t = t;
			this.tri = tri;
			this.nt = nt;
			this.edg = e;
		}
				
		private Double[] getAngles() { return angles;} 
		private Tet getT()  { return t; }
		private Tet getNT() { return nt; }
		private Tri getTri() { return tri; }
		private Edge getEdg() { return edg; }
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
	
	public KineticACD(double alpha, List<Point> points) {
		if (testingScreen || screenAlpha)
			scene = null;
			//scene  = J3DScene.createJ3DSceneInFrame();

		Tetrahedron bigT = Tetrahedron.regularTetrahedron();
		bigT.blowUp(1000);
		lastTet = new Tet(bigT);
		for(int i=0;i<4;i++) 
			vertices.add(lastTet.getCorner(i));
		
		
		for (Point p: points) {
			insertPoint(p);	
		}
		lastTet = null;
//		System.out.println("Initial coords for vertex indexed at 10 = "+vertices.get(10).getCoord(0)+" , "+vertices.get(10).getCoord(1)+" , "+vertices.get(10).getCoord(2));
		this.setAlpha(alpha);
		
		/** Each vertex gets a pointer to one of its faces */
		for (Tet tet : tets) {
			for (int i = 0; i < 4; i++) tet.getCorner(i).setTet(tet);
		}
		setDirection(KineticACD.Direction.CCW);
		initializeAlphaComplex();
	}
	
	/** Adds new event to the heap */
	private void addToHeap(Double[] angles, Tet t, Tet nt) {
		heap.insert(new HeapItem(angles, t, nt, null, null));
	}
	private void addToHeap(Double[] angles, Tet t) {
		HeapItem newItem =  new HeapItem(angles, t, null, null, null);
		heap.insert(newItem);
	}
	private void addToHeap(Double[] angles, Tri t) {
		heap.insert(new HeapItem(angles, null, null, t, null));
	}
	private void addToHeap(Double[] angles, Edge e) {
		heap.insert(new HeapItem(angles, null, null, null, e));
	}
	
	/** Rotates and draws kinetic DT */
	public void animate(J3DScene scene, double alpha) {
		int steps = 1;
		double angleStep = alpha/steps;
		for (int k = 0; k < steps; k++) {
			for (int i = 0; i < rotIndx.size(); i++) {
				Vertex v = vertices.get(rotIndx.get(i));
				v.rotationCW(rotationAxis, angleStep);
				if (screenAlpha) v.toScene(scene, Color.RED, 0.1);
			}
			
			if ( (testingScreen && sphereAnimation) || screenAlpha){ 
				try { Thread.sleep(30); } catch (InterruptedException e) {}
				animateSpheres();
			}
			//if (k==0) isAlpha(); 
//			System.out.println("AlphaComplex? "+isAlpha()+" : "+alphaTets.toString()+" , "+alphaTris.toString()+", "+alphaEdges.toString());

		}
	}
	private void isAlpha() {
		boolean ret = true;
		Set<Tri> trisHandl = new HashSet<Tri>();
		Set<Edge> edgHandl = new HashSet<Edge>();
		
		for (Tet t: tets) {
			if (t.toString().equals("[35, 319, 324, 326]")) {
				System.out.println(t+" is looked at. "+t.getAlph()+" r = "+t.getCircumSphereRadius()+" in AC? "+alphaTets.contains(t));
			}
			if (alphaTets.contains(t)) {
				if ((!t.isAlpha(alphaVal))) {
					System.out.println("Tet "+t+" wrong in alpha at "+angleTotal+" with radius of "+t.getCircumSphereRadius());
				}
			} else {
				if (t.isAlpha(alphaVal)) {
					System.out.println("Tet "+t+" miss from alpha at "+angleTotal+" with radius of "+t.getCircumSphereRadius());
					continue;
				}
				for (Tri tri : t.getTris()) {
					if (trisHandl.contains(tri)) continue;
					if (alphaTris.contains(tri)) {
						edgHandl.add(new Edge(tri.getCorner(0), tri.getCorner(1)));
						edgHandl.add(new Edge(tri.getCorner(1), tri.getCorner(2)));
						edgHandl.add(new Edge(tri.getCorner(0), tri.getCorner(2)));
						if ((!isGabriel(tri) || !tri.isAlpha(alphaVal))) {
							System.out.println("Triangle "+tri+" wrong in alpha at "+angleTotal+" with radius of "+tri.getCircumRadius()+" gab = "+isGabriel(tri)+" alph = "+tri.isAlpha(alphaVal));
						}
					} else {
						if (isGabriel(tri) && tri.isAlpha(alphaVal)) {
							System.out.println("Triangle "+tri+" miss from alpha at "+angleTotal+" with radius of "+tri.getCircumRadius()+" gab = "+isGabriel(tri)+" alph = "+tri.isAlpha(alphaVal));
						}
					}
					trisHandl.add(tri);
				}
				for (Edge e : t.getEdges()) {
					if (edgHandl.contains(e)) continue;
					if (alphaEdges.contains(e) && (!isGabriel(e) || !e.isAlpha(alphaVal))) {
						System.out.println("Edge "+e+" wrong in alpha at "+angleTotal+" with radius of "+e.getCircumRadius()+" count = "+e.getCount()+" e.alph = "+e.getAlph());
					} else if (!alphaEdges.contains(e) && isGabriel(e) && e.isAlpha(alphaVal)){
						System.out.println("Edge "+e+" miss from alpha at "+angleTotal+" with radius of "+e.getCircumRadius()+" count = "+e.getCount()+" e.alph = "+e.getAlph());
					}
					edgHandl.add(e);
				}
			}
		}
		
		
		
		
		
		
		
		
		
		
		
		
		
		
/*		
		for (Tet t : tets) {
			if (!t.isAlive()) {
				System.out.println("Tet "+t.toString()+" is not alive!");
				break;
			}
			if ((t.isAlpha(alphaVal) && (!(alphaTets.contains(t)))) || ((!(t.isAlpha(alphaVal))) && alphaTets.contains(t))) {
				if ((t.getCircumSphereRadius()<(alphaVal-Constants.EPSILON) || (t.getCircumSphereRadius()>(alphaVal+Constants.EPSILON)))) {
					System.out.println("Wrong tet in alpha : "+t.toString()+" radius = "+t.getCircumSphereRadius()+" count = "+t.getCount()+" at angle "+angleTotal);
					double e = Math.abs(t.getCircumSphereRadius()-alphaVal);
					error += e;
					if (e>maxError) maxError=e;
					nrErrors++;
				}
				ret = false;
			}
			for (Tri tri : t.getTris()) {
				if (trisHandl.contains(tri)) continue;
				if (!tri.isAlive()) {
					System.out.println("Tri "+tri.toString()+" is not alive!");
					break;
				}
				if ((tri.isAlpha(alphaVal) && (!(alphaTris.contains(tri)))) || ((!(tri.isAlpha(alphaVal))) && alphaTris.contains(tri))) {
					if ((tri.getCircumRadius()<(alphaVal-Constants.EPSILON) || (tri.getCircumRadius()>(alphaVal+Constants.EPSILON)))) {
						System.out.println("Wrong tri in alpha : "+tri.toString()+" radius = "+tri.getCircumRadius()+" count = "+tri.getCount()+" at angle "+angleTotal);
					}
					ret = false;
				}
				trisHandl.add(tri);
			}
			for (Edge e : t.getEdges()) {
				if (edgHandl.contains(e)) continue;
				if (!e.isAlive()) {
					System.out.println("Edge "+e.toString()+" is not alive!");
					break;
				}
				if ((e.isAlpha(alphaVal) && (!(alphaEdges.contains(e)))) || ((!(e.isAlpha(alphaVal))) && alphaEdges.contains(e))) {
					if ((e.getCircumRadius()<(alphaVal-Constants.EPSILON) || (e.getCircumRadius()>(alphaVal+Constants.EPSILON)))) {
						System.out.println("Wrong edge in alpha : "+e.toString()+" radius = "+e.getCircumRadius()+" count = "+e.getCount()+" at angle "+angleTotal);
					}
					return false;
				}
				edgHandl.add(e);
			}
		}
		return ret;*/
	}
	
	/** Returns the rotation direction (clockwise or counterclockwise */
	public Direction getDirection() { return rotDir; }

	/** Reads or creates a problem instance */
	private List<Point> getInstance(ProblemInstanceType instanceType, int size) {
		
		if ( getInstanceType() == ProblemInstanceType.cro) {
			f = new PDBFile("/home/daisy/Desktop/KinAlpha/VoidRemover/Experiments/2cro.pdb", false);
			return f.getAtomCoords();
		}
		if ( getInstanceType() == ProblemInstanceType.chain2) {
			f = new PDBFile("/home/daisy/Downloads/2oed_cs_244_samples/sample_000267200000_2_91.576430.pdb", true);
			return f.getAtomCoords();
		}
		if ( getInstanceType() == ProblemInstanceType.chain) {
			f = new PDBFile("/home/daisy/Downloads/2oed_cs_244_samples/sample_000000200000_2_284.178215.pdb", true);
			return f.getAtomCoords();
		}
		if (instanceType == ProblemInstanceType.pdb2) {
			f = new PDBFile("/home/daisy/Downloads/2oed_cs_244_samples/sample_000000200000_2_284.178215.pdb", false);
			return f.getAtomCoords();
		}
		if (instanceType == ProblemInstanceType.pdb1 || instanceType == ProblemInstanceType.pdb1_AA) {
			f = new PDBFile("/home/daisy/Downloads/1X0O.pdb", false);
			return f.getAtomCoords();
		}
		if (instanceType == ProblemInstanceType.pdbNCC) {
			// reads a pdb-file - at the moment reads the backbone atoms only
			PDBFile f = new PDBFile("/Users/pawel/Downloads/1X0O.pdb", false);
			return f.getAtomCoords("N,Ca,C");	
		}
		if (instanceType == ProblemInstanceType.random) {
		// Generates a set of uniformly distributed points inside a unit cube
			Randomization.seed(17);
			return PointList.generatePointsInCube(size, -2.5, 2.5, -2.5, 2.5, -2.5, 2.5);
		}
		if (instanceType == ProblemInstanceType.toy) {
		// Generates a toy example
			List<Point> points = new java.util.LinkedList<Point>();
			points.add(new Point(0,0,0));
			points.add(new Point(1.5,0.5,0));
			points.add(new Point(0,1.5,0.5));
			points.add(new Point(0,0.5,1.5));
			points.add(new Point(1.1,1.1,1.1));
			points.add(new Point(0.1,1.5,0.5));
			points.add(new Point(1.1,1.5,0.5));
			points.add(new Point(0.1,1.3,0.2));
			points.add(new Point(0.4,0.5,1.5));
			points.add(new Point(1.0,1.2,1.2));
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
//	public Set<Tet> getTets() { return tets; }

	/** Returns i-th vertex */
	public Vertex getVertex(int i) { return vertices.get(i+4); }
	
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
	
	private void initializeRadiusEvents() {
		initDT();
		
		//Test each tetrahedron
		for (Tet T : tets) {
			if (T.isAlive()) {
				angles = getRoot(T.getCorner(0), T.getCorner(1), T.getCorner(2), T.getCorner(3), T.getCount());
				if ((angles != null) && (angles[0] < angleLimit)) {
					addToHeap(angles, T);
				}
			}
		}
			
		//Test each triangle
		for (Tri t : initDTtriangles) {
			angles = getRoot(t.getCorner(0), t.getCorner(1), t.getCorner(2), t.getCount());
			if ((angles != null) && (angles[0] < angleLimit)) {
				addToHeap(angles, t);
			}
		}
		
		//Test each edge
		for (Edge e : initDTedges) {
			angles = getRoot(e.getCorner(0), e.getCorner(1), e.getCount());
			if ((angles != null) && (angles[0] < angleLimit)) {
				addToHeap(angles, e);
			}
		}
	}
	
	// Daisy
	private boolean isGabriel(Edge e) {
		Sphere s = new Sphere(Point.getMidpoint(e.getCorner(0), e.getCorner(1)), e.getCircumRadius());
		for (Vertex p : vertices) {
//			if (e.hasVertex(new Vertex(p))) continue;
			if (!e.hasVertex(p) && s.contains(p)) return false;
		}
		return true;
	}
	private boolean isGabriel(Tri t) {
		Sphere s = new Sphere(new Circle(t.getCorner(0), t.getCorner(1), t.getCorner(2)).getCenter(), t.getCircumRadius());
		for (Vertex p : vertices) {
//			if (t.hasVertex(new Vertex(p))) continue;
//			if (s.contains(p)) return false;
			if( !t.hasVertex(p) && s.contains(p) ) return false;
		}
		return true;
	}
	private void initializeAlphaComplex() {
		for (Tet t : tets) {
			if (!t.isAlive()) continue;
//			System.out.println("Tet = "+t.toString());
			Double[] angles;
			Vertex v0 = t.getCorner(0);
			Vertex v1 = t.getCorner(1);
			Vertex v2 = t.getCorner(2);
			Vertex v3 = t.getCorner(3);
			
			if (t.isAlpha(alphaVal)) {
				alphaTets.add(t);
				t.setAlph(1);
			} else t.setAlph(0);
			Tri tri;
			if (!mapTris.containsKey(new TrianglePoints(v0, v1, v2))) {
				tri = new Tri(v0, v1, v2);
				mapTris.put(new TrianglePoints(v0, v1, v2), tri);
				t.setTri(tri, t.indexOf(v3));
//				tri.setTet(t);
				if (tri.isAlpha(alphaVal)) {
					tri.setAlph(1);
					if (isGabriel(tri)) {
						alphaTris.add(tri);
					}
				} else tri.setAlph(0);
			} else {
				tri = mapTris.get(new TrianglePoints(v0, v1, v2));
//				tri.setTet(t);
				t.setTri(tri, t.indexOf(v3));
			}
			if (!mapTris.containsKey(new TrianglePoints(v0, v1, v3))) {
				tri = new Tri(v0, v1, v3);
				mapTris.put(new TrianglePoints(v0, v1, v3), tri);
				t.setTri(tri, t.indexOf(v2));
//				tri.setTet(t);
				if (tri.isAlpha(alphaVal)) {
					tri.setAlph(1);
					if (isGabriel(tri)) {
						alphaTris.add(tri);
					}
				} else tri.setAlph(0);
			} else {
				tri = mapTris.get(new TrianglePoints(v0, v1, v3));
//				tri.setTet(t);
				t.setTri(tri, t.indexOf(v2));
			}
			if (!mapTris.containsKey(new TrianglePoints(v1, v2, v3))) {
				tri = new Tri(v1, v2, v3);
				mapTris.put(new TrianglePoints(v1, v2, v3), tri);
				t.setTri(tri, t.indexOf(v0));
//				tri.setTet(t);
				if (tri.isAlpha(alphaVal)) {
					tri.setAlph(1);
					if (isGabriel(tri)) {
						alphaTris.add(tri);
					}
				} else tri.setAlph(0);
			} else {
				tri = mapTris.get(new TrianglePoints(v1, v2, v3));
//				tri.setTet(t);
				t.setTri(tri, t.indexOf(v0));
			}
			if (!mapTris.containsKey(new TrianglePoints(v0, v2, v3))) {
				tri = new Tri(v0, v2, v3);
				mapTris.put(new TrianglePoints(v0, v2, v3), tri);
				t.setTri(tri, t.indexOf(v1));
//				tri.setTet(t);
				if (tri.isAlpha(alphaVal)) {
					tri.setAlph(1);
					if (isGabriel(tri)) {
						alphaTris.add(tri);
					}
				} else tri.setAlph(0);
			} else {
				tri = mapTris.get(new TrianglePoints(v0, v2, v3));
//				tri.setTet(t);
				t.setTri(tri, t.indexOf(v1));
			}
			
			//Test each edge
			Edge tmp;
			if (!mapEdges.containsKey(new EdgePoints(v0, v1))) {
				tmp =  new Edge(v0, v1);
				if (tmp.getLength()<shortestEdge) shortestEdge = tmp.getLength();
				mapEdges.put(new EdgePoints(v0, v1), tmp);
				//edges.add(tmp);
//				tmp.setTet(t);
				t.setEdge(tmp);
				if (tmp.isAlpha(alphaVal)) {
					tmp.setAlph(true);
					if (isGabriel(tmp)) {
						alphaEdges.add(tmp);
					}
				}
				angles = getRoot(v0, v1, tmp.getCount());
				if ((angles != null) && (angles[0] < angleLimit)) {
					addToHeap(angles, tmp);
				}
			} else t.setEdge(mapEdges.get(new EdgePoints(v0, v1)));
			
			if (!mapEdges.containsKey(new EdgePoints(v0, v2))) {
				tmp =  new Edge(v0, v2);
				if (tmp.getLength()<shortestEdge) shortestEdge = tmp.getLength();
				mapEdges.put(new EdgePoints(v0, v2), tmp);
				//edges.add(tmp);
//				tmp.setTet(t);
				t.setEdge(tmp);
				if (tmp.isAlpha(alphaVal)) {
					tmp.setAlph(true);
					if (isGabriel(tmp)) {
						alphaEdges.add(tmp);
					}
				}
				angles = getRoot(v0, v2, tmp.getCount());
				if ((angles != null) && (angles[0] < angleLimit)) {
					addToHeap(angles, tmp);
				}
			} else t.setEdge(mapEdges.get(new EdgePoints(v0, v2)));
			
			if (!mapEdges.containsKey(new EdgePoints(v0, v3))) {
				tmp =  new Edge(v0, v3);
				if (tmp.getLength()<shortestEdge) shortestEdge = tmp.getLength();
				mapEdges.put(new EdgePoints(v0, v3), tmp);
//				tmp.setTet(t);
				t.setEdge(tmp);
				if (tmp.isAlpha(alphaVal)) {
					tmp.setAlph(true);
					if (isGabriel(tmp)) {
						alphaEdges.add(tmp);
					}
				}
				angles = getRoot(v0, v3, tmp.getCount());
				if ((angles != null) && (angles[0] < angleLimit)) {
					addToHeap(angles, tmp);
				}
			} else t.setEdge(mapEdges.get(new EdgePoints(v0, v3)));
			
			if (!mapEdges.containsKey(new EdgePoints(v1, v2))) {
				tmp =  new Edge(v1, v2);
				if (tmp.getLength()<shortestEdge) shortestEdge = tmp.getLength();
				mapEdges.put(new EdgePoints(v1, v2), tmp);
//				tmp.setTet(t);
				t.setEdge(tmp);
				if (tmp.isAlpha(alphaVal)) {
					tmp.setAlph(true);
					if (isGabriel(tmp)) {
						alphaEdges.add(tmp);
					}
				}
				angles = getRoot(v1, v2, tmp.getCount());
				if ((angles != null) && (angles[0] < angleLimit)) {
					addToHeap(angles, tmp);
				}
			} else t.setEdge(mapEdges.get(new EdgePoints(v1, v2)));
			
			if (!mapEdges.containsKey(new EdgePoints(v1, v3))) {
				tmp =  new Edge(v1, v3);
				if (tmp.getLength()<shortestEdge) shortestEdge = tmp.getLength();
				mapEdges.put(new EdgePoints(v1, v3), tmp);
//				tmp.setTet(t);
				t.setEdge(tmp);
				if (tmp.isAlpha(alphaVal)) {
					tmp.setAlph(true);
					if (isGabriel(tmp)) {
						alphaEdges.add(tmp);
					}
				}
				angles = getRoot(v1, v3, tmp.getCount());
				if ((angles != null) && (angles[0] < angleLimit)) {
					addToHeap(angles, tmp);
				}
			} else t.setEdge(mapEdges.get(new EdgePoints(v1, v3)));
			
			if (!mapEdges.containsKey(new EdgePoints(v2, v3))) {
				tmp =  new Edge(v2, v3);
				if (tmp.getLength()<shortestEdge) shortestEdge = tmp.getLength();
				mapEdges.put(new EdgePoints(v2, v3), tmp);
//				tmp.setTet(t);
				t.setEdge(tmp);
				if (tmp.isAlpha(alphaVal)) {
					tmp.setAlph(true);
					if (isGabriel(tmp)) {
						alphaEdges.add(tmp);
					}
				}
				angles = getRoot(v2, v3, tmp.getCount());
				if ((angles != null) && (angles[0] < angleLimit)) {
					addToHeap(angles, tmp);
				}
			} else t.setEdge(mapEdges.get(new EdgePoints(v2, v3)));
		}
	}
	
	private void noRotate(Tet t, Vertex s0, Vertex s1, Vertex s2, Vertex s3) {
		if (t.isAlpha(alphaVal)) {
			alphaTets.add(t);
			t.setAlph(1);
		}
		Tri tri;
		if (!mapTris.containsKey(new TrianglePoints(s0, s1, s2))) {
			tri = new Tri(s0, s1, s2);
			mapTris.put(new TrianglePoints(s0, s1, s2), tri);
			t.setTri(tri, t.indexOf(s3));
//			tri.setTet(t);
			if (tri.isAlpha(alphaVal)) {
				tri.setAlph(1);
				if (isGabriel(tri)) {
					alphaTris.add(tri);
				}
			}
		} else {
			tri = mapTris.get(new TrianglePoints(s0, s1, s2));
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(s3));
		}
		if (!mapTris.containsKey(new TrianglePoints(s0, s1, s3))) {
			tri = new Tri(s0, s1, s3);
			mapTris.put(new TrianglePoints(s0, s1, s3), tri);
			t.setTri(tri, t.indexOf(s2));
//			tri.setTet(t);
			if (tri.isAlpha(alphaVal)) {
				tri.setAlph(1);
				if (isGabriel(tri)) {
					alphaTris.add(tri);
				}
			}
		} else {
			tri = mapTris.get(new TrianglePoints(s0, s1, s3));
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(s2));
		}
		if (!mapTris.containsKey(new TrianglePoints(s1, s2, s3))) {
			tri = new Tri(s1, s2, s3);
			mapTris.put(new TrianglePoints(s1, s2, s3), tri);
			t.setTri(tri, t.indexOf(s0));
//			tri.setTet(t);
			if (tri.isAlpha(alphaVal)) {
				tri.setAlph(1);
				if (isGabriel(tri)) {
					alphaTris.add(tri);
				}
			}
		} else {
			tri = mapTris.get(new TrianglePoints(s1, s2, s3));
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(s0));
		}
		if (!mapTris.containsKey(new TrianglePoints(s0, s2, s3))) {
			tri = new Tri(s0, s2, s3);
			mapTris.put(new TrianglePoints(s0, s2, s3), tri);
			t.setTri(tri, t.indexOf(s1));
//			tri.setTet(t);
			if (tri.isAlpha(alphaVal)) {
				tri.setAlph(1);
				if (isGabriel(tri)) {
					alphaTris.add(tri);
				}
			}
		} else {
			tri = mapTris.get(new TrianglePoints(s0, s2, s3));
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(s1));
		}
	}
	private void oneRotate(Tet t, Vertex s0, Vertex s1, Vertex s2, Vertex r) {
		double d0 = s0.distanceSquared(s1);
		double d1 = s1.distanceSquared(s2);
		double d2 = s2.distanceSquared(s0);
		double a0 = r.distanceSquared(s0);
		double a1 = r.distanceSquared(s1);
		double a2 = r.distanceSquared(s2);
		if (t.isAlpha(alphaVal)) {
			alphaTets.add(t);
			t.setAlph(1);
		} else if (a0 < d0 && a1 < d0 && a1 < d1 && a2 < d2 && a2 < d2 && a0 < d2) {
			t.setAlph(2);
		} else t.setAlph(0);
		angles = getRoot(t.getCorner(0), t.getCorner(1), t.getCorner(2), t.getCorner(3), t.getCount());
		if ((angles != null) && (angles[0] < angleLimit)) {
			addToHeap(angles, t);
		}
		Tri tri;
		if (!mapTris.containsKey(new TrianglePoints(s0, s1, s2))) {
			tri = new Tri(s0, s1, s2);
			mapTris.put(new TrianglePoints(s0, s1, s2), tri);
//			tri.setTet(t);
			t.setTri(tri, t.indexOf_slow(r));
			if (tri.isAlpha(alphaVal)) {
				tri.setAlph(1);
				if (isGabriel(tri)) {
					alphaTris.add(tri);
				}
			}
		} else {
			tri = mapTris.get(new TrianglePoints(s0, s1, s2));
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(r));
		}
		if (!mapTris.containsKey(new TrianglePoints(s0, s1, r))) {
			tri = new Tri(s0, s1, r);
			mapTris.put(new TrianglePoints(s0, s1, r), tri);
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(s2));
			if (tri.isAlpha(alphaVal)) {
				tri.setAlph(1);
				if (isGabriel(tri)) {
					alphaTris.add(tri);
				}
			} else {
				Point center = Point.getMidpoint(s0, s1);
				Circle C = new Circle(center, alphaVal, center.vectorTo(s0));
				double dist = C.getDistanceSquared(r);
				if (dist < (alphaVal*alphaVal)+Constants.EPSILON) {
					tri.setAlph(2);
				} else tri.setAlph(0);
			}
		} else {
			tri = mapTris.get(new TrianglePoints(s0, s1, r));
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(s2));
		}
		if (!mapTris.containsKey(new TrianglePoints(s0, s2, r))) {
			tri = new Tri(s0, s2, r);
			mapTris.put(new TrianglePoints(s0, s2, r), tri);
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(s1));
			if (tri.isAlpha(alphaVal)) {
				tri.setAlph(1);
				if (isGabriel(tri)) {
					alphaTris.add(tri);
				}
			} else {
				Point center = Point.getMidpoint(s0, s2);
				Circle C = new Circle(center, alphaVal, center.vectorTo(s0));
				double dist = C.getDistanceSquared(r);
				if (dist < (alphaVal*alphaVal)+Constants.EPSILON) {
					tri.setAlph(2);
				} else tri.setAlph(0);
			}
		} else {
			tri = mapTris.get(new TrianglePoints(s0, s2, r));
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(s1));
		}
		if (!mapTris.containsKey(new TrianglePoints(s1, s2, r))) {
			tri = new Tri(s1, s2, r);
			mapTris.put(new TrianglePoints(s1, s2, r), tri);
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(s0));
			if (tri.isAlpha(alphaVal)) {
				tri.setAlph(1);
				if (isGabriel(tri)) {
					alphaTris.add(tri);
				}
			} else {
				Point center = Point.getMidpoint(s1, s2);
				Circle C = new Circle(center, alphaVal, center.vectorTo(s1));
				double dist = C.getDistanceSquared(r);
				if (dist < (alphaVal*alphaVal)+Constants.EPSILON) {
					tri.setAlph(2);
				} else tri.setAlph(0);
			}
		} else {
			tri = mapTris.get(new TrianglePoints(s1, s2, r));
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(s0));
		}
	}
	
	
	private void twoRotate(Tet t, Vertex s0, Vertex s1, Vertex r0, Vertex r1) {
		if (t.isAlpha(alphaVal)) {
			alphaTets.add(t);
			t.setAlph(1);
		}
		angles = getRoot(t.getCorner(0), t.getCorner(1), t.getCorner(2), t.getCorner(3), t.getCount());
		if ((angles != null) && (angles[0] < angleLimit)) {
			addToHeap(angles, t);
		}
		Tri tri;
		if (!mapTris.containsKey(new TrianglePoints(s0, s1, r0))) {
			tri = new Tri(s0, s1, r0);
			mapTris.put(new TrianglePoints(s0, s1, r0), tri);
//			tri.setTet(t);
			t.setTri(tri, t.indexOf_slow(r1));
			if (tri.isAlpha(alphaVal)) {
				tri.setAlph(1);
				if (isGabriel(tri)) {
					alphaTris.add(tri);
				}
			} else {
				Point center = Point.getMidpoint(s0, s1);
				Circle C = new Circle(center, alphaVal, center.vectorTo(s0));
				double dist = C.getDistanceSquared(r0);
				if (dist < (alphaVal*alphaVal)+Constants.EPSILON) {
					tri.setAlph(2);
				} else tri.setAlph(0);
			}
		} else {
			tri = mapTris.get(new TrianglePoints(s0, s1, r0));
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(r1));
		}
		if (!mapTris.containsKey(new TrianglePoints(s0, s1, r1))) {
			tri = new Tri(s0, s1, r1);
			mapTris.put(new TrianglePoints(s0, s1, r1), tri);
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(r0));
			if (tri.isAlpha(alphaVal)) {
				tri.setAlph(1);
				if (isGabriel(tri)) {
					alphaTris.add(tri);
				}
			} else {
				Point center = Point.getMidpoint(s0, s1);
				Circle C = new Circle(center, alphaVal, center.vectorTo(s0));
				double dist = C.getDistanceSquared(r1);
				if (dist < (alphaVal*alphaVal)+Constants.EPSILON) {
					tri.setAlph(2);
				} else tri.setAlph(0);
			}
		} else {
			tri = mapTris.get(new TrianglePoints(s0, s1, r1));
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(r0));
		}
		if (!mapTris.containsKey(new TrianglePoints(s0, r0, r1))) {
			tri = new Tri(s0, r0, r1);
			mapTris.put(new TrianglePoints(s0, r0, r1), tri);
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(s1));
			if (tri.isAlpha(alphaVal)) {
				tri.setAlph(1);
				if (isGabriel(tri)) {
					alphaTris.add(tri);
				}
			} else {
				Point center = Point.getMidpoint(r0, r1);
				Circle C = new Circle(center, alphaVal, center.vectorTo(r0));
				double dist = C.getDistanceSquared(s0);
				if (dist < (alphaVal*alphaVal)+Constants.EPSILON) {
					tri.setAlph(2);
				} else tri.setAlph(0);
			}
		} else {
			tri = mapTris.get(new TrianglePoints(s0, r0, r1));
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(s1));
		}
		if (!mapTris.containsKey(new TrianglePoints(s1, r0, r1))) {
			tri = new Tri(s1, r0, r1);
			mapTris.put(new TrianglePoints(s1, r0, r1), tri);
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(s0));
			if (tri.isAlpha(alphaVal)) {
				tri.setAlph(1);
				if (isGabriel(tri)) {
					alphaTris.add(tri);
				}
			} else {
				Point center = Point.getMidpoint(r0, r1);
				Circle C = new Circle(center, alphaVal, center.vectorTo(s1));
				double dist = C.getDistanceSquared(s1);
				if (dist < (alphaVal*alphaVal)+Constants.EPSILON) {
					tri.setAlph(2);
				} else tri.setAlph(0);
			}
		} else {
			tri = mapTris.get(new TrianglePoints(s1, r0, r1));
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(s0));
		}
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
	public void setRotationAxis(Vector v)  { rotationAxis = v.normalize(); }
	
	public void setRotationAxis(Point a, Point b)  {
		Vector translate = new Vector(-a.get(0),-a.get(1),-a.get(2));
		for (Point v : vertices) v.addThis(translate);

		Vector dir = a.vectorTo(b);
		rotationAxis = new Vector(0,0,1);
		if (!dir.isParallel(rotationAxis)) { //Rotate
			Vector rotAxis = dir.cross(rotationAxis);
			rotAxis.normalizeThis();
			double angle = dir.angle(rotationAxis); 
			Matrix rotMatrix = Matrix.createRotationMatrix(angle, rotAxis);
			for (Point v : vertices) {
				rotMatrix.multiplyIn(v);
			}
		}
/*		Vector dir = a.vectorTo(b);
		Vector e3 = new Vector(0,0,1);
		Point rotPoint = a;
		Point rotPoint1 = b;
		if (!dir.isParallel(e3)) { //Rotate
			Vector rotAxis = dir.cross(e3);
			rotAxis.normalizeThis();
			double angle = dir.angle(e3);
			Matrix rotMatrix = Matrix.createRotationMatrix(angle, rotAxis);
			for (Point v : vertices) {
				v = (Point)rotMatrix.multiplyIn(v);
			}
			rotPoint = rotMatrix.multiply(a);
			rotPoint1 = rotMatrix.multiply(b);
		}
		
		//Translate
		Line axis = new Line(dir);
		Point linePoint = axis.orthogonalProjection(rotPoint);
		Vector translate = rotPoint.vectorTo(linePoint);
		for (Point v : vertices) {
			v.addThis(translate);
		}
		rotPoint.addThis(translate);
		rotPoint1.addThis(translate);
		rotationAxis = new Vector(0,0,1);*/
	}

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

	//Daisy
	private double angleCW(Vector v1, Vector v2) {
		double ret = Math.atan2(v1.y(), v1.x()) - Math.atan2(v2.y(), v2.x());
		if (ret<0) { ret = Constants.TAU+ret; return ret; }
		if (ret>0) { /*ret = (-1)*ret;*/ return ret; }
		return ret;
	}
	private double angleCCW(Vector v1, Vector v2) {
		double ret = Math.atan2(v1.y(), v1.x()) - Math.atan2(v2.y(), v2.x());
		if (ret<0) { ret = (-1)*ret; return ret; }
		if (ret>0) { ret = Constants.TAU-ret; return ret; }
		return ret;
	}
	/** Create a sphere with center A and radius 2*alphaVal, check for intersection
	 *  with the circle path of B 												  **/
	private Double[] getRootSR(Vertex A, Vertex B, int dir) {
//		if (A.distanceSquared(B) > 4*alphaVal*alphaVal) { System.out.println("Distance A to B is too big : "+A.distance(B)+" squared : "+A.distanceSquared(B)); return null; }
		Double[] angles = new Double[4];
		Sphere S = new Sphere(A, 2*alphaVal);
		Circle C = new Circle(new Point(0.0,0.0,B.getCoord(2)), B.distanceXY(), rotationAxis);
		Point[] intersectionPoints = S.getIntersections(C);
		if (!(intersectionPoints==null) && intersectionPoints.length>1) {
			int i = 0;
			for (Point p : intersectionPoints) {
				Point z = new Point(0.0, 0.0, B.getCoord(2));
				if (dir == 1) {
					angles[i] = (angleTotal+angleCCW(z.vectorTo(p), z.vectorTo((Point)B)));
				} else	angles[i] = (angleTotal+angleCW(z.vectorTo(p), z.vectorTo((Point)B)));
				i += 1;
			}
		}
		
		angles = getRotAngle(angles, 0);
		return angles;
	}
	//Daisy : Create torus from static points and intersect it with the moving
	private Double[] getRootSSR(Vertex A, Vertex B, Vertex C, int dir) {
/*		if (A.getId()==118 && B.getId()==145 && C.getId()== 18) {
			System.out.println("NOW!");
		}*/
		Double[] angles = new Double[4];
		if (A.distanceSquared(B) >= 4*alphaVal*alphaVal) return null;
		Point torusCenter = Point.getMidpoint(A, B);
		//Find major radius using pythagoras:
		double R = Math.pow(Math.abs(alphaVal*alphaVal-A.distanceSquared(torusCenter)), 0.5);
		//Torus
		Vector torusNormal = A.vectorTo(torusCenter).normalize();
		Torus torus = new Torus(torusCenter, torusNormal, R, alphaVal);
		//Circle
		double circleRadius = C.distanceXY();
		Circle circle = new Circle(new Point(0,0,C.getCoord(2)), circleRadius, rotationAxis);
		
		Point[] intersections = torus.getIntersectionCircle(circle);
		int i = 0;
		if (!(intersections==null)) {
			for (Point p : intersections) {
				if (p==null) break;
				Point z = new Point(0.0, 0.0, C.getCoord(2));
				if (dir==1) {
					angles[i] = (angleTotal+angleCCW(z.vectorTo(p), z.vectorTo((Point)C)));
				} else angles[i] = (angleTotal+angleCW(z.vectorTo(p), z.vectorTo((Point)C)));
				i += 1;
			}
		}
		return getRotAngle(angles, 0);
	}
	// Daisy : Create two spheres from three static points and radius alphaVal. Check for intersections between these and the rotation circle
	private Double[] getRootSSSR(Vertex A, Vertex B, Vertex C, Vertex D, int dir) {
		Double[] angles = new Double[4];
		Circle c = new Circle(A, B, C);
		Point midpoint = c.getCenter();
		// Pythagoras
		double cMinusa = alphaVal*alphaVal-A.distanceSquared(midpoint);
		if (cMinusa<0) return null;
		double distToCenter = Math.sqrt(Math.abs(cMinusa));
		Vector vecToCenter = c.getNormal().clone();
		vecToCenter.scaleToLengthThis(distToCenter);
		
		Point center0 = midpoint.clone().addThis(vecToCenter);
		Point center1 = midpoint.clone().addThis(vecToCenter.multiply(-1.0));
		Circle rotationC = new Circle(new Point(0.0, 0.0, D.getCoord(2)), D.distanceXY(), rotationAxis);
//		rotationC.toScene(scene, 0.1, java.awt.Color.RED);
		Sphere S0 = new Sphere(center0, alphaVal);
		Point[] intersectionPoints1 = S0.getIntersections(rotationC);
		int i;
		if (!(intersectionPoints1==null)) {
			i = 0;
			for (Point p : intersectionPoints1) {
				Point z = new Point(0.0, 0.0, D.getCoord(2));
				if (dir==1) {
					angles[i] = (angleTotal+angleCCW(z.vectorTo(p), z.vectorTo((Point)D)));
				} else angles[i] = (angleTotal+angleCW(z.vectorTo(p), z.vectorTo((Point)D)));
				i += 1;
			}
		} else i = 0;
		Sphere S1 = new Sphere(center1, alphaVal);
		Point[] intersectionPoints2 = S1.getIntersections(rotationC);
		if (!(intersectionPoints2==null)) {
//			int i = 2;
			for (Point p : intersectionPoints2) {
				Point z = new Point(0.0, 0.0, D.getCoord(2));
				if (dir==1) {
					angles[i] = (angleTotal+angleCCW(z.vectorTo(p), z.vectorTo((Point)D)));
				} else angles[i] = (angleTotal+angleCW(z.vectorTo(p), z.vectorTo((Point)D)));
				i += 1;
			}
		}
		return getRotAngle(angles, 0);
	}
	
	// Daisy : numerical method
	private Double[] getRootSSRR(Vertex A, Vertex B, Vertex C, Vertex D, int dir) {
		tet2Rot++;
//		System.out.println("Called with ["+A.getId()+", "+B.getId()+", "+C.getId()+", "+D.getId()+"]");
		double precision = 0.05;
		Double[] angles = new Double[4];
		int numOfAngles = 0;
//		double offRadius = new Sphere(A, B, C, D).getRadius()-alphaVal;
		double offRadius = Sphere.computeSphere_fast(A, B, C, D, null)-alphaVal;
		double sign = Math.signum(offRadius);
		Point Cnew = C.clone();
		Point Dnew = D.clone();
//		if (dir == 1) k = -1;
		double angleSoFar = Math.floor(Math.toDegrees(angleTotal));
		double max = Math.ceil(Math.toDegrees(angleLimit)-angleSoFar);
		for (double i = precision ; i<=max ; i+=precision) {
//			Cnew = C.clone();
//			Dnew = D.clone();
			Cnew.set(C);
			Dnew.set(D);
			if (dir == 0) {
				Cnew.rotationCW(rotationAxis, Math.toRadians(i));
				Dnew.rotationCW(rotationAxis, Math.toRadians(i));
			} else {
//				Cnew.rotationCCW(Z, Math.toRadians(i));
//				Dnew.rotationCCW(Z, Math.toRadians(i));
			}
//			double radius = new Sphere(A, B, Cnew, Dnew).getRadius();
			double radius = Sphere.computeSphere_fast(A, B, Cnew, Dnew, null);
/*			if (A.getId()==5 && B.getId()==38 && C.getId()== 11 && D.getId()==13) {
				System.out.println("Radius = "+radius+" angletotal = "+angleTotal);
				System.out.println("Radius-alphaVal = "+(radius-alphaVal));
			}*/
			offRadius = radius-alphaVal;
			if (Math.signum(offRadius) != sign) {
/*				if (A.getId()==5 && B.getId()==38 && C.getId()== 11 && D.getId()==13) {
					System.out.println("sign changed");
				}*/
//				System.out.println("Change in sign at "+i);
				sign = Math.signum(offRadius);
				Double angle = findAngle(i-precision, i, A, B, C, D, 0);
				if (angle==null) {
/*					if (A.getId()==5 && B.getId()==38 && C.getId()== 11 && D.getId()==13) {
						System.out.println("Angles was null for this Tet");
					}*/
//					System.out.println("findAngle returned null");
					return null;
				}
				
				angles[numOfAngles] = Math.toRadians(angle)+angleTotal;
				
//				System.out.println("Angle in radians = "+angles[numOfAngles]);
				numOfAngles++;
			}
		}
		if (numOfAngles==0) {
			return null;
		} else return angles;
	}
	
	private Double findAngle(double start, double end, Vertex A, Vertex B, Vertex C, Vertex D, int dir) {
		Point Cnew = C.clone();
		Point Dnew = D.clone();
		if (dir == 0) {
			Cnew.rotationCW(rotationAxis, Math.toRadians(start));
			Dnew.rotationCW(rotationAxis, Math.toRadians(start));
		}
//		double radius = new Sphere(A, B, Cnew, Dnew).getRadius();
		double radius = Sphere.computeSphere_fast(A, B, Cnew, Dnew, null);
//		if (Math.abs(radius-alphaVal)<Math.pow(10, -9)) return start;
		if (Math.abs(radius-alphaVal)<0.000000001) return start;
		double sign = Math.signum(radius-alphaVal);
		double newStart = start;
		double newEnd = end;
		double angle = 0.0;
		double newSign;
		
		for (int i = 0; i<99 ; i++) {
			angle = (newEnd+newStart)/2.0;
			Cnew.set(C);
			Dnew.set(D);
//			Cnew = C.clone();
//			Dnew = D.clone();
			if (dir == 0) {
				Cnew.rotationCW(rotationAxis, Math.toRadians(angle));
				Dnew.rotationCW(rotationAxis, Math.toRadians(angle));
			}
//			radius = new Sphere(A, B, Cnew, Dnew).getRadius();
			radius = Sphere.computeSphere_fast(A, B, Cnew, Dnew, null);
//			if (Math.abs(radius-alphaVal)<Math.pow(10, -9)) {
			if (Math.abs(radius-alphaVal)<0.000000001) {
				return angle;
			}
			newSign = Math.signum(radius-alphaVal);
			if (newSign!=sign) {
				newEnd = angle;
			} else newStart = angle;
		}
		System.out.println("All 99 iterations used!!!! ");
		return angle;
	}
	
	private Double[] getRootSSSSR(Vertex A, Vertex B, Vertex C, Vertex D, Vertex E, int dir) {
		Double[] angles = new Double[4];
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
	    Double[] angles_tmp = Trigonometry.solveAsinXPlusBcosXplusC(coefSin, coefCos, coef);
	    if (angles_tmp!=null) {
		    for (int i = 0 ; i < angles_tmp.length ; i++) {
		    	angles[i] = angles_tmp[i];
		    }
	    }
	    return getRotAngle(angles, dir);
	}

	private Double[] getRootSSSRR(Vertex A, Vertex B, Vertex C, Vertex D, Vertex E, int dir) {
		Double[] angles = new Double[4];
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
	    Double[] angles_tmp = Trigonometry.solveAsinXPlusBcosXplusC(coefSin, coefCos, coef);
	    if (angles_tmp!=null) {
		    for (int i = 0 ; i < angles_tmp.length ; i++) {
		    	angles[i] = angles_tmp[i];
		    }
	    }
	    angles = getRotAngle(angles, dir);
	    return angles;
	}

	public Double[] getRotAngle(Double[] angles, int dir) {
		if (angles == null || (angles[1]==null && angles[3]==null)) return null;
		Double[] oldAngles = new Double[4];
		int i = 0;
		for (Double angle : angles) {
			if (angle == null) continue;
			oldAngles[i] = angle;
			i += 1;
		}
		Double[] newAngles = new Double[i];
		int j = 0;
		for (Double angle : oldAngles) {
			if (angle == null) break;
			if (dir == 1){
				angle = Constants.TAU-angle+angleTotal;
			}
			if (rotDir == Direction.CCW) {
				if (angle < angleTotal + Constants.EPSILON) angle = angleLimit;//Constants.TAU;
			}
			else {
				if (Constants.TAU - angle < angleTotal + Constants.EPSILON) angle = angleLimit;//Constants.TAU;
			}
			newAngles[j] = angle;
			j += 1;
		}
		Arrays.sort(newAngles);
		return newAngles;
	}

	// Edge
	public Double[] getRoot(Vertex v0, Vertex v1, int count) {
		if (count == 0 || count == 3) return null;
		if (count == 1) {
			return getRootSR(v1, v0, 0);
		}
		return getRootSR(v0, v1, 0);
	}
	
	// Triangle
	public Double[] getRoot(Vertex v0, Vertex v1, Vertex v2, int count) {
		if (count == 0 || count == 7) return null;
		if (count <=3 ) {
			if (count == 1) return getRootSSR(v1, v2, v0, 0);
			if (count == 2) return getRootSSR(v0, v2, v1, 0);
			return getRootSSR(v0, v1, v2, 1);
		}
		if (count <= 5) {
			if (count == 4) return getRootSSR(v0, v1, v2, 0);
			return getRootSSR(v0, v2, v1, 1);
		}
		return getRootSSR(v1, v2, v0, 1);
	}
	
	// Tetrahedron
	public Double[] getRoot(Vertex v0, Vertex v1, Vertex v2, Vertex v3, int count) {
		if ( count == 0 || count == 15 ) { return null; }
		if ( count <= 7 ) {
			if (count <=3 ) {
				if (count == 1) return getRootSSSR(v1, v2, v3, v0, 0);
				if (count == 2) return getRootSSSR(v0, v2, v3, v1, 0);
				long start = System.nanoTime();
				Double[] a = getRootSSRR(v2, v3, v0, v1, 0); 
				numericalTime += (System.nanoTime() - start)/1000000.0;
//				getRootSSSR(v2, v3, v0, v1, 0);
				return a;
			}
			if (count <= 5) {
				if (count == 4) return getRootSSSR(v0, v1, v3, v2, 0);
				long start = System.nanoTime();
				Double[] a = getRootSSRR(v1, v3, v0, v2, 0);
				numericalTime += (System.nanoTime() - start)/1000000.0;
//				getRootSSR(v1, v3, v0, 0);
				return a;
			}
			if (count == 6) {
				long start = System.nanoTime();
				Double[] a = getRootSSRR(v0, v3, v1, v2, 0);
				numericalTime += (System.nanoTime() - start)/1000000.0;
//				getRootSSR(v0, v3, v1, 0);
				return a;
			}
			return getRootSSSR(v0, v1, v2, v3, 1);
		}
		if (count <= 11) {
			if (count <= 9) {
				if (count == 8) return getRootSSSR(v0, v1, v2, v3, 0);
				long start = System.nanoTime();
				Double[] a = getRootSSRR(v1, v2, v0, v3, 0);
				numericalTime += (System.nanoTime() - start)/1000000.0;
//				getRootSSR(v1, v2, v0, 0);
				return a;
			}
			if (count == 10) {
				long start = System.nanoTime();
				Double[] a = getRootSSRR(v0, v2, v1, v3, 0);
				numericalTime += (System.nanoTime() - start)/1000000.0;
//				getRootSSR(v0, v2, v1, 0);
				return a;
			}
			return getRootSSSR(v0, v1, v3, v2, 1);
		}
		if (count <= 13) {
			if (count == 12) {
				long start = System.nanoTime();
				Double[] a = getRootSSRR(v0, v1, v2, v3, 0);
				numericalTime += (System.nanoTime() - start)/1000000.0;
//				getRootSSR(v0, v1, v2, 0);
				return a;
			}
			return getRootSSSR(v0, v2, v3, v1, 1);
		}
		return getRootSSSR(v1, v2, v3, v0, 1);
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
				return getRootSSSRR(t.getCorner(0), t.getCorner(2), oppV, t.getCorner(1), t.getCorner(3), 1);
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

	private void initDT() {
		initDTtriangles.clear();
		initDTedges.clear();
		for (Tet t : tets) {
			if (t.isAlive()) {
				initDTtriangles.addAll(Arrays.asList(t.getTris()));
				initDTedges.addAll(t.getEdges());
			}
		}
	}
	
	public void initializeRotation(List<Integer> rotIndices, int a, int b) {
		heap.clear();
		angleTotal = 0;
		
		if(rotIndices!=null){//Don't reuse existing axis
			setRotationAxis(vertices.get(a+4), vertices.get(b+4));
			setRotVertices(rotIndices);

			// set constants associated with vertices - does not include vertices of big points
			for (Vertex v : vertices) {
				if (testingPrint) System.out.print(v.getId() + ": " + v.toString(2));
				double zDist = Math.sqrt(v.get(0)*v.get(0) + v.get(1)*v.get(1));
				v.setSquaredPolarRadius(v.distanceSquared());
				if(zDist<Constants.EPSILON){
//				if (v.distanceSquared() < Constants.EPSILON) { //takes care of the special case when the rotation center overlaps with one of the given points
					
//					v.setSquaredPolarRadius(0.0);                            
					//				v.setPolarRadius(0.0);                        
					//				v.setPolarAngle(0.0);			
					if (testingPrint) System.out.println(", polar angle: 0.0"); 
					v.setCosAngle(1.0);
					v.setSinAngle(0.0);
					
					v.setType(Vertex.VertexType.S);
					rotIndx.remove((Integer)v.getId());
				} else {
//					v.setSquaredPolarRadius(v.distanceSquared());       // remains unchanged when rotating around the same point
					//				v.setPolarRadius(v.distance());                     // remains unchanged when rotating around the same point
					//				v.setPolarAngle(v.polarAngleXY());     			
					//				if (testingPrint) System.out.println(", initial polar angle: " + Functions.toDeg(v.getPolarAngle())); 
					v.setCosAngle(v.polarAngleCosXY());
					v.setSinAngle(v.polarAngleSinXY());
					
				}
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
		
		initializeRadiusEvents();
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
		long start = System.nanoTime();
		nrFlips++;
		Vertex vns0 = t0.getCorner(t0.apex(t1));
		Vertex vns1 = t1.getCorner(t1.apex(t0));
		
		Edge commonEdge = new Edge(vns0, vns1); //To be created in complex
		if (commonEdge.isAlpha(alphaVal)) { 
			commonEdge.setAlph(true);
			if (isGabriel(commonEdge)) {
				alphaEdges.add(commonEdge);
				if (!commonEdge.isBig() && screenAlpha) commonEdge.toSceneEdge(scene, Color.BLACK, 0.005);
			}
		}
		angles = getRoot(commonEdge.getCorner(0), commonEdge.getCorner(1), commonEdge.getCount());
		if (angles!=null && angles[0]< angleLimit) {
			addToHeap(angles, commonEdge);
		}
		
		Tri commonTri = t0.getTri(vns0); //To be deleted in complex
		if (commonTri.getAlph()==1) {
			alphaTris.remove(commonTri);
		}
		commonTri.setAlive(false);
		if (screenAlpha) commonTri.fromScene(scene);
		
		// identify three shared vertices
		Vertex[] vs = new Vertex[3];
		int k = 0;
		for (int i = 0; i < 4; i++) if (t0.getCorner(i) != vns0) vs[k++] = t0.getCorner(i);
		
		Tet[] nt = new Tet[3];
		Tet nt0 = new Tet(vns0, vns1, vs[0], vs[1]);
		nt0.setEdge(commonEdge);
		nt[0] = nt0;
		Tet nt1 = new Tet(vns0, vns1, vs[1], vs[2]);
		nt1.setEdge(commonEdge);
		nt[1] = nt1;
		Tet nt2 = new Tet(vns0, vns1, vs[2], vs[0]);
		nt2.setEdge(commonEdge);
		nt[2] = nt2;
		
		// New triangles
		for (int i = 0 ; i<3 ; i++) {
			Tet t = nt[i];
			Tri tri = new Tri(vns0, vns1, vs[i]);
//			tri.setTet(t);
			t.setTri(tri, t.indexOf(vs[(i+1)%3]));
			Tet n = nt[(i+2)%3];
//			tri.setTet(n);
			n.setTri(tri, n.indexOf(vs[(i+2)%3]));
			if (tri.isAlpha(alphaVal)) {
				tri.setAlph(1);
				if (isGabriel(tri)) {
					alphaTris.add(tri);
					if (!tri.isBig() && screenAlpha) tri.toScene(scene, new Color(0,0,200, 100));
				}
			} else {
				tri.setAlph(0);
			}
			angles = getRoot(tri.getCorner(0), tri.getCorner(1), tri.getCorner(2), tri.getCount());
			if ((angles != null) && (angles[0] < angleLimit)) {
				addToHeap(angles, tri);
			}
		}
		for (Edge e : t0.getEdges()) {
			if (e.equals(new Edge(vns0, vs[0]))) {
				nt0.setEdge(e);
				nt2.setEdge(e);
			} else if (e.equals(new Edge(vns0, vs[1]))) {
				nt0.setEdge(e);
				nt1.setEdge(e);
			} else if (e.equals(new Edge(vns0, vs[2]))) {
				nt2.setEdge(e);
				nt1.setEdge(e);
			} else if (e.equals(new Edge(vs[0], vs[1]))) {
				nt0.setEdge(e);
			} else if (e.equals(new Edge(vs[1], vs[2]))) {
				nt1.setEdge(e);
			} else if (e.equals(new Edge(vs[0], vs[2]))) {
				nt2.setEdge(e);
			}
		}
		for (Edge e : t1.getEdges()) {
			if (e.equals(new Edge(vns1, vs[0]))) {
				nt0.setEdge(e);
				nt2.setEdge(e);
			} else if (e.equals(new Edge(vns1, vs[1]))) {
				nt0.setEdge(e);
				nt1.setEdge(e);
			} else if (e.equals(new Edge(vns1, vs[2]))) {
				nt2.setEdge(e);
				nt1.setEdge(e);
			}
		}
		
//		commonEdge.setTet(nt0);
		nt0.setEdge(commonEdge);
		nt1.setEdge(commonEdge);
		nt2.setEdge(commonEdge);
		// Events for new tetrahedra
		for (Tet t : nt) {
			if (t.isAlpha(alphaVal)) { 
				t.setAlph(1);
				alphaTets.add(t);
			} else t.setAlph(0);
			angles = getRoot(t.getCorner(0), t.getCorner(1), t.getCorner(2), t.getCorner(3), t.getCount());
			if (angles!=null && angles[0]< angleLimit) {
				addToHeap(angles, t);
			}
		}
		alphaTets.remove(t0);
		alphaTets.remove(t1);
		
		
		// vertices get new tetrahedron pointers
		vns0.setTet(nt0);
		vns1.setTet(nt0);
		vs[0].setTet(nt0);
		vs[1].setTet(nt0);
		vs[2].setTet(nt1);
		
		// nt0 contain vs[0] and vs[1]
		if (!nt1.hasVertex(vs[0])) nt0.neighbors[nt0.indexOf(vs[0])] = nt1; else nt0.neighbors[nt0.indexOf(vs[0])] = nt2;
		if (!nt1.hasVertex(vs[1])) nt0.neighbors[nt0.indexOf(vs[1])] = nt1; else nt0.neighbors[nt0.indexOf(vs[1])] = nt2;
		Tet t = t1.neighbors[t1.indexOf(vs[2])];
		Tri tri = t1.getTri(t1.indexOf(vs[2]));
//		tri.setTet(nt0);
		nt0.neighbors[nt0.indexOf(vns0)] = t;
		nt0.setTri(tri, nt0.indexOf(vns0));
		if (t != null) t.neighbors[t.apex(t1)] = nt0;
		t = t0.neighbors[t0.indexOf(vs[2])];
		tri = t0.getTri(t0.indexOf(vs[2]));
//		tri.setTet(nt0);
		nt0.setTri(tri, nt0.indexOf(vns1));
		nt0.neighbors[nt0.indexOf(vns1)] = t;
		if (t != null) t.neighbors[t.apex(t0)] = nt0;
		
		// nt1 contains vs[1] and vs[2]
		if (!nt0.hasVertex(vs[1])) nt1.neighbors[nt1.indexOf(vs[1])] = nt0; else nt1.neighbors[nt1.indexOf(vs[1])] = nt2;
		if (!nt0.hasVertex(vs[2])) nt1.neighbors[nt1.indexOf(vs[2])] = nt0; else nt1.neighbors[nt1.indexOf(vs[2])] = nt2;
		t = t1.neighbors[t1.indexOf(vs[0])];
		tri = t1.getTri(t1.indexOf(vs[0]));
//		tri.setTet(nt1);
		nt1.setTri(tri, nt1.indexOf(vns0));
		nt1.neighbors[nt1.indexOf(vns0)] = t;
		if (t != null) t.neighbors[t.apex(t1)] = nt1;
		t = t0.neighbors[t0.indexOf(vs[0])];
		tri = t0.getTri(t0.indexOf(vs[0]));
//		tri.setTet(nt1);
		nt1.setTri(tri, nt1.indexOf(vns1));
		nt1.neighbors[nt1.indexOf(vns1)] = t;
		if (t != null) t.neighbors[t.apex(t0)] = nt1;
		
		// nt2 contains vs[2] and vs[0]
		if (!nt0.hasVertex(vs[2])) nt2.neighbors[nt2.indexOf(vs[2])] = nt0; else nt2.neighbors[nt2.indexOf(vs[2])] = nt1;
		if (!nt0.hasVertex(vs[0])) nt2.neighbors[nt2.indexOf(vs[0])] = nt0; else nt2.neighbors[nt2.indexOf(vs[0])] = nt1;
		t = t1.neighbors[t1.indexOf(vs[1])];
		tri = t1.getTri(t1.indexOf(vs[1]));
//		tri.setTet(nt2);
		nt2.setTri(tri, nt2.indexOf(vns0));
		nt2.neighbors[nt2.indexOf(vns0)] = t;
		if (t != null) t.neighbors[t.apex(t1)] = nt2;
		t = t0.neighbors[t0.indexOf(vs[1])];
		tri = t0.getTri(t0.indexOf(vs[1]));
//		tri.setTet(nt2);
		nt2.setTri(tri, nt2.indexOf(vns1));
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
		if (testingScreen) {
			for (int i = 0; i < 3; i++) newTets[i].toSceneEdges(scene, Color.black, 0.001, 0.0001);
		}
		tets.add(nt0);
		tets.add(nt1);
		tets.add(nt2);
		//Flip event
		angles = getRoot(nt0, nt1);
		if (angles!=null && angles[0] < angleLimit) {
			addToHeap(angles, nt0, nt1);
			if (testingPrint) System.out.println(nt0 + " " + nt1 + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
		}
		angles = getRoot(nt0, nt2);
		if (angles!=null && angles[0] < angleLimit) {
			addToHeap(angles, nt0, nt2);
			if (testingPrint) System.out.println(nt0 + " " + nt2 + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
		}
		angles = getRoot(nt1, nt2);
		if (angles!=null && angles[0] < angleLimit) {
			addToHeap(angles, nt1, nt2);
			if (testingPrint) System.out.println(nt1 + " " + nt2 + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
		}
				
		// check for new events for new tetrahedra neighbours
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
		flipTime += (System.nanoTime() - start)/1000000.0;
/*		boolean isDT = isDelaunay();
		if (!isDT) {
			throw new RuntimeException("Flip failed between "+t0+" and "+t1);
		}*/
		
//		assert(nt0.getEdges().size()==6) : "flip23, nt0";
//		assert(nt1.getEdges().size()==6) : "flip23, nt1";
//		assert(nt2.getEdges().size()==6) : "flip23, nt2";
		assert(!tets.contains(t0));
		assert(!tets.contains(t1));
		
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
		long start = System.nanoTime();
		nrFlips++;
		//Locate three non-shared vertices
		Vertex[] vns = new Vertex[3];
		vns[0] = t0.getCorner(t0.apex(t1));
		vns[1] = t1.getCorner(t1.apex(t2));
		vns[2] = t2.getCorner(t2.apex(t0));
		
		//kill common edge
		boolean foundEdge = false;
		for (Edge e : t0.getEdges()) {
			if (t1.hasEdge(e) && t2.hasEdge(e)) {
				foundEdge = true;
				e.setAlive(false);
				if (e.getAlph()) { 
					e.setAlph(false); 
					alphaEdges.remove(e);
					if (screenAlpha) e.fromSceneEdge(scene);
				}
				break;
			}
		}
		if (!foundEdge)	{
			System.out.println("Did not find common edge of edges ");
			for (Edge e : t0.getEdges()) {
				System.out.println("   "+e.toString());
			}
			System.out.println();
			for (Edge e : t1.getEdges()) {
				System.out.println("   "+e.toString());
			}
			System.out.println();
			for (Edge e : t2.getEdges()) {
				System.out.println("   "+e.toString());
			}
		}
		//remove common triangles
		Tri tri0 = t0.getTri(vns[0]);
		if (tri0.getAlph()==1) {
			alphaTris.remove(tri0);
		}
		tri0.setAlive(false);
		if (screenAlpha) tri0.fromScene(scene);
		Tri tri1 = t1.getTri(vns[1]);
		if (tri1.getAlph()==1) {
			alphaTris.remove(tri1);
		}
		tri1.setAlive(false);
		if (screenAlpha) tri1.fromScene(scene);
		Tri tri2 = t2.getTri(vns[2]);
		if (tri2.getAlph()==1) {
			alphaTris.remove(tri2);
		}
		tri2.setAlive(false);
		if (screenAlpha) tri2.fromScene(scene);

		//Locate two shared vertices
		Vertex[] vs = new Vertex[2];
		int nrSharedFound = 0;
		for (int i = 0; (i < 4) && (nrSharedFound < 2); i++) {
			Vertex v = t0.getCorner(i);
			if (t1.hasVertex(v) && t2.hasVertex(v)) vs[nrSharedFound++] = v;
		}
		
		Tet[] nt = new Tet[2];
		Tet nt0 = new Tet(vns[0], vns[1], vns[2], vs[0]);
		nt[0] = nt0;
		Tet nt1 = new Tet(vns[0], vns[1], vns[2], vs[1]);
		nt[1] = nt1;
		
		// New common triangle
		Tri tri = new Tri(vns[0], vns[1], vns[2]);
//		tri.setTet(nt0);
//		tri.setTet(nt1);
		nt0.setTri(tri, nt0.indexOf(vs[0]));
		nt1.setTri(tri, nt1.indexOf(vs[1]));
		if (tri.isAlpha(alphaVal)) {
			tri.setAlph(1);
			if (isGabriel(tri) || t0.getAlph()==1) {
				alphaTris.add(tri);
				if (!tri.isBig() && screenAlpha) tri.toScene(scene, new Color(0,0,200, 100));
			}
		} else tri.setAlph(0);
		angles = getRoot(tri.getCorner(0), tri.getCorner(1), tri.getCorner(2), tri.getCount());
		if ((angles != null) && (angles[0] < angleLimit)) {
			addToHeap(angles, tri);
		}
		
		for (Edge e : t0.getEdges()) {
			if (e.equals(new Edge(vs[0], vns[0]))) {
				nt0.setEdge(e);
			} else if (e.equals(new Edge(vs[0], vns[1]))) {
				nt0.setEdge(e);
			} else if (e.equals(new Edge(vs[1], vns[0]))) {
				nt1.setEdge(e);
			} else if (e.equals(new Edge(vs[1], vns[1]))) {
				nt1.setEdge(e);
			} else if (e.equals(new Edge(vns[1], vns[0]))) {
				nt1.setEdge(e);
				nt0.setEdge(e);
			}
		}
		for (Edge e : t1.getEdges()) {
			if (e.equals(new Edge(vs[0], vns[2]))) {
				nt0.setEdge(e);
			} else if (e.equals(new Edge(vs[0], vns[1]))) {
				nt0.setEdge(e);
			} else if (e.equals(new Edge(vs[1], vns[2]))) {
				nt1.setEdge(e);
			} else if (e.equals(new Edge(vs[1], vns[1]))) {
				nt1.setEdge(e);
			} else if (e.equals(new Edge(vns[1], vns[2]))) {
				nt1.setEdge(e);
				nt0.setEdge(e);
			}
		}
		for (Edge e : t2.getEdges()) {
			if (e.equals(new Edge(vs[0], vns[2]))) {
				nt0.setEdge(e);
			} else if (e.equals(new Edge(vs[0], vns[0]))) {
				nt0.setEdge(e);
			} else if (e.equals(new Edge(vs[1], vns[2]))) {
				nt1.setEdge(e);
			} else if (e.equals(new Edge(vs[1], vns[0]))) {
				nt1.setEdge(e);
			} else if (e.equals(new Edge(vns[0], vns[2]))) {
				nt1.setEdge(e);
				nt0.setEdge(e);
			}
		}
		
		// Events for new tetrahdra
		for (Tet t : nt) {
			if (t.isAlpha(alphaVal)) {
				t.setAlph(1);
				alphaTets.add(t);
			} else  t.setAlph(0);
			angles = getRoot(t.getCorner(0), t.getCorner(1), t.getCorner(2), t.getCorner(3), t.getCount());
			if (angles!=null && angles[0]< angleLimit) {
				addToHeap(angles, t);
			}
		}
		
		alphaTets.remove(t0); alphaTets.remove(t1); alphaTets.remove(t2);

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
			Tri triangle;
			Tet n0 = tn0.neighbors[i];
			triangle = tn0.getTri(i);
			int idx;
			if (n0 == null) { 
				if (tn0.getCorner(i)==vs[1]) {
					nt0.setTri(triangle, nt0.indexOf_slow(vns[0]));
				} else nt1.setTri(triangle, nt1.indexOf_slow(vns[0]));
			} else if((n0 != tn1) && (n0 != tn2)) {
				idx = n0.apex(tn0);
				if (n0.hasVertex(vs[0])) {
					n0.neighbors[idx] = nt0;
					triangle = n0.getTri(idx);
//					triangle.setTet(nt0);
					nt0.setTri(triangle, nt0.indexOf_slow(vns[0]));
					nt0.neighbors[nt0.indexOf_slow(vns[0])] = n0;
				}
				else {
					n0.neighbors[idx] = nt1;
					triangle = n0.getTri(idx);
//					triangle.setTet(nt1);
					nt1.setTri(triangle, nt1.indexOf_slow(vns[0]));
					nt1.neighbors[nt1.indexOf_slow(vns[0])] = n0;
				}
			}
		}

		for (int i = 0; i < 4; i++) {
			Tri triangle;
			Tet n1 = tn1.neighbors[i];
			triangle = tn1.getTri(i);
			int idx;
			if (n1 == null) { 
				if (tn1.getCorner(i)==vs[1]) {
					nt0.setTri(triangle, nt0.indexOf_slow(vns[1]));
				} else nt1.setTri(triangle, nt1.indexOf_slow(vns[1]));
			} else if((n1 != tn0) && (n1 != tn2)) {
				idx = n1.apex(tn1);
				if (n1.hasVertex(vs[0])) {
					n1.neighbors[idx] = nt0;
					triangle = n1.getTri(idx);
//					triangle.setTet(nt0);
					nt0.setTri(triangle, nt0.indexOf_slow(vns[1]));
					nt0.neighbors[nt0.indexOf_slow(vns[1])] = n1;
				}
				else {
					n1.neighbors[idx] = nt1;
					triangle = n1.getTri(idx);
//					triangle.setTet(nt1);
					nt1.setTri(triangle, nt1.indexOf_slow(vns[1]));
					nt1.neighbors[nt1.indexOf_slow(vns[1])] = n1;
				}
			}
		}

		for (int i = 0; i < 4; i++) {
			Tri triangle;
			Tet n2 = tn2.neighbors[i];
			triangle = tn2.getTri(i);
			int idx;
			if (n2 == null) { 
				if (tn2.getCorner(i)==vs[1]) {
					nt0.setTri(triangle, nt0.indexOf_slow(vns[2]));
				} else nt1.setTri(triangle, nt1.indexOf_slow(vns[2]));
			} else if((n2 != tn0) && (n2 != tn1)) {
				idx = n2.apex(tn2);
				if (n2.hasVertex(vs[0])) {
					n2.neighbors[idx] = nt0;
					triangle = n2.getTri(idx);
//					triangle.setTet(nt0);
					nt0.setTri(triangle, nt0.indexOf_slow(vns[2]));
					nt0.neighbors[nt0.indexOf_slow(vns[2])] = n2;
				}
				else {
					n2.neighbors[idx] = nt1;
					triangle = n2.getTri(idx);
//					triangle.setTet(nt1);
					nt1.setTri(triangle, nt1.indexOf_slow(vns[2]));
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
		//Flip event
		angles = getRoot(nt0, nt1);
		if (angles!=null && angles[0] < angleLimit) {
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
					if ((nt1.toString().equals("[408, 410, 491, 972]")) && nti.toString().equals("[408, 415, 491, 972]")) {
						System.out.println("Tets checked and results in: "+angles[0]);
					}
					addToHeap(angles, nt1, nti);
					if (testingPrint) System.out.println(nt1 + " " + nti + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
				}
			}
		}
		flipTime += (System.nanoTime() - start)/1000000.0;
/*		boolean isDT = isDelaunay();
		if (!isDT) {
			throw new RuntimeException("Flip failed between "+t0+" and "+t1+" and "+t2);
		}*/
		
//		assert(nt0.getEdges().size()==6) : "flip32, nt0";
//		assert(nt1.getEdges().size()==6) : "flip32, nt1";
		assert(!tets.contains(t0));
		assert(!tets.contains(t1));
		assert(!tets.contains(t2));
		
		return newTets;
	}
	
	/** Returns TRUE id the tetrahedron is Delaunay (circumscribing sphere contains no point) */
	public boolean isDelaunay() {
		boolean cont = true;
		for (Tet t : tets) {
			if (t.isBig()) continue;
			t.setCircumSphere();
//			List<Vertex> noBigpoints = new ArrayList<Vertex>(vertices.size());
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
	
	private void showAlpha() {
		for (Edge ae : alphaEdges) {
			if (ae.isBig()) continue;
			Color clr = Color.BLACK;
//			if (ae.getCount()==3) clr = Color.RED;
			ae.toSceneEdge(scene, clr, 0.005);
		}
		for (Tri at : alphaTris) {
			if (at.isBig()) continue;
			at.toScene(scene, new Color(0,0,200, 100));
		}
	}
	
	public void rotateTo(double rotateTo) {
		if (rotateTo<angleTotal) {
			throw new RuntimeException("Cannot rotate clockwise");
		}
		if (rotateTo>angleLimit) {
			rotateTo(angleLimit);
//			System.out.println("Initial coords for vertex indexed at 10 = "+vertices.get(10).getCoord(0)+" , "+vertices.get(10).getCoord(1)+" , "+vertices.get(10).getCoord(2));
//			initializeRotation(null, null);
			initializeRotation(null, 0,1);
			rotateTo(rotateTo-angleLimit);
			return;
			//throw new RuntimeException("Angle is bigger than limit=360");
		}
		HeapItem heapItem;
		Tet t, nt, tt;
		Edge edg;
		Tri tri;
		int nrFlips = 0;
		double rotAngle;
		if (screenAlpha) {
			showAlpha();
			for (int i = 0; i < vertices.size(); i++) {
				Vertex v = vertices.get(i);
//				scene.addShape(new TextShape(Integer.toString(i), v.returnAsPoint(), 0.3));
				Color clr = Color.BLACK;
				if (v.getType()==VertexType.R) clr = Color.RED;
				v.toScene(scene, clr, 0.1);
			}
		}
		clashFirst = alphaEdges.isEmpty();
		clash = clashFirst;
		while (!heap.isEmpty() && angleTotal < rotateTo) {
			HeapItem angleCheck = (HeapItem) heap.peek();
			if (angleCheck.angles[0]>rotateTo) break;
			heapItem = (HeapItem) heap.extract();
			t = heapItem.getT();
			tri = heapItem.getTri();
			if (t!=null) {
				nt = heapItem.getNT();
				if (nt!=null) {
					if (testingPrint) {
						System.out.print(t.toString());  if (!t.isAlive())  System.out.print("");
						System.out.print(nt.toString()); if (!nt.isAlive()) System.out.print("");
					}
					angles = heapItem.getAngles();
					rotAngle = angles[0];
					if (rotAngle > rotateTo) break;
					if (t.isAlive() && nt.isAlive()) {
						
						animate(scene, rotAngle-angleTotal);
						angleTotal = rotAngle;
//						System.out.println("AngleTotal = "+angleTotal);
						if (!t.isConvex(nt)) {
							tt = KineticToolbox.getThirdTet(t, nt);
							if (tt != null) {
//								System.out.println("Flip event with "+t.toString()+" "+(t.getAlph()==1)+" and "+nt.toString()+" "+(nt.getAlph()==1)+" and "+tt.toString()+" "+(tt.getAlph()==1));
								if (testingPrint) {
									System.out.print(tt.toString()); 
									if (!tt.isAlive())  System.out.print("�");
									System.out.print(" rotated to angle = " + Functions.toDeg(heapItem.getAngles()[0]));						
									System.out.println();
								}
								flip32(t, nt, tt);
								
								if (testingPrint) System.out.println(++nrFlips + ". flip.");
							}
							else { 
								System.out.println(" not convex but tt = null");
/*								System.out.println("t = "+t+" count = "+t.getCount()+" nt = "+nt+" count = "+nt.getCount());
								System.out.println("radius t = "+t.getCircumSphereRadius()+" radius nt = "+nt.getCircumSphereRadius());
								J3DScene scene = J3DScene.createJ3DSceneInFrame();
								t.toSceneFaces(scene, new Color(200,0,0,100));
								nt.toSceneFaces(scene, new Color(0,0,200,100));
								System.out.println("t's neighbours :");
								for (int i = 0 ; i<4 ; i++) {
									System.out.println(i+" = "+t.getNeighbour(i));
								}
								System.out.println("nt's neighbours :");
								for (int i = 0 ; i<4 ; i++) {
									System.out.println(i+" = "+nt.getNeighbour(i));
								}
								vertices.get(972).toScene(scene, 0.1, new Color(0,0,0,255));
								System.out.println("point id : "+vertices.get(972).toString()+" type = "+vertices.get(972).getType());
								throw new RuntimeException("WuaaaaH");*/
//								break;
							}
						}
						else {
//							System.out.println("Flip event with "+t.toString()+" "+(t.getAlph()==1)+" and "+nt.toString()+" "+(nt.getAlph()==1));
							if (testingPrint) {
								System.out.print(" rotated to angle = " + Functions.toDeg(heapItem.getAngles()[0]));
								System.out.println();
							}
							if (testingScreen) {
								t.fromSceneFaces(scene);
								nt.fromSceneFaces(scene);
							}
							flip23(t, nt);
							if (testingPrint) System.out.println(++nrFlips + ". flip.");
						}
						if (clash!=alphaEdges.isEmpty()) {
							clash = alphaEdges.isEmpty();
							clashes.add(Math.toDegrees(angleTotal));
						}
						if (clash!=alphaEdges.isEmpty()) {
							clash = alphaEdges.isEmpty();
							clashes.add(Math.toDegrees(angleTotal));
						}
//						isAlpha();
//						isDelaunay();
					}
				}
				else { //radiusevent tetrahedron
//					System.out.println("Radius event : tet "+t.toString());
					angles = heapItem.getAngles();
					rotAngle = angles[0];
					if (t.isAlive() && rotAngle < angleLimit) {
						if (t.nrRotating()==1) { // one corner is rotating
							animate(scene, rotAngle-angleTotal);
							angleTotal = rotAngle;
//							System.out.println("AngleTotal = "+angleTotal);
							int alph = t.getAlph();
							if (alph == 0) {
								//if (a0 <= Constants.EPSILON || a1 <= Constants.EPSILON || a2 <= Constants.EPSILON) {
								t.setAlph(1);
								alphaTets.add(t);
							}
							else if (alph == 1) {
								//if (a0 <= Constants.EPSILON || a1 <= Constants.EPSILON || a2 <= Constants.EPSILON) {
								t.setAlph(0);
								alphaTets.remove(t);
							}
							Double[] tmp = new Double[4];
							int i = -1;
							for (Double angle : angles){
								if (i==-1) { i += 1; continue; }
								if (angle==null) break;
								tmp[i] = angle;
								i += 1;
							}
							if (i>0) {
								Double[] newangles = new Double[i];
								for (int j = 0 ; j<i; j++) {
									newangles[j] = tmp[j];
								}
								addToHeap(newangles, t);
							}
						} else { // two corners rotating	
							animate(scene, rotAngle-angleTotal);
							angleTotal = rotAngle;
//								System.out.println("AngleTotal = "+angleTotal);
							if (t.getAlph()==1) {
								t.setAlph(0);
								alphaTets.remove(t);
							}
							else {
								t.setAlph(1);
								alphaTets.add(t);
							}
							Double[] tmp = new Double[4];
							int i = -1;
							for (Double angle : angles){
								if (i==-1) { i += 1; continue; }
								if (angle==null) break;
								tmp[i] = angle;
								i += 1;
							}
							if (i>0) {
								Double[] newangles = new Double[i];
								for (int j = 0 ; j<i; j++) {
									newangles[j] = tmp[j];
								}
								addToHeap(newangles, t);
							}
						}
						if (clash!=alphaEdges.isEmpty()) {
							clash = alphaEdges.isEmpty();
							clashes.add(Math.toDegrees(angleTotal));
						}
//						isAlpha();
					}
				}
			} else if (tri!=null) {
//				System.out.println("Radius event : tri "+tri.toString());
				angles = heapItem.getAngles();
				rotAngle = angles[0];
				if (tri.isAlive() && rotAngle < angleLimit) {
					animate(scene, rotAngle-angleTotal);
					angleTotal = rotAngle;
//					System.out.println("AngleTotal = "+angleTotal);
					int alph = tri.getAlph();
					if (alph == 1) {
						alphaTris.remove(tri);
						if (screenAlpha) tri.fromScene(scene);
						tri.setAlph(0);
					} else if(alph == 0) {
						tri.setAlph(1);
						if (isGabriel(tri)) {
							alphaTris.add(tri);
							if (!tri.isBig() && screenAlpha) tri.toScene(scene, new Color(0,0,200, 100));
						}
					}
					Double[] tmp = new Double[4];
					int i = -1;
					for (Double angle : angles){
						if (i==-1) { i += 1; continue; }
						if (angle==null) break;
						tmp[i] = angle;
						i += 1;
					}
					if (i>0) {
						Double[] newangles = new Double[i];
						for (int j = 0 ; j<i; j++) {
							newangles[j] = tmp[j];
						}
						addToHeap(newangles, tri);
					}
					if (clash!=alphaEdges.isEmpty()) {
						clash = alphaEdges.isEmpty();
						clashes.add(Math.toDegrees(angleTotal));
					}
//					isAlpha();
				}
			} else {
				edg = heapItem.getEdg();
//				System.out.println("Radius event : edg "+edg.toString());
				angles = heapItem.getAngles();
				rotAngle = angles[0];
				if (edg.isAlive() && edg!=null && rotAngle < angleLimit) {
					animate(scene, rotAngle-angleTotal);
					angleTotal = rotAngle;
//					System.out.println("AngleTotal = "+angleTotal);
					if (edg.getAlph()) {
						edg.setAlph(false);
						if (screenAlpha) edg.fromSceneEdge(scene);
						alphaEdges.remove(edg);
					}
					else {
						edg.setAlph(true);
						if (isGabriel(edg)) {
							if (!edg.isBig() && screenAlpha) edg.toSceneEdge(scene, Color.BLACK, 0.005);
							alphaEdges.add(edg);
						}
					}
					Double[] newangles = new Double[4];
					if (angles[1]!=null) {
						newangles[0] = angles[1];
						addToHeap(newangles, edg);
					}
					if (clash!=alphaEdges.isEmpty()) {
						clash = alphaEdges.isEmpty();
						clashes.add(Math.toDegrees(angleTotal));
					}
//					isAlpha();
				}
			}
			
			
		}
		
		animate(scene, rotateTo-angleTotal);
		angleTotal = rotateTo;
//		System.out.println("AlphaComplex = "+alphaTets.toString()+" and "+alphaTris.toString()+" and "+alphaEdges.toString());
	}
	
	/* specify what is rotating - 3 methods to do it */
	public void setRotVertices() { for (Vertex v : vertices) v.setType(Vertex.VertexType.S); }
	public void setRotVertices(int a) { setRotVertices(a, a); }
	public void setRotVertices(int a, int b) { 
		for (int i = a; i <= b; i++) rotIndx.add(i); 
		for (Vertex v : vertices) v.setType(Vertex.VertexType.S); 
		for (int i : rotIndx) {
			Vertex v = vertices.get(i);
			v.setType(Vertex.VertexType.R);
		}
	}
	
	public void setRotVertices(List<Integer> rotList) { 
		rotIndx.clear();
		for (Integer  i : rotList) {
			rotIndx.add(i+4);
		}
		for (Vertex v : vertices) v.setType(Vertex.VertexType.S); 
		for (int i : rotIndx) vertices.get(i).setType(Vertex.VertexType.R);
	}

	/** returns the scene of this kDT */
	public J3DScene getScene() { return scene; }
	
	/** returns the list of tetrahedra in this kinetic Dalaunay tessellation */
	public Set<Tet> getTetrahedra(){ return new HashSet<Tet>(tets); }
	
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
//		Color blue_tr = new Color(0,0,255,50);
//		Color green_tr = Color.green;
//		Color green_tr = new Color(0,128,0,50);
//		Color red_tr = new Color(255,0,0,50);
		Vector tr = new Vector (0.02, 0.02, 0.02);
//		int nrTetrahedra = 0;
		for (Tet tet : tets) { 
//			nrTetrahedra++;
//			System.out.println(nrTetrahedra);
//			tet.getCircumSphere();
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
//		for (Vertex v : vertices) {
//			if (v.getType() == VertexType.S) v.toScene(scene, 0.01, Color.BLACK, 5);
//			else {
//				v.toScene(scene, 0.01, Color.RED, 5);
////				Circle c = new Circle(new Point(0, 0, v.z()), v, getRotationAxis());
//			}
//			if (this.vertices.size() < 20) scene.addText(String.valueOf(v.getId()), v);
//		}
	}
	
	//Daisy:
	public void toScene(Tet t) {
		scene.addShape(new Tetrahedron(t.getCorner(0), t.getCorner(1), t.getCorner(2), t.getCorner(3)), new Color(200, 0, 0, 100));
	}
	public void toScene(Tri t) {
		scene.addShape(new Triangle(t.getCorner(0), t.getCorner(1), t.getCorner(2)), new Color(0, 2000, 0, 100));
	}
	public void toScene(Edge e) {
		scene.addShape(new LSS(e.getCorner(0), e.getCorner(1), 0.01), new Color(0, 0, 200, 100));
	}
	
	public static void main(String[] args){
		System.out.println(java.lang.Runtime.getRuntime().maxMemory());
		long start = System.nanoTime();
		List<Point> points = new PDBFile("/home/daisy/Desktop/KinAlpha/VoidRemover/Experiments/1KT3/start.pdb", false).getAtomCoords();
		KineticACD kDT = new KineticACD(2.5, points);
//		vertices.get(28+4), vertices.get(29+4));
		List<Integer> rotIndices = new ArrayList<Integer>(Arrays.asList(30, 31));
		kDT.initializeRotation(rotIndices, 28, 29);
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
		
//		System.out.printf(" time in miliseconds %.2f\n", (System.nanoTime() - start)/1000000.0);
//		start = System.nanoTime();
		kDT.rotateTo(Math.toRadians(360));
//		System.out.println("Flips : "+kDT.nrFlips);
		System.out.printf(" time in miliseconds %.2f\n", (System.nanoTime() - start)/1000000.0);
//		start = System.nanoTime();
		//kDT.rotate(Math.toRadians(360));
		System.out.println("Flips : "+kDT.nrFlips);
//		System.out.printf(" time in miliseconds %.2f\n", (System.nanoTime() - start)/1000000.0);
		
		boolean curClash = kDT.clashFirst;
		System.out.println("length of shortest edge is "+kDT.shortestEdge);
		System.out.println("Clash at first ? "+!curClash);
		System.out.println("Alpha edges empty? "+kDT.alphaEdges.isEmpty());
		for (Double clash : kDT.clashes) {
			System.out.println(!curClash?"No clash from "+clash:"Clash from "+clash);
			if (curClash) {
				curClash=false;
			} else curClash=true;
		}
	}

	public ProblemInstanceType getInstanceType() {
		return instanceType;
	}

	public void setInstanceType(ProblemInstanceType instanceType) {
		this.instanceType = instanceType;
	}

	public double getAlpha() {
		return alphaVal;
	}

	public void setAlpha(double alpha) {
		this.alphaVal = alpha;
	}

	/** Returns the set of tetrahedra that are currently part of the alpha complex */
	public Set<Tet> getAlphaTetrahedra(){ return alphaTets; }
	
}
