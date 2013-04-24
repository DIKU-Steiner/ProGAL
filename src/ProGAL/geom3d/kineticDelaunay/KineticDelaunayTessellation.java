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

public class KineticDelaunayTessellation {
	private static enum Direction { CW, CCW }

	private List<Vertex> vertices = new ArrayList<Vertex>();

	private List<Tet> tets = new LinkedList<Tet>();
	private Tet lastTet;
	private Double[] angles = new Double[2];
	private double angleTotal = 0.0;
	
	private Point rotationPoint;   // rotation center
	private Vector rotationAxis; // rotation vector
	private Direction rotDir;
	private List<Integer> rotIndx = new ArrayList<Integer>();     // indicies of vertices that move

	private Heap heap = new Heap(this.vertices.size(), new SortToolHeapItems());
	private boolean testing = true;
	private boolean sphereAnimation = false;
	J3DScene scene  = J3DScene.createJ3DSceneInFrame();

	Triangle commonTriangle = null;
	Triangle commonTriangle02 = null;
	Triangle commonTriangle12 = null;
	Tetrahedron commonTriangleShape = null;
	Tetrahedron commonTriangleShape02 = null;
	Tetrahedron commonTriangleShape12 = null;
	LineSegment apexSegment = null;
	LSS apexLSS = null;

	double alpha;
	
	private class HeapItem {
		private Double[] angles;
		private Tet t;
		private Tet nt;
		
		private HeapItem(Double[] angles, Tet t, Tet nt) {
			this.angles = angles;
			this.t = t;
			this.nt = nt;
		}
				
		private Double getAngle() { return angles[0];} 
		private Tet getT()  { return t; }
		private Tet getNT() { return nt; }
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

	
	public KineticDelaunayTessellation(List<Point> points) {
		Tetrahedron bigT = Tetrahedron.regularTetrahedron();
		bigT.blowUp(5);
		lastTet = new BigTet(bigT, vertices);
		for (Point p: points) insertPoint(p);		
	}

	public Direction getDirection() { return rotDir; }
	public void setDirection(Direction rotDir) { this.rotDir = rotDir; }

	public Point getRotationPoint() { return rotationPoint; }
	/* Sets the rotation point and translates all vertices so that the rotation point is in the origo.  */
	public void setRotationPoint(Point rotationPoint) { 
		if (this.rotationPoint != null) for (Vertex v : vertices) v.addThis(this.rotationPoint);
		this.rotationPoint = rotationPoint; 
		for (Vertex v : vertices) v.subtractThis(rotationPoint);	
	}

	public void setRotationAxis(Vector v)  { rotationAxis = v; }
	public Vector getRotationAxis() { return rotationAxis; }
	
	private void addToHeap(Double[] angles, Tet t, Tet nt) {
		heap.insert(new HeapItem(angles, t, nt));
	}

	
	
	public void insertPoint(Point p){
		Vertex v = new Vertex(p);
		System.out.println("Inserting "+v);
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
					for(int i=1;i<4;i++) corners[i-1] = c.corners[(f+i)%4];
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

	private void restoreNeighborhood(List<Tet> newCells){
		for(Tet c1: newCells){
			cellLoop: for(Tet c2: newCells){
				if(c1==c2) break;
				
				//Check if c1 and c2 share a face based on corner vertices
				//Both will contain the last inserted vertex (highest index)
				//i and j run over the arrays, I and J record the location of differences
				int i=0,j=0, I=-1, J=-1;
				while(i<=3&&j<=3){
					if(c1.corners[i]==c2.corners[j]){
						i++;j++;
					}else if(i<3 && c1.corners[i+1]==c2.corners[j]){
						if(I>=0) continue cellLoop;
						I=i;i++;
					}else {
						if(J>=0) continue cellLoop;
						J=j;j++;
					}
				}
				c1.neighbors[I] = c2;
				c2.neighbors[J] = c1;
			}
		}
	}
	

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
	
	public void animate(J3DScene scene, double alpha) {
		int steps = 1;
		double angleStep = alpha/steps;
		for (int k = 0; k < steps; k++) {
			for (int i = 0; i < rotIndx.size(); i++) {
				Vertex v = vertices.get(rotIndx.get(i));
				v.rotation(rotationAxis, angleStep);
			}
			try { Thread.sleep(10); } catch (InterruptedException e) {}
			if (sphereAnimation) animateSpheres();			
			scene.repaint();
		}
	}

	public void animateSpheres() {
		Sphere d;
		for (Tet t : tets) {
			if (!t.isFlat()) {
				Sphere c = t.getCircumSphere();
				if (c != null) {
					d = new Sphere(t.corners[0], t.corners[1], t.corners[2], t.corners[3]);
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
	    
	    double M12 = A.z()*(bb-cc)                   + C.z()*(aa-bb)                   + B.z()*(cc-aa);
	    double M13 = A.y()*(bb-cc)                   + C.y()*(aa-bb)                   + B.y()*(cc-aa);
	    double M14 = A.y()*(B.z()-C.z())             + C.y()*(A.z()-B.z())             + B.y()*(C.z()-A.z());
	    double M15 = A.y()*(B.z()*cc-C.z()*bb)       + C.y()*(A.z()*bb-B.z()*aa)       + B.y()*(C.z()*aa-A.z()*cc);
	    double M23 = A.x()*(bb-cc)                   + C.x()*(aa-bb)                   + B.x()*(cc-aa);
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

	public Double[] getRoot(Tet t, Tet oppT) {
		int count = t.getCount();
		Vertex oppV = oppT.corners[oppT.apex(t)];
		if (oppV.getType() == Vertex.VertexType.R) count = count + 16;
		return getRoot(t, oppV, count);
	}
	
	public Double[] getRoot(Tet t, Vertex oppV, int count) {		
		if ((count == 0) || (count == 31)) return null;
		if (count <= 15) {	
			if (count <= 7) {
				if (count <= 3) {
					if (count == 1) return getRootSSSSR(t.corners[1], t.corners[2], t.corners[3], oppV, t.corners[0], 0);
					if (count == 2) return getRootSSSSR(t.corners[0], t.corners[2], t.corners[3], oppV, t.corners[1], 0);
					return getRootSSSRR(t.corners[2], t.corners[3], oppV, t.corners[0], t.corners[1], 0);
				}
				if (count <= 5) {
					if (count == 4) return getRootSSSSR(t.corners[0], t.corners[1], t.corners[3], oppV, t.corners[2], 0);
					return getRootSSSRR(t.corners[1], t.corners[3], oppV, t.corners[0], t.corners[2], 0);
				}
				if (count == 6) return getRootSSSRR(t.corners[0], t.corners[3], oppV, t.corners[1], t.corners[2], 0);					
				return getRootSSSRR(t.corners[0], t.corners[1], t.corners[2], t.corners[3], oppV, 1);
			}
			if (count <= 11) {
				if (count <= 9) {
					if (count == 8) return getRootSSSSR(t.corners[0], t.corners[1], t.corners[2], oppV, t.corners[3], 0);
					return getRootSSSRR(t.corners[1], t.corners[2], oppV, t.corners[0], t.corners[3], 0);
				}
				if (count == 10) return getRootSSSRR(t.corners[0], t.corners[2], oppV, t.corners[1], t.corners[3], 0);
				return getRootSSSRR(t.corners[0], t.corners[1], t.corners[3], t.corners[2], oppV, 1);
			}
			if (count <= 13) {
				if (count == 12) return getRootSSSRR(t.corners[0], t.corners[1], oppV, t.corners[2], t.corners[3], 0);
				return getRootSSSRR(t.corners[0], t.corners[2], t.corners[3], t.corners[1], oppV, 1);
			}
			if (count == 14) return getRootSSSRR(t.corners[1], t.corners[2], t.corners[3], t.corners[0], oppV, 1);
			return getRootSSSSR(t.corners[0], t.corners[1], t.corners[2], t.corners[3], oppV, 1);
		}
		if (count <= 23) {
			if (count <= 19) {
				if (count <= 17) {
					if (count == 16) return getRootSSSSR(t.corners[0], t.corners[1], t.corners[2], t.corners[3], oppV, 0);
					return getRootSSSRR(t.corners[1], t.corners[2], t.corners[3], t.corners[0], oppV, 0);
				}
				if (count == 18) return getRootSSSRR(t.corners[0], t.corners[2], t.corners[3], t.corners[1], oppV, 0);
				return getRootSSSRR(t.corners[0], t.corners[1], oppV, t.corners[2], t.corners[3], 1);
			}
			if (count <= 21) {
				if (count == 20) return getRootSSSRR(t.corners[0], t.corners[1], t.corners[3], t.corners[2], oppV, 0);
				return getRootSSSRR(t.corners[1], t.corners[3], oppV, t.corners[0], t.corners[2], 0);
			}
			if (count == 22) return getRootSSSRR(t.corners[1], t.corners[2], oppV, t.corners[0], t.corners[3], 1);					
			return getRootSSSSR(t.corners[0], t.corners[1], t.corners[2], oppV, t.corners[3], 1);
		}
		if (count <= 27) {
			if (count <= 25) {
				if (count == 24) return getRootSSSRR(t.corners[0], t.corners[1], t.corners[2], t.corners[3], oppV, 0);
				return getRootSSSRR(t.corners[0], t.corners[3], oppV, t.corners[1], t.corners[2], 1);
			}
			if (count == 26) return getRootSSSRR(t.corners[1], t.corners[3], oppV, t.corners[0], t.corners[2], 1);
			return getRootSSSSR(t.corners[0], t.corners[1], t.corners[3], oppV, t.corners[2], 1);
		}
		if (count <= 29) {
			if (count == 28) return getRootSSSRR(t.corners[2], t.corners[3], oppV, t.corners[0], t.corners[1], 1);
			return getRootSSSSR(t.corners[0], t.corners[2], t.corners[3], oppV, t.corners[1], 1);
		}
		return getRootSSSSR(t.corners[1], t.corners[2], t.corners[3], oppV, t.corners[0], 1);
	}
	
	private void testInfo(Tet t, Color clr, boolean second) {
		if (second) System.out.print(" ");
		System.out.println(t.toString()); 
		t.circumSphere = new Sphere(t.corners[0], t.corners[1], t.corners[2], t.corners[3]);
		scene.addShape(t.circumSphere, new Color(255, 255, 0, 100), 32);
	}

	
	private void initializeRotation() {
		heap.clear();
		// set constants associated with vertices - does not include vertices of big points
		for (Vertex v : vertices) {
			if (testing) System.out.print(v.getId() + ": " + v.toString(2));
			if (v.distanceSquared() < Constants.EPSILON) { //takes care of the special case when the rotation center overlaps with one of the given points
				v.setSquaredPolarRadius(0.0);                            
				v.setPolarRadius(0.0);                        
				v.setPolarAngle(0.0);			
				if (testing) System.out.println(", polar angle: 0.0"); 
				v.setCosAngle(1.0);
				v.setSinAngle(0.0);
				v.setType(Vertex.VertexType.S);
				rotIndx.remove((Integer)v.getId());
			}
			else {
				v.setSquaredPolarRadius(v.distanceSquared());       // remains unchanged when rotating around the same point
				v.setPolarRadius(v.distance());                     // remains unchanged when rotating around the same point
				v.setPolarAngle(v.polarAngleXY());     			
				if (testing) System.out.println(", initial polar angle: " + Functions.toDeg(v.getPolarAngle())); 
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
					oppV = nt.corners[nt.apex(t)];  
					if (t.corners[t.apex(nt)].getId() < oppV.getId()) {
						count = count4;
						if (oppV.getType() == Vertex.VertexType.R) count = count + 16;	
						angles = getRoot(t, oppV, count);
						if (angles != null) {
							addToHeap(angles, t, nt);
							if (testing) System.out.println(t + " " + nt + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
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

		// identify 2 not shared vertices
		int a0 = t0.apex(t1);
		int a1 = t1.apex(t0);
		Vertex[] vns = new Vertex[2];
		vns[0] = t0.corners[a0];
		vns[1] = t1.corners[a1];
		
		// identify three shared vertices
		Vertex[] vs = new Vertex[3];
		int k = 0;
		for (int i = 0; i < 4; i++) if (t0.corners[i] != t0.corners[a0]) vs[k++] = t0.corners[i];
		
		Tet nt0 = new Tet(vns[0], vns[1], vs[0], vs[1]);
		Tet nt1 = new Tet(vns[0], vns[1], vs[1], vs[2]);
		Tet nt2 = new Tet(vns[0], vns[1], vs[2], vs[0]);
		nt0.sortCorners();
		nt1.sortCorners();
		nt2.sortCorners();
		
		// nt0 contain vs[0]Êand vs[1]
		
		if (!nt1.hasVertex(vs[0])) nt0.neighbors[nt0.indexOf(vs[0])] = nt1; else nt0.neighbors[nt0.indexOf(vs[0])] = nt2;
		if (!nt1.hasVertex(vs[1])) nt0.neighbors[nt0.indexOf(vs[1])] = nt1; else nt0.neighbors[nt0.indexOf(vs[1])] = nt2;
		Tet t = t1.neighbors[t1.indexOf(vs[2])];
		nt0.neighbors[nt0.indexOf(vns[0])] = t;
		if (t != null) t.neighbors[t.apex(t1)] = nt0;
		t = t0.neighbors[t0.indexOf(vs[2])];
		nt0.neighbors[nt0.indexOf(vns[1])] = t;
		if (t != null) t.neighbors[t.apex(t0)] = nt0;
		
		// nt1 contains vs[1] and vs[2]
		
		if (!nt0.hasVertex(vs[1])) nt1.neighbors[nt1.indexOf(vs[1])] = nt0; else nt1.neighbors[nt1.indexOf(vs[1])] = nt2;
		if (!nt0.hasVertex(vs[2])) nt1.neighbors[nt1.indexOf(vs[2])] = nt0; else nt1.neighbors[nt1.indexOf(vs[2])] = nt2;
		t = t1.neighbors[t1.indexOf(vs[0])];
		nt1.neighbors[nt1.indexOf(vns[0])] = t;
		if (t != null) t.neighbors[t.apex(t1)] = nt1;
		t = t0.neighbors[t0.indexOf(vs[0])];
		nt1.neighbors[nt1.indexOf(vns[1])] = t;
		if (t != null) t.neighbors[t.apex(t0)] = nt1;
		
		// nt2 contains vs[2] and vs[0]
		
		if (!nt0.hasVertex(vs[2])) nt2.neighbors[nt2.indexOf(vs[2])] = nt0; else nt2.neighbors[nt2.indexOf(vs[2])] = nt1;
		if (!nt0.hasVertex(vs[0])) nt2.neighbors[nt2.indexOf(vs[0])] = nt0; else nt2.neighbors[nt2.indexOf(vs[0])] = nt1;
		t = t1.neighbors[t1.indexOf(vs[1])];
		nt2.neighbors[nt2.indexOf(vns[0])] = t;
		if (t != null) t.neighbors[t.apex(t1)] = nt2;
		t = t0.neighbors[t0.indexOf(vs[1])];
		nt2.neighbors[nt2.indexOf(vns[1])] = t;
		if (t != null) t.neighbors[t.apex(t0)] = nt2;
		
		tets.remove(t0);
		t0.setAlive(false);
		if (testing) t0.fromSceneEdges(scene);
		tets.remove(t1);
		t1.setAlive(false);		
		if (testing) t1.fromSceneEdges(scene);
		Tet[] newTets = new Tet[3];
		newTets[0] = nt0;
		newTets[1] = nt1;
		newTets[2] = nt2;
		if (testing) {
			newTets[0].toSceneEdges(scene, Color.black, 0.005, 0.0001);
			newTets[1].toSceneEdges(scene, Color.black, 0.005, 0.0001);
			newTets[2].toSceneEdges(scene, Color.black, 0.005, 0.0001);
		}
		tets.add(nt0);
		tets.add(nt1);
		tets.add(nt2);
		angles = getRoot(nt0, nt1);
		if (angles != null) {
			addToHeap(angles, nt0, nt1);
			if (testing) System.out.println(nt0 + " " + nt1 + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
		}
		angles = getRoot(nt0, nt2);
		if (angles != null) {
			addToHeap(angles, nt0, nt2);
			if (testing) System.out.println(nt0 + " " + nt2 + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
		}
		angles = getRoot(nt1, nt2);
		if (angles != null) {
			addToHeap(angles, nt1, nt2);
			if (testing) System.out.println(nt1 + " " + nt2 + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
		}
		
		// check for new events
		Tet nti;
		for (int i = 0; i < 4; i++) {
			nti = nt0.neighbors[i];
			if ((nti != null) && (nti != nt1) && (nti != nt2)) {
				angles = getRoot(nt0, nti);
				if (angles != null) {
					addToHeap(angles, nt0, nti);
					if (testing) System.out.println(nt0 + " " + nti + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
				}
			}
		}
		
		for (int i = 0; i < 4; i++) {
			nti = nt1.neighbors[i];
			if ((nti != null) && (nti != nt0) && (nti != nt2)) {
				angles = getRoot(nt1, nti);
				if (angles != null) {
					addToHeap(angles, nt1, nti);
					if (testing) System.out.println(nt1 + " " + nti + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
				}
			}
		}

		for (int i = 0; i < 4; i++) {
			nti = nt2.neighbors[i];
			if ((nti != null) && (nti != nt0) && (nti != nt1)) {
				angles = getRoot(nt2, nti);
				if (angles != null) {
					addToHeap(angles, nt2, nti);
					if (testing) System.out.println(nt2 + " " + nti + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
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
		vns[0] = t0.corners[t0.apex(t1)];
		vns[1] = t1.corners[t1.apex(t2)];
		vns[2] = t2.corners[t2.apex(t0)];

		//Locate two shared vertices
		Vertex[] vs = new Vertex[2];
		int nrSharedFound = 0;
		for (int i = 0; (i < 4) && (nrSharedFound < 2); i++) {
			Vertex v = t0.corners[i];
			if (t1.hasVertex(v) && t2.hasVertex(v)) vs[nrSharedFound++] = v;
		}
		
		Tet nt0 = new Tet(vns[0], vns[1], vns[2], vs[0]);
		Tet nt1 = new Tet(vns[0], vns[1], vns[2], vs[1]);
		nt0.sortCorners();
		nt1.sortCorners();

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
		if (testing) t0.fromSceneEdges(scene);
		tets.remove(t1);
		t1.setAlive(false);		
		if (testing) t1.fromSceneEdges(scene);
		tets.remove(t2);
		t2.setAlive(false);
		if (testing) t2.fromSceneEdges(scene);
		Tet[] newTets = new Tet[2];
		newTets[0] = nt0;
		newTets[1] = nt1;
		newTets[0].toSceneEdges(scene, Color.black, 0.005, 0.0001);
		newTets[1].toSceneEdges(scene, Color.black, 0.005, 0.0001);
		tets.add(nt0);
		tets.add(nt1);
		angles = getRoot(nt0, nt1);
		if (angles != null) {
			addToHeap(angles, nt0, nt1);
			if (testing) System.out.println(nt0 + " " + nt1 + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
		}

		// check for new events
		
		Tet nti;
		for (int i = 0; i < 4; i++) {
			nti = nt0.neighbors[i];
			if ((nti != null) && (nti != nt1)) {
				angles = getRoot(nt0, nti);
				if (angles != null) {
					addToHeap(angles, nt0, nti);
					if (testing) System.out.println(nt0 + " " + nti + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
				}
			}
		}

		for (int i = 0; i < 4; i++) {
			nti = nt1.neighbors[i];
			if ((nti != null) && (nti != nt0)) {
				angles = getRoot(nt1, nti);
				if (angles != null) {
					addToHeap(angles, nt1, nti);
					if (testing) System.out.println(nt1 + " " + nti + " " + Functions.toDeg(angles[0]) + " " + Functions.toDeg(angles[1]));
				}
			}
		}	
		return newTets;
	}
	
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
		apexSegment = new LineSegment(t.corners[t.apex(nt)], nt.corners[nt.apex(t)]);
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
			if (testing) {
				System.out.print(t.toString());  if (!t.isAlive()) System.out.print(" ");
				System.out.print(nt.toString()); if (!nt.isAlive()) System.out.print(" ");
			}
			rotAngle = heapItem.getAngle();
			if (t.isAlive() && nt.isAlive() && (rotAngle < Constants.TAU)) {
				if (testing) {
//					t.toSceneFaces(scene, Color.blue);
//					nt.toSceneFaces(scene, Color.red);
					animate(scene, (rotAngle-angleTotal)/2.0);
					isDelaunay();
					animate(scene, (rotAngle-angleTotal)/2.0);	
				}
				angleTotal = rotAngle;
				
				if (!t.isConvex(nt)) {
					tt = KineticToolbox.getThirdTet(t, nt);
					if (testing) {
						System.out.print(tt.toString()); 
						if (!tt.isAlive())  System.out.print(" ");
						System.out.print(" rotated to angle = " + Functions.toDeg(heapItem.getAngle()));						
						System.out.println();
					}
					newTets = flip32(t, nt, tt);
					if (testing) System.out.println(++nrFlips + ". flip.");
				}
				else {
					if (testing) {
						System.out.print(" rotated to angle = " + Functions.toDeg(heapItem.getAngle()));
						System.out.println();
					}
					if (testing) {
						t.fromSceneFaces(scene);
						nt.fromSceneFaces(scene);
					}
					newTets = flip23(t, nt);
					if (testing) {
//						newTets[0].toSceneFaces(scene, Color.blue);
//						newTets[1].toSceneFaces(scene, Color.red);
//						newTets[2].toSceneFaces(scene, Color.green);
//						newTets[0].fromSceneFaces(scene);
//						newTets[1].fromSceneFaces(scene);
//						newTets[2].fromSceneFaces(scene);
					}
					if (testing) System.out.println(++nrFlips + ". flip.");
				}
			}
			else System.out.println();
		}
		animate(scene, Constants.TAU-angleTotal);

	}
	/* specify what is rotating - 3 methods to do it */
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

	
	public List<Tet> getTetrahedra(){ return tets; }
	
	public void toScene(J3DScene scene) { toScene(scene, Constants.bigDouble); }
	
	public void toScene(J3DScene scene, double alpha) {
		scene.removeAllShapes();
		Vector tr = new Vector (0.02, 0.02, 0.02);
		for (Tet tet : tets) { 
			if (tet.circumRadius() <= alpha) tet.toSceneEdges(scene, Color.black, 0.005, 0.0005); 
			}
		for (Vertex v : vertices) {
			if (v.getType() == VertexType.S) {
				v.toScene(scene, 0.01, Color.red);
				scene.addText(String.valueOf(v.getId()), v.add(tr));
			}
			else {
				v.toScene(scene, 0.01, Color.blue);
				Circle c = new Circle(new Point(0, 0, v.z()), v, getRotationAxis());
				c.toScene(scene, 0.002, 32);
				scene.addText(String.valueOf(v.getId()), v);
			}
		}
	}
	
	public static void main(String[] args){
//		List<Point> points = new java.util.LinkedList<Point>();
//		points.add(new Point(0,0,0));
//		points.add(new Point(1,0,0));
//		points.add(new Point(0,1,0));
//		points.add(new Point(0,0,1));
//		points.add(new Point(1.1,1.1,1.1));

//		Vertex v0 = new Vertex(new Point(0,0,0));
//		Vertex v1 = new Vertex(new Point(1,0,0));
//		Vertex v2 = new Vertex(new Point(0,1,0));
//		Vertex v3 = new Vertex(new Point(0,0,1));
//		Tet t = new Tet(new Vertex[]{v0,v1,v2,v3});
//		System.out.println(t.inSphere(new Point(0,0,0)));
//		System.out.println();
//		System.out.println(t.inSphere(new Point(-0.1,0,0)));
//		System.out.println();
//		System.out.println(t.inSphere(new Point(0.1,0,0)));
		
		Randomization.seed(3);
		List<Point> points = PointList.generatePointsInCube(100, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5);
		KineticDelaunayTessellation kDT = new KineticDelaunayTessellation(points);
		kDT.setRotationPoint(new Point(0.5,0.5,0.5));
		kDT.setRotationAxis(new Vector(0, 0, 1));
		kDT.setDirection(KineticDelaunayTessellation.Direction.CCW);
		kDT.setRotVertices(4, 14);

		for(Tet t: kDT.getTetrahedra()) System.out.println(t);
		
		
		kDT.alpha = 0.425;
		
		kDT.toScene(kDT.scene);
/*		for (Vertex v : kDT.vertices) {
			if (v.getType() == VertexType.S) v.toScene(kDT.scene, 0.05, Color.red);
			else {
				v.toScene(kDT.scene, 0.05, Color.blue);
				Circle c = new Circle(new Point(0, 0, v.z()), v, kDT.getRotationAxis());
				c.toScene(kDT.scene, 0.005, 32);
			}
		}
*/	
		kDT.rotate();
/*		Vector u = new Vector(1, 0, 0);	
		Vector v = new Vector(0, 1, 0);
		Vector c = u.cross(v).normalize();
		u.toScene(scene, Color.yellow, 0.05);
		v.toScene(scene, Color.pink, 0.05);
		c.toScene(scene, Color.green, 0.05);
		double alpha = Vector.getAngle(u, v);
		System.out.println("alpha = " + Functions.toDeg(alpha));
		for (Vertex p : kDT.vertices) { 
			Vector pOld = new Vector(p.x(), p.y(), p.z());
			p.rotation(c, alpha);
			Vector pNew = new Vector(p.x(), p.y(), p.z());
			System.out.println(Functions.toDeg(Vector.getAngle(pOld, pNew)));
			pOld.toScene(scene, Color.yellow, 0.02);
			pNew.toScene(scene, Color.pink, 0.02);
		}
		System.out.println(kDT.getTetrahedra().size());
		for(Tet t: kDT.getTetrahedra()){
			System.out.println(t);
			t.toScene(scene, Color.red);
		}
*/
	}
	
}
