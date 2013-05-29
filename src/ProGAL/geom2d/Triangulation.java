package ProGAL.geom2d;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import ProGAL.dataStructures.Heap;
import ProGAL.dataStructures.Set;
import ProGAL.dataStructures.SortTool;
import ProGAL.dataStructures.SortToolLineSegment2dAroundCommonPoint;
import ProGAL.dataStructures.SortToolLineSegment2dByLength;
import ProGAL.dataStructures.SortToolPoint2dXY;
import ProGAL.dataStructures.Sorter;
import ProGAL.dataStructures.SorterQuick;
import ProGAL.geom2d.Point;
import ProGAL.geom2d.TriangulationVertex;
import ProGAL.geom2d.TriangulationFace;
import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.geom2d.viewer.TextShape;
import ProGAL.geom3d.kineticDelaunay.Tet;
import ProGAL.geom3d.predicates.ExactJavaPredicates;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.math.Constants;

public class Triangulation {
	
	public List<TriangulationVertex> vertices = new ArrayList<TriangulationVertex>();
	public List<TriangulationFace> triangulationFaces = new ArrayList<TriangulationFace>();

	public static enum TriangulationAlgorithm  { Greedy, Delaunay };
	private final ExactJavaPredicates pred = new ExactJavaPredicates();
	Sorter sort = new SorterQuick();
//	J2DScene scene = J2DScene.createJ2DSceneInFrame();

	/** creates a triangulation of a set of points 
	 * Two triangulation algorithms are so far implemented: Greedy, Delaunay */
	public Triangulation(PointSet points, TriangulationAlgorithm algorithm) {
		switch (algorithm) {
		case Delaunay:
			sort.Sort(points, new SortToolPoint2dXY());
			for (Point p: points) vertices.add(new TriangulationVertex(p.x(), p.y()));
			for (int i = 0; i < vertices.size(); i++) { vertices.get(i).setId(i); }
			TriangulationFace newTriangulationFace;
			if (Point.leftTurn(vertices.get(0), vertices.get(1), vertices.get(2))) newTriangulationFace = new TriangulationFace(vertices.get(0), vertices.get(1), vertices.get(2));
			else newTriangulationFace = new TriangulationFace(vertices.get(0), vertices.get(2), vertices.get(1));
			for (int i = 0; i < 3; i++) vertices.get(i).face = newTriangulationFace;
			triangulationFaces.add(newTriangulationFace);
			newTriangulationFace.id = 0;

			for (int i=3; i<vertices.size(); i++) {
				TriangulationVertex p = vertices.get(i);
				TriangulationVertex q = vertices.get(i-1);
				TriangulationFace oppTriangulationFace = q.face;
				TriangulationVertex r = oppTriangulationFace.corners[(oppTriangulationFace.getIndex(q)+1)%3];
				while (Point.rightTurn(q, r, p)) {
					newTriangulationFace = new TriangulationFace(p, r, q);//TODO optimize the rightturn in triangle constructor
					newTriangulationFace.neighbors[0] = oppTriangulationFace;
					if (p.face == null)  q.face = p.face = newTriangulationFace; 
					else {
						newTriangulationFace.neighbors[1] = p.face;
						p.face.neighbors[2] = newTriangulationFace;
						p.face = newTriangulationFace;
					}
					oppTriangulationFace.neighbors[(oppTriangulationFace.getIndex(q)+2)%3] = newTriangulationFace;
					triangulationFaces.add(newTriangulationFace);
					newTriangulationFace.id = triangulationFaces.size()-1;
					legalizeEdge(newTriangulationFace, 0, true);
					q = r;
					oppTriangulationFace = q.face;
					r = oppTriangulationFace.corners[(oppTriangulationFace.getIndex(q)+1)%3];
				}
				
				q = vertices.get(i-1);
				oppTriangulationFace = findLastTriangulationFace(q);
				r = oppTriangulationFace.corners[(oppTriangulationFace.getIndex(q)+2)%3];
				while(Point.leftTurn(q, r, p)){ 
					newTriangulationFace = new TriangulationFace(p,q,r);//TODO optimize the rightturn in triangle constructor
					newTriangulationFace.neighbors[0] = oppTriangulationFace;
					oppTriangulationFace.neighbors[(oppTriangulationFace.getIndex(q)+1)%3] = newTriangulationFace;
					if (p.face == null) r.face = p.face = newTriangulationFace;
					else {
						TriangulationFace lastTriangulationFace = findLastTriangulationFace(p);
						newTriangulationFace.neighbors[2] = lastTriangulationFace;
						lastTriangulationFace.neighbors[1] = newTriangulationFace;
						r.face = newTriangulationFace;
					}
					triangulationFaces.add(newTriangulationFace);
					newTriangulationFace.id = triangulationFaces.size()-1;
					legalizeEdge(newTriangulationFace, 0, true);
					q = r;
					oppTriangulationFace = findLastTriangulationFace(q);
					r = oppTriangulationFace.corners[(oppTriangulationFace.getIndex(q)+2)%3];
				}
			}
			break;

		case Greedy:
			for (Point p : points) vertices.add(new TriangulationVertex(p.x(), p.y()));
			Set<LineSegment> segments = new Set<LineSegment>();
			for (int i = 0; i < vertices.size(); i++) {
				TriangulationVertex u = vertices.get(i);
				u.setId(i);
				for (int j = i+1; j < vertices.size(); j++) {
					TriangulationVertex v = vertices.get(j);
					segments.insert(new LineSegment(u, v));
				}
			}
			sort.Sort(segments, new SortToolLineSegment2dByLength());
			List<LineSegment> acceptedSegments = new ArrayList<LineSegment>();
			for (LineSegment seg : segments) {
				boolean cont = true;
				for (LineSegment segA : acceptedSegments) {
					cont = !seg.intersects(segA);
					if (!cont) break;
				}
				if (cont) {
					acceptedSegments.add(seg);
				}
			}
			Set<LineSegment> incidentSegments = new Set<LineSegment>();
			for (TriangulationVertex v : vertices) {
				for (LineSegment seg : acceptedSegments) {
					if ((TriangulationVertex)seg.a == v) incidentSegments.insert(seg);
					else {
						if ((TriangulationVertex)seg.b == v) incidentSegments.insert(seg.reverse());
					}
				}
				sort.Sort(incidentSegments, new SortToolLineSegment2dAroundCommonPoint((Point)v));
				int sz = incidentSegments.getSize();
				TriangulationVertex p; 
				TriangulationVertex q = (TriangulationVertex)incidentSegments.get(0).getB();
				for (int j = 0; j < sz; j++) {
					p = q;
					q = (TriangulationVertex)incidentSegments.get((j+1)%sz).getB();
					System.out.println(p.id + "," + v.id + "," + q.id);
					if (Point.rightTurn(p, v, q)) {
						TriangulationFace TriangulationFace = findTriangulationFace(v, p, q);
//						TriangulationFace.draw(scene);		
					}
				}
				incidentSegments.clear();
			}
			for (TriangulationFace TriangulationFace : triangulationFaces) {
				for (int i = 0; i < 3; i++) {
					TriangulationVertex v = TriangulationFace.corners[i];
					TriangulationVertex w = TriangulationFace.corners[(i+1)%3];
					TriangulationVertex z = TriangulationFace.corners[(i+2)%3];
					TriangulationFace prevTriangulationFace = findPrevTriangulationFace(v, w);
					if ((v.face == null) || (prevTriangulationFace == null)) v.face = TriangulationFace;
					TriangulationFace.neighbors[(i+2)%3] = prevTriangulationFace;
					TriangulationFace.neighbors[(i+1)%3] = findNextTriangulationFace(v,z);
					TriangulationFace.neighbors[i] = findPrevTriangulationFace(w, z);
				}
			}
			break;
		}
	}

	/** returns TriangulationFace with specified vertices. If not found, creates such a TriangulationFace */
	public TriangulationFace findTriangulationFace(TriangulationVertex u, TriangulationVertex v, TriangulationVertex w) {
		for (TriangulationFace TriangulationFace : triangulationFaces) {
			if (TriangulationFace.hasVertex(u) & TriangulationFace.hasVertex(v) & TriangulationFace.hasVertex(w)) return TriangulationFace;
		}
		TriangulationFace newTriangulationFace = new TriangulationFace(u, v, w);
		triangulationFaces.add(newTriangulationFace);
		newTriangulationFace.id = triangulationFaces.size()-1;
		return newTriangulationFace;
	}

	/** returns next (left when looking from vertex u to vertex v) TriangulationFace. */
	public TriangulationFace findNextTriangulationFace(TriangulationVertex u, TriangulationVertex v) {
		for (TriangulationFace TriangulationFace : triangulationFaces) {
			if (TriangulationFace.hasVertex(u) && TriangulationFace.hasVertex(v) && 
				(TriangulationFace.corners[(TriangulationFace.getIndex(u)+1)%3] == v)) return TriangulationFace;
		}
		return null;
	}
	
	/** returns the previous (right when looking from vertex u to vertex v) TriangulationFace. */
	public TriangulationFace findPrevTriangulationFace(TriangulationVertex u, TriangulationVertex v) {
		for (TriangulationFace TriangulationFace : triangulationFaces) {
			if (TriangulationFace.hasVertex(u) && TriangulationFace.hasVertex(v) && 
				(TriangulationFace.corners[(TriangulationFace.getIndex(v)+1)%3] == u)) return TriangulationFace;
		}
		return null;
	}
	
	/** returns the last TriangulationFace incident with vertex v (null if v is not on the convex hull 
	 * and therefore has no last TriangulationFace */
	public TriangulationFace findLastTriangulationFace(TriangulationVertex v) {
		TriangulationFace TriangulationFace = v.face;
		int indx = TriangulationFace.getIndex(v);
		if (TriangulationFace.neighbors[(indx+2)%3] != null) return null;
		TriangulationFace nextTriangulationFace = TriangulationFace.neighbors[(indx+1)%3];
		while (nextTriangulationFace != null) {
			TriangulationFace = nextTriangulationFace;
			indx = TriangulationFace.getIndex(v);
			nextTriangulationFace = TriangulationFace.neighbors[(indx+1)%3];
		} 
		return TriangulationFace;
	}

	/** returns the first TriangulationFace incident with vertex v (null if v is not on the convex hull
	 * and therefore has no last TriangulationFace */
	public TriangulationFace findFirstTriangulationFace(TriangulationVertex v) {
		int indx = v.face.getIndex(v);
		if (v.face.neighbors[(indx+2)%3] != null) return null; else return v.face;
	}
	 
	public void legalizeEdge(TriangulationFace triangulationFace, int indx, boolean recursive) {
		legalizeEdge(triangulationFace, indx, recursive, null);
	}
	
	public void legalizeEdge(TriangulationFace triangulationFace, int indx, boolean recursive, J2DScene scene){
		TriangulationFace oppTriangulationFace = triangulationFace.neighbors[indx];
		int oppIndx = (oppTriangulationFace.getIndex(triangulationFace.corners[(indx+1)%3])+1)%3;
		double inc = pred.incircle(triangulationFace.corners[0].getCoords(), triangulationFace.corners[1].getCoords(), triangulationFace.corners[2].getCoords(), 
								   oppTriangulationFace.corners[oppIndx].getCoords());
		if (inc > 0) flip(triangulationFace, oppTriangulationFace, recursive, scene, false);
	}

	public TriangulationFace[] flip(TriangulationFace t013, TriangulationFace t123, boolean recursive) {
		return flip(t013, t123, recursive, null, false);
	}

	
	
	public TriangulationFace[] flip(TriangulationFace t013, TriangulationFace t123, boolean recursive, 
					 J2DScene scene, boolean testing){
		t013.setAlive(false);
		t013.hide(scene, testing);
		triangulationFaces.remove(t013);
		t123.setAlive(false);
		t123.hide(scene, testing);
		triangulationFaces.remove(t123);
		int p0Indx = t013.getIndex(t123);
		int p2Indx = t123.getIndex(t013);
		
		TriangulationVertex p0 = t013.corners[p0Indx];
		TriangulationVertex p1 = t013.corners[(p0Indx+1)%3];
		TriangulationVertex p2 = t123.corners[p2Indx];
		TriangulationVertex p3 = t123.corners[(p2Indx+1)%3];
		
		TriangulationFace t012 = new TriangulationFace(p0, p1, p2); 
		triangulationFaces.add(t012);
		t012.draw(scene, testing);
		TriangulationFace t023 = new TriangulationFace(p0, p2, p3); 
		triangulationFaces.add(t023);
		t023.draw(scene, testing);

		TriangulationFace t01 = t013.neighbors[(p0Indx+2)%3];
		TriangulationFace t12 = t123.neighbors[(p2Indx+1)%3];
		TriangulationFace t23 = t123.neighbors[(p2Indx+2)%3];
		TriangulationFace t30 = t013.neighbors[(p0Indx+1)%3];

		t012.neighbors[0] = t12;
		t012.neighbors[1] = t023;
		t012.neighbors[2] = t01;
		t023.neighbors[0] = t23;
		t023.neighbors[1] = t30;
		t023.neighbors[2] = t012;
		if (t01 != null) t01.neighbors[t01.getIndex(t013)] = t012;
		if (t12 != null) t12.neighbors[t12.getIndex(t123)] = t012;
		if (t23 != null) t23.neighbors[t23.getIndex(t123)] = t023;
		if (t30 != null) t30.neighbors[t30.getIndex(t013)] = t023;

		if (p0.face == t013) p0.face = t012;
		if ((p1.face == t123) || (p1.face == t013)) p1.face = t012;
		if (p2.face == t123) p2.face = t023;
		if ((p3.face == t013) || (p3.face == t123)) p3.face = t023;
		
		t012.setId(t013.id);
		t023.setId(t123.id);
		t013 = t012;
		t123 = t023;
		if (recursive) {
			if (t12 != null) legalizeEdge(t012, 0, true, scene);
			if (t23 != null) legalizeEdge(t023, 0, true, scene);
		}
		TriangulationFace[] newFaces = new TriangulationFace[2];
		newFaces[0] = t012;
		newFaces[1] = t023;
		return newFaces;
	}
	
	/* v becomes a new vertex on the convex hull */
	public void boundaryFlipOut(TriangulationVertex v, TriangulationFace face, J2DScene scene, boolean testing) {	
		int indx = face.getIndex(v);
		face.setAlive(false);
		face.hide(scene, testing);
		triangulationFaces.remove(face);
		TriangulationVertex a = face.getCorner((indx+2)%3);
		TriangulationVertex b = face.getCorner((indx+1)%3);
		if (triangulationFaces.isEmpty()) {
			TriangulationFace newFace = new TriangulationFace(v, a, b);
			triangulationFaces.add(newFace);
			newFace.setId(triangulationFaces.size());
			if (testing) newFace.draw(scene, testing);
		}
		else {
			TriangulationFace firstFace = v.face;
			TriangulationFace lastFace = v.getLastFace();
			TriangulationFace nextFace = face.getNeighbor((indx+1)%3);
			TriangulationFace prevFace  = face.getNeighbor((indx+2)%3);
			if (firstFace == lastFace) {
				v.face = nextFace;
				b.face = prevFace;
				nextFace.setNeighbor((nextFace.getIndex(v)+2)%3, null);
				prevFace.setNeighbor((prevFace.getIndex(v)+1)%3, null);
			}
			else {
				if (prevFace != null) b = firstFace.getCorner((firstFace.getIndex(v)+1)%3);
				if (nextFace != null) a = lastFace.getCorner((lastFace.getIndex(v)+2)%3);
				TriangulationFace newFace = new TriangulationFace(v, a, b);
				newFace.setNeighbor(0, null);
				newFace.setNeighbor(1, prevFace);
				newFace.setNeighbor(2, nextFace);
				if (prevFace != null) prevFace.setNeighbor((prevFace.getIndex(v)+2)%3, newFace);
				if (nextFace != null) nextFace.setNeighbor((nextFace.getIndex(v)+1)%3, newFace); 
				triangulationFaces.add(newFace);
				newFace.setId(triangulationFaces.size());
				a.face = newFace;
				v.face = newFace;
				if (prevFace != null) prevFace.setNeighbor((prevFace.getIndex(v)+2)%3, newFace); 
				else face.setNeighbor((face.getIndex(v)+2)%3, newFace);
				if (nextFace != null) nextFace.setNeighbor((nextFace.getIndex(v)+1)%3, newFace);
				else face.setNeighbor((face.getIndex(v)+1)%3, newFace);
			}
		}
	}
	
	public void boundaryFlipIn(TriangulationVertex a, TriangulationVertex b, TriangulationVertex c) {
		boundaryFlipIn(a, b, c, null, false);
	}
	
	public void boundaryFlipIn(TriangulationVertex a, TriangulationVertex b, TriangulationVertex c, 
							   J2DScene scene, boolean testing) {
		TriangulationFace bcFace = b.getFace();
		if (b.getNextFace(bcFace) != null) {
			TriangulationFace newFace = new TriangulationFace(a, c, b);
			newFace.neighbors[0] = bcFace;
			TriangulationFace abFace = a.getFace();
			newFace.neighbors[1] = abFace;
			newFace.neighbors[2] = null;
			bcFace.setNeighbor((bcFace.getIndex(b)+2)%3, newFace);
			abFace.setNeighbor((abFace.getIndex(a)+2)%3, newFace);
			a.face = newFace;
			triangulationFaces.add(newFace);
			newFace.id = triangulationFaces.size()-1;
			newFace.draw(scene, testing);
		}
		else {
			if (bcFace.getOppFace(b) == null) {   // this case applies only if DT has one triangle
				b.getFace().setAlive(false);
				TriangulationFace newFace = new TriangulationFace(a, c, b, scene, testing);
				for (int i = 0; i < 3; i++) newFace.neighbors[i] = null;
				a.face = newFace;
				b.face = newFace;
				c.face = newFace;
				triangulationFaces.add(newFace);
				newFace.id = triangulationFaces.size()-1;
				newFace.draw(scene, testing);
			}
			else {  // this case probably never applies 
				System.out.println("This boundaryFlipIn case should not occur");
			}
		}
	}
	
	private class HeapItem<T> {
		private double power;
		private TriangulationFace face;
		private TriangulationFace nextFace;
		
		private HeapItem(TriangulationFace face, TriangulationFace nextFace, double power) {
			this.power = power;
			this.face = face;
			this.nextFace = nextFace;
		}
				
		private double getPower() { return power;} 
		private TriangulationFace getFace() { return face; }
		private TriangulationFace getNextFace() { return nextFace; }
	}
	private class SortToolHeapItems implements SortTool {
		public int compare(Object x1, Object x2) {
			if ((x1 instanceof HeapItem) && (x2 instanceof HeapItem)) {
				double d1 = ((HeapItem)x1).getPower();
				double d2 = ((HeapItem)x2).getPower();
				if (d1 < d2) return COMP_LESS;
				else { if (d1 > d2) return COMP_GRTR; else return COMP_EQUAL; }
			}
			else throw SortTool.err1;
		}
	}

	private double getPower(Point a, Point b, Point c, Point p) {
		double area = Point.area(a, b, c);
		if (area <= 0.0) return Constants.bigDouble;
		return -Point.inCircle(b, a, c, p)/area;
	}
	
	/* deletes vertex v and all incident faces. Determines Delauney triangulation of the remaining points. */
	public void delete(TriangulationVertex u, J2DScene scene) {
		boolean testing = true;
		int indx, indxN;
		TriangulationFace face, face1, face2, face3, prevFace, nextFace, newFace, newFace1, newFace2, oppFace1, oppFace2, oppFace3;
		TriangulationVertex a, b, c;
		
		// set up the heap of consecutive pairs of faces around vertex u
		Heap heap = new Heap(vertices.size(), new SortToolHeapItems());
		double power;
		int sz = u.getFaces().size();
		face = u.getFaces().get(sz-1);
		indx = face.getIndex(u);
		for (int i = 0; i < sz; i++) {
			nextFace = u.getFaces().get(i);
			indxN = nextFace.getIndex(u);
			power = getPower(face.corners[(indx+1)%3], face.corners[(indx+2)%3], nextFace.corners[(indxN+2)%3],   u);
			if (power < Constants.bigDouble) heap.insert(new HeapItem(face, nextFace, power));
			if (testing) {
				Circle cir = new Circle(face.corners[(indx+1)%3], face.corners[(indx+2)%3], nextFace.corners[(indxN+2)%3]);
				cir.toScene(scene, Color.blue);
				System.out.println(power);
				scene.removeShape(cir);
			}
			face = nextFace;
			indx = indxN;
		}
				
		// loop
		TriangulationFace[] newFaces = new TriangulationFace[2];
		while (sz != 3) {
			HeapItem item = (HeapItem)heap.extract();
			face1 = item.getFace();
			face2 = item.getNextFace();
			if (face1.isAlive() && face2.isAlive()) {
			
				// add new face to the triangulation
				int indx1 = face1.getIndex(u);
				int indx2 = face2.getIndex(u);
				a = face1.corners[(indx1+1)%3];
				b = face1.corners[(indx1+2)%3];
				c = face2.corners[(indx2+2)%3];
				oppFace1 = face1.getNeighbor(indx1);
				oppFace2 = face2.getNeighbor(indx2);
			
				newFace1 = new TriangulationFace(a, b, c);
				newFace2 = new TriangulationFace(a, c, u);

				newFace1.setNeighbor(0, oppFace2);
				if (oppFace2 != null) oppFace2.setNeighbor((oppFace2.getIndex(b)+1)%3, newFace1);
				newFace1.setNeighbor(1, newFace2);
				newFace1.setNeighbor(2, oppFace1);
				if (oppFace1 != null) oppFace1.setNeighbor((oppFace1.getIndex(a)+1)%3, newFace1);
				if (a.face == face1) a.face = newFace1;
				if ((b.face == face1) || (b.face == face2)) b.face = newFace1;
				if (c.face == face2) c.face = newFace1;
			    triangulationFaces.add(newFace1);
				if (testing) {
					Circle cir = new Circle(a, b, c);
					cir.toScene(scene, Color.blue);
					scene.removeShape(cir);
				}

			    // update star of u
				prevFace = u.getPrevFace(face1);
				nextFace = u.getNextFace(face2);
				
			    newFace2.setNeighbor(0, nextFace);
			    newFace2.setNeighbor(1, prevFace);
			    newFace2.setNeighbor(2, newFace1);
			    
			    nextFace.setNeighbor((nextFace.getIndex(c)+1)%3, newFace2);
			    prevFace.setNeighbor((prevFace.getIndex(a)+2)%3, newFace2);
                if ((u.face == face1) || (u.face == face2)) u.face = newFace2;
                triangulationFaces.add(newFace2);
                
			    power = getPower(prevFace.corners[(prevFace.getIndex(a)+2)%3], a, c, u);
			    if (power < Constants.bigDouble) heap.insert(new HeapItem(prevFace, newFace2, power));
			    power = getPower(a, c, nextFace.corners[(nextFace.getIndex(c)+1)%3], u);
			    if (power < Constants.bigDouble) heap.insert(new HeapItem(newFace2, nextFace, power));
			    
			    face1.setAlive(false);
				triangulationFaces.remove(face1);
			    face2.setAlive(false);
				triangulationFaces.remove(face2);
			    sz--;
       		}   
		}		
			
       	// remove u and its 3 faces
		face1 = u.getFace();
		face2 = u.getNextFace(face1);
		face3 = u.getNextFace(face2);
		a = face1.corners[(face1.getIndex(u)+1)%3];
		b = face2.corners[(face2.getIndex(u)+1)%3];
		c = face3.corners[(face3.getIndex(u)+1)%3];
		oppFace1 = face1.getNeighbor(face1.getIndex(u));
		oppFace2 = face2.getNeighbor(face2.getIndex(u));
		oppFace3 = face3.getNeighbor(face3.getIndex(u));
		 
		newFace = new TriangulationFace(a, b, c);
		newFace.setNeighbor(0, oppFace2);
		newFace.setNeighbor(1, oppFace3);
		newFace.setNeighbor(2, oppFace1);
		
		if (oppFace1 != null) oppFace1.setNeighbor((oppFace1.getIndex(a)+1)%3, newFace);
		if (oppFace2 != null) oppFace2.setNeighbor((oppFace2.getIndex(b)+1)%3, newFace);
		if (oppFace3 != null) oppFace3.setNeighbor((oppFace3.getIndex(c)+1)%3, newFace);

		face1.setAlive(false);
		triangulationFaces.remove(face1);
		face2.setAlive(false);
		triangulationFaces.remove(face2);
		face3.setAlive(false);
		triangulationFaces.remove(face3);
		if ((a.face == face1) || (a.face == face2)) a.face = newFace;
		if ((b.face == face2) || (b.face == face3)) b.face = newFace;
		if ((c.face == face3) || (c.face == face1)) c.face = newFace;
		u.face = null;
		
       	triangulationFaces.add(newFace);
       	if (testing) {
       		scene.removeAllShapes();
       		draw(scene);
       	}
	}
	
	public void print() {
		for (TriangulationVertex v : vertices) {
			System.out.print("Vertex " + v.id);  
		    if (v.face != null) System.out.println(" has first TriangulationFace " + v.face.toString()); else System.out.println(" has no TriangulationFaces.");
		}
		for (TriangulationFace TriangulationFace : triangulationFaces) {
			System.out.print("TriangulationFace " + TriangulationFace.id + " " + TriangulationFace.toString() + " has neighbors: ");
			for (int i = 0; i < 3; i++) {
				if (TriangulationFace.neighbors[i] != null) System.out.print(TriangulationFace.neighbors[i].toString() + " ");
			}
			System.out.println();
		}
	}
	public void draw(J2DScene scene) {
		scene.removeAllShapes();
		for (TriangulationFace t : triangulationFaces) {
			t.draw(scene);
		}
		for (TriangulationVertex v: vertices) scene.addShape(new TextShape(String.valueOf(v.id), v, 0.1));
	}
	public static void main(String[] args) {
		PointSet points = new PointSet(150);
		Triangulation tr = new Triangulation(points, TriangulationAlgorithm.Delaunay);
		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		tr.draw(scene);
		tr.print();
		
	}
}