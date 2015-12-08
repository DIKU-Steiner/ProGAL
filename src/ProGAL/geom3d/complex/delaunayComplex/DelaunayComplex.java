package ProGAL.geom3d.complex.delaunayComplex;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.complex.CEdge;
import ProGAL.geom3d.complex.CTriangle;
import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.complex.SimplicialComplex;
import ProGAL.geom3d.predicates.*;
import ProGAL.geom3d.predicates.Predicates.SphereConfig;
import ProGAL.math.Randomization;

/** <p>
 *  A Delaunay complex for a set of d-dimensional points is a tesselation of the points such that no point is inside 
 *  the circumscribing hypersphere of the d-simplices (for the 3D case: Tetrahedra). 
 *  </p>
 *  
 *  <p>
 *  This class builds a three-dimensional Delaunay complex in the constructor and accesses it using e.g. the 
 *  <code>getTetrahedra</code> method. The following example displays the Delaunay complex of ten random points. 
 *  <pre>
 *  {@code
 *  //Generate the complex
 *  List<Point> pl = PointList.generatePointsInCube(10);
 *  DelaunayComplex dc = new DelaunayComplex(pl);
 *    
 *  //Display the complex
 *  J3DScene scene = J3DScene.createJ3DSceneInFrame();
 *  for(CTetrahedron t: dc.getTetrahedra()){
 *    scene.addShape(t, new Color(200,100,100,100));
 *  }
 *  }
 *  </pre>     
 *  </p>
 *  <p>The original point-set is left unaltered and non-referenced by this class. A new set of vertices is 
 *  allocated using the CVertex class. These are randomly perturbed to avoid degeneracies. If one wishes to 
 *  associate the original points with a vertex in the complex it would be sufficient to test if the distance 
 *  between the point and the vertex is less than 0.0001.</p>
 *  
 *  <p>The complex is bounded by a big tetrahedron whose corner-points are located sufficiently far from any of 
 *  the vertices of the complex. The simplices that have one of these 'big points' as corners can not be accessed 
 *  directly via the getter-methods, but they will be neighbors of other normal simplices. For instance:
 *  <pre>
 *  DelaunayComplex dc = new DelaunayComplex( PointList.generatePointsInCube(4) );
 *  for(CTetrahedron t: dc.getTetrahedra()){
 *  	System.out.println( t.containsBigPoint() );
 *  	System.out.println( t.getNeighbor(0).containsBigPoint() );
 *  	System.out.println( t.getNeighbor(1).containsBigPoint() );
 *  	System.out.println( t.getNeighbor(2).containsBigPoint() );
 *  	System.out.println( t.getNeighbor(3).containsBigPoint() );
 *  }
 *  </pre> 
 *  Will print false, true, true, true and true.</p>
 *  
 *  @author R.Fonseca
 */
public class DelaunayComplex implements SimplicialComplex{
	private final List<CVertex> points;
	private final List<CEdge> edges;
	private final List<CTriangle> triangles;
	private final List<CTetrahedron> tetrahedra;
	private final Predicates predicates;
	private final Walk walk;
	private final Flip14 f14;
	private final Flips flips;

	/** Builds the Delaunay complex of the specified point-set */
	public DelaunayComplex(List<Point> points) {
		this.points = new ArrayList<CVertex>(points.size());
		
		int i=0;
		for(Point p: points) this.points.add(new CVertex(p, i++));
		this.edges = new ArrayList<CEdge>(points.size()*6);
		this.triangles = new ArrayList<CTriangle>(points.size()*6);
		this.tetrahedra = new ArrayList<CTetrahedron>(points.size()*6);
		
		this.predicates = new ExactJavaPredicates();
		this.walk = new Walk(predicates);
		this.flips = new Flips(predicates);
		this.f14 = new Flip14(flips);
		
		compute();
		completeComplex();
	}
	
	/** TODO: Finish */
	public boolean isDelaunay() {
		for (CTetrahedron tetr : tetrahedra) {
			Sphere sphere = new Sphere(tetr);
		}
		return true;
	}
	
	/** Get the tetrahedra in the complex. The tetrahedra that has 'big points' as corners are not returned */
	public List<CTetrahedron> getTetrahedra() {
		return new ArrayList<CTetrahedron>(tetrahedra);
	}
	
	
	/** returns all tetrahedra (including tetrahedra with bigPoint */
	public List<CTetrahedron> getAllTetrahedra() {
		List<CTetrahedron> allTetrahedra = new ArrayList<CTetrahedron>();
		for (CVertex v : getVertices() ) {
			List<CTetrahedron> adjacentTetrahedra = v.getAllAdjacentTetrahedra();
			for (CTetrahedron tetr : adjacentTetrahedra) {
				if (!allTetrahedra.contains(tetr)) allTetrahedra.add(tetr);
			}
		}
		return allTetrahedra;
	}
	
	/** returns all big tetrahedra */
	public List<CTetrahedron> getBigTetrahedra() {
		List<CTetrahedron> bigTetrahedra = new ArrayList<CTetrahedron>();
		for (CVertex v : getVertices() ) {
			List<CTetrahedron> adjacentTetrahedra = v.getAllAdjacentTetrahedra();
			for (CTetrahedron tetr : adjacentTetrahedra) {
				if (tetr.containsBigPoint() && !bigTetrahedra.contains(tetr)) bigTetrahedra.add(tetr);
			}
		}
		return bigTetrahedra;
	}


	/** Get the triangles in the complex. The triangles that has 'big points' as corners are not returned */
	public List<CTriangle> getTriangles(){
		return new ArrayList<CTriangle>(triangles);
	}
	/** Get the edges in the complex. The edges that has 'big points' as end-points are not returned */
	public List<CEdge> getEdges(){
		return new ArrayList<CEdge>(edges);
	}
	/** Get the vertices in the complex. The 'big points' are not returned */
	public List<CVertex> getVertices(){
		return new ArrayList<CVertex>(points);
	}

	public CVertex getVertex(int i) { return points.get(i); }
	
	protected void compute() {
		double max = 1000;//TODO find a more meaningful max
		
		//TODO: Take care of degeneracies in a better way than perturbation
		for(CVertex v: points){
			v.addThis(new Vector(
					Randomization.randBetween(-0.00001, 0.00001),
					Randomization.randBetween(-0.00001, 0.00001),
					Randomization.randBetween(-0.00001, 0.00001)
					)
			);
		}

		//Find the enclosing tetrahedron
		CTetrahedron next_t = new FirstTetrahedron(max);
		flips.addTetrahedron(next_t);

		//Iterer over punkterne
		for(CVertex p: points){
			next_t = walk.walk(next_t, p);
			next_t = f14.flip14(next_t, p);
			CTetrahedron tmp = flips.fixDelaunay();

			if (tmp != null) 
				next_t = tmp;
		}
		
	}

	/** Add edges and triangles and remove auxiliary tetrahedra. */
	protected void completeComplex() {
		tetrahedra.clear();
		triangles.clear();
		edges.clear();
	
		//Add the non-modified tetrahedra that doesnt contain one of the big points
		for(CTetrahedron t: flips.getTetrahedrastack()){
			if (!t.isModified() && !t.containsBigPoint()){
				tetrahedra.add(t);
			}
		}
//		flips.setTetrahedrastack(null);//Should free up some memory after garbage collection
//		flips.getTetrahedrastack().clear();
		
		class VertexPair{
			CVertex p1, p2;
			VertexPair(CVertex p1, CVertex p2){
				this.p1 = p1;
				this.p2 = p2;
			}
			public boolean equals(Object o){
				return (((VertexPair)o).p1==p1 && ((VertexPair)o).p2==p2)||(((VertexPair)o).p1==p2 && ((VertexPair)o).p2==p1); 
			}
			public int hashCode(){
				return p1.hashCode()^p2.hashCode();
			}
		}
		
		//Construct edges
		Map<VertexPair, CEdge> edgeMap = new HashMap<VertexPair, CEdge>();
		for(CTetrahedron t: tetrahedra){
			edgeMap.put( new VertexPair(t.getPoint(0),t.getPoint(1)),new CEdge(t.getPoint(0),t.getPoint(1)) );
			edgeMap.put( new VertexPair(t.getPoint(0),t.getPoint(2)),new CEdge(t.getPoint(0),t.getPoint(2)) );
			edgeMap.put( new VertexPair(t.getPoint(0),t.getPoint(3)),new CEdge(t.getPoint(0),t.getPoint(3)) );
			edgeMap.put( new VertexPair(t.getPoint(1),t.getPoint(2)),new CEdge(t.getPoint(1),t.getPoint(2)) );
			edgeMap.put( new VertexPair(t.getPoint(1),t.getPoint(3)),new CEdge(t.getPoint(1),t.getPoint(3)) );
			edgeMap.put( new VertexPair(t.getPoint(2),t.getPoint(3)),new CEdge(t.getPoint(2),t.getPoint(3)) );
		}
		edges.addAll(edgeMap.values());
		
		for(CEdge e: edges){
			((CVertex)e.getA()).addAdjacentEdge(e);
			((CVertex)e.getB()).addAdjacentEdge(e);
		}
		
		//Construct triangles
		Set<CTriangle> triangleSet = new HashSet<CTriangle>();
		for(CTetrahedron t: tetrahedra){
			triangleSet.add(new CTriangle(t.getPoint(1),t.getPoint(2),t.getPoint(3), t,t.getNeighbour(0)));
			triangleSet.add(new CTriangle(t.getPoint(0),t.getPoint(2),t.getPoint(3), t,t.getNeighbour(1)));
			triangleSet.add(new CTriangle(t.getPoint(0),t.getPoint(1),t.getPoint(3), t,t.getNeighbour(2)));
			triangleSet.add(new CTriangle(t.getPoint(0),t.getPoint(1),t.getPoint(2), t,t.getNeighbour(3)));
		}
		triangles.addAll(triangleSet);
		for(CTriangle t: triangles){
			CEdge e1 = edgeMap.get(new VertexPair((CVertex)t.getPoint(0),(CVertex)t.getPoint(1)));
			CEdge e2 = edgeMap.get(new VertexPair((CVertex)t.getPoint(1),(CVertex)t.getPoint(2)));
			CEdge e3 = edgeMap.get(new VertexPair((CVertex)t.getPoint(2),(CVertex)t.getPoint(0)));
			t.setEdge(0, e1);
			t.setEdge(1, e2);
			t.setEdge(2, e3);
			e1.addTriangle(t);
			e2.addTriangle(t);
			e3.addTriangle(t);
		}
		
		//Set faces of tetrahedra
		for(CTriangle t: triangles){
			CTetrahedron tet = t.getAdjacentTetrahedron(0);
			if(tet.getNeighbour(0).containsTriangle(t)) tet.setTriangle(0, t);
			else if(tet.getNeighbour(1).containsTriangle(t)) tet.setTriangle(1, t);
			else if(tet.getNeighbour(2).containsTriangle(t)) tet.setTriangle(2, t);
			else if(tet.getNeighbour(3).containsTriangle(t)) tet.setTriangle(3, t);

			tet = t.getAdjacentTetrahedron(1);
			if(tet.getNeighbour(0).containsTriangle(t)) tet.setTriangle(0, t);
			else if(tet.getNeighbour(1).containsTriangle(t)) tet.setTriangle(1, t);
			else if(tet.getNeighbour(2).containsTriangle(t)) tet.setTriangle(2, t);
			else if(tet.getNeighbour(3).containsTriangle(t)) tet.setTriangle(3, t);
		}
	}
	
	/** The vertex-hull of v is the set of all tetrahedrons that has v as a corner-point */
	public Set<CTetrahedron> getVertexHull(CVertex v){
		Set<CTetrahedron> hull = new HashSet<CTetrahedron>();
		for(CEdge e: v.getAdjacentEdges()){
			for(CTriangle tri: e.getAdjacentTriangles()){
				hull.add(tri.getAdjacentTetrahedron(0));
				hull.add(tri.getAdjacentTetrahedron(1));
			}
		}
		return hull;
	}
	

	/** Checks that all tetrahedra comply with the Delaunay-criteria. */
	public boolean checkTetrahedra() {
		for(CTetrahedron t: tetrahedra){
			for(CVertex p: points){
				if(		t.getPoint(0)!=p && 
						t.getPoint(1)!=p && 
						t.getPoint(2)!=p && 
						t.getPoint(3)!=p &&
						predicates.insphere(t, p).equals(SphereConfig.INSIDE) )
				{
					return false;
				}
			}
		}
		return true;
	}
}
