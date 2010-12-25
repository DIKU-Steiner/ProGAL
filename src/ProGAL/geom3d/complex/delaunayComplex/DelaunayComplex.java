package ProGAL.geom3d.complex.delaunayComplex;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.complex.CEdge;
import ProGAL.geom3d.complex.CTriangle;
import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.complex.SimplicialComplex;
import ProGAL.geom3d.predicates.*;
import ProGAL.geom3d.predicates.Predicates.SphereConfig;

public class DelaunayComplex implements SimplicialComplex{

	private final List<CVertex> points;
	private final List<CEdge> edges;
	private final List<CTriangle> triangles;
	private final List<CTetrahedron> tetrahedra;
	private final Predicates predicates;
	private final Walk walk;
	private final Flip14 f14;
	private final Flips flips;

	public DelaunayComplex(List<Point> points) {
		this.points = new ArrayList<CVertex>(points.size());
		for(Point p: points) this.points.add(new CVertex(p));
		this.edges = new ArrayList<CEdge>(points.size()*6);
		this.triangles = new ArrayList<CTriangle>(points.size()*6);
		this.tetrahedra = new ArrayList<CTetrahedron>(points.size()*6);
		this.predicates = new ExactJavaPredicates();
		this.walk = new Walk(predicates);
		this.flips = new Flips(predicates);
		this.f14 = new Flip14(flips);
		
		compute();
	}

	public List<CTetrahedron> getTetrahedra() {
		return new ArrayList<CTetrahedron>(tetrahedra);
	}
	public List<CTriangle> getTriangles(){
		return new ArrayList<CTriangle>(triangles);
	}
	public List<CEdge> getEdges(){
		return new ArrayList<CEdge>(edges);
	}
	public List<CVertex> getVertices(){
		return new ArrayList<CVertex>(points);
	}

	private void compute() {
		double max = 1000;//TODO find a more meaningful max

		//Find the enclosing tetrahedron
		CTetrahedron next_t = new FirstTetrahedron(max);
		flips.addTetrahedron(next_t);

		//Iterér over punkterne
		for(CVertex p: points){
			next_t = walk.walk(next_t, p);
			next_t = f14.flip14(next_t, p);
			CTetrahedron tmp = flips.fixDelaunay();

			if (tmp != null) 
				next_t = tmp;
		}
		
		//Add edges and triangles
		completeComplex();
	}

	/** Add edges and triangles and remove auxiliary tetrahedra. */
	private void completeComplex() {
		//Add the non-modified tetrahedra that doesnt contain one of the big points
		for(CTetrahedron t: flips.getTetrahedrastack()){
			if (!t.isModified() &&!t.containsBigPoint()){
				tetrahedra.add(t);
			}
		}
		
		
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
			e1.addTriangle(t);
			e2.addTriangle(t);
			e3.addTriangle(t);
		}
		
		//Set faces of tetrahedra
		for(CTriangle t: triangles){
			CTetrahedron tet = t.getNeighbour(0);
			if(tet.getNeighbour(0).containsTriangle(t)) tet.setTriangle(0, t);
			else if(tet.getNeighbour(1).containsTriangle(t)) tet.setTriangle(1, t);
			else if(tet.getNeighbour(2).containsTriangle(t)) tet.setTriangle(2, t);
			else if(tet.getNeighbour(3).containsTriangle(t)) tet.setTriangle(3, t);
		}
	}

	/** Checks that the tetrahedra comply with the Delaunay-criteria. */
	public boolean checkTetrahedra() {
		for(CTetrahedron t: flips.getTetrahedrastack()){
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


