package ProGAL.geom3d.complex.alphaComplex;

import java.util.ArrayList; 
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import java.util.HashMap;

import ProGAL.geom3d.PointList;
import ProGAL.geom3d.Simplex;
import ProGAL.geom3d.complex.*;
import ProGAL.geom3d.complex.delaunayComplex.DelaunayComplex;
import ProGAL.geom3d.predicates.*;
import ProGAL.geom3d.predicates.Predicates.SphereConfig;

public class AlphaFiltration {

	private final DelaunayComplex del3d;
	private Predicates p = new ExactJavaPredicates();
	private final AlphaComparator alphaOrdering = new AlphaComparator();
	private final Map<Simplex, SimplexAlphaProperties> propertyMap = new HashMap<Simplex, SimplexAlphaProperties>();
	
	private ArrayList<CTetrahedron> tetrahedra = new ArrayList<CTetrahedron>();
	private ArrayList<CTriangle> triangles = new ArrayList<CTriangle>();
	private ArrayList<CEdge> edges = new ArrayList<CEdge>();
	private ArrayList<CVertex> vertices = new ArrayList<CVertex>();
	private ArrayList<Simplex> simplices = new ArrayList<Simplex>();
	

	public AlphaFiltration(PointList pl){
		this.del3d = new DelaunayComplex(pl);
		compute();
	}

	public AlphaFiltration(DelaunayComplex d3d){
		this.del3d = d3d;
		compute();
	}

	private void compute(){
		tetrahedra.addAll(del3d.getTetrahedra());
		triangles.addAll(del3d.getTriangles());
		edges.addAll(del3d.getEdges());
		vertices.addAll(del3d.getVertices());
		computeTetrahedraIntervals();
		computeTriangleIntervals();
		computeEdgeIntervals();
		computeVertexIntervals();
		simplices.addAll(vertices);
		simplices.addAll(edges);
		simplices.addAll(triangles);
		simplices.addAll(tetrahedra);
		Collections.sort(simplices,alphaOrdering);
	}


	public List<CTetrahedron> getTetrahedra(double alpha){
		List<CTetrahedron> ret = new ArrayList<CTetrahedron>();
		for(CTetrahedron t: tetrahedra)	if(getInAlpha(t)<alpha) ret.add(t);
		return ret;
	}
	public List<CTriangle> getTriangles(double alpha){
		List<CTriangle> ret = new ArrayList<CTriangle>();
		for(CTriangle t: triangles)	if(getInAlpha(t)<alpha) ret.add(t);
		return ret;
	}
	public List<CEdge> getEdges(double alpha){
		List<CEdge> ret = new ArrayList<CEdge>();
		for(CEdge t: edges)	if(getInAlpha(t)<alpha) ret.add(t);
		return ret;
	}
	public List<CVertex> getVertices(){
		return new ArrayList<CVertex>(vertices);
	}
	public List<Simplex> getSimplices(double alpha) {
		List<Simplex> ret = new ArrayList<Simplex>();
		for(Simplex s: simplices)	if(getInAlpha(s)<alpha) ret.add(s);
		return ret;
	}


	

	private double getInAlpha(Simplex s){
		return propertyMap.get(s).getInAlphaComplex();
	}
	private boolean getOnCH(Simplex s){
		SimplexAlphaProperties prop = propertyMap.get(s);
		switch(prop.getSimplexType()){
		case 0: return ((VertexAlphaProperties)prop).getOnConvexHull();
		case 1: return ((EdgeAlphaProperties)prop).getOnConvexHull();
		case 2: return ((TriangleAlphaProperties)prop).getOnConvexHull();
		case 3: throw new Error("Tetrahedrons are never completely on convex hull");
		}
		throw new Error("Undefined simplex type");
	}

	private void computeTetrahedraIntervals(){	//this method does the actual interval computation
		for(CTetrahedron t: del3d.getTetrahedra()){
			tetrahedra.add(t);
			
			double r = p.circumradius(t);
			TetrahedronAlphaProperties properties = new TetrahedronAlphaProperties(r);
			propertyMap.put(t, properties);
		}
		Collections.sort(tetrahedra,alphaOrdering);
	}

	public ArrayList<CTriangle> getCTriangles(){
		return triangles;
	}	

	private void computeTriangleIntervals(){
		for(CTriangle tri: del3d.getTriangles()){
			triangles.add(tri);
			
			boolean ch = tri.getNeighbour(0).containsBigPoint()||tri.getNeighbour(1).containsBigPoint();
			boolean att = 
				p.insphere(tri, tri.getNeighbour(0).oppositeVertex(tri))==SphereConfig.INSIDE ||
				p.insphere(tri, tri.getNeighbour(1).oppositeVertex(tri))==SphereConfig.INSIDE ;
			double minmu = triminmu(tri, ch);
			double maxmu = trimaxmu(tri, ch);
			double rho = p.circumradius(tri);
			TriangleAlphaProperties prop = new TriangleAlphaProperties(minmu, maxmu, rho, ch, att);
			propertyMap.put(tri, prop);
		}	
		Collections.sort(triangles,alphaOrdering);
	}

	private double triminmu(CTriangle tri, boolean ch){		//computes minmu 
		if(ch)	return getInAlpha(tri.getNeighbour(0));     //always put in place 0
		else	return Math.min(getInAlpha(tri.getNeighbour(0)), getInAlpha(tri.getNeighbour(1)));
	}
	private double trimaxmu(CTriangle tri, boolean ch){		//computes maxmu 
		if(ch)	return getInAlpha(tri.getNeighbour(0));     //always put in place 0
		else	return Math.max(getInAlpha(tri.getNeighbour(0)), getInAlpha(tri.getNeighbour(1)));
	}

	

	private void computeEdgeIntervals(){
		for(CEdge e: del3d.getEdges()){
			boolean ch = false;
			for(CTriangle t: e.getAdjacentTriangles()) ch|=getOnCH(t);
			boolean att = false;
			for(CTriangle t: e.getAdjacentTriangles()) 
				att|=p.edgeinsphere(e, t.oppositeVertex(e))==SphereConfig.INSIDE;
			double rho = p.edgecircumradius(e);
			double minmu = edgeminmu(e);
			double maxmu = edgemaxmu(e);
			EdgeAlphaProperties prop = new EdgeAlphaProperties(minmu,maxmu,rho,ch,att);
			propertyMap.put(e, prop);
		}
		Collections.sort(edges,alphaOrdering);
	}

	private double edgeminmu(CEdge e){		 		
		double min = p.edgecircumradius(e);
		for(CTriangle tri: e.getAdjacentTriangles()){
			TriangleAlphaProperties triProps = (TriangleAlphaProperties)propertyMap.get(tri);

			//TODO: The following line read:
			//min = Math.min(min,Math.min(triProps.getSingularInterval().getLeft(), triProps.getRegularInterval().getLeft()));
			//but that made no sense to me so i changed it
			if(triProps.isAttached()) 	min = Math.min(min, triProps.getRegularInterval().getLeft());
			else 						min = Math.min(min, triProps.getSingularInterval().getLeft());
		}
		return min;
	}

	private double edgemaxmu(CEdge e){
		double max =0;
		for(CTriangle tri: e.getAdjacentTriangles()){
			TriangleAlphaProperties triProps = (TriangleAlphaProperties)propertyMap.get(tri);
			if(!triProps.getOnConvexHull()) 	
				max = Math.max(max,triProps.getRegularInterval().getRight());
		}
		return max;
	}


	private void computeVertexIntervals(){
		for(CVertex v: vertices){
			VertexAlphaProperties prop = new VertexAlphaProperties(1,1, false);
			propertyMap.put(v, prop);
		}
	}


	
	private class AlphaComparator implements Comparator<Simplex>{
	    public int compare(Simplex s1, Simplex s2) {
	    	SimplexAlphaProperties p1 = propertyMap.get(s1);
	    	SimplexAlphaProperties p2 = propertyMap.get(s2); 
	    	
	        int c = Double.compare(p1.getInAlphaComplex(), p2.getInAlphaComplex());
	        if(c!=0) 	return c;
	        else		return p1.getSimplexType()-p2.getSimplexType();
	    }
	}
}
