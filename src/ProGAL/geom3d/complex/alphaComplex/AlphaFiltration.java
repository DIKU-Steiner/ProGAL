package ProGAL.geom3d.complex.alphaComplex;

import java.util.ArrayList; 
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashMap;

import ProGAL.dataStructures.UnionFind;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Simplex;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.complex.*;
import ProGAL.geom3d.complex.delaunayComplex.DelaunayComplex;
import ProGAL.geom3d.predicates.*;
import ProGAL.geom3d.predicates.Predicates.SphereConfig;
import ProGAL.geom3d.volumes.Tetrahedron;

/**	<p>
 *  An alpha complex for a set of d-dimensional points and a real number alpha is a subset of the Delaunay complex 
 *  where all simplices that can be enclosed by an alpha-probe (a hypersphere of radius alpha), without the probe 
 *  enclosing any points, are removed. The alpha-filtration is a generalization of the alpha complex that associates
 *  with each simplex the alpha value at which it enters the alpha-complex. 
 *  </p>
 *  
 *  <p>
 *  This class builds the three-dimensional alpha filtration. The alpha complex is easily extracted by using the 
 *  <code>List<Simplex> getSimplices(double alpha)</code> method. For convenience there are also methods to retrieve 
 *  only the edges, triangles and tetrahedra of the alpha complex. The alpha complex is built in the constructor of 
 *  <code>AlphaComplex</code>. Lists of simplices are always returned in non-decreasing order of alpha.
 *  
 *  
 *  The following example displays the alpha-complex of 20 random points using an alpha probe of radius 0.2.   
 *  <pre>
 *  {@code
 *  //Generate the filtration
 *  List<Point> pl = PointList.generatePointsInCube(20);
 *  AlphaFiltration af = new AlphaFiltration(pl);
 *  
 *  //Display the alpha complex
 *  J3DScene scene = J3DScene.createJ3DSceneInFrame();
 *  for(CTetrahedron t: af.getTetrahedra(0.2)){
 *     scene.addShape(t, new Color(200,100,100,100));
 *  }
 *  }
 *  </pre>     
 *  </p>
 *  
 *  <p>
 *  It should be noted that the simplices returned by this class are all connected as if they were in a 
 *  <code>DelaunayComplex</code> and if one wishes e.g. to find a component in the alpha-complex then the 
 *  breadth-first-search should be modified to only traverse simplices that enter the alpha complex 'before' the 
 *  desired alpha-value. To determine at which time a simplex enters the alpha complex the  
 *  <code>double getInAlpha(Simplex s)</code>-method can be used. If the result is less than the probes radius then 
 *  the simplex is in the complex.
 *  
 *  <it>TODO: Describe the use of SimplexAlphaProperties and write a good example</it>
 *  </p>
 * @author R. Fonseca
 */
public class AlphaFiltration {

	private final DelaunayComplex del3d;
	private Predicates p = new ExactJavaPredicates();
	private final AlphaComparator alphaOrdering = new AlphaComparator();
	private final Map<Simplex, SimplexAlphaProperties> propertyMap = new HashMap<Simplex, SimplexAlphaProperties>();
	private int[][] bettiNumbers = null;
	
	protected ArrayList<CTetrahedron> tetrahedra = new ArrayList<CTetrahedron>();
	protected ArrayList<CTriangle> triangles = new ArrayList<CTriangle>();
	protected ArrayList<CEdge> edges = new ArrayList<CEdge>();
	protected ArrayList<CVertex> vertices = new ArrayList<CVertex>();
	protected ArrayList<Simplex> simplices = new ArrayList<Simplex>();
	
	/** 
	 * Build the alpha-filtration of the specified point-list. Note that an entire Delaunay complex 
	 * is built as part of this constructor.
	 */
	public AlphaFiltration(List<Point> pl){
		this(new DelaunayComplex(pl));
	}

	/** Build the alpha-filtration of the specified Delaunay complex. */
	public AlphaFiltration(DelaunayComplex d3d){
		this.del3d = d3d;
		compute();
	}

	private void compute(){
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


	public List<Simplex> getSimplices() {		return simplices;	}
	public List<CTetrahedron> getTetrahedra() {	return tetrahedra;	}
	public List<CTriangle> getTriangles() {		return triangles;	}
	public List<CEdge> getEdges() {				return edges;	}
	
	/** Get a list of tetrahedra that are part of the alpha-complex with the specified probe radius. */
	public List<CTetrahedron> getTetrahedra(double alpha){
		List<CTetrahedron> ret = new ArrayList<CTetrahedron>();
		for(CTetrahedron t: tetrahedra)	
			if (getInAlpha(t) < alpha) ret.add(t);
		return ret;
	}
	
	/** Get a list of tetrahedra that are part of the alpha-complex with the probe radius in a specified interval.  Added by PW*/
	public List<CTetrahedron> getTetrahedra(double alphaLow, double alphaHigh){
		List<CTetrahedron> ret = new ArrayList<CTetrahedron>();
		for(CTetrahedron t: tetrahedra)	
			if ((getInAlpha(t)>=alphaLow) && (getInAlpha(t)<alphaHigh)) ret.add(t);
		return ret;
	}

	/** Get a list of triangles that are part of the alpha-complex with the specified probe radius. */
	public List<CTriangle> getTriangles(double alpha){
		List<CTriangle> ret = new ArrayList<CTriangle>();
		for(CTriangle t: triangles)	if(getInAlpha(t)<alpha) ret.add(t);
		return ret;
	}
	
	/** Returns triangles in one tetrahedron only */
	public List<Triangle> getSurfaceTriangles(double alpha) {
		List<CTetrahedron> tetrahedra = getTetrahedra(alpha);
		List<Triangle> triangles = new ArrayList<Triangle>();
		for (CTetrahedron tetrahedron : tetrahedra) {
			for (int i = 0; i < 4; i++) {
				CTetrahedron neighbor = tetrahedron.getNeighbour(i);
				if ((neighbor == null) || (getInAlpha(neighbor) >= alpha)) triangles.add(tetrahedron.getTriangle(i));
			}
		}
		return triangles;
	}

	/** Get a list of edges that are part of the alpha-complex with the specified probe radius. */
	public List<CEdge> getEdges(double alpha){
		List<CEdge> ret = new ArrayList<CEdge>();
		for(CEdge t: edges)	if(getInAlpha(t)<alpha) ret.add(t);
		return ret;
	}

	/** Get a list of the vertices of the complex. */
	public List<CVertex> getVertices(){
		return new ArrayList<CVertex>(vertices);
	}

	/** Get a list of simplices that are part of the alpha-complex with the specified probe radius. */
	public List<Simplex> getSimplices(double alpha) {
		List<Simplex> ret = new ArrayList<Simplex>();
		for(Simplex s: simplices)	if(getInAlpha(s)<alpha) ret.add(s);
		return ret;
	}

	public List<CTriangle> getAlphaShape(double alpha) {
		List<CTriangle> ret = new ArrayList<CTriangle>();
		for (CTriangle t: getTriangles(alpha)) {
			try{
				double a0 = getInAlpha(t.getAdjacentTetrahedron(0));
				double a1 = getInAlpha(t.getAdjacentTetrahedron(1));
				if( (a0>alpha)^(a1>alpha) )
					ret.add(t);
			}catch(NullPointerException exc){
				ret.add(t);
			}
		}
		return ret;
	}

	public int getDim(Simplex s){
		SimplexAlphaProperties prop = propertyMap.get(s);
		return prop.getSimplexType();
	}
	
	/** Return the probe-radius at which the simplex <code>s</code> enters the alpha complex. */
	public double getInAlpha(Simplex s){
		return propertyMap.get(s).getInAlphaComplex();
	}
	
	public boolean getAttached(Simplex s){
		return propertyMap.get(s).isAttached();
	}
	
	/** 
	 * Return true iff the simplex is on the convex hull. Calling this method with a CTetrahedra 
	 * will throw an error. 
	 */
	public boolean getOnCH(Simplex s){
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

	private void computeTriangleIntervals(){
		for(CTriangle tri: del3d.getTriangles()){
			triangles.add(tri);
			
			boolean ch = tri.getAdjacentTetrahedron(0).containsBigPoint()||tri.getAdjacentTetrahedron(1).containsBigPoint();
			boolean att = 
				p.insphere(tri, tri.getAdjacentTetrahedron(0).oppositeVertex(tri))==SphereConfig.INSIDE ||
				p.insphere(tri, tri.getAdjacentTetrahedron(1).oppositeVertex(tri))==SphereConfig.INSIDE ;
			double minmu = triminmu(tri, ch);
			double maxmu = trimaxmu(tri, ch);
			double rho = p.circumradius(tri);
			TriangleAlphaProperties prop = new TriangleAlphaProperties(minmu, maxmu, rho, ch, att);
			propertyMap.put(tri, prop);
		}	
		Collections.sort(triangles,alphaOrdering);
	}

	private double triminmu(CTriangle tri, boolean ch){		//computes minmu 
		if(tri.getAdjacentTetrahedron(0).containsBigPoint())	return getInAlpha(tri.getAdjacentTetrahedron(1));
		if(tri.getAdjacentTetrahedron(1).containsBigPoint())	return getInAlpha(tri.getAdjacentTetrahedron(0));
		else	return Math.min(getInAlpha(tri.getAdjacentTetrahedron(0)), getInAlpha(tri.getAdjacentTetrahedron(1)));
	}
	private double trimaxmu(CTriangle tri, boolean ch){		//computes maxmu 
		if(tri.getAdjacentTetrahedron(0).containsBigPoint())	return getInAlpha(tri.getAdjacentTetrahedron(1));
		if(tri.getAdjacentTetrahedron(1).containsBigPoint())	return getInAlpha(tri.getAdjacentTetrahedron(0));
//		if(ch)	return getInAlpha(tri.getNeighbour(0));     //always put in place 0
		else	return Math.max(getInAlpha(tri.getAdjacentTetrahedron(0)), getInAlpha(tri.getAdjacentTetrahedron(1)));
	}

	

	private void computeEdgeIntervals(){
		for(CEdge e: del3d.getEdges()){
			edges.add(e);
			
			boolean ch = false;
			for(CTriangle t: e.getAdjacentTriangles()) ch |= getOnCH(t);
			boolean att = false;
			for(CTriangle t: e.getAdjacentTriangles()) 
				att |= p.edgeinsphere(e, t.oppositeVertex(e))==SphereConfig.INSIDE;
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
		for(CVertex v: del3d.getVertices()){
			vertices.add(v);
			
			VertexAlphaProperties prop = new VertexAlphaProperties(1,1, false);
			propertyMap.put(v, prop);
		}
	}
	
	
	/** 
	 * The vertex-hull of v is the set of all tetrahedrons that has v as a corner-point. Notice that 
	 * that this method operates on the entire Delaunay complex and is not affected by any probe size. 
	 */
	public Set<CTetrahedron> getVertexHull(CVertex v){
		return del3d.getVertexHull(v);
	}
	

	/** 
	 * Get the Betti-numbers of the alpha-filtration. The implementation follows the specification  
	 * of Delfinado and Edelsbrunner 1995. This algorithm runs in O(n alpha(n)) time and O(n) space 
	 * where alpha(n) is the inverse Ackerman function.
	 * 
	 * <b>Reference</b>: An incremental algorithm for Betti numbers of simplicial 
	 * complexes on the 3-sphere. C.J.A. Delfinado and H. Edelsbrunner. Computer Aided Geometric 
	 * Design, vol. 12, p. 771-784, 1995. 
	 * 
	 * @return an array of int-arrays. The zero'th entry holds the zero'th Betti numbers for each 
	 * simplex, the next entry the first Betti numbers etc. The fifth entry holds 1 if the 
	 * corresponding simplex was marked and 0 otherwise. The length of each entry is n, where n 
	 * is the number of simplices in the AlphaFiltration.  
	 */
	public int[][] getBettiNumbers(){
		if(bettiNumbers==null){
			int n = simplices.size();
			//The first phase of the algorithm marks every simplex that belongs to a cycle of the same 
			//dimension in K_i
			boolean[] marked = new boolean[n];
			
			//Each vertex belongs to a 0-cycle, so all vertices get marked.
			for(CVertex v: vertices) {
				int i = simplices.indexOf(v);
				marked[i] = true;
			}
			
			//To mark the appropriate edges we process the simplices in forward direction and maintain a 
			//union-find structure for K_i^(1). An edge is marked iff it does not cause a UNION operation
			DisjointSet ds = new DisjointSet();
			for(CVertex v: vertices) ds.makeSet(v); 
			for(CEdge e: edges) {
				int i = simplices.indexOf(e);
				if(  ds.find(e.getA())==ds.find(e.getB())  )
					marked[i] = true;
				else
					ds.union(e.getA(), e.getB());
			}
			
			//For marking the appropriate triangles we process the simplices in backward direction, from 
			//sigma_n down to sigma_1. A union-find structure representing the dual graph, G_i of  
			//K_i^tilda = T - K_i is maintained, and a triangle is marked iff it causes a UNION operation.
			ds = new DisjointSet();
			UnionFind<CTetrahedron> uf = new UnionFind<CTetrahedron>();
			CTetrahedron bigTet = new CTetrahedron(null,null,null,null);
			for(int iT=triangles.size()-1;iT>=0;iT--){
				CTriangle t = triangles.get(iT);
				CTetrahedron n0 = t.getAdjacentTetrahedron(0);
				CTetrahedron n1 = t.getAdjacentTetrahedron(1);
				if(n0.containsBigPoint()) n0 = bigTet;
				if(n1.containsBigPoint()) n1 = bigTet;
				CTetrahedron s0 = uf.find(n0);
				CTetrahedron s1 = uf.find(n1);
				
				if(  s0!=s1  ){
					marked[simplices.indexOf(t)] = true;
					uf.union(n0, n1);
				}
			}
			
			//The only tetrahedron that belongs to a 3-cycle at the time it is processed is sigma_n. 
			//This is the only tetrahedron that gets marked.
			int lastTetIdx = simplices.indexOf(tetrahedra.get(tetrahedra.size()-1));
			marked[lastTetIdx] = true;
			
			//The second phase counts the marked and unmarked simplices and derives the Betti numbers as 
			//simple sums of these numbers.
			bettiNumbers = new int[6][n];
			int[] b = new int[4];
			for(int i=0;i<n;i++){
				Simplex s = simplices.get(i);
				int k = s.getDimension();//getDim(s);
				if(marked[i]) 	b[k]++;
				else			b[k-1]--;

				bettiNumbers[0][i] = b[0];
				bettiNumbers[1][i] = b[1];
				bettiNumbers[2][i] = b[2];
				bettiNumbers[3][i] = b[3];
				bettiNumbers[4][i] = marked[i]?1:0;//For testing purposes the marks are stored as well.
				bettiNumbers[5][i] = k;
			}
		}
		return bettiNumbers;
	}
	
	/**
	 * TODO: Comment
	 */
	public int[][][] getBettiPersistence(){
		List<int[]>[] pairs = getPairSimplices();
		int n = simplices.size();
		int maxPersistence = 0;
		int[][][] ret = new int[3][n][n];
		for(int d=0;d<3;d++){
			for(int[] pair: pairs[d]){
				int persistence = pair[1]-pair[0];
				maxPersistence = Math.max(maxPersistence, persistence);
				for(int i=pair[0];i<pair[1];i++){
					for(int p=0;p<persistence;p++){
						ret[d][i][p]++;
					}
					persistence--;
				}
				
			}
		}
		
		int[][][] compressedRet = new int[3][maxPersistence+1][n];
		for(int d=0;d<3;d++){
			for(int i=0;i<n;i++){
				for(int p=0;p<maxPersistence+1;p++){
					compressedRet[d][i][p] = ret[d][i][p];
				}
			}
		}
		return compressedRet;
	}
	
	/**
	 * TODO: Comment
	 */
	@SuppressWarnings("unchecked")
	public List<int[]>[] getPairSimplices(){
		List<int[]>[] ret = new List[3];
		ret[0] = new LinkedList<int[]>();
		ret[1] = new LinkedList<int[]>();
		ret[2] = new LinkedList<int[]>();

		int n = simplices.size();
		int[][] bettiNumbers = getBettiNumbers();
		List<Integer>[] T = new List[n];
		
		for(int j=0;j<n;j++){
			int k = getDim( simplices.get(j) );
			if( bettiNumbers[4][j]==0 ){
				int i=youngest(j, bettiNumbers[4], T);
				ret[k-1].add( new int[]{i, j} );
			}
		}
		
		return ret;
	}
	
	private int youngest(int j, int[] marked, List<Integer>[] T){
		List<Integer> Lambda = positiveD(j,marked);
		
		while(true){
			int i=Lambda.get(Lambda.size()-1);
			if(T[i]==null){
				T[i] = new LinkedList<Integer>();
				T[i].add(j);
				T[i].addAll(Lambda);
				return i;
			}
			Lambda = addLists(  Lambda , T[i].subList(1, T[i].size())  );
		}
	}
	
	private List<Integer> positiveD(int j, int[] marked){
		List<Integer> ret = new LinkedList<Integer>();
		switch(getDim(simplices.get(j))){
		case 0: return ret;
		case 1:
			CEdge e = (CEdge)simplices.get(j);
			addPositive(simplices.indexOf(e.getA()), marked, ret);
			addPositive(simplices.indexOf(e.getB()), marked, ret);
			break;
		case 2:
			CTriangle t = (CTriangle)simplices.get(j);
			addPositive(simplices.indexOf(t.getEdge(0)), marked, ret);
			addPositive(simplices.indexOf(t.getEdge(1)), marked, ret);
			addPositive(simplices.indexOf(t.getEdge(2)), marked, ret);
			break;
		case 3:
			CTetrahedron tet = (CTetrahedron)simplices.get(j);
			addPositive(simplices.indexOf(tet.getTriangle(0)), marked, ret);
			addPositive(simplices.indexOf(tet.getTriangle(1)), marked, ret);
			addPositive(simplices.indexOf(tet.getTriangle(2)), marked, ret);
			addPositive(simplices.indexOf(tet.getTriangle(3)), marked, ret);
			break;
		}
		Collections.sort(ret);
		return ret;
	}
	private static void addPositive(int i, int[] marked, List<Integer> ret){
		if(marked[i]==1) ret.add(i);
	}
	private static List<Integer> addLists(List<Integer> a, List<Integer> b){
		List<Integer> ret = new LinkedList<Integer>();
		for(Integer i: a)
			if(!b.contains(i)) ret.add(i);
		for(Integer i: b)
			if(!a.contains(i)) ret.add(i);
		Collections.sort(ret);
		return ret;
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
