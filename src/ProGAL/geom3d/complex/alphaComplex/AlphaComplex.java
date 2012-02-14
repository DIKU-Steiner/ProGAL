package ProGAL.geom3d.complex.alphaComplex;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Simplex;
import ProGAL.geom3d.complex.*;
import ProGAL.geom3d.complex.delaunayComplex.RegularComplex;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.proteins.PDBFile;

/**	<p>
 *  An alpha complex for a set of d-dimensional points and a real number alpha is a subset of the Delaunay complex 
 *  where all simplices that can be enclosed by an alpha-probe (a hypersphere of radius alpha), without the probe 
 *  enclosing any points, are removed. 
 *  </p>
 * @author R. Fonseca
 */
public class AlphaComplex extends AlphaFiltration implements SimplicialComplex{
	protected double alpha;

	private Map<Simplex, Integer> depthMap = null;

	/** 
	 * Build the alpha-complex of the specified point-list. Note that an entire Delaunay complex 
	 * is built as part of this constructor.
	 */
	public AlphaComplex(List<Point> pl, double alpha){
		super(pl);
		this.alpha = alpha;
	}

	/** Build the alpha-complex of the specified Delaunay complex. */
	public AlphaComplex(RegularComplex d3d, double alpha){
		super(d3d);
		this.alpha = alpha;
	}

	public double getAlpha() {
		return alpha;
	}
	
	public List<CEdge> getEdges() {
		return super.getEdges(alpha);
	}
	
	public List<CTriangle> getTriangles() {
		return getTriangles(alpha);
	}

	public List<CTetrahedron> getTetrahedra() {
		return super.getTetrahedra(alpha);
	}

	public List<Simplex> getSimplices() {
		return super.getSimplices(alpha);
	}
	public List<Simplex> getAllSimplices(){
		return super.getSimplices();
	}
	public List<CEdge> getAllEdges(){
		return super.edges;
	}

	/** 
	 * Return the depth of simplex from the surface of the alpha complex. Note that a depth 
	 * is returned for all simplices in the Delaunay complex. The depth of exposed simplices 
	 * not in the alpha complex is -1. The depth of all other k-simplices (including simplices 
	 * in voids) is the smallest depth of adjacent k-simplices plus one.
	 */
	public int getDepth(Simplex s){
		if(depthMap==null)
			calculateDepthsOld();
		return depthMap.get(s);
	}
	
	public void setAlpha(double alpha){
		this.alpha = alpha;
		calculateDepths();
	}

	private boolean isBoundaryTetrahedron(CTetrahedron t) {
		return  t.getNeighbour(0).containsBigPoint() || 
				t.getNeighbour(1).containsBigPoint() || 
				t.getNeighbour(2).containsBigPoint() || 
				t.getNeighbour(3).containsBigPoint();		
	}
	
	private void calculateDepths() {
		calculateDepthsInit();

		boolean update = true;
		int dist = 0;
		while (update) {
			update = false;
			for (CTetrahedron t: super.tetrahedra) {
				if (depthMap.get(t) == dist) {
					boolean inAC = getInAlpha(t) <= alpha;
					for (int i = 0; i < 4; i++) {
						CTetrahedron tx = t.getNeighbour(i);
						if (!tx.containsBigPoint()) {
							double txAlpha = getInAlpha(tx);
							if ((inAC &&  (txAlpha <= alpha)) || (!inAC && (txAlpha > alpha))) {
								if (depthMap.get(tx) == Integer.MAX_VALUE) {
									depthMap.put(tx, dist+1);
									update = true;
								}
							}
						}
					}
				}
			}
			dist++;
		}
		for (CTetrahedron t: super.tetrahedra) {
			System.out.println("depth = " + depthMap.get(t) + ", alpha = " + getInAlpha(t));
		}
	}
	
	private void calculateDepthsInit() {
		depthMap = new HashMap<Simplex,Integer>();
		for(CTetrahedron t: super.tetrahedra){
			depthMap.put(t, Integer.MAX_VALUE);
			for (int n = 0;n < 4; n++) {
				CTetrahedron neighbor = t.getNeighbour(n);
				if(neighbor.containsBigPoint()) {
					depthMap.put(neighbor, -1);
					depthMap.put(t, 0);
				}
			}
		}
	}
	
	private void calculateDepthsOld(){
		calculateDepthsInit();

		//Now iterate through all tetrahedra and update depths as long as they havent converged 
		//This could probably be done much more efficiently (e.g. single source shortest path) 
		boolean update = true;
		while(update){
			update = false;
			for(CTetrahedron t: super.tetrahedra){
				int oldDepth = depthMap.get(t);
				int newDepth = Math.min(
						Math.min(depth(t.getNeighbour(0)),depth(t.getNeighbour(1))), 
						Math.min(depth(t.getNeighbour(2)),depth(t.getNeighbour(3))) 
				); 
				if(getInAlpha(t)>alpha && newDepth==-1) 
					newDepth = -1;
				else if(newDepth!=Integer.MAX_VALUE)
					newDepth++;
				
				if(oldDepth!=newDepth){
					depthMap.put(t, newDepth);
					update = true;
				}
			}
		}
	}
	private int depth(Simplex s){
		Integer d = depthMap.get(s);
		if(d==null) return Integer.MAX_VALUE;
		return d;
	}
	
	public void getVoids(J3DScene scene) {
		for(CTetrahedron tetr : super.tetrahedra) { 
			if (depthMap.get(tetr) == Integer.MAX_VALUE) {
				System.out.println("MAX_VALUE");
				scene.addShape(tetr, java.awt.Color.YELLOW);
			}
		}
	}
	
	public CTetrahedron getDeepestCavityTetrahedron() {
		double alpha = getAlpha();
		CTetrahedron deep = null;
		int deepDepth = -1;
		for(CTetrahedron t: super.tetrahedra) {
			int depth = getDepth(t);
			if ((getInAlpha(t) > alpha) && (depth > deepDepth) && (depth != Integer.MAX_VALUE)) {
				deepDepth = depth;
				deep = t;
			}
		}
		return deep;
	}
	
	public ArrayList<CTetrahedron> getAllDeepestCavityTetrahedra(int depthBound) {
		double alpha = getAlpha();
		ArrayList<CTetrahedron> tetras = new ArrayList<CTetrahedron>();

		for (CTetrahedron t: super.tetrahedra) {
			double tDepth = getDepth(t);
			CTetrahedron n0 = t.getNeighbour(0);
			CTetrahedron n1 = t.getNeighbour(1);
			CTetrahedron n2 = t.getNeighbour(2);
			CTetrahedron n3 = t.getNeighbour(3);
			if ((tDepth > depthBound) && (getInAlpha(t) > alpha) && (tDepth != Integer.MAX_VALUE) &&
				(n0.containsBigPoint() || (getInAlpha(n0) <= alpha) || (getDepth(n0) <= tDepth)) &&
				(n1.containsBigPoint() || (getInAlpha(n1) <= alpha) || (getDepth(n1) <= tDepth)) &&
				(n2.containsBigPoint() || (getInAlpha(n2) <= alpha) || (getDepth(n2) <= tDepth)) &&
				(n3.containsBigPoint() || (getInAlpha(n3) <= alpha) || (getDepth(n3) <= tDepth))) {
				tetras.add(t);
				System.out.println("Deepest tetrahedron found. Depth = " + getDepth(t) + ", alpha = " + getInAlpha(t));
			}
		}
		return tetras;	
	}
	
	public void getCavity(CTetrahedron t, int lowerBound, J3DScene scene) {
		double alpha = getAlpha();
		Stack<CTetrahedron> stack = new Stack<CTetrahedron>();
		if (getDepth(t) > lowerBound) stack.push(t);
		scene.addShape(t, java.awt.Color.YELLOW);
		t.setModified(true);
		while (!stack.isEmpty()) {
			CTetrahedron ct = stack.pop();
			scene.addShape(ct, java.awt.Color.RED);
			System.out.println("Tetrahedron of depth " + getDepth(ct) + " added to the cavity. Its alpha is " + getInAlpha(ct));
			for (int i = 0; i < 4; i++) {
				CTetrahedron nt = ct.getNeighbour(i);
				if (!nt.isModified() && (getInAlpha(nt) > alpha) && (getDepth(nt) > lowerBound)) {
					stack.push(nt);
					nt.setModified(true);
				}
			}
			
		}
	}
	
	public void getAllCavities(ArrayList<CTetrahedron> tetras, int lowerBound, J3DScene scene) {
		for (CTetrahedron t: tetras) getCavity(t, lowerBound, scene);
	}
	
	public void getAllCavityPaths(ArrayList<CTetrahedron> tetras, J3DScene scene) {
		for (CTetrahedron t: tetras) getCavityPath(t, scene);
	}
	
	public void getCavityPath(CTetrahedron deep, J3DScene scene) {
		double alpha = getAlpha();
		scene.addShape(deep,java.awt.Color.YELLOW);
		deep.setModified(true);
		int depth = getDepth(deep);
		while (depth != 0) {
			CTetrahedron t = deep.getNeighbour(0);
			if ((getDepth(t) == depth-1)  && (getInAlpha(t) > alpha)) deep = t;
			else {
				t = deep.getNeighbour(1);
				if ((getDepth(t) == depth-1)  && (getInAlpha(t) > alpha)) deep = t;
				else {
					t = deep.getNeighbour(2);
					if ((getDepth(t) == depth-1)   && (getInAlpha(t) > alpha)) deep = t;
					else {
						t = deep.getNeighbour(3);
						if ((getDepth(t) == depth-1)  && (getInAlpha(t) > alpha)) deep = t;
					}
				}
			}
			if (t.isModified()) depth = 0;
			else {
				System.out.println("red tetrahedron with depth = " + getDepth(deep) + " and alpha = " + getInAlpha(deep));
				scene.addShape(deep,java.awt.Color.RED);
				deep.setModified(true);
				depth--;
			}
		}
		
	}
	
	public static void main(String[] args) {
//		String[] pdbs = new String[]{"1X5RA","1X0OA","1XDXA","1AKPA","1Y6DA", "1CTFA"};
//		List<Point> points = new PDBFile("/Users/pawel/BioRepo/MotherOfAllProjects/pdb_files/1F94.pdb").getAtomCoords();
		List<Point> points = new PDBFile("/Users/rfonseca/Downloads/3SQF.pdb").getAtomCoords();
//		ProteinComplex pc = new ProteinComplex(f);
//		List<Point> points = new PDBFile(PDBWebReader.downloadPDBFile("1XDXA")).getAtomCoords();
		AlphaComplex ac = new AlphaComplex(points, 2.8);
		
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		for (CTetrahedron tetr : ac.getTetrahedra(0, 2.8)) {
			scene.addShape(tetr, new java.awt.Color(0,0,255,10));
		}

//		ac.calculateDepths();
//		ac.getVoids(scene);
//		CTetrahedron t = ac.getDeepestCavityTetrahedron();
//		ac.getCavityPath(t, scene);
//		ac.getAllCavityPaths(ac.getAllDeepestCavityTetrahedra(9), scene);
//		ac.getAllCavities(ac.getAllDeepestCavityTetrahedra(9), 5, scene);
//		for(CTetrahedron tetr : ac.getTetrahedra(2.8, 100)) { 
//			if (ac.getDepth(tetr) > -1) 	scene.addShape(tetr, java.awt.Color.YELLOW);
//		}

	
/*		if (true) return;
		int depth  = ac.getDepth(deepestTetra);
		Stack<CTetrahedron> tetras = new Stack<CTetrahedron>(); 
		tetras.push(deepestTetra);
		deepestTetra.setModified(true);
		int counter = 0;
		System.out.println("highest depth =" + depth);
		while (!tetras.empty()) {
			CTetrahedron tetra = tetras.pop();
			counter++;
			scene.addShape(tetra, java.awt.Color.YELLOW);
			tetra.setModified(true);
			System.out.println(counter + ". depth = " + ac.getDepth(tetra) + " alpha = " + ac.getInAlpha(tetra));
			for (int i = 0; i < 4; i++) {
				CTetrahedron nextTetra = tetra.getNeighbour(i);
				if (!nextTetra.isModified() && (ac.getDepth(nextTetra) > -2) && (ac.getInAlpha(nextTetra) > 2.8)) {
					tetras.push(nextTetra);
					nextTetra.setModified(true);
				}
			}
		}
		
		if (true) return;
		
		CTetrahedron currentTetra = deepestTetra;
		int currentDepth = depth;
		while ((currentTetra != null) && (ac.getDepth(currentTetra) >= 0)) {
			scene.addShape(currentTetra, Color.RED);
			System.out.println("current depth = " + currentDepth);
			int i = 0;
			CTetrahedron nextTetra = null;
			while (i < 4) {
				nextTetra = currentTetra.getNeighbour(i);
				System.out.println("  next tetra depth = " + ac.getDepth(nextTetra) + " alpha = " + ac.getInAlpha(nextTetra));
				if ((currentDepth - ac.getDepth(nextTetra) <= 1) && (ac.getInAlpha(nextTetra) > 2.8)) {
					currentTetra = nextTetra;
					currentDepth--;
					i = 4;
				} 
				else {
					nextTetra = null;
					i++;
				}
			}
			if (nextTetra == null) currentTetra = null;
		}
*/
	}
	
}
