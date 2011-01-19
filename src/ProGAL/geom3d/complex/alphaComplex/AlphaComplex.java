package ProGAL.geom3d.complex.alphaComplex;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Simplex;
import ProGAL.geom3d.complex.*;
import ProGAL.geom3d.complex.delaunayComplex.DelaunayComplex;

/**	<p>
 *  An alpha complex for a set of d-dimensional points and a real number alpha is a subset of the Delaunay complex 
 *  where all simplices that can be enclosed by an alpha-probe (a hypersphere of radius alpha), without the probe 
 *  enclosing any points, are removed. 
 *  </p>
 * @author R. Fonseca
 */
public class AlphaComplex extends AlphaFiltration implements SimplicialComplex{
	private double alpha;

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
	public AlphaComplex(DelaunayComplex d3d, double alpha){
		super(d3d);
		this.alpha = alpha;
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

	/** 
	 * Return the depth of simplex from the surface of the alpha complex. Note that a depth is 
	 * is returned for all simplices in the Delaunay complex. The depth of exposed simplices 
	 * not in the alpha complex is -1. The depth of all other k-simplices (including simplices 
	 * in voids) is the smallest depth of adjacent k-simplices plus one.
	 */
	public int getDepth(Simplex s){
		if(depthMap==null)
			calculateDepths();
		return depthMap.get(s);
	}

	private void calculateDepths(){
		depthMap = new HashMap<Simplex,Integer>();

		//Set the depth of all tetrahedra touching the big-points to -1 
		for(CTetrahedron t: super.tetrahedra){
			for(int n=0;n<4;n++) {
				CTetrahedron neighbor = t.getNeighbour(n);
				if(neighbor.containsBigPoint()) {
					depthMap.put(neighbor, -1);
				}
			}
		}

		//Now iterate through all tetrahedra and update depths as long as they havent converged 
		//This could probably be done much more efficiently (e.g. shortest path) 
		boolean update = true;
		while(update){
			update = false;
			for(CTetrahedron t: super.tetrahedra){
				int oldDepth = depth(t);
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
}
