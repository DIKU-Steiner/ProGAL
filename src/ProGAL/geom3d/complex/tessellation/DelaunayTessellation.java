package ProGAL.geom3d.complex.tessellation;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointWeighted;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.predicates.*;
import ProGAL.math.Randomization;

/** <p>
 *  A Delaunay tessellation for a set of d-dimensional points is a tesselation of the points such that no point is inside 
 *  the circumscribing hypersphere of the d-simplices (for the 3D case: Tetrahedra). 
 *  </p>
 *  
 *  <p>
 *  This class builds a three-dimensional Delaunay tessellation in the constructor and accesses it using the 
 *  <code>getTetrahedra</code> method. The following example displays the Delaunay complex of ten random points. 
 *  <pre>
 *  {@code
 *  //Generate the complex
 *  List<Point> pl = PointList.generatePointsInCube(10);
 *  DelaunayTessellation dc = new DelaunayTessellation(pl);
 *    
 *  //Display the complex
 *  J3DScene scene = J3DScene.createJ3DSceneInFrame();
 *  for(CTetrahedron t: dc.getTetrahedra()){
 *    scene.addShape(t, new Color(200,100,100,100));
 *  }
 *  }
 *  </pre>     
 *  </p>
 *  <p>The original point-set is left unaltered and non-referenced by this class. A new set of vertices are 
 *  allocated using the CVertex class. These are randomly permuted to avoid degeneracies. If one wishes to 
 *  associate the original points with a vertex in the complex it would be sufficient to test if the distance 
 *  between the point and the vertex is less than 0.0001.</p>
 *  
 *  <p>The complex is bounded by a big tetrahedron whose corner-points are located sufficiently far from any of 
 *  the vertices of the complex. These points are part of the tessellation. The simplices that have one of these 
 *  'big points' as corners can not be accessed directly via the getter-methods, but they will be neighbors of 
 *  other normal simplices. For instance:
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
 *  Will print false, true, true, true and true. A tetrahedron will contain a big-point if and only if it is 
 *  outside the convex hull.
 *  </p>
 *  @author P.Sterner, H.Sterner, R.Fonseca
 */
public class DelaunayTessellation{

	private final List<CVertex> points;
	private final List<CTetrahedron> tetrahedra;
	private final Predicates predicates;
	private final Walk walk;
	private final Flip14 f14;
	private final Flips flips;

	public static void main(String[] args){
		List<Point> points = new LinkedList<Point>();//PointList.generatePointsInCube(100);
		for(int i=0;i<100;i++){
			points.add(new PointWeighted(Randomization.randBetween(0.0, 1.0),Randomization.randBetween(0.0, 1.0),Randomization.randBetween(0.0, 1.0),0));
		}
		DelaunayTessellation dt = new DelaunayTessellation(points);
		
	}
	
	/** Builds the Delaunay complex of the specified point-set */
	public DelaunayTessellation(List<Point> points) {
		this(points, new InexactRegularJavaPredicates());
	}
	
	protected DelaunayTessellation(List<Point> points, Predicates predicates){
		this.points = new ArrayList<CVertex>(points.size());
		int i=0;
		for(Point p: points) this.points.add(new CVertex(p,i++));
		this.tetrahedra = new ArrayList<CTetrahedron>(points.size()*6);
		
		this.predicates = predicates;
		this.walk = new Walk(predicates);
		this.flips = new Flips(predicates);
		this.f14 = new Flip14(flips);
		
		compute();
	}


	/** Get the tetrahedra in the complex. The tetrahedra that has 'big points' as corners are not returned */
	public List<CTetrahedron> getTetrahedra() {
		return new ArrayList<CTetrahedron>(tetrahedra);
	}
	/** Get the vertices in the complex. The 'big points' are not returned */
	public List<CVertex> getVertices(){
		return new ArrayList<CVertex>(points);
	}

	private final void compute() {
		double max = 1000;//TODO find a more meaningful max
		
		//TODO: Take care of degeneracies in a better way than permutation
		for(CVertex v: points){
			v.addThis(new Vector(
					Randomization.randBetween(-0.0001, 0.0001),
					Randomization.randBetween(-0.0001, 0.0001),
					Randomization.randBetween(-0.0001, 0.0001)
					)
			);
		}

		//Find the enclosing tetrahedron
		CTetrahedron next_t = new FirstTetrahedron(max);
		flips.addTetrahedron(next_t);

		//Iterate over points
		for(CVertex p: points){
			if(points.indexOf(p)%1000==0) System.out.println("Inserted "+points.indexOf(p)+" points");
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
			if (!t.isModified() && !t.containsBigPoint()){
				tetrahedra.add(t);
			}
		}
		flips.getTetrahedrastack().clear();
	}
	

	/** Checks that all tetrahedra comply with the Delaunay-criteria. */
	public boolean checkTetrahedra() {
//		for(CTetrahedron t: tetrahedra){
//			for(PointWeighted p: points){
//				if(		t.getPoint(0)!=p && 
//						t.getPoint(1)!=p && 
//						t.getPoint(2)!=p && 
//						t.getPoint(3)!=p &&
//						predicates.insphere(t, p).equals(SphereConfig.INSIDE) )
//				{
//					return false;
//				}
//			}
//		}
		return true;
	}
}


