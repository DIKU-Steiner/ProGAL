package ProGAL.geomNd.complex.delaunayComplex;

import java.util.ArrayList;
import java.util.List;

import ProGAL.geomNd.Point;
import ProGAL.geomNd.Vector;
import ProGAL.geomNd.complex.Vertex;
import ProGAL.geomNd.complex.Tessel;
import ProGAL.math.Matrix;
import ProGAL.math.Randomization;

/** <p>
 *  A Delaunay complex for a set of d-dimensional points is a tesselation of the points such that no point is inside 
 *  the circumscribing hypersphere of the d-simplices (for the 3D case: Tetrahedra). 
 *  
 *  This class builds a d-dimensional Delaunay complex in the constructor and enables access to it using the public 
 *  methods <code>getTessels</code> and <code>getAllTessels</code>. </p>
 *  
 *  <p>The original point-set is left unaltered and non-referenced by this class. A new set of vertices is 
 *  allocated using the CVertex class. The position of these are randomly perturbed to avoid degeneracies. If one wishes to 
 *  associate the original points with a vertex in the complex it would be sufficient to test if the distance 
 *  between the point and the vertex is less than 0.0001.</p>
 *  
 *  <p>The complex is bounded by a big tetrahedron whose corner-points are located sufficiently far from any of 
 *  the vertices of the complex. The simplices that have one of these 'big points' as corners can not be accessed 
 *  directly via the <code>tessels</code> field, but they will be neighbors of other normal simplices.</p>
 *  
 *  @author Desiree M. S. Joergensen
 *  @author Annie J. Pinder
 */
public class BowyerWatson{
	/** The vertices of the tessellation. These are slightly perturbed compared to the originally added points. */
	private final List<Vertex> points = new ArrayList<Vertex>();
	
	/** The tessels of the tessellation. Includes tessels with 'big points' */
	private final List<Tessel> tessels = new ArrayList<Tessel>();
	
	private List<Tessel> newTets = new ArrayList<Tessel>(); // New tessels created after retesselation
	
	/** The dimension of the tessellation */
	private final int dimension;

	/** Builds the Delaunay complex of the specified point-set using Bowyer-Watson algorithm */
	public BowyerWatson(List<Point> points) {
		this.dimension= points.get(0).getCoords().length;
		for (Point p: points) {
			this.points.add(new Vertex(p));
			if (p.getCoords().length!=dimension) {
				throw new RuntimeException("Mismatch in dimensions of points");
			}
		}
		compute();
	}

	/** Get the tesselts in the complex. The tessels that has 'big points' as corners are not returned */
	public List<Tessel> getTessels() {
		List<Tessel> tess = new ArrayList<Tessel>();
		for (int i=0;i<tessels.size();i++) {
			if (!(tessels.get(i).containsBigPoint())) {
				tess.add(tessels.get(i));
			}
		}
		return tess;
	}


	/** Returns all tessels (including tessels with 'big points') */
	public List<Tessel> getAllTessels() {
		return tessels;
	}

	/** The dimension of the tessellation */
	public int getDim() {
		return dimension;
	}

	/** Get the vertices in the complex. The 'big points' are not returned */
	public List<Vertex> getVertices(){
		return new ArrayList<Vertex>(points);
	}

	/** Get the i'th vertex */
	public Vertex getVertex(int i) { return points.get(i); }

	/** Constructs the corners of a regular simplex all with distance 1 to origo */ 
	private Vertex[] regularSimplex(){
		double angle = -1.0/dimension;
		Vertex[] ret = new Vertex[dimension+1];

		for(int s=0;s<dimension+1;s++) ret[s] = new Vertex(new Point(new double[dimension]));

		double tmp = 0;
		for(int i=0;i<dimension;i++){
			double nCoord = Math.sqrt(1-tmp); 
			ret[i].set(i, nCoord);
			double rCoord = (angle - tmp)/nCoord;

			for(int s=i+1;s<dimension+1;s++)
				ret[s].set(i, rCoord);

			tmp += rCoord*rCoord;
		}
		return ret;
	}

	/** Walk from tessel t to the tessel containing p */
	private Tessel walk(Tessel t, Vertex p){
		int i;
		boolean next = false;

		while(true){
			for(i=0;i<dimension+1;i++){
				//Get points
				List<Point> opPoints = new ArrayList<Point>();
				for (int l=1;l<=dimension;l++) {
					opPoints.add(t.getPoint((i+l)%(dimension+1)));
				}

				//Construct orientation matrix
				double[][] coords = new double[dimension+1][dimension+1];
				for (int l=0;l<dimension;l++) {
					double[] coord = opPoints.get(l).getCoords();
					for (int j=0;j<coord.length;j++) {
						coords[l][j] = coord[j]; 
					}
					coords[l][dimension] = 1;
				}
				double[] pCoords = p.getCoords();
				for (int k=0;k<pCoords.length;k++) {
					coords[dimension][k] = pCoords[k];
				}
				coords[dimension][dimension]=1;

				Matrix matrix = new Matrix(coords);
				double detP = matrix.determinant();

				double[] iCoords = t.getPoint(i).getCoords();
				for (int m=0;m<iCoords.length;m++) {
					coords[dimension][m] = iCoords[m];
				}

				matrix = new Matrix(coords);
				double detI = matrix.determinant();

				if (Math.signum(detI)!=Math.signum(detP)) { //Same sign indicates that p is inside
					t = t.getNeighbour(i);
					next=true;
					break;
				}
			}
			if(next) {
				next=false;
			}
			else return t;
		}
	}

	protected void compute() {

		//TODO: Take care of degeneracies in a better way than perturbation
		for(Vertex v: points){
			double[] coords = new double[dimension];
			for (int j=0;j<dimension;j++) {
				coords[j] = Randomization.randBetween(-0.00001, 0.00001);
			}
			v.addThis(new Vector(coords));
		}

		//Find the enclosing tetrahedron
		double max = 1000;//TODO find a more meaningful max
		Vertex[] firstTess = regularSimplex();
		for (int j=0;j<dimension+1;j++) {
			double[] tmpCoords = new double[dimension]; 
			for (int k=0;k<dimension;k++) {
				tmpCoords[k] = firstTess[j].getCoord(k)*max;
			}
			firstTess[j]=new Vertex(new Point(tmpCoords), true);
		}
		Tessel next_t = new Tessel(firstTess);
		Point inP = Point.getMidpoint(next_t.getPoint(0), next_t.getPoint(1));
		if (inSphere(next_t, new Point(inP))==false) {
			Vertex[] tPoints = next_t.getPoints();
			Vertex l = tPoints[dimension];
			Vertex s = tPoints[dimension-1];
			tPoints[dimension-1]=l;
			tPoints[dimension]=s;
			next_t = new Tessel(tPoints);
		}
		tessels.add(next_t);

		//Iterer over punkterne
		for(Vertex p: points){
			next_t = walk(next_t, p);

			//Collect conflicting tessels, create new ones and merge them non-conclicting tessels (removing the old)
			findNTes(next_t, p, -1);

			//Connect new tessels to each other
			int n = newTets.size();
			next_t=newTets.get(0);
			for (int j=0;j<n;j++) {
				Tessel tess = newTets.get(j);
				for (int k=0;k<n;k++) {
					if (k!=j) {
						Tessel neigh = newTets.get(k); 
						List<Vertex> comPoints = tess.getCommonVertices(neigh);
						if (comPoints.size()==dimension) {
							int index = tess.findpoint(tess.findVertex(neigh));
							tess.setNeighbour(index, neigh);
						}
					}
				}
			}
			newTets.clear();
		}
	}


	/** Find star shaped polytope and retriangulate */
	public void findNTes(Tessel tess, Vertex p, int apexID) {
		tess.setModified(true);
		tessels.remove(tess);
		Tessel newTes;
		for (int i=0;i<dimension+1;i++) {
			Tessel neigh = tess.getNeighbour(i);
			if (neigh==null) {
				List<Vertex> neighPoints = tess.oppositeVertices(i);
				Vertex[] newTess = new Vertex[dimension+1];
				for (int j=0;j<neighPoints.size();j++) {
					newTess[j]=neighPoints.get(j);
				}
				newTess[dimension] = p;
				newTes = new Tessel(newTess);
				Point inP = Point.getMidpoint(newTes.getPoint(0), newTes.getPoint(1));
				if (inSphere(newTes, new Point(inP))==false) {
					Vertex l = newTess[dimension];
					Vertex s = newTess[dimension-1];
					newTess[dimension-1]=l;
					newTess[dimension]=s;
					newTes = new Tessel(newTess);
				}
				tessels.add(newTes);
				newTets.add(newTes);
				continue;
			}
			if (neigh.isModified()) {
				continue;
			}
			else {
				if(inSphere(neigh,p)==false) {
					List<Vertex> points = tess.getCommonVertices(neigh);
					Vertex[] newTess = new Vertex[dimension+1];
					for (int j=0;j<points.size();j++) {
						newTess[j]=points.get(j);
					}
					newTess[dimension] = p;
					newTes = new Tessel(newTess);
					Point inP = Point.getMidpoint(newTes.getPoint(0), newTes.getPoint(1));
					if (inSphere(newTes, new Point(inP))==false) {
						Vertex l = newTess[dimension];
						Vertex s = newTess[dimension-1];
						newTess[dimension-1]=l;
						newTess[dimension]=s;
						newTes = new Tessel(newTess);
					}
					int pid= newTes.getID(p);
					int nid = -1;
					for (int j=0; j<dimension+1;j++) {
						if (!(points.contains(neigh.getPoint(j)))) {
							nid = j;
						}
					}
					newTes.setNeighbour(pid, neigh);
					neigh.setNeighbour(nid, newTes);
					tessels.add(newTes);
					newTets.add(newTes);
				}
				else {
					int apex = tess.apexid(i);
					findNTes(neigh, p, apex);	
				}
			}
		}
	}

	/** Return true if and only if the point p is inside or on the boundary of the circumsphere of t */
	public static boolean inSphere(Tessel t, Point p) {
		int dimension = p.getDimensions();
		//Create inSphere matrix (orientation of lifted coordinates)
		double[][] coords = new double[dimension+2][dimension+2];
		for (int i=0;i<dimension+1;i++) {
			double[] coord = t.getPoint(i).getCoords();
			for (int j=0;j<coord.length;j++) {
				coords[i][j] = coord[j]; 
			}
			coords[i][dimension] = t.getPoint(i).dot(t.getPoint(i));
			coords[i][dimension+1] = 1;
		}
		double[] pCoords = p.getCoords();
		for (int k=0;k<pCoords.length;k++) {
			coords[dimension+1][k] = pCoords[k];
		}
		coords[dimension+1][dimension]=p.dot(p);
		coords[dimension+1][dimension+1]=1;

		Matrix matrix = new Matrix(coords);

		double detP = matrix.determinant();
		return detP>0;			// Positive means p is inside circumsphere
	}


	/** Checks that all tetrahedra comply with the Delaunay-criteria. */
	public boolean checkTessels() {

		for(Tessel t: tessels){
			for(Vertex p: points){
				if(		!(t.containsPoint(p)) && inSphere(t, p)==true) {
					return false;
				}
			}
		}
		return true;
	}

	/** Convenient static method for creating delaunay tessellation of a point-set */
	public static List<Tessel> createDelaunay(List<Point> points) {
		if (!(points.isEmpty())) {
			BowyerWatson bw = new BowyerWatson(points);
			return bw.getTessels();
		}
		return new ArrayList<Tessel>();
	}

	/** Convenient static method for creating delaunay tessellation of a point-set 
	 * and if its 2D or 3D also displaying it. */
	public static List<Tessel> createAndDrawDelaunay(List<Point> points) {
		if (!(points.isEmpty())) {
			BowyerWatson bw = new BowyerWatson(points);
			BWViewer.displayTessel(bw);
			return bw.getTessels();
		}
		return new ArrayList<Tessel>();
	}

	public static void main(String[] args) {
		int setSize = 0;
		int dim = 0;
		if(args.length>0){
			try {
				setSize = Integer.parseInt(args[0]);
				dim = Integer.parseInt(args[1]);
			} catch (NumberFormatException e) {
				System.err.println("Argument" + " must be an integer");
				System.exit(1);
			}
		}else{
			setSize = 50;
			dim = 3;
		}
		List<Point> points = ProGAL.geomNd.PointList.generatePointsInCube(setSize, dim);
		BowyerWatson.createAndDrawDelaunay(points);
	}

}
