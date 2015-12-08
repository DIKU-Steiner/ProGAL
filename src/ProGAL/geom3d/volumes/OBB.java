package ProGAL.geom3d.volumes;

import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.Vector;
import ProGAL.math.Matrix.EigenvalueDecomposition;
import ProGAL.math.Matrix3x3;

public class OBB implements Volume{
	protected Point anchor;  // midpoint of the box
	protected Vector[] bases;   // unit vectors
	public double[] extents;

	public OBB(Point anchor, Vector xdir, Vector ydir, Vector zdir) {
		this.anchor = anchor;
		extents = new double[]{xdir.length(), ydir.length(), zdir.length()}; 
		bases = new Vector[]{xdir.scaleToLength(1), ydir.scaleToLength(1), zdir.scaleToLength(1)};
	}
	public OBB(Point anchor, Vector[] bases, double[] extents){
		this.anchor = anchor;
		this.bases = bases;
		this.extents = extents;
	}

	public Point getAnchor() { return anchor; }
	public Vector getXDir() { return bases[0]; }
	public Vector getYDir() { return bases[1]; }
	public Vector getZDir() { return bases[2]; }
	public Vector[] getBases() { return bases; }
	public void setAnchor(Point pos) { anchor = pos; }

	public double cutArealYZ() {
		return 4*extents[1]*extents[2];
	}

//	public static Box3d createBoundingBox(List<Point> points, int i, int j) {
//		List<Point> pointSubset = new List<Point>();
//		for (int k = i; k <= j; k++) pointSubset.insert(points.get(k));
//		return createBoundingBox_Covariance(pointSubset);
//	}

	/* 458HOps */
	public static OBB createBoundingBox_Covariance(PointList points){
		Matrix3x3 m = points.getCovariance();
		EigenvalueDecomposition ed = m.getEigenvalueDecomposition();
		ProGAL.geomNd.Vector r = ed.getV().getColumn(2);
		ProGAL.geomNd.Vector s = ed.getV().getColumn(1);
		ProGAL.geomNd.Vector t = ed.getV().getColumn(0);
		OBB ret = createBoxFromBases(new Vector[]{new Vector(r),new Vector(s),new Vector(t)}, points);
//		double[][] coords = new double[3][3];
//		for(int i=0;i<3;i++) for(int j=0;j<3;j++) coords[i][j] = m.get(i, j);
//		Jama.Matrix A = new Jama.Matrix(coords);
//		EigenvalueDecomposition ed = new EigenvalueDecomposition(A);//84HOps
//		Jama.Matrix eV = ed.getV();

//		Vector r = new Vector(eV.get(0, 2), eV.get(1, 2), eV.get(2, 2));
//		Vector s = new Vector(eV.get(0, 1), eV.get(1, 1), eV.get(2, 1));
//		Vector t = r.cross(s);//6HOps
////		Vector t = new Vector(eV.get(0, 0), eV.get(1, 0), eV.get(2, 0));
//		//		Vector[] eigens = m.getEigenvectors();
//		Box3d ret = createBoxFromBases( new Vector[]{r,s,t}, points);//16*18+71
		
		return ret;
	}

//	public static Box3d createBoundingBox_Covariance(Box3d b1, Box3d b2){
//		List<Point> corners = new List<Point>();
//		corners.append(b1.getVertices());//72HOps
//		corners.append(b2.getVertices());//72HOps
//		return createBoundingBox_Covariance(corners);
//	}


	/**
	 * Creates a bounding box enclosing the points in <code>points</code> using the supplied 
	 * bases.
	 * @hops n*9+15
	 */
	public static OBB createBoxFromBases(Vector[] bases, List<Point> points){
		//Find point farthest back along each base
//		Point[] fPoints = new Point[3];
		double[][] extremeDots = {
				{Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY},
				{Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY},
				{Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY} };
		for(int b=0;b<3;b++){
			for(Point p: points){
				double dot = bases[b].dot(p.toVector());//3HOps
				if(dot<extremeDots[b][0]) 
					extremeDots[b][0] = dot;
				if(dot>extremeDots[b][1])
					extremeDots[b][1] = dot;
			}
		}//n*3*3

		Point center = bases[0].multiply((extremeDots[0][0]+extremeDots[0][1])/2).toPoint();//4HOps
		center.addThis(bases[1].multiply((extremeDots[1][0]+extremeDots[1][1])/2));//4HOps
		center.addThis(bases[2].multiply((extremeDots[2][0]+extremeDots[2][1])/2));//4HOps
		
		double[] extents = new double[3];
		extents[0] = (extremeDots[0][1]-extremeDots[0][0])/2;//1HOps
		extents[1] = (extremeDots[1][1]-extremeDots[1][0])/2;//1HOps
		extents[2] = (extremeDots[2][1]-extremeDots[2][0])/2;//1HOps
		
		return new OBB(center,bases,extents);
		
//		//Solve Ax=b
//		double[][] array = {
//				{bases[0].x,bases[0].y,bases[0].z},
//				{bases[1].x,bases[1].y,bases[1].z},
//				{bases[2].x,bases[2].y,bases[2].z}};
//		Jama.Matrix A = new Jama.Matrix(array);
//		Jama.Matrix b = new Jama.Matrix(new double[][]{
//				{Vector.dotProduct(fPoints[0],bases[0])},
//				{Vector.dotProduct(fPoints[1],bases[1])},
//				{Vector.dotProduct(fPoints[2],bases[2])}
//		});//9HOps
//		Jama.Matrix x = A.solve(b);//Est. 50HOps
//
//		Point corner = new Point(x.get(0, 0),x.get(1, 0),x.get(2, 0));
//
//		//Find extents
//		double[] extents = new double[3];
//		for(int base=0;base<3;base++){
//			double maxDot = -1;
//			for(Point point: points){
//				double dot = corner.vectorTo(point).dotProduct(bases[base]);
//				if(dot>maxDot) maxDot = dot;
//			}
//			extents[base] = Math.max(maxDot,Constants.epsilon);
//
//		}//n*3*3
//		bases[0].scaleToLength(extents[0]/2);//4HOps
//		bases[1].scaleToLength(extents[1]/2);//4HOps
//		bases[2].scaleToLength(extents[2]/2);//Bases now have the right dimension//4HOps
//		corner.addIn(bases[0]);
//		corner.addIn(bases[1]);
//		corner.addIn(bases[2]);//Corner is now center
		
//		return new Box3d(corner, bases[0], bases[1], bases[2]);
	}


	/** Returns true iff this box has an overlap with b. Uses the separating axis implementation 
	 * described in "Realtime collision detection" by Christer Ericsson.  30 -> 27+9+9+18+54=117HOps*/
	public boolean overlaps(OBB b){
		double[][] R = new double[3][3], AbsR = new double[3][3];
		double[] ae = extents, be = b.extents;

		// Compute rotation matrix expressing b in a�s coordinate frame
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				R[i][j] = bases[i].dot(b.bases[j]);	
		//Total: 3*3*3 = 27HOps

		// Compute translation vector t
		Vector tmp = this.anchor.vectorTo(b.anchor);

		// Bring translation into a�s coordinate frame
		tmp = new Vector(tmp.dot(bases[0]), tmp.dot(bases[1]), tmp.dot(bases[2]));//9HOps
		double[] t = new double[]{tmp.x(), tmp.y(), tmp.z()};

		// Compute common subexpressions. Add in an epsilon term to
		// counteract arithmetic errors when two edges are parallel and
		// their cross product is (near) null (see text for details)
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				AbsR[i][j] = Math.abs(R[i][j])+0.00001;

		// Test axes L = A0, L = A1, L = A2
		for (int i = 0; i < 3; i++) 
			if (Math.abs(t[i]) > ae[i] + be[0]*AbsR[i][0]+be[1]*AbsR[i][1]+be[2]*AbsR[i][2]) return false; //9HOps
		// Test axes L = B0, L = B1, L = B2
		for (int i = 0; i < 3; i++) 
			if (Math.abs(t[0]*R[0][i] + t[1]*R[1][i] + t[2]*R[2][i]) > ae[0]*AbsR[0][i] + ae[1]*AbsR[1][i] + ae[2]*AbsR[2][i] + be[i]) return false;  //18HOps

		//For all the rest: 9*6 = 54HOps
		// Test axis L = A0 x B0
		if (Math.abs(t[2] * R[1][0] - t[1] * R[2][0]) > ae[1] * AbsR[2][0] + ae[2] * AbsR[1][0] + be[1] * AbsR[0][2] + be[2] * AbsR[0][1]) return false;
		// Test axis L = A0 x B1
		if (Math.abs(t[2] * R[1][1] - t[1] * R[2][1]) > ae[1] * AbsR[2][1] + ae[2] * AbsR[1][1] + be[0] * AbsR[0][2] + be[2] * AbsR[0][0]) return false;
		// Test axis L = A0 x B2
		if (Math.abs(t[2] * R[1][2] - t[1] * R[2][2]) > ae[1] * AbsR[2][2] + ae[2] * AbsR[1][2] + be[0] * AbsR[0][1] + be[1] * AbsR[0][0]) return false;
		// Test axis L = A1 x B0
		if (Math.abs(t[0] * R[2][0] - t[2] * R[0][0]) > ae[0] * AbsR[2][0] + ae[2] * AbsR[0][0] + be[1] * AbsR[1][2] + be[2] * AbsR[1][1]) return false;
		// Test axis L = A1 x B1
		if (Math.abs(t[0] * R[2][1] - t[2] * R[0][1]) > ae[0] * AbsR[2][1] + ae[2] * AbsR[0][1] + be[0] * AbsR[1][2] + be[2] * AbsR[1][0]) return false;
		// Test axis L = A1 x B2
		if (Math.abs(t[0] * R[2][2] - t[2] * R[0][2]) > ae[0] * AbsR[2][2] + ae[2] * AbsR[0][2] + be[0] * AbsR[1][1] + be[1] * AbsR[1][0]) return false;
		// Test axis L = A2 x B0
		if (Math.abs(t[1] * R[0][0] - t[0] * R[1][0]) > ae[0] * AbsR[1][0] + ae[1] * AbsR[0][0] + be[1] * AbsR[2][2] + be[2] * AbsR[2][1]) return false;
		// Test axis L = A2 x B1
		if (Math.abs(t[1] * R[0][1] - t[0] * R[1][1]) > ae[0] * AbsR[1][1] + ae[1] * AbsR[0][1] + be[0] * AbsR[2][2] + be[2] * AbsR[2][0]) return false;
		// Test axis L = A2 x B2
		if (Math.abs(t[1] * R[0][2] - t[0] * R[1][2]) > ae[0] * AbsR[1][2] + ae[1] * AbsR[0][2] + be[0] * AbsR[2][1] + be[1] * AbsR[2][0]) return false;
		// Since no separating axis is found, the OBBs must be intersecting
		return true;
	}


	public Point getCenter() {
		return this.getAnchor();
	}


	/* 3*3*8 = 72HOps */
	public Point[] getVertices(){
		Point[] ret = new Point[8];
		ret[0] = anchor.subtract(bases[0].multiply(extents[0])).subtractThis(bases[1].multiply(extents[1])).subtractThis(bases[2].multiply(extents[2]));
		ret[1] = anchor.subtract(bases[0].multiply(extents[0])).subtractThis(bases[1].multiply(extents[1])).addThis(bases[2].multiply(extents[2]));
		ret[2] = anchor.subtract(bases[0].multiply(extents[0])).addThis(bases[1].multiply(extents[1])).subtractThis(bases[2].multiply(extents[2]));
		ret[3] = anchor.subtract(bases[0].multiply(extents[0])).addThis(bases[1].multiply(extents[1])).addThis(bases[2].multiply(extents[2]));
		ret[4] = anchor.add(bases[0].multiply(extents[0])).subtractThis(bases[1].multiply(extents[1])).subtractThis(bases[2].multiply(extents[2]));
		ret[5] = anchor.add(bases[0].multiply(extents[0])).subtractThis(bases[1].multiply(extents[1])).addThis(bases[2].multiply(extents[2]));
		ret[6] = anchor.add(bases[0].multiply(extents[0])).addThis(bases[1].multiply(extents[1])).subtractThis(bases[2].multiply(extents[2]));
		ret[7] = anchor.add(bases[0].multiply(extents[0])).addThis(bases[1].multiply(extents[1])).addThis(bases[2].multiply(extents[2]));
		return ret;
	}


	public boolean overlaps(Volume s) {
		if(s instanceof OBB) return overlaps((OBB)s);
		throw new Error("Unknown volume");
	}

	public double volume() {
		return extents[0]*extents[1]*extents[2]*8;
	}

	public OBB clone(){
		Point aClone = anchor.clone();
		Vector[] bClone = {bases[0].clone(),bases[1].clone(),bases[2].clone()};
		double[] eClone = {extents[0],extents[1],extents[2]};
		return new OBB(aClone,bClone,eClone);

	}

	public String toString(){
		return "Box3d["+getAnchor()+","+bases[0]+","+bases[1]+","+bases[2]+"]";
	}

	@Override
	public double getVolume() {
		// TODO Auto-generated method stub
		return 0;
	}

}
