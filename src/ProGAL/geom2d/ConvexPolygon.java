package ProGAL.geom2d;

import ProGAL.dataStructures.*;
import ProGAL.geom2d.viewer.J2DScene;


public class ConvexPolygon extends Polygon{
	private static final long serialVersionUID = 1L;

	public ConvexPolygon() { super(); }

	/** creates a convex hull of a simple polygon in O(n) time */
	public ConvexPolygon(Polygon pol) {
		int i = pol.leftExtremePointIndx();
		int n = pol.size();
		int k = (i+n-1)%n;
		add(pol.get(k));
		add(pol.get(i));
		add(pol.get((i+1)%n));
		int j = (i+2)%n;
		while (j != k) {
			Point p = pol.get(j);
			if (contains(p)) j++;
			else while (!Point.rightTurn(p, get(size()-1), get(size()-2))) deleteLast();
			add(p);
			j = (j+1)%n;
		}
	}

		public static enum ConvexHullAlgorithm  { JarvisMarch, GrahamsScan };
	
	/** creates a convex polygon with given three points as corners (counterclockwise) */
	public ConvexPolygon(Point p0, Point p1, Point p2) { super(p0, p1, p2); }
	
	/** creates a convex hull of a set of points 
	 * Jarvis March O(hn)
	 * Grahams Scan O(nlogn)  */
	public ConvexPolygon(PointSet points, ConvexHullAlgorithm algorithm) {
		switch (algorithm) {
		case JarvisMarch:
			Point p = points.leftExtremePoint();
			Point r = p;
			do add(r = points.getNextExtremePoint(r)); while (r != p);
			break;
		case GrahamsScan:
			Sorter sort = new SorterQuick();
			sort.Sort(points, new SortToolPoint2dAroundPoint(points.getCentroid()), false);
			Polygon pol = new Polygon(points);   
			int indx = pol.leftExtremePointIndx();
			int k1 = (indx + pol.size() - 1)%pol.size();
			int k2 = indx;
			int k3 = (indx + 1)%pol.size();
			boolean back = true;
			while ((k2 != indx) || back) {
				if (Point.leftTurn(pol.get(k1), pol.get(k2), pol.get(k3))) {
					k1 = k2;
					k2 = k3;
					k3 = (k2 + 1)%pol.size();
					back = false;
				}
				else {
					pol.remove(k2); 
					if (k2 < indx) indx--;  // 012  123   234  345  450 501       401 012 123 234 340 340
					if (k1 == 0) k1 = pol.size()-1; else { if (k1 == pol.size()) k1 = k1-2; else k1--; }
					k2 = (k1 + 1)%pol.size();
					k3 = (k2 + 1)%pol.size();
					back = true;
				}
			}
			for (Point q : pol) add(q);
			break;
		}
	}
	

	
	/** returns true if point p is inside the convex polygon */
	public boolean contains(Point p) {
		for (int i = 0; i < size()-1; i++) if (Point.rightTurn(get(i), get(i+1), p)) return false; 
		return Point.leftTurn(get(size()-1), get(0), p);
	}
	
	
	/** returns the segment between most distant pair of points. */
	public LineSegment getDiameter() {
		int n = size();
		int i = leftExtremePointIndx();  
		int iStart = i;
		int iNext = (i + 1 == n)? 0 : i+1;
		Vector vi = new Vector(get(i), get(iNext));
		int j = rightExtremePointIndx(); 
		int jStart = j;
		int jNext = (j + 1 == n)? 0 : j+1;
		Vector vj = new Vector(get(jNext), get(j));
		LineSegment seg = new LineSegment(get(i), get(j));
		LineSegment diameter = seg.clone();
		double diamLength = seg.getSquaredLength();
		double lng;
		do {
//			System.out.println(i + ": " + (get(i)).toString() + ",  " + j +": " + (get(j)).toString() + ",  " + seg.getSquaredLength());
//			System.out.println("Vectors vi and vj: " + vi.toString() + "  " + vj.toString());
			if  (Vector.rightTurn(vi, vj)) {
				j = jNext;
				jNext = (j + 1 == n)? 0 : j+1;
				vj = new Vector(get(jNext), get(j));
				seg.setB(get(j));
			}
			else {
				i = iNext;
				iNext = (i + 1 == n)? 0 : i+1;
				vi = new Vector(get(i), get(iNext));
				seg.setA(get(i));
			}
			lng = seg.getSquaredLength();
			if (lng > diamLength) { diameter = seg.clone(); diamLength = lng; }
		} while  ((i != jStart) || (j != iStart)) ;
		return diameter;
	}

	public double[][] beamDetector() {
		int n = size();
		double[][] bp = new double[n][n];
		int[][] p = new int[n][n];
		for (int k = 0; k < n; k++) { bp[k][(k+1)%n] = 0.0; p[k][(k+1)%n] = k; }
		for (int l = 2; l < n; l++) {
			for (int i = 0; i < n; i++) {
				int j = (i+l)%n;
				LineSegment seg = new LineSegment(get(i), get(j));
				bp[i][j] = 99999999.9;
				int k = i;
				while (k != j) {
					k = (k+1)%n;
					double length = bp[i][k] + bp[k][j] + get(k).getDistance(seg);
					if (bp[i][j] > length) {
						p[i][j] = k;
						bp[i][j] = length;
					}
				}
			}
		}
		return bp;
	}
	
	public int farthestVertex(int p, int r) {
		int n = size();
		LineSegment seg = new LineSegment(get(p), get(r));
		double dist = 0.0;
		int q = -1;
		int i = (p+1)%n; 
		while (i != r) {
			double d = get(i).getSquaredDistance(seg);
			if (d > dist) {
				dist = d;
				q = i;
			}
			i = (i+1)%n;
		}
		return q;
	}
	
	public static void main(String[] args) {
		PointSet points = new PointSet(200);
		ConvexPolygon cPol = new ConvexPolygon(points, ConvexPolygon.ConvexHullAlgorithm.GrahamsScan);
		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		cPol.draw(scene);
	}

}
