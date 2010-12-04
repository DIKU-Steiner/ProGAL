package ProGAL.geom3d.volumes;

import ProGAL.geom3d.Point;

/** 
 *  
 */
public class Cone implements Volume{
	
	/**  */
	public Cone() {
	}

	
	
	/**
	 * Returns true if the specified point p is inside cone with the specified apex.
	 */
	public static boolean isInCone(Point p0, Point[] face, Point p) {
//		if (Point.isBehind(p0, face[0], face[2], face[1])) {
//			Point temp = face[1];
//			face[1]= face[2];
//			face[2] = temp;
//		}
//		//		System.out.println("p0 = " + p0.toString(3));
//		//		System.out.println("Face: " + face[0] + " " + face[1] + " " + face[2]);
//		for (int i = 0; i < 3 ; i++) {
//			//			System.out.println("p = " + p.toString(3));
//			if (!Point.isBehind(p, face[i], p0, face[(i+1)%3])) return false;
//		}
//		return true;
		// TODO Adapt to Cone. Moved from Point
		return false;
	}
	
	
	public double getVolume() {
		// TODO Auto-generated method stub
		return 0;
	}

	public boolean overlaps(Volume vol) {
		// TODO Auto-generated method stub
		return false;
	}

	public Point getCenter() {
		// TODO Auto-generated method stub
		return null;
	}

	public Cone clone(){
		return null;
	}
}

