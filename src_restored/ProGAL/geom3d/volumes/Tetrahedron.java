package ProGAL.geom3d.volumes;

import ProGAL.geom3d.Point3d;

public class Tetrahedron implements Volume3d {
	private Point3d[] corners = new Point3d[4];
	
	public Tetrahedron(Point3d p1, Point3d p2, Point3d p3, Point3d p4){
		corners[0] = p1;
		corners[1] = p2;
		corners[2] = p3;
		corners[3] = p4;
	}
	
	public double getVolume() {
		// TODO Auto-generated method stub
		return 0;
	}

	public boolean overlaps(Volume3d vol) {
		// TODO Auto-generated method stub
		return false;
	}

	public Point3d getCenter() {
		// TODO Auto-generated method stub
		return null;
	}

	public Volume3d clone(){
		return null;
	}
}
