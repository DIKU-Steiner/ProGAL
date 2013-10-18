package ProGAL.geom2d;
/**
The standard equation of an ellipse is

(x-h)^2/rx^2 + (y-k)^2/ry^2 = 1

where (h,k) is the center of the ellipse, rx is the distance from the
center in x-direction and ry is the distance from the center in the 
y-direction

An ellipse has two focal points in which the sum of the length of 
both focal points to any given point on the ellipse is always the same.

The foci of an ellipse are on the major axis at the distance c form the c
enter. The value of c is rx^2-ry^2 if 
the major axis is horizontal and ry^2-rx^2 if the major axis is vertical.
*/


public class Ellipse {

	Point center;
	double rx;
	double ry;
	
	public Ellipse(Point center, double rx, double ry) {
		this.center = center;
		this.rx = rx;
		this.ry = ry;
	}
	
	public Point getCenter() { return center; }
	public double getRx() { return rx; }
	public double getRy() { return ry; }
	
	public Point[] getFoci() {
		Point[] foci = new Point[2];
		double c;
		if (rx >= ry) {
			c = Math.sqrt(rx*rx - ry*ry);
			foci[0] = center.add(-c,0);
			foci[1] = center.add(c,0);
		}
		else {
			c = Math.sqrt(ry*ry - rx*rx);
			foci[0] = center.add(0,-c);
			foci[1] = center.add(0,c);
		}
		return foci;
	}
	
}
