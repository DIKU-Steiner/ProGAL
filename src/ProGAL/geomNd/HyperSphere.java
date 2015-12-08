package ProGAL.geomNd;

public class HyperSphere {
	private Point center;
	private double radius;
	
	public HyperSphere(Point center, double radius){
		this.center = center;
		this.radius = radius;
	}
	
	public Point getCenter(){	return center;	}
	public double getRadius(){	return radius;	}

	public void setCenter(Point p){ center = p; }
	public void setRadius(double r){ radius = r; }
	
	
	public String toString(){ return toString(2); }
	
	public String toString(int dec){
		return String.format("HyperSphere[%s, %."+dec+"f]", center.toString(), radius);
	}
}
