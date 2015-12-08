package ProGAL.geom2d;

/**
 * An interface for 2d shapes such as circles, rectangles and text. 
 */
public interface Shape {
	
	/** Get the geometric center of the shape. The center of a shape can be interpreted 
	 * in many ways (center of mass, circumcenter, inscribed center etc.). No strict 
	 * requirement is given here, but typically the circumcenter should be supplied.
	 * 
	 * This method is most prominently used to find the average center-position of a 
	 * collection of shapes such that scene-viewers can be centered on the scene.
	 */
	public Point getCenter();
	
	public boolean contains(Point p);
}
