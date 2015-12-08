package ProGAL.geom3d.volumes;

import ProGAL.geom3d.Shape;

/**
 * An interface for 3d volumes such as spheres and boxes. Implementing classes should consider 
 * a volume as all points within and on the boundary of the volume. This implies, for instance, 
 * that two spheres overlap when the distance between their center-points is <b>equal to or  
 * less</b> than their combined radii. 
 * 
 * If only the shell should be considered then subclasses with adequate overriding methods should 
 * be created.
 * @author R.Fonseca
 */
public interface Volume extends Shape{
	
	/** Determine if this volume overlaps vol. Two volumes overlap if their surfaces touch 
	 * or if the union of interiors is non-empty.
	 */
	public boolean overlaps(Volume vol);
	
	/** Get the volume. */ 
	public double getVolume();

	/** Make a deep clone this volume. */
	public Volume clone();
}
