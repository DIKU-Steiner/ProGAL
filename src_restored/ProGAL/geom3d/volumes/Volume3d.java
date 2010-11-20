package ProGAL.geom3d.volumes;

import ProGAL.geom3d.Shape3d;

/**
 * An interface for 3d volumes such as spheres and boxes. Implementing classes should consider 
 * a volume as all points within and on the boundary of the volume. This implies, for instance, 
 * that two spheres overlap when the distance between their center-points is <b>equal to or  
 * less</b> than their combined radii. 
 * 
 * If only the shell should be considered use the VolumeShell3d interface.
 */
public interface Volume3d extends Shape3d{
	public boolean overlaps(Volume3d vol);
	public double getVolume();
	public Volume3d clone();
}
