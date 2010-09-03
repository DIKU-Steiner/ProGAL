package ProGAL.geom3d.volumes;

import ProGAL.geom3d.Shape3d;

/**
 * An interface for 3d volumes such as spheres and boxes. 
 */
public interface Volume3d extends Shape3d{
	public boolean overlaps(Volume3d vol);
	public double getVolume();
	public Volume3d clone();
}
