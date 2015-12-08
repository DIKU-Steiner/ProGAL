package ProGAL.geom3d.complex.betaComplex;

import ProGAL.geom3d.volumes.Sphere;

public class CVer extends Sphere{
	CTet tet;

	public CVer(Sphere s) {
		super(s.getCenter(), s.getRadius());
	}

	
}
