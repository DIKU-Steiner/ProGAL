package ProGAL.geom3d.complex.betaComplex;

import ProGAL.geom3d.volumes.Tetrahedron;

public class CTet extends Tetrahedron{
	CVer[] corners = new CVer[4];
	
	public CTet(CVer v0, CVer v1, CVer v2, CVer v3) {
		super(v0.getCenter(), v1.getCenter(), v2.getCenter(), v3.getCenter());
		corners[0] = v0;
		corners[1] = v1;
		corners[2] = v2;
		corners[3] = v3;
	}
	
	
}
