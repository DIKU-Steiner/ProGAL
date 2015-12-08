package ProGAL.geom3d.complex.betaComplex;

import java.util.ArrayList;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.complex.*;
import ProGAL.geom3d.complex.alphaComplex.AlphaFiltration;
import ProGAL.geom3d.volumes.Sphere;

public class BetaComplex {
	private List<CVer> vertices = new ArrayList<CVer>();
	private List<CTet> tetrahedra = new ArrayList<CTet>();
	
	
	public BetaComplex(List<Sphere> spheres){
		
		List<Point> centers = new ArrayList<Point>(spheres.size());
		AlphaFiltration af = new AlphaFiltration(centers);
		
		
		
	}
	
	
	public static void main(String[] args) {

	}

}
