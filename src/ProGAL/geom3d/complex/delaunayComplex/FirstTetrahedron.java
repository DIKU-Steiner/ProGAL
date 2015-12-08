package ProGAL.geom3d.complex.delaunayComplex;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.complex.CTetrahedron;

class FirstTetrahedron extends CTetrahedron{

	FirstTetrahedron(double max){
		super();
		double k = 200;
		double m=k*max;
		
//		m = 30;
		
		setPoint(new CVertex(new Point( m, 0,-m),true, -4), 0 );
		setPoint(new CVertex(new Point( 0, m,-m),true, -3), 1 );
		setPoint(new CVertex(new Point(-m,-m,-m),true, -2), 2 );
		setPoint(new CVertex(new Point( 0, 0, m),true, -1), 3 );
	}

}