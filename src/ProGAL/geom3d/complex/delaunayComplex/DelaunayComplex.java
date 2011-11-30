package ProGAL.geom3d.complex.delaunayComplex;

import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.predicates.*;

public class DelaunayComplex extends RegularComplex {

	public DelaunayComplex(List<Point> points) {
		super(points, new InexactJavaPredicates());
	}

}
