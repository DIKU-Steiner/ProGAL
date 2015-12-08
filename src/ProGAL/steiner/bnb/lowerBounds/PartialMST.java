package ProGAL.steiner.bnb.lowerBounds;

import ProGAL.geomNd.Point;
import ProGAL.steiner.bnb.LowerBound;
import ProGAL.steiner.bnb.Node;
import ProGAL.steiner.bnb.PointPlacementWS;
import ProGAL.steiner.bnb.Topology;

public class PartialMST implements LowerBound{
	private final Topology top;
//	private final PointPlacement placement;
	private final PointPlacementWS placement;
	private final double tolerance;
	
	public PartialMST(Point[] sites, double tolerance){
		int N = sites.length;
		top = new Topology(N);
//		placement = new PointPlacement(sites);
		placement = new PointPlacementWS(sites);
		this.tolerance = tolerance;
	}
	
	@Override
	public double lowerBound(Node n) {
		top.updateFromNode(n);
		double len = placement.updateSteinerPoints(top, tolerance);
		return len;
	}

}
