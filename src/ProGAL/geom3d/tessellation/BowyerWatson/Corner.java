package ProGAL.geom3d.tessellation.BowyerWatson;

import ProGAL.geom3d.PointWeighted;

public class Corner {
	final Tet tet;
	final PointWeighted point;
	
	
	public Corner(Tet tet, PointWeighted point){
		this.tet = tet;
		this.point = point;
	}
}
