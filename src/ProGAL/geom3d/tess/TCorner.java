package ProGAL.geom3d.tess;

public class TCorner {
	TPoint v; // index of vertex
	int opp; // pointer to opposite corner in neighboring tet
	TSphere plane; // plane equation of separating plane (v * plane > 0)

}
