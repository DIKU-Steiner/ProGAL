package ProGAL.geom3d.complex;

import java.util.List;

public interface SimplicialComplex {
	public List<CVertex> getVertices();
	public List<CEdge> getEdges();
	public List<CTriangle> getTriangles();
	public List<CTetrahedron> getTetrahedra();
}
