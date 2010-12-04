package ProGAL.geom3d.viewer;

import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import javax.media.j3d.*;

import ProGAL.geom3d.*;

/** 
 * @author R.Fonseca
 */
public class TriangleSet3D extends Shape3D {

	public TriangleSet3D(Collection<Triangle> triangles, Appearance app) {
		super();

		
		List<Point> verts = new LinkedList<Point>();
		List<Vector> normals = new LinkedList<Vector>();

		for(Triangle t: triangles){
			Vector n = t.getP1().vectorTo(t.getP2()).crossThis(t.getP1().vectorTo(t.getP3()));
			if(n.getLength()>0.000001){
				verts.add(t.getP1().clone());
				verts.add(t.getP2().clone());
				verts.add(t.getP3().clone());
				normals.add(n);
				normals.add(n);
				normals.add(n);
			}
		}

		int i=0;
		float[] vertArr = new float[verts.size()*3];
		for(Point v: verts){ vertArr[i++] = (float)v.getX(); vertArr[i++] = (float)v.getY(); vertArr[i++] = (float)v.getZ(); }
		i=0;
		float[] normArr = new float[normals.size()*3];
		for(Vector v: normals){ normArr[i++] = (float)v.getX(); normArr[i++] = (float)v.getY(); normArr[i++] = (float)v.getZ(); }

		TriangleArray caps = new TriangleArray(vertArr.length/3, TriangleArray.COORDINATES | TriangleArray.NORMALS);

		caps.setCoordinates(0, vertArr);
		caps.setNormals(0, normArr);

		caps.setCapability(Geometry.ALLOW_INTERSECT);
		setGeometry(caps);
		if(app==null)
			setAppearance(new Appearance());
		else
			setAppearance(app);
	}



}