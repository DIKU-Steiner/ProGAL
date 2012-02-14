package ProGAL.geom3d.viewer;

import javax.media.j3d.*;

/** 
 * @author R.Fonseca
 */
class Triangle3D extends Shape3D {

	public Triangle3D(float size, Appearance app) {
		super();

		float[] vertArr = {0,0,0,  size,0,0,  0,size,0};//new float[verts.size()*3];
		float[] normArr = {0,0,1,  0,0,1,  0,0,1};//new float[normals.size()*3];

		TriangleArray caps = new TriangleArray(3, TriangleArray.COORDINATES | TriangleArray.NORMALS);

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