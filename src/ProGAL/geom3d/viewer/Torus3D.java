package ProGAL.geom3d.viewer;

import java.util.LinkedList;
import java.util.List;

import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.PI;

import javax.media.j3d.*;
import javax.vecmath.*;


/** Java3D torus-shape centered in origo and a normal along the z-axis 
 * It is rendered using a TriangleArray 
 * @author R.Fonseca
 */
class Torus3D extends Shape3D {

	/** Construct the torus shape. 
	 * @param height The distance between the defining points.
	 * @param radius The radius of the cylinder. 
	 * @param app Appearance of the cylinder. 
	 */
	public Torus3D(double majorRad, double minorRad, Appearance app, int divisions) {
		super();

		List<Point3d> verts = new LinkedList<Point3d>();
		List<Vector3d> normals = new LinkedList<Vector3d>();

		double d = PI/divisions;

		for(double phi=0;phi<PI*2;phi+=d){
			for(double psi=0;psi<PI*2;psi+=d){
				verts.add(new Point3d(  majorRad*cos(phi  ) + minorRad*cos(phi  )*cos(psi  ), majorRad*sin(phi  ) + minorRad*sin(phi  )*cos(psi  ), sin(psi  )  ));
				verts.add(new Point3d(  majorRad*cos(phi+d) + minorRad*cos(phi+d)*cos(psi  ), majorRad*sin(phi+d) + minorRad*sin(phi+d)*cos(psi  ), sin(psi  )  ));
				verts.add(new Point3d(  majorRad*cos(phi  ) + minorRad*cos(phi  )*cos(psi+d), majorRad*sin(phi  ) + minorRad*sin(phi  )*cos(psi+d), sin(psi+d)  ));
				verts.add(new Point3d(  majorRad*cos(phi+d) + minorRad*cos(phi+d)*cos(psi+d), majorRad*sin(phi+d) + minorRad*sin(phi+d)*cos(psi+d), sin(psi+d)  ));

				normals.add(new Vector3d(  cos(phi)*cos(psi)  , sin(phi)*cos(psi),   sin(psi)  ));
			}
		}

		int i=0;
		float[] vertArr = new float[verts.size()*4];
		for(Point3d v: verts){ vertArr[i++] = (float)v.x; vertArr[i++] = (float)v.y; vertArr[i++] = (float)v.z; }
		i=0;
		float[] normArr = new float[normals.size()*4];
		for(Vector3d v: normals){ normArr[i++] = (float)v.x; normArr[i++] = (float)v.y; normArr[i++] = (float)v.z; }

		QuadArray quads = new QuadArray(vertArr.length/4, QuadArray.COORDINATES | QuadArray.NORMALS );
		quads.setCoordinates(0, vertArr);
		quads.setNormals(0, normArr);

		quads.setCapability(Geometry.ALLOW_INTERSECT);
		setGeometry(quads);
		if(app==null)
			setAppearance(new Appearance());
		else
			setAppearance(app);
	}



}