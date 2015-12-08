package ProGAL.geom3d.viewer;

import java.util.LinkedList;
import java.util.List;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.PI;

import javax.media.j3d.*;
import javax.vecmath.*;

/** Java3D capsule-shape (also known as line-swept-sphere) centered in origo and extending 
 * along the y-axis. It is rendered using a TriangleArray and by default uses 210 triangles. 
 * @author R.Fonseca
 */
class QuarterSphere3D extends Shape3D {

	/** Construct the hemisphere shape. 
	 * @param radius The radius of the hemisphere. 
	 * @param app Appearance of the hemisphere. 
	 */
	public QuarterSphere3D(float radius, Appearance app, int divisions) {
		super();

		List<Point3d> verts = new LinkedList<Point3d>();
		List<Vector3d> normals = new LinkedList<Vector3d>();

		float r = radius;
		double d = PI/(2*divisions);

		for(double theta=0;theta<PI;theta+=d){
			for(double phi=0;phi<PI/2;phi+=d){
				verts.add(new Point3d(  r*sin(theta  )*cos(phi  )  ,  r*sin(theta  )*sin(phi  )  ,  r*cos(theta  )  ));
				verts.add(new Point3d(  r*sin(theta+d)*cos(phi  )  ,  r*sin(theta+d)*sin(phi  )  ,  r*cos(theta+d)  ));
				verts.add(new Point3d(  r*sin(theta  )*cos(phi+d)  ,  r*sin(theta  )*sin(phi+d)  ,  r*cos(theta  )  ));

				verts.add(new Point3d(  r*sin(theta+d)*cos(phi+d)  ,  r*sin(theta+d)*sin(phi+d)  ,  r*cos(theta+d)  ));
				verts.add(new Point3d(  r*sin(theta  )*cos(phi+d)  ,  r*sin(theta  )*sin(phi+d)  ,  r*cos(theta  )  ));
				verts.add(new Point3d(  r*sin(theta+d)*cos(phi  )  ,  r*sin(theta+d)*sin(phi  )  ,  r*cos(theta+d)  ));
				
				normals.add(new Vector3d(  sin(theta  )*cos(phi  )  ,  sin(theta  )*sin(phi  )  ,  cos(theta  )  ));
				normals.add(new Vector3d(  sin(theta+d)*cos(phi  )  ,  sin(theta+d)*sin(phi  )  ,  cos(theta+d)  ));
				normals.add(new Vector3d(  sin(theta  )*cos(phi+d)  ,  sin(theta  )*sin(phi+d)  ,  cos(theta  )  ));

				normals.add(new Vector3d(  sin(theta+d)*cos(phi+d)  ,  sin(theta+d)*sin(phi+d)  ,  cos(theta+d)  ));
				normals.add(new Vector3d(  sin(theta  )*cos(phi+d)  ,  sin(theta  )*sin(phi+d)  ,  cos(theta  )  ));
				normals.add(new Vector3d(  sin(theta+d)*cos(phi  )  ,  sin(theta+d)*sin(phi  )  ,  cos(theta+d)  ));
				//normals.add(c.vectorTo(new Vector(sin(theta+d/2)*cos(phi+d/2)  ,  c.y() + r*sin(theta+d/2)*sin(phi+d/2)  ,  c.z() + r*cos(theta+d/2))));
			}
		}

		int i=0;
		float[] vertArr = new float[verts.size()*3];
		for(Point3d v: verts){ vertArr[i++] = (float)v.x; vertArr[i++] = (float)v.y; vertArr[i++] = (float)v.z; }
		i=0;
		float[] normArr = new float[normals.size()*3];
		for(Vector3d v: normals){ normArr[i++] = (float)v.x; normArr[i++] = (float)v.y; normArr[i++] = (float)v.z; }

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