package ProGAL.geom3d.viewer;

import java.util.LinkedList;
import java.util.List;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.PI;

import javax.media.j3d.*;
import javax.vecmath.*;



/** Java3D cylinder-shape centered in origo and extending 
 * along the y-axis. It is rendered using a TriangleArray 
 * and by default uses 14 triangles. 
 * @author R.Fonseca
 */
class HalfCylinder3D extends Shape3D {

	/** Construct the cylinder shape. 
	 * @param height The distance between the defining points.
	 * @param radius The radius of the cylinder. 
	 * @param app Appearance of the cylinder. 
	 */
	public HalfCylinder3D(float height, float radius, Appearance app, int divisions) {
		super();

		List<Point3d> verts = new LinkedList<Point3d>();
		List<Vector3d> normals = new LinkedList<Vector3d>();

		float r = radius;
		double d = PI/divisions;

		double dH = height;
		for(double theta=0;theta<PI;theta+=d){
			double y=-height/2;
			verts.add(new Point3d(  r*sin(theta  )  , y   , r*cos(theta  )  ));
			verts.add(new Point3d(  r*sin(theta+d)  , y   , r*cos(theta+d)  ));
			verts.add(new Point3d(  r*sin(theta  )  , y+dH, r*cos(theta  )  ));
			
			verts.add(new Point3d(  r*sin(theta+d)  , y+dH, r*cos(theta+d)  ));
			verts.add(new Point3d(  r*sin(theta  )  , y+dH, r*cos(theta  )  ));
			verts.add(new Point3d(  r*sin(theta+d)  , y   , r*cos(theta+d)  ));
			normals.add(new Vector3d(  sin(theta  )  , 0 , cos(theta  )  ));
			normals.add(new Vector3d(  sin(theta+d)  , 0 , cos(theta+d)  ));
			normals.add(new Vector3d(  sin(theta  )  , 0 , cos(theta  )  ));
			
			normals.add(new Vector3d(  sin(theta+d)  , 0 , cos(theta+d)  ));
			normals.add(new Vector3d(  sin(theta  )  , 0 , cos(theta  )  ));
			normals.add(new Vector3d(  sin(theta+d)  , 0 , cos(theta+d)  ));
		}

		int i=0;
		float[] vertArr = new float[verts.size()*3];
		for(Point3d v: verts){ vertArr[i++] = (float)v.x; vertArr[i++] = (float)v.y; vertArr[i++] = (float)v.z; }
		i=0;
		float[] normArr = new float[normals.size()*3];
		for(Vector3d v: normals){ normArr[i++] = (float)v.x; normArr[i++] = (float)v.y; normArr[i++] = (float)v.z; }

		TriangleArray caps = new TriangleArray(vertArr.length/3, TriangleArray.COORDINATES | TriangleArray.NORMALS );

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