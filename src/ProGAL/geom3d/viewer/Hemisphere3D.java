package ProGAL.geom3d.viewer;

import java.util.LinkedList; 
import java.util.List;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.PI;

import javax.media.j3d.*;

import ProGAL.geom3d.*;

/** Java3D capsule-shape (also known as line-swept-sphere) centered in origo and extending 
 * along the y-axis. It is rendered using a TriangleArray and by default uses 210 triangles. 
 * @author R.Fonseca
 */
class Hemisphere3D extends Shape3D {

	/** Construct the hemisphere shape. 
	 * @param radius The radius of the hemisphere. 
	 * @param app Appearance of the hemisphere. 
	 */
	public Hemisphere3D(float radius, Appearance app, int divisions) {
		super();
		
		List<Point> verts = new LinkedList<Point>();
		List<Vector> normals = new LinkedList<Vector>();

		float r = radius;
		double d = PI/divisions;

		for(double theta=0;theta<PI;theta+=d){
			double sinTh = sin(theta);
			double sinThD = sin(theta+d);
			double cosTh = cos(theta);
			double cosThD = cos(theta+d);
			for(double phi=0;phi<PI;phi+=d){
				double cosPh = cos(phi);
				double cosPhD = cos(phi+d);
				double sinPh = sin(phi);
				double sinPhD = sin(phi+d);
				
				verts.add(new Point(  r*sinTh*cosPh  ,  r*sinTh*sinPh  ,  r*cosTh  ));
				verts.add(new Point(  r*sinThD*cosPh  ,  r*sinThD*sinPh  ,  r*cosThD  ));
				verts.add(new Point(  r*sinTh*cosPhD  ,  r*sinTh*sinPhD  ,  r*cosTh  ));
				verts.add(new Point(  r*sinThD*cosPhD  ,  r*sinThD*sinPhD  ,  r*cosThD  ));
				verts.add(new Point(  r*sinTh*cosPhD  ,  r*sinTh*sinPhD  ,  r*cosTh  ));
				verts.add(new Point(  r*sinThD*cosPh  ,  r*sinThD*sinPh  ,  r*cosThD  ));
				
				normals.add(new Vector(  sinTh*cosPh  ,  sinTh*sinPh  ,  cosTh  ));
				normals.add(new Vector(  sinThD*cosPh  ,  sinThD*sinPh  ,  cosThD  ));
				normals.add(new Vector(  sinTh*cosPhD  ,  sinTh*sinPhD  ,  cosTh  ));
				normals.add(new Vector(  sinThD*cosPhD  ,  sinThD*sinPhD  ,  cosThD  ));
				normals.add(new Vector(  sinTh*cosPhD  ,  sinTh*sinPhD  ,  cosTh  ));
				normals.add(new Vector(  sinThD*cosPh  ,  sinThD*sinPh  ,  cosThD  ));
			}
		}

		int i=0;
		float[] vertArr = new float[verts.size()*3];
		for(Point v: verts){ vertArr[i++] = (float)v.x(); vertArr[i++] = (float)v.y(); vertArr[i++] = (float)v.z(); }
		i=0;
		float[] normArr = new float[normals.size()*3];
		for(Vector v: normals){ normArr[i++] = (float)v.x(); normArr[i++] = (float)v.y(); normArr[i++] = (float)v.z(); }

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