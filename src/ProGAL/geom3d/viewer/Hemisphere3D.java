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
	 * @param topAngle Indicates the starting angle of the phi parameter. <code>topAngle==0</code> will 
	 *  result in a complete sphere and <code>topAngle==Math.PI/2</code> is a half-sphere.  
	 * @param app Appearance of the hemisphere. 
	 */
	public Hemisphere3D(float radius, double topAngle, Appearance app, int divisions) {
		super();
		
		List<Point> verts = new LinkedList<Point>();
		List<Vector> normals = new LinkedList<Vector>();

		float r = radius;
		double d = PI/divisions;
		double cosPh0 = cos(topAngle);
		double sinPh0 = sin(topAngle);

		for(double theta=0;theta<2*PI;theta+=d){
			double sinTh = sin(theta);
			double sinThD = sin(theta+d);
			double cosTh = cos(theta);
			double cosThD = cos(theta+d);
			for(double phi=topAngle;phi<PI;phi+=d){
				double cosPh = cos(phi);
				double cosPhD = cos(phi+d);
				double sinPh = sin(phi);
				double sinPhD = sin(phi+d);
				
				verts.add(new Point(  r*cosTh *sinPh /sinPh0  ,  (-r*cosPh +cosPh0)/(1+cosPh0)  ,  r*sinTh *sinPh /sinPh0  ));
				verts.add(new Point(  r*cosTh *sinPhD/sinPh0  ,  (-r*cosPhD+cosPh0)/(1+cosPh0)  ,  r*sinTh *sinPhD/sinPh0  ));
				verts.add(new Point(  r*cosThD*sinPh /sinPh0  ,  (-r*cosPh +cosPh0)/(1+cosPh0)  ,  r*sinThD*sinPh /sinPh0  ));
				verts.add(new Point(  r*cosThD*sinPh /sinPh0  ,  (-r*cosPh +cosPh0)/(1+cosPh0)  ,  r*sinThD*sinPh /sinPh0  ));
				verts.add(new Point(  r*cosTh *sinPhD/sinPh0  ,  (-r*cosPhD+cosPh0)/(1+cosPh0)  ,  r*sinTh *sinPhD/sinPh0  ));
				verts.add(new Point(  r*cosThD*sinPhD/sinPh0  ,  (-r*cosPhD+cosPh0)/(1+cosPh0)  ,  r*sinThD*sinPhD/sinPh0  ));
				

				normals.add(new Vector(  cosTh *sinPh /sinPh0  ,  (-cosPh +cosPh0)/(1+cosPh0)  ,  sinTh *sinPh /sinPh0  ));
				normals.add(new Vector(  cosTh *sinPhD/sinPh0  ,  (-cosPhD+cosPh0)/(1+cosPh0)  ,  sinTh *sinPhD/sinPh0  ));
				normals.add(new Vector(  cosThD*sinPh /sinPh0  ,  (-cosPh +cosPh0)/(1+cosPh0)  ,  sinThD*sinPh /sinPh0  ));
				normals.add(new Vector(  cosThD*sinPh /sinPh0  ,  (-cosPh +cosPh0)/(1+cosPh0)  ,  sinThD*sinPh /sinPh0  ));
				normals.add(new Vector(  cosTh *sinPhD/sinPh0  ,  (-cosPhD+cosPh0)/(1+cosPh0)  ,  sinTh *sinPhD/sinPh0  ));
				normals.add(new Vector(  cosThD*sinPhD/sinPh0  ,  (-cosPhD+cosPh0)/(1+cosPh0)  ,  sinThD*sinPhD/sinPh0  ));
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