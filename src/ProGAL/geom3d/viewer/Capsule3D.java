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
public class Capsule3D extends Shape3D {
	private int divisions = 200;

	/** Construct the capsule shape. 
	 * @param height The distance between the defining points. Notice that the actual height 
	 * of the capsule will be <code>height+2*radius</code>
	 * @param radius The radius of the capsule. 
	 * @param app Appearance of the capsule. 
	 */
	public Capsule3D(float height, float radius, Appearance app) {
		super();

		List<Point> verts = new LinkedList<Point>();
		List<Vector> normals = new LinkedList<Vector>();

		//Upper half-sphere
		Point c = new Point(0,height/2,0);
		float r = radius;
		double d = PI/divisions;

		for(double theta=0;theta<PI;theta+=d){
			for(double phi=0;phi<PI;phi+=d){
				verts.add(new Point(  c.getX() + r*sin(theta  )*cos(phi  )  ,  c.getY() + r*sin(theta  )*sin(phi  )  ,  c.getZ() + r*cos(theta  )  ));
				verts.add(new Point(  c.getX() + r*sin(theta+d)*cos(phi  )  ,  c.getY() + r*sin(theta+d)*sin(phi  )  ,  c.getZ() + r*cos(theta+d)  ));
				verts.add(new Point(  c.getX() + r*sin(theta  )*cos(phi+d)  ,  c.getY() + r*sin(theta  )*sin(phi+d)  ,  c.getZ() + r*cos(theta  )  ));
                                                                                   
				verts.add(new Point(  c.getX() + r*sin(theta+d)*cos(phi+d)  ,  c.getY() + r*sin(theta+d)*sin(phi+d)  ,  c.getZ() + r*cos(theta+d)  ));
				verts.add(new Point(  c.getX() + r*sin(theta  )*cos(phi+d)  ,  c.getY() + r*sin(theta  )*sin(phi+d)  ,  c.getZ() + r*cos(theta  )  ));
				verts.add(new Point(  c.getX() + r*sin(theta+d)*cos(phi  )  ,  c.getY() + r*sin(theta+d)*sin(phi  )  ,  c.getZ() + r*cos(theta+d)  ));
				
				normals.add(new Vector(  sin(theta  )*cos(phi  )  ,  sin(theta  )*sin(phi  )  ,  cos(theta  )  ));
				normals.add(new Vector(  sin(theta+d)*cos(phi  )  ,  sin(theta+d)*sin(phi  )  ,  cos(theta+d)  ));
				normals.add(new Vector(  sin(theta  )*cos(phi+d)  ,  sin(theta  )*sin(phi+d)  ,  cos(theta  )  ));

				normals.add(new Vector(  sin(theta+d)*cos(phi+d)  ,  sin(theta+d)*sin(phi+d)  ,  cos(theta+d)  ));
				normals.add(new Vector(  sin(theta  )*cos(phi+d)  ,  sin(theta  )*sin(phi+d)  ,  cos(theta  )  ));
				normals.add(new Vector(  sin(theta+d)*cos(phi  )  ,  sin(theta+d)*sin(phi  )  ,  cos(theta+d)  ));
				//normals.add(c.vectorTo(new Vector(sin(theta+d/2)*cos(phi+d/2)  ,  c.y() + r*sin(theta+d/2)*sin(phi+d/2)  ,  c.z() + r*cos(theta+d/2))));
			}
		}
		//Lower half-sphere
		c = new Point(0,-height/2,0);
		for(double theta=0;theta<PI;theta+=d){
			for(double phi=PI;phi<PI*2;phi+=d){
				verts.add(new Point(  c.getX() + r*sin(theta  )*cos(phi  )  ,  c.getY() + r*sin(theta  )*sin(phi  )  ,  c.getZ() + r*cos(theta  )  ));
				verts.add(new Point(  c.getX() + r*sin(theta+d)*cos(phi  )  ,  c.getY() + r*sin(theta+d)*sin(phi  )  ,  c.getZ() + r*cos(theta+d)  ));
				verts.add(new Point(  c.getX() + r*sin(theta  )*cos(phi+d)  ,  c.getY() + r*sin(theta  )*sin(phi+d)  ,  c.getZ() + r*cos(theta  )  ));

				verts.add(new Point(  c.getX() + r*sin(theta+d)*cos(phi+d)  ,  c.getY() + r*sin(theta+d)*sin(phi+d)  ,  c.getZ() + r*cos(theta+d)  ));
				verts.add(new Point(  c.getX() + r*sin(theta  )*cos(phi+d)  ,  c.getY() + r*sin(theta  )*sin(phi+d)  ,  c.getZ() + r*cos(theta  )  ));
				verts.add(new Point(  c.getX() + r*sin(theta+d)*cos(phi  )  ,  c.getY() + r*sin(theta+d)*sin(phi  )  ,  c.getZ() + r*cos(theta+d)  ));
				
				normals.add(new Vector(  sin(theta+d)*cos(phi  )  ,  sin(theta+d)*sin(phi  )  ,  cos(theta+d)  ));
				normals.add(new Vector(  sin(theta  )*cos(phi  )  ,  sin(theta  )*sin(phi  )  ,  cos(theta  )  ));
				normals.add(new Vector(  sin(theta  )*cos(phi+d)  ,  sin(theta  )*sin(phi+d)  ,  cos(theta  )  ));

				normals.add(new Vector(  sin(theta+d)*cos(phi+d)  ,  sin(theta+d)*sin(phi+d)  ,  cos(theta+d)  ));
				normals.add(new Vector(  sin(theta  )*cos(phi+d)  ,  sin(theta  )*sin(phi+d)  ,  cos(theta  )  ));
				normals.add(new Vector(  sin(theta+d)*cos(phi  )  ,  sin(theta+d)*sin(phi  )  ,  cos(theta+d)  ));
				//normals.add(c.vectorTo(new Vector(sin(theta+d/2)*cos(phi+d/2)  ,  c.y() + r*sin(theta+d/2)*sin(phi+d/2)  ,  c.z() + r*cos(theta+d/2))));
			}
		}

		//Cylinder
		double dH = height;
		for(double theta=0;theta<PI*2;theta+=d){
			double y=-height/2;
			verts.add(new Point(  r*sin(theta  )  , y   , r*cos(theta  )  ));
			verts.add(new Point(  r*sin(theta+d)  , y   , r*cos(theta+d)  ));
			verts.add(new Point(  r*sin(theta  )  , y+dH, r*cos(theta  )  ));
			verts.add(new Point(  r*sin(theta+d)  , y+dH, r*cos(theta+d)  ));
			verts.add(new Point(  r*sin(theta  )  , y+dH, r*cos(theta  )  ));
			verts.add(new Point(  r*sin(theta+d)  , y   , r*cos(theta+d)  ));
			normals.add(new Vector(  sin(theta  )  , 0 , cos(theta  )  ));
			normals.add(new Vector(  sin(theta+d)  , 0 , cos(theta+d)  ));
			normals.add(new Vector(  sin(theta  )  , 0 , cos(theta  )  ));
			normals.add(new Vector(  sin(theta+d)  , 0 , cos(theta+d)  ));
			normals.add(new Vector(  sin(theta  )  , 0 , cos(theta  )  ));
			normals.add(new Vector(  sin(theta+d)  , 0 , cos(theta+d)  ));
		}

		int i=0;
		double[] vertArr = new double[verts.size()*3];
		for(Point v: verts){ vertArr[i++] = v.getX(); vertArr[i++] = v.getY(); vertArr[i++] = v.getZ(); }
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