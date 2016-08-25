# ProGAL: Geometric Algorithm Library

ProGAL is a Java-library containing data structures and implementations of a range of geometric algorithms. Many of the implementations are based on research projects at Department of Computer Science, University of Copenhagen.

[ProGAL.jar](ProGAL.jar)

## Usage
To use ProGAL include the jar in the class path of your project. For example, let ProGALTest.java have the following contents:
```
import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.geom2d.*;

class ProGALTest{
   public static void main(String[] args){
      Circle c = new Circle(new Point(0,0), 1);
      J2DScene scene = J2DScene.createJ2DSceneInFrame();
      scene.addShape(c, java.awt.Color.BLUE);
   }
}
```
From the command-line run:
```
javac -cp ProGAL.jar:. ProGALTest.java
java -cp ProGAL.jar:. ProGALTest
```
to run a small test. Using eclipse, simply create a project and go to the Project-menu -> Properties -> Java Build Path -> Libraries and locate the ProGAL.jar file after pressing the 'Add external jar' button. 

If you wish to use the 3D visualization classes (J3DScene) you need to have Java3D installed. Java3D is included in Java 6 but disappeared from Java 7 and up. To temporarily change the java-version, type:
```
export JAVA_HOME=/System/Library/Frameworks/JavaVM.framework/Versions/1.6.0/Home/ #Mac
export JAVA_HOME=/usr/lib/jvm/java-6-sun/ #Ubuntu
```
we're working on an updated 3D viewer using JavaFX.
     
## Documentation

Documentation of ProGAL is maintained partly through a [JavaDoc homepage](link) and partly through a [collection of examples](link) that demonstrate the use of the library.

Projects that contributed to ProGAL:
* Computational Geometry and Bioinformatics. Peter & Henrik Sterner. Masters thesis. 2008
* Rectangular Swept Spheres. Mikkel Kjær Jensen. 2010
* Fitting an All-atom Protein Model to a Ca-trace. Martin Dybdal, Anders Boesen Lindbo Larsen & Esben Skaarup. 2010
* Applying Inverse Kinematics to Adjustable ChainTrees. Hans-Kristian Bjerregaard. Bachelor project. 2011
* Void-trees. Desirée Malene Schreyer Jørgensen. Bachelor project. 2012
* Bowyer-Watson any-dimensional Delaunay tessellations. Desirée Malene Schreyer Jørgensen & Anni Jane Pinder. 2012
* Range-trees. Søren Lynnerup. 2012
* Kinetic alpha-complexs. Desirée Malene Schreyer Jørgensen. Masters thesis. 2014
