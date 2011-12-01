package ProGAL.geom3d.complex.alphaComplex;


public class TetrahedronAlphaProperties implements SimplexAlphaProperties {
	private Interval interior;
    private double inAlphaComplex;

    public TetrahedronAlphaProperties(double rho){
    	this.interior = new Interval(rho,Double.POSITIVE_INFINITY);
    	this.inAlphaComplex = rho;
    }
    
    
    public double getInAlphaComplex(){ 		return inAlphaComplex; 		}
    public boolean isAttached(){ 			return false; 				}

    public Interval getInteriorInterval(){	return interior; 			}

	public int getSimplexType() {	return 3;	}
}
