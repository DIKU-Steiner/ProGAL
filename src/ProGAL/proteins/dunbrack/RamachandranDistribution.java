package ProGAL.proteins.dunbrack;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

import ProGAL.io.IOToolbox;
import ProGAL.math.Randomization;

/**
 * An aminoacid dependent distribution of Ramachandran angles based on the Dunbrack TCBIG dataset. 
 *  
 * @author R.Fonseca
 *
 */
public class RamachandranDistribution {
	private static RamachandranDistribution instance;//Singleton pattern.

	public static synchronized RamachandranDistribution getDistribution(){
		if(instance==null) instance = new RamachandranDistribution();

		return instance;
	}


//	public double getProbSumLeftAll(int aaType, double phi, double psi){
//		int phiBin = (int)((phi+180)/5.0);
//		int psiBin = (int)((psi+180)/5.0);
//		if(phiBin>=72) phiBin=71;
//		if(psiBin>=72) psiBin=71;
//		return probSumsLeftAll[aaType][phiBin*72+psiBin];
//	}
	public double[] samplePhiPsi(int aaType){
		return binToTorsions(sample(probSumsLeftAll[aaType]));
	}
	private static int sample(double[] arr){
		double rand = Randomization.randBetween(0.0, 1.0);
		int low = 0;
        int high = arr.length - 1;
        int mid=0;

        while( low <= high ){
            mid = ( low + high ) / 2;
            if( arr[mid]<rand )			low = mid + 1;
            else if( arr[mid]>rand )	high = mid - 1;
            else						return mid;
        }
        return low;
	}

	private static double[] binToTorsions(int bin){
		int phiBin = bin/72;
		int psiBin = bin%72;
		double phi = phiBin*5-180+Randomization.randBetween(0.0, 5.0);
		double psi = psiBin*5-180+Randomization.randBetween(0.0, 5.0);
		return new double[]{phi,psi};
	}
		

	private File dunbrackDir = null;
	
	private List<String> types = new ArrayList<String>(20);
	//left all
	private double[][] probSumsLeftAll;
	
	private RamachandranDistribution(){
		String[] configLines;
		try{
			configLines = IOToolbox.readFromFile(System.getProperty("user.home")+File.separatorChar+".progal").split("\n");
		}catch(Exception exc){
			throw new RuntimeException("Error creating distribution: The file '~/.progal' must specify a full directory " +
					"with the Dunbrack datasets, e.g.:\ndunbrackDir = /Users/me/Datasets/DunbrackBackbone/\n" +
					"where DunbrackBackbone contains the file 'NDRD_TCBIG.txt.gz'.");
		}

		//Read dunbrackDir
		for(String l: configLines){
			String line = l.trim();
			if(line.startsWith("dunbrackDir")) {
				String dirString = line.substring(line.indexOf('=')+1).trim();
				dunbrackDir = new File(dirString);
				if(!dunbrackDir.exists() || !dunbrackDir.isDirectory()){
					throw new RuntimeException("Error creating distribution: Directory '"+dunbrackDir+"' does not exist");
				}
			}
		}

		String fName = dunbrackDir.getAbsolutePath()+"/NDRD_TCBIG.txt.gz";
		BufferedReader in;
		try {
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream( fName ))));

			String line;
			probSumsLeftAll = new double[21][72*72];
			int aa = 0,c=0;
			while( (line=in.readLine())!=null ){
				if(line.startsWith("#") || line.trim().isEmpty()) continue;
				if(		line.substring(4,8).equalsIgnoreCase("left") && 
						line.substring(10,13).equalsIgnoreCase("ALL")){
					double probSum = Double.parseDouble(line.substring(52,64));
					probSumsLeftAll[aa][c++] = probSum; 
					if(c==72*72) {
						types.add(line.substring(0,3));
						c=0;
						aa++;
					}
				}
				
			}

			in.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	public List<String> getTypes(){ return types; }


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		RamachandranDistribution distr = RamachandranDistribution.getDistribution();
		for(int i=0;i<1000;i++){
			double[] sample = distr.samplePhiPsi(20);
			System.out.printf("%.2f %.2f\n",sample[0],sample[1]);
		}
	}

}
