package Correlescence;

import java.util.ArrayList;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class KymographGenerator implements PlugIn {

	/**kymo width, px**/
	int imWidth;
	/**kymo height, frames**/
	int imHeight;
	/** SD of PSF**/
	double dPSF;
	/** noise offset **/
	double dNoiseOffset;
	/** noise standard deviation **/
	double dNoiseSD;
	/** number of particles populations **/
	int nPartGroupsN;
	/** array containing particles groups**/
	ArrayList<ParticleGroup> particlegroups;
	/** final kymograph image **/
	FloatProcessor finalip;
	
	@Override
	public void run(String arg) {
		
		int i,j;
		int nPGN;
		
		ParticleGroup temppG;
		double [] nVelParams;
		particlegroups = new ArrayList<ParticleGroup>();
		//General parameters
		if(!kgGeneralDialog())
			return;
		//each group parameters
		for (i=0;i<nPartGroupsN;i++)
		{
			//TODO something with return of value
			temppG=kgGroupDialog(i+1);
			temppG.imWidth=imWidth;
			if(temppG.nParticleN==0)
				return;
			else
			{
				if(temppG.nVComponentsN>0)
				{
					//ini vel parameters
					temppG.dSpeedM=new double [temppG.nVComponentsN];
					temppG.dSpeedSD=new double [temppG.nVComponentsN];
					temppG.dSpeedDurM=new double [temppG.nVComponentsN];
					temppG.dSpeedDurSD=new double [temppG.nVComponentsN];
					temppG.dComponentsP=new double [temppG.nVComponentsN];
					for(j=0;j<temppG.nVComponentsN;j++)
					{
						nVelParams=kgVelocityComponentsDialog(i+1,j+1);
						if(nVelParams==null)
							return;
						temppG.dSpeedM[j]=nVelParams[0];
						temppG.dSpeedSD[j]=nVelParams[1];
						temppG.dSpeedDurM[j]=nVelParams[2];
						temppG.dSpeedDurSD[j]=nVelParams[3];
						//TODO add interface
						//now all components probability is equal
						temppG.dComponentsP[j]=1.0/(double)temppG.nVComponentsN;
						
					}
				}
			}
			particlegroups.add(temppG);
		}
		finalip = new FloatProcessor(imWidth,imHeight);
		finalip.noise(dNoiseSD);
		finalip.add(dNoiseOffset);
		
		//initialize groups
		for(nPGN=0;nPGN<nPartGroupsN;nPGN++)
		{
			particlegroups.get(nPGN).init();
		}
		//render kymos
		for (j=0;j<imHeight;j++)
		{
			for(nPGN=0;nPGN<nPartGroupsN;nPGN++)
			{
				
				temppG=particlegroups.get(nPGN);
				//update positions
				temppG.update();
				//update row
				addLine(temppG.getCurrentLine(dPSF),j);
			}
			
		}
		
		new ImagePlus("art_kymo", finalip).show();
		IJ.run("Enhance Contrast", "saturated=0.35");
	}
	
	void addLine(float [] patInt, int nRow)
	{

		for(int i=0;i<imWidth;i++)
		{
			finalip.setf(i, nRow, finalip.getf(i,nRow)+patInt[i]);
		}
		
	}
	
	/** 
	 * Dialog displaying general parameters
	 * **/
	public boolean kgGeneralDialog() 
	{
		GenericDialog kgGDialog = new GenericDialog("General artificial kymograph parameters");
		kgGDialog.addNumericField("Image width", Prefs.get("KymographGenerator.imWidth", 100), 0, 5, " pixels");
		kgGDialog.addNumericField("Image height", Prefs.get("KymographGenerator.imHeight", 100), 0, 5, " frames");
		kgGDialog.addNumericField("SD of PSF", Prefs.get("KymographGenerator.dPSF", 1.5), 2, 5, " pixels");
		kgGDialog.addNumericField("Background noise offset", Prefs.get("KymographGenerator.dNoiseOffset", 10), 0, 5, " a.u. (intensity)");
		kgGDialog.addNumericField("Noise SD", Prefs.get("KymographGenerator.dNoiseSD", 3), 0, 5, " a.u. (intensity)");
		kgGDialog.addNumericField("Number of particles groups", Prefs.get("KymographGenerator.nPartGroupsN", 1), 0, 1, " ");
		kgGDialog.setResizable(false);
		kgGDialog.showDialog();
		if (kgGDialog.wasCanceled())
	        return false;
		
		imWidth = (int)kgGDialog.getNextNumber();
		Prefs.set("KymographGenerator.imWidth", imWidth);
		imHeight = (int)kgGDialog.getNextNumber();
		Prefs.set("KymographGenerator.imHeight", imHeight);
		dPSF = kgGDialog.getNextNumber();
		Prefs.set("KymographGenerator.dPSF", dPSF);
		dNoiseOffset = kgGDialog.getNextNumber();
		Prefs.set("KymographGenerator.dNoiseOffset", dNoiseOffset);
		dNoiseSD = kgGDialog.getNextNumber();
		Prefs.set("KymographGenerator.dNoiseSD", dNoiseSD);	
		
		nPartGroupsN = (int)kgGDialog.getNextNumber();
		Prefs.set("KymographGenerator.nPartGroupsN", nPartGroupsN);	
		return true;
	}
	
	/** 
	 * Dialog displaying parameters for one particle group
	 * and in case of "OK"
	 * adding Group to the main Group array
	 * **/
	public ParticleGroup kgGroupDialog(int nGroup) 
	{
		//int nPatN;
		//double dDiff;
		ParticleGroup pG=new ParticleGroup();
		
		String [] sIniPos = new String [] {
				"Random","at x specified below"};
		String [] sBoundary = new String [] {
				"Ignore","Reappears from another side","Reappears from the same side","Reappears at random location"};
		String [] sEndOfLifeTime = new String [] {
				"Disappears","Reappears from the left side","Reappears from the right side", "Reappears at random location"};
		GenericDialog kgGrDialog = new GenericDialog("Parameters for Group "+Integer.toString(nGroup));
		kgGrDialog.addNumericField("Initial number of particles", Prefs.get("KymographGenerator.nParticleN"+Integer.toString(nGroup), 10), 0, 5, " ");
		kgGrDialog.addChoice("Particles initial distribution:", sIniPos, Prefs.get("KymographGenerator.iniDistr"+Integer.toString(nGroup), "Random"));
		kgGrDialog.addNumericField("X coord initial distribution", Prefs.get("KymographGenerator.nXiniDistr"+Integer.toString(nGroup), 0), 1, 5, " px");
		kgGrDialog.addNumericField("Intensity mean", Prefs.get("KymographGenerator.dIntensityM"+Integer.toString(nGroup), 0), 2, 5, " a.u.");
		kgGrDialog.addNumericField("SD of intensity", Prefs.get("KymographGenerator.dIntensitySD"+Integer.toString(nGroup), 0), 2, 5, " a.u.");
		kgGrDialog.addNumericField("Diffusion coeff (0=no diffusion)", Prefs.get("KymographGenerator.dDiffusion"+Integer.toString(nGroup), 0), 2, 5, " px^2/frame");
		kgGrDialog.addNumericField("Average lifetime (0=infinite)", Prefs.get("KymographGenerator.dLifeTimeMean"+Integer.toString(nGroup), 0), 0, 5, " frames");		
		//kgGrDialog.addNumericField("Lifetime SD", Prefs.get("KymographGenerator.dLifeTimeSD"+Integer.toString(nGroup), 0), 0, 5, " frames");
		kgGrDialog.addChoice("At the end of lifetime, particle:", sEndOfLifeTime, Prefs.get("KymographGenerator.sEndOfLifeTime"+Integer.toString(nGroup), "Disappears"));
		kgGrDialog.addNumericField("Mean pause duration (0=no pauses)", Prefs.get("KymographGenerator.dPauseM"+Integer.toString(nGroup), 0), 0, 5, " frames");
		kgGrDialog.addNumericField("Pause duration SD", Prefs.get("KymographGenerator.dPauseSD"+Integer.toString(nGroup), 0), 0, 5, " frames");
		kgGrDialog.addNumericField("Probability of pause", Prefs.get("KymographGenerator.dPauseP"+Integer.toString(nGroup), 0), 2, 5, " ");
		kgGrDialog.addNumericField("Velocities components number (0=absent)", Prefs.get("KymographGenerator.nVComponentsN"+Integer.toString(nGroup), 0), 0, 5, " ");
		kgGrDialog.addChoice("Boundary condition for particle:", sBoundary, Prefs.get("KymographGenerator.sBoundary"+Integer.toString(nGroup), "Ignore"));
		kgGrDialog.setResizable(false);
		kgGrDialog.showDialog();
		if (kgGrDialog.wasCanceled())
		{
			pG.nParticleN=0;
	        return pG;
		}
		pG.nParticleN=(int)kgGrDialog.getNextNumber();
		Prefs.set("KymographGenerator.nParticleN"+Integer.toString(nGroup), pG.nParticleN);
		pG.iniDistr= kgGrDialog.getNextChoiceIndex();
		Prefs.set("KymographGenerator.iniDistr"+Integer.toString(nGroup), sIniPos[pG.iniDistr]);
		pG.nXiniDistr=kgGrDialog.getNextNumber();
		if(pG.nXiniDistr<0)
			pG.nXiniDistr=0;
		if(pG.nXiniDistr>(imWidth-1))
			pG.nXiniDistr=(imWidth-1);
		Prefs.set("KymographGenerator.nXiniDistr"+Integer.toString(nGroup), pG.nXiniDistr);
		pG.dIntensityM=kgGrDialog.getNextNumber();
		Prefs.set("KymographGenerator.dIntensityM"+Integer.toString(nGroup), pG.dIntensityM);
		pG.dIntensitySD=kgGrDialog.getNextNumber();
		Prefs.set("KymographGenerator.dIntensitySD"+Integer.toString(nGroup), pG.dIntensitySD);
		pG.dDiffusion=kgGrDialog.getNextNumber();
		Prefs.set("KymographGenerator.dDiffusion"+Integer.toString(nGroup), pG.dDiffusion);
		pG.dLifeTimeMean=kgGrDialog.getNextNumber();
		Prefs.set("KymographGenerator.dLifeTimeMean"+Integer.toString(nGroup), pG.dLifeTimeMean);
		//pG.dLifeTimeSD=kgGrDialog.getNextNumber();
		//Prefs.set("KymographGenerator.dLifeTimeSD"+Integer.toString(nGroup), pG.dLifeTimeSD);
		pG.nEndOfLifeTime= kgGrDialog.getNextChoiceIndex();
		Prefs.set("KymographGenerator.sEndOfLifeTime"+Integer.toString(nGroup), sBoundary[pG.nEndOfLifeTime]);
		pG.dPauseM=kgGrDialog.getNextNumber();
		Prefs.set("KymographGenerator.dPauseM"+Integer.toString(nGroup), pG.dPauseM);
		pG.dPauseSD=kgGrDialog.getNextNumber();
		Prefs.set("KymographGenerator.dPauseSD"+Integer.toString(nGroup), pG.dPauseSD);
		pG.dPauseP=kgGrDialog.getNextNumber();
		Prefs.set("KymographGenerator.dPauseP"+Integer.toString(nGroup), pG.dPauseP);
		pG.nVComponentsN=(int)kgGrDialog.getNextNumber();
		Prefs.set("KymographGenerator.nVComponentsN"+Integer.toString(nGroup), pG.nVComponentsN);	
		pG.nBoundary= kgGrDialog.getNextChoiceIndex();
		Prefs.set("KymographGenerator.sBoundary"+Integer.toString(nGroup), sBoundary[pG.nBoundary]);
		return pG;
	}
	/** 
	 * Dialog displaying parameters for one particle group/velocity components
	 * returns array of values;
	 * 	 * **/
	public double [] kgVelocityComponentsDialog(int nGroup, int nComponent) 
	{
		double [] nVelParam=null;
		GenericDialog kgVelDialog = new GenericDialog("Parameters for Group "+Integer.toString(nGroup)+" Velocity component "+Integer.toString(nComponent));
		kgVelDialog.addMessage("Velocity component "+Integer.toString(nComponent)+" distribution");
		kgVelDialog.addNumericField("Average speed", Prefs.get("KymographGenerator.dSpeedM"+Integer.toString(nGroup)+Integer.toString(nComponent), 0), 2, 5, " px/frame");
		kgVelDialog.addNumericField("Speed SD", Prefs.get("KymographGenerator.dSpeedSD"+Integer.toString(nGroup)+Integer.toString(nComponent), 0), 2, 5, " px/frame");
		kgVelDialog.addMessage("Run's duration distribution for component "+Integer.toString(nComponent));
		kgVelDialog.addNumericField("Duration mean", Prefs.get("KymographGenerator.dPauseM"+Integer.toString(nGroup)+Integer.toString(nComponent), 0), 2, 5, " frames");
		kgVelDialog.addNumericField("SD of duration", Prefs.get("KymographGenerator.dPauseSD"+Integer.toString(nGroup)+Integer.toString(nComponent), 0), 2, 5, " frames");
		kgVelDialog.setResizable(false);
		kgVelDialog.showDialog();

		
		if (kgVelDialog.wasCanceled())
		{
	        return nVelParam;
		}
		
		nVelParam=new double [4];
		nVelParam[0]=kgVelDialog.getNextNumber();
		Prefs.set("KymographGenerator.dSpeedM"+Integer.toString(nGroup)+Integer.toString(nComponent), nVelParam[0]);
		nVelParam[1]=kgVelDialog.getNextNumber();
		Prefs.set("KymographGenerator.dSpeedSD"+Integer.toString(nGroup)+Integer.toString(nComponent), nVelParam[1]);
		nVelParam[2]=kgVelDialog.getNextNumber();
		Prefs.set("KymographGenerator.dPauseM"+Integer.toString(nGroup)+Integer.toString(nComponent), nVelParam[2]);
		nVelParam[3]=kgVelDialog.getNextNumber();
		Prefs.set("KymographGenerator.dPauseSD"+Integer.toString(nGroup)+Integer.toString(nComponent), nVelParam[3]);
		return nVelParam;
		
	}
}
