package Correlescence;

import java.util.ArrayList;
import java.util.Random;

public class ParticleGroup {
	/** total current number of particles **/
	int nParticleN;
	
	/** Mean Intensity**/
	double dIntensityM;
	
	/** Intensity SD**/
	double dIntensitySD;
	
	/** diffusion coefficient (px^2)/frame**/
	double dDiffusion;
	
	/** sqrt(2*diffusion coefficient)**/
	double dDiffusionSD;
	
	/**number of velocity components**/
	int nVComponentsN;
	
	/** velocities  mean**/
	double [] dSpeedM;
	
	/** velocities sd**/
	double [] dSpeedSD;
	
	/** velocity duration mean**/
	double [] dSpeedDurM;
	
	/** velocities duration SD**/
	double [] dSpeedDurSD;
	
	/**probability to pick speed from one or another component**/
	double [] dComponentsP;
	
	/** pauses duration mean**/
	double dPauseM;
	
	/** pauses duration SD**/
	double dPauseSD;
	
	/** pauses to movement probability**/
	double dPauseP;
	
	/**average particle lifetime 0=infinite**/
	double dLifeTimeMean;
	
	/** SD of particle lifetime **/
	//double dLifeTimeSD;
	
	/** At the end of lifetime particle
	 * 0=Disappears
	 * 1=Reappears from the left side
	 * 2=Reappears from the right side
	 * 3=Reappears at random location}; 
	 * **/
	int nEndOfLifeTime;
			
	/** initial particles distribution 
	 * 0=random, 
	 * 1=all at nXiniDistr**/
	int iniDistr;
	
	/**initial position of particles in case of non-random distribution**/
	double nXiniDistr;
	
	/**how to handle boundary: 
	 * 0=Ignore
	 * 1=Reappears from another side
	 * 2=Reappears from the same side
	 * 3=Reappears at random location**/
	int nBoundary;
	
	/** image width**/
	int imWidth;
	
	/** particles **/
	ArrayList<Particle1D> particles;
	
	public ParticleGroup() 
	{
		nParticleN=0;
	}
	
	/** initializes new group**/
	public void init()
	{
		int i;
		Particle1D tempP;
		particles=new ArrayList<Particle1D>();
		Random rnd=new Random();
		//some intermediate parameters
		dDiffusionSD = Math.sqrt(2.0*dDiffusion);
		
		for(i=0;i<nParticleN;i++)
		{
			tempP=new Particle1D();
			//initial position
			if(iniDistr==0)
				tempP.x=rnd.nextDouble()*(imWidth-1);
			else
				tempP.x=nXiniDistr;
			//intensity
			if(dIntensitySD==0.0)
				tempP.dInt=dIntensityM;
			else
				tempP.dInt=dIntensityM+(rnd.nextGaussian()*dIntensitySD);
			
			tempP.bVisible=true;
			//state
			assignState(tempP,true);
			assignLifeTime(tempP,true);
			
			particles.add(tempP);
		}	
	}
	
	/** updates current state of particles **/
	public void update()
	{
		int i;
		Particle1D tempP;
		Random rnd=new Random();
		boolean bLeft;
		for(i=0;i<nParticleN;i++)
		{
			tempP=particles.get(i);
			
			//increase time
			//=decrease lifetime
			tempP.dStateTimer--;
			
			if(dLifeTimeMean>0)
			{
				tempP.dLTTimer--;
			}
			
			//state update
			if(tempP.dStateTimer<0.0)
			{

				//state
				assignState(tempP,false);
				
			}
			
			//diffusion
			if(dDiffusion>0.0)
			{
				tempP.x+=((rnd.nextInt(2)*2)-1)*(rnd.nextGaussian()*dDiffusionSD);
			}
			
			//speed
			if(nVComponentsN>0 && tempP.nState==1)
				tempP.x+=tempP.v;
			
			bLeft = false;
			//verify boundary conditions
			if(tempP.x<0 || tempP.x>(imWidth-1))
			{
				if(tempP.x<0)
					bLeft=true;
				switch (nBoundary)
				{
					case 0:
						//do nothing
						break;
						//appear from opposite side
					case 1:
						if(bLeft)
							tempP.x=(imWidth-1);
						else
							tempP.x=0;
						assignState(tempP,false);
						assignLifeTime(tempP,false);
						break;
						//appear from the same side
					case 2:
						if(bLeft)
							tempP.x=0;
						else
							tempP.x=(imWidth-1);
						assignState(tempP,false);
						assignLifeTime(tempP,false);
						break;

					case 3:
						//appear at random location
						tempP.x=rnd.nextDouble()*(imWidth-1);
						assignState(tempP,false);
						assignLifeTime(tempP,false);
						break;
				}
			}
			
			
			//check for lifetime
			if(dLifeTimeMean>0 && tempP.dLTTimer<0)
			{
			   switch(nEndOfLifeTime)
			   {
			   case 0:
				//do nothing
					break;
				//appear from left side
				case 1:
					tempP.x=0;
					assignState(tempP,false);
					assignLifeTime(tempP,false);
					break;
				//appear from right side
				case 2:
					tempP.x=(imWidth-1);
					assignState(tempP,false);
					assignLifeTime(tempP,false);
					break;
				//appear at random location
				case 3:
					tempP.x=rnd.nextDouble()*(imWidth-1);
					assignState(tempP,false);
					assignLifeTime(tempP,false);
					break;
			   }
			}
		}
	}
	
	/** returns convolution of current particles states with provided PSF**/
	public float [] getCurrentLine(double dSDPSF)
	{
		float [] finline = new float[imWidth];
		int nPatN;
		int i;
		double dSDdenum = 1.0/(2.0*dSDPSF*dSDPSF);
		Particle1D tempP;
		for(nPatN=0;nPatN<nParticleN;nPatN++)
		{
			tempP=particles.get(nPatN);
			if(tempP.bVisible)
			{
				//very simple case for now
				for(i=0;i<imWidth;i++)
				{
					finline[i]+=tempP.dInt*Math.exp(-Math.pow(tempP.x-i+0.5, 2)*dSDdenum);
				}
			}
			
		}
		
		return finline;
		
	}
	
	
	public void assignLifeTime(Particle1D mP, boolean bIni)
	{
		Random rnd=new Random();
		//lifetime is exponentially distributed
		if(dLifeTimeMean>0)
		{
			if(bIni)
				mP.dLTTimer=(-1.0)*dLifeTimeMean*Math.log(rnd.nextDouble())*rnd.nextDouble();
			
			else
				mP.dLTTimer=(-1.0)*dLifeTimeMean*Math.log(rnd.nextDouble());
		}
	}
	
	/** sets  particles state for particle
	 * boolean parameter defines whether it is done
	 * during initialization = shift in lifetime/state time **/
	public void assignState(Particle1D mP, boolean bIni)
	{
		Random rnd=new Random();
		if(dPauseM>0.0)
		{
			//pause state
			if(rnd.nextDouble()<dPauseP)
				setPauseState(mP,bIni);
			//run state
			else
				setMovingState(mP,bIni);
		}
		else
			setMovingState(mP,bIni);
		
	
	}
	
	/** sets pause state for particle **/
	public void setPauseState(Particle1D mP, boolean bIni)
	{
		Random rnd=new Random();

		mP.nState=0;	
		mP.dStateTimer=dPauseM+rnd.nextGaussian()*dPauseSD;
		//if initial moment, add random time shift;
		if(bIni)
			mP.dStateTimer=mP.dStateTimer*rnd.nextDouble();
		
	}
	/** sets moving state for particle **/
	public void setMovingState(Particle1D mP, boolean bIni)
	{
		int nCompN;
		Random rnd=new Random();
		float dRand;
		float dUpper;
		boolean bYes=true;
		
		mP.nState=1;
		
		dRand=rnd.nextFloat();
		dUpper=0;
		nCompN=0;
		
		//if there are components
		
		if(dComponentsP!=null)
		{
			while (bYes)
			{
				dUpper+=dComponentsP[nCompN];
				if(dRand<=dUpper)
					bYes=false;
				else
				{
					nCompN++;
				}
			}
			
			if(dSpeedDurSD[nCompN]>0.0)
				mP.dStateTimer=dSpeedDurM[nCompN]+rnd.nextGaussian()*dSpeedDurSD[nCompN];
			else 
				mP.dStateTimer=dSpeedDurM[nCompN];
			
			if(dSpeedSD[nCompN]>0.0)
				mP.v=dSpeedM[nCompN]+rnd.nextGaussian()*dSpeedSD[nCompN];
			else
				mP.v=dSpeedM[nCompN];
		}		
		
		//if initial moment, add random time shift;
		if(bIni)
			mP.dStateTimer=mP.dStateTimer*rnd.nextDouble();
		
	}
}
