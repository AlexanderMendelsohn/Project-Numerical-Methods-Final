package projectTest;

import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;
import cern.jet.random.engine.MersenneTwister64;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.AbstractRandomVariableFactory;

/**
 * This class implements a Compound Poisson Process with intensity \lambda, where the jump-height has distribution 
 * \[
 * \exp{\mu + \sigma  Z) - 1
 * \]
 * and Z is normally distributed.
 * 
 * @author A V L
 * @see PointProcessInterface
 * @version 1.0
 */

public class CompoundPoissonProcess implements PointProcessInterface {
	
	private AbstractRandomVariableFactory randomVariableFactory = new RandomVariableFactory();
	
	private RandomVariableInterface[] compoundPoissonProcess;
	private RandomVariableInterface[] compoundPoissonProcessIncrements;
	
	private double intensity;
	private double mean;
	private double variance;
	private TimeDiscretizationInterface timeDiscretization;	
	private int numberOfPaths;
	private int seed;

	
	private final		Object						compoundPoissonProcessIncrementsLazyInitLock = new Object();
	
	/**
	 * @param intensity The intensity of the Poisson Process
	 * @param mean The mean and variance of the jumps
	 * @param variance exp (mean + variance * Z) - 1,  Z normally distributed
	 */
	
	public CompoundPoissonProcess(double intensity, double mean, double variance,
			TimeDiscretizationInterface timeDiscretization, int numberOfPaths,
			int seed) {
		super();
		this.intensity = intensity;
		this.mean = mean;
		this.variance = variance;
		this.timeDiscretization = timeDiscretization;
		this.numberOfPaths = numberOfPaths;
		this.seed = seed;
		this.compoundPoissonProcess = null;
		this.compoundPoissonProcessIncrements = null;
	}

	@Override
	public RandomVariableInterface getProcess(int timeIndex) {
		synchronized(compoundPoissonProcessIncrementsLazyInitLock) {
			if (compoundPoissonProcess == null) doGenerateCompoundPoissonProcess();
		}
		return compoundPoissonProcess[timeIndex];
	}

	@Override
	public RandomVariableInterface getProcessIncrements(int timeIndex) {
		synchronized(compoundPoissonProcessIncrementsLazyInitLock) {
			if (compoundPoissonProcess == null) doGenerateCompoundPoissonProcess();
		}
		return compoundPoissonProcessIncrements[timeIndex];
	}
	
	private void doGenerateCompoundPoissonProcess(){
		
		if (compoundPoissonProcess != null) return;
		
		MersenneTwister64		mersenneTwister		= new MersenneTwister64(seed);
		PointProcessInterface poissonProcess = new PoissonProcess(intensity, timeDiscretization, numberOfPaths, seed);
		
		double[][] compoundPoissonProcessArray = new double[timeDiscretization.getNumberOfTimes()][numberOfPaths];
		double[][] compoundPoissonProcessIncrementsArray = new double[timeDiscretization.getNumberOfTimes()][numberOfPaths];
		
		for (int path = 0; path < numberOfPaths; path++){
			compoundPoissonProcessArray[0][path] = 0.0;
			for (int timeIndex = 1; timeIndex < timeDiscretization.getNumberOfTimes(); timeIndex ++){
				compoundPoissonProcessArray[timeIndex][path] = compoundPoissonProcessArray[timeIndex-1][path];
				int difference = (int)(poissonProcess.getProcess(timeIndex).get(path) - poissonProcess.getProcess(timeIndex-1).get(path));
				if (difference > 0){
					for (int i = 1; i<=difference; i++){
						double uniformNumber = mersenneTwister.nextDouble();
						compoundPoissonProcessArray[timeIndex][path] += Math.exp(mean + variance * net.finmath.functions.NormalDistribution.inverseCumulativeDistribution(uniformNumber))-1;
					}
				}
			}
		}
		
		for (int path = 0; path<numberOfPaths;path++){
			for (int timeIncrements = 0; timeIncrements < timeDiscretization.getNumberOfTimeSteps();timeIncrements++){
				compoundPoissonProcessIncrementsArray[timeIncrements][path] = 
						compoundPoissonProcessArray[timeIncrements+1][path]
								- compoundPoissonProcessArray[timeIncrements][path];
			}	
		}
		
		compoundPoissonProcess = new RandomVariableInterface[timeDiscretization.getNumberOfTimes()];
		compoundPoissonProcessIncrements = new RandomVariableInterface[timeDiscretization.getNumberOfTimeSteps()];
		
		for(int timeIndex=0; timeIndex<timeDiscretization.getNumberOfTimeSteps(); timeIndex++) {
			double time = timeDiscretization.getTime(timeIndex+1);	
			compoundPoissonProcessIncrements[timeIndex] =
					randomVariableFactory.createRandomVariable(time, compoundPoissonProcessIncrementsArray[timeIndex]);		
		}
		
		for(int timeIndex=0; timeIndex<timeDiscretization.getNumberOfTimes();timeIndex++){
			double time = timeDiscretization.getTime(timeIndex);
			compoundPoissonProcess[timeIndex]=
					randomVariableFactory.createRandomVariable(time, compoundPoissonProcessArray[timeIndex]);
					
		}
		
		
		
	}
	
	

	
	
}
