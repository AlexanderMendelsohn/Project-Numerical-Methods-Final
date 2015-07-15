package projectTest;


import net.finmath.time.TimeDiscretizationInterface;
import net.finmath.stochastic.RandomVariableInterface;
import cern.jet.random.engine.MersenneTwister64;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.AbstractRandomVariableFactory;

public class PoissonProcess implements PointProcessInterface{
	
	private AbstractRandomVariableFactory randomVariableFactory = new RandomVariableFactory();
	
	private double intensity;
	private TimeDiscretizationInterface timeDiscretization;
	private RandomVariableInterface[] poissonProcess;
	private RandomVariableInterface[] poissonProcessIncrements;
	private int numberOfPaths;
	private int seed;
	
	private final		Object						poissonProcessIncrementsLazyInitLock = new Object();

	
	public PoissonProcess(double intensity,
			TimeDiscretizationInterface timeDiscretization,
			int numberOfPaths,
			int seed) {
		super();
		this.intensity = intensity;
		this.timeDiscretization = timeDiscretization;
		this.poissonProcess = null;
		this.poissonProcessIncrements = null;
		this.numberOfPaths = numberOfPaths;
		this.seed = seed;
	}

	public RandomVariableInterface getProcessIncrements(int timeIndex){
		synchronized(poissonProcessIncrementsLazyInitLock) {
			if (poissonProcessIncrements == null) doGeneratePoissonProcess();
		}
		return poissonProcessIncrements[timeIndex];
	}
	
	public RandomVariableInterface getProcess(int timeIndex){
		synchronized(poissonProcessIncrementsLazyInitLock) {
			if (poissonProcess == null) doGeneratePoissonProcess();
		}
		return poissonProcess[timeIndex];
	}
	
	private void doGeneratePoissonProcess(){
		if (poissonProcess != null) return;
		
		MersenneTwister64		mersenneTwister		= new MersenneTwister64(seed);
		
		double PoissonProcessArray[][] = new double[timeDiscretization.getNumberOfTimes()][numberOfPaths];
		double PoissonProcessIncrementsArray[][] = new double[timeDiscretization.getNumberOfTimeSteps()][numberOfPaths];


		for (int path = 0; path<numberOfPaths; path++){
			double sum = inverseExponentialFunction(mersenneTwister.nextDouble());
			PoissonProcessArray[0][path] = 0.0;
			for (int timeIndex = 1; timeIndex < timeDiscretization.getNumberOfTimes();timeIndex++){
				PoissonProcessArray[timeIndex][path] = PoissonProcessArray[timeIndex-1][path];
				while (timeDiscretization.getTime(timeIndex) >= sum) {
					//System.out.println("sum "+ sum + "    time " + timeDiscretization.getTime(timeIndex) + "  " + timeIndex + "  " + path);
					double uniformIncrement = mersenneTwister.nextDouble();
					double exponentialIncrement = inverseExponentialFunction(uniformIncrement);
					sum += exponentialIncrement;
					PoissonProcessArray[timeIndex][path] +=1;
					//System.out.println(PoissonProcessArray[timeIndex][path]);
					
				} 
			}
		}
		
		for (int path = 0; path<numberOfPaths;path++){
			for (int timeIncrements = 0; timeIncrements < timeDiscretization.getNumberOfTimeSteps();timeIncrements++){
				PoissonProcessIncrementsArray[timeIncrements][path] = 
						PoissonProcessArray[timeIncrements+1][path]
								- PoissonProcessArray[timeIncrements][path];
			}	
		}
		
		poissonProcess = new RandomVariableInterface[timeDiscretization.getNumberOfTimes()];
		poissonProcessIncrements = new RandomVariableInterface[timeDiscretization.getNumberOfTimeSteps()];
		
		for(int timeIndex=0; timeIndex<timeDiscretization.getNumberOfTimeSteps(); timeIndex++) {
			double time = timeDiscretization.getTime(timeIndex+1);	
			poissonProcessIncrements[timeIndex] =
					randomVariableFactory.createRandomVariable(time, PoissonProcessIncrementsArray[timeIndex]);		
		}
		
		for(int timeIndex=0; timeIndex<timeDiscretization.getNumberOfTimes();timeIndex++){
			double time = timeDiscretization.getTime(timeIndex);
			poissonProcess[timeIndex]=
					randomVariableFactory.createRandomVariable(time, PoissonProcessArray[timeIndex]);
					
		}
	}
	
	private double inverseExponentialFunction(double X){
		return -1.0 / intensity * Math.log(1.0 - X);
	}
	
	

}
