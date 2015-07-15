/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 20.01.2004
 */
package projectTest;

import java.util.ArrayList;
import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.model.AbstractModel;
import net.finmath.montecarlo.process.AbstractProcess;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;
import net.finmath.montecarlo.assetderivativevaluation.AssetModelMonteCarloSimulationInterface;
import projectTest.JumpProcessEulerScheme;


/**
 * This class glues together a <code>BlackScholeModel</code> with random jumps via <code>CompoundPoissonProcess</code> and a Monte-Carlo implementation of a <code>AbstractProcess</code>
 * and forms a Monte-Carlo implementation of the Black-Scholes-Jump-Diffusion Model by implementing <code>AssetModelMonteCarloSimulationInterface</code>.
 *
 * The model is
 * \[
 * 	dS = r S dt + \sigma S dW + SdJ, \quad S(0) = S_{0},
 * \]
 * \[
 * 	dN = r N dt, \quad N(0) = N_{0},
 * \]
 * 
 * where J defines a Compound Poisson Process with jumps 
 * \[
 *  \exp{\mu + \sigma  Z) - 1,
 * \]
 * where Z is normally distributed.
 * 
 * The stateSpaceTransform f is set to id.
 * 
 * 
 */
public class MonteCarloMertonJumpDiffusionModel extends AbstractModel implements AssetModelMonteCarloSimulationInterface {

	private final double initialValue;
	private final double riskFreeRate;		// Actually the same as the drift (which is not stochastic)
	private final double volatility;
	private final double poissonIntensity;
	private final double jumpMean;
	private final double jumpVariance;
	
	private final int seed = 3141;

	private final RandomVariableInterface[]	initialValueVector	= new RandomVariableInterface[1];
	private final RandomVariableInterface	drift;
	private final RandomVariableInterface	volatilityOnPaths;

	/**
	 * Create a Monte-Carlo simulation using given time discretization.
	 * 
	 * @param timeDiscretization The time discretization.
	 * @param numberOfPaths The number of Monte-Carlo path to be used.
	 * @param initialValue Spot value.
	 * @param riskFreeRate The risk free rate.
	 * @param volatility The log volatility.
	 * @param poissonIntensity The intensity of the Poisson Process
	 * @param jumpMean The mean and variance of the jumps
	 * @param jumpVariance exp (mean + variance * Z) - 1,  Z normally distributed
	 */
	public MonteCarloMertonJumpDiffusionModel(
			TimeDiscretizationInterface timeDiscretization,
			int numberOfPaths,
			double initialValue,
			double riskFreeRate,
			double volatility,
			double poissonIntensity,
			double jumpMean,
			double jumpVariance) {
		super();

		this.initialValue	= initialValue;
		this.riskFreeRate	= riskFreeRate;
		this.volatility		= volatility;
		this.poissonIntensity = poissonIntensity;
		this.jumpMean = jumpMean;
		this.jumpVariance = jumpVariance;
		// Create a corresponding MC process
		AbstractProcess process = new JumpProcessEulerScheme(new BrownianMotion(timeDiscretization, 1 /* numberOfFactors */, numberOfPaths, seed),
															new CompoundPoissonProcess(poissonIntensity, jumpMean, jumpVariance, timeDiscretization, numberOfPaths, seed + 300));
		

		
		this.initialValueVector[0]	= process.getBrownianMotion().getRandomVariableForConstant(initialValue);
		this.drift					= process.getBrownianMotion().getRandomVariableForConstant(riskFreeRate);
		this.volatilityOnPaths		= process.getBrownianMotion().getRandomVariableForConstant(volatility);

		// Link model and process for delegation
		process.setModel(this);
		this.setProcess(process);
	}

	/**
	 * Create a Monte-Carlo simulation using given process discretization scheme.
	 * 
	 * @param initialValue Spot value
	 * @param riskFreeRate The risk free rate
	 * @param volatility The log volatility
	 * @param process The process discretization scheme which should be used for the simulation.
	 */

	 public MonteCarloMertonJumpDiffusionModel(
	 
			double initialValue,
			double riskFreeRate,
			double volatility,
			AbstractProcess process,
			double poissonIntensity,
			double jumpMean,
			double jumpVariance) {
		super();

		this.initialValue	= initialValue;
		this.riskFreeRate	= riskFreeRate;
		this.volatility		= volatility;
		this.poissonIntensity = poissonIntensity;
		this.jumpMean = jumpMean;
		this.jumpVariance = jumpVariance;

		/*
		 * The interface definition requires that we provide the drift and the volatility in terms of random variables.
		 * We construct the corresponding random variables here and will return (immutable) references to them.
		 */
		this.initialValueVector[0]	= process.getBrownianMotion().getRandomVariableForConstant(initialValue);
		this.drift					= process.getBrownianMotion().getRandomVariableForConstant(riskFreeRate);
		this.volatilityOnPaths		= process.getBrownianMotion().getRandomVariableForConstant(volatility);
		
		// Link model and process for delegation
		process.setModel(this);
		this.setProcess(process);
	}


	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.model.AbstractModelInterface#getInitialState()
	 */
	@Override
	public RandomVariableInterface[] getInitialState() {
		return initialValueVector;
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.model.AbstractModelInterface#getDrift(int, net.finmath.stochastic.RandomVariableInterface[], net.finmath.stochastic.RandomVariableInterface[])
	 */
	@Override
	public RandomVariableInterface[] getDrift(int timeIndex, RandomVariableInterface[] realizationAtTimeIndex, RandomVariableInterface[] realizationPredictor) {
		return new RandomVariableInterface[] { drift.mult(realizationAtTimeIndex[0])};
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.model.AbstractModelInterface#getFactorLoading(int, int, net.finmath.stochastic.RandomVariableInterface[])
	 */
	@Override
	public RandomVariableInterface[] getFactorLoading(int timeIndex, int component, RandomVariableInterface[] realizationAtTimeIndex) {
		return new RandomVariableInterface[] { volatilityOnPaths.mult(realizationAtTimeIndex[0]) };
	}

	public double getPoissonIntensity() {
		return poissonIntensity;
	}

	public double getJumpMean() {
		return jumpMean;
	}

	public double getJumpVariance() {
		return jumpVariance;
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.model.AbstractModelInterface#applyStateSpaceTransform(int, net.finmath.stochastic.RandomVariableInterface)
	 * StateSpaceTransformation modified to identity. We are not consiedering the log-Euler scheme in the Diffusion Process
	 */
	@Override
	public RandomVariableInterface applyStateSpaceTransform(int componentIndex, RandomVariableInterface randomVariable) {
		return randomVariable;
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.assetderivativevaluation.AssetModelMonteCarloSimulationInterface#getAssetValue(double, int)
	 */
	@Override
	public RandomVariableInterface getAssetValue(double time, int assetIndex) throws CalculationException {
		return getAssetValue(getTimeIndex(time), assetIndex);
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.assetderivativevaluation.AssetModelMonteCarloSimulationInterface#getAssetValue(int, int)
	 */
	@Override
	public RandomVariableInterface getAssetValue(int timeIndex, int assetIndex) throws CalculationException {
		return getProcessValue(timeIndex, assetIndex);
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.MonteCarloSimulationInterface#getMonteCarloWeights(double)
	 */
	@Override
	public RandomVariableInterface getMonteCarloWeights(double time) throws CalculationException {
		return getMonteCarloWeights(getTimeIndex(time));
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.assetderivativevaluation.AssetModelMonteCarloSimulationInterface#getNumeraire(int)
	 */
	@Override
	public RandomVariableInterface getNumeraire(int timeIndex) {
		double time = getTime(timeIndex);

		return getNumeraire(time);
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.model.AbstractModelInterface#getNumeraire(double)
	 */
	@Override
	public RandomVariableInterface getNumeraire(double time) {
		double numeraireValue = Math.exp(riskFreeRate * time);

		return getRandomVariableForConstant(numeraireValue);
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.assetderivativevaluation.AssetModelMonteCarloSimulationInterface#getRandomVariableForConstant(double)
	 */
	@Override
	public RandomVariableInterface getRandomVariableForConstant(double value) {
		return getProcess().getBrownianMotion().getRandomVariableForConstant(value);
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.model.AbstractModelInterface#getNumberOfComponents()
	 */
	@Override
	public int getNumberOfComponents() {
		return 1;
	}

	/* (non-Javadoc)
	 * @see net.finmath.montecarlo.assetderivativevaluation.AssetModelMonteCarloSimulationInterface#getNumberOfAssets()
	 */
	@Override
	public int getNumberOfAssets() {
		return 1;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return super.toString() + "\n" +
				"MonteCarloBlackScholesModel:\n" +
				"  initial value...:" + initialValue + "\n" +
				"  risk free rate..:" + riskFreeRate + "\n" +
				"  volatiliy.......:" + volatility;
	}

	/**
	 * Returns the risk free rate parameter of this model.
	 *
	 * @return Returns the riskFreeRate.
	 */
	public double getRiskFreeRate() {
		return riskFreeRate;
	}

	/**
	 * Returns the volatility parameter of this model.
	 * 
	 * @return Returns the volatility.
	 */
	public double getVolatility() {
		return volatility;
	}

	@Override
	public AssetModelMonteCarloSimulationInterface getCloneWithModifiedData(Map<String, Object> dataModified) {
		/*
		 * Determine the new model parameters from the provided parameter map.
		 */
		double	newInitialTime	= dataModified.get("initialTime") != null	? ((Number)dataModified.get("initialTime")).doubleValue() : getTime(0);
		double	newInitialValue	= dataModified.get("initialValue") != null	? ((Number)dataModified.get("initialValue")).doubleValue() : initialValue;
		double	newRiskFreeRate	= dataModified.get("riskFreeRate") != null	? ((Number)dataModified.get("riskFreeRate")).doubleValue() : riskFreeRate;
		double	newVolatility	= dataModified.get("volatility") != null	? ((Number)dataModified.get("volatility")).doubleValue()	: volatility;
		int		newSeed			= dataModified.get("seed") != null			? ((Number)dataModified.get("seed")).intValue()				: seed;
		
		/*selfmade*/
		double newPoissonIntensity 	= dataModified.get("poissonIntensity") != null ? ((Number)dataModified.get("poissonIntensity")).doubleValue() : poissonIntensity;
		double newJumpMean 			= dataModified.get("jumpMean") != null ? ((Number)dataModified.get("jumpMean")).doubleValue() : jumpMean;
		double newJumpVariance		= dataModified.get("jumpVariance") != null ? ((Number)dataModified.get("jumpVariance")).doubleValue() : jumpVariance;

		/*
		 * Create a new model with the new model parameters
		 */
		BrownianMotionInterface brownianMotion;
		if(dataModified.get("seed") != null) {
			// The seed has changed. Hence we have to create a new BrownianMotion.
			brownianMotion = new BrownianMotion(this.getTimeDiscretization(), 1, this.getNumberOfPaths(), newSeed);
		}
		else
		{
			// The seed has not changed. We may reuse the random numbers (Brownian motion) of the original model
			brownianMotion = this.getProcess().getBrownianMotion();
		}

		double timeShift = newInitialTime - getTime(0);
		if(timeShift != 0) {
			ArrayList<Double> newTimes = new ArrayList<Double>();
			newTimes.add(newInitialTime);
			for(Double time : getProcess().getBrownianMotion().getTimeDiscretization()) if(time > newInitialTime) newTimes.add(time);
			TimeDiscretizationInterface newTimeDiscretization = new TimeDiscretization(newTimes);
			brownianMotion = brownianMotion.getCloneWithModifiedTimeDiscretization(newTimeDiscretization);
		}
		
		/* selfmade : new parameters in constructor*/
		AbstractProcess process = new ProcessEulerScheme(brownianMotion);
		return new MonteCarloMertonJumpDiffusionModel(newInitialValue, newRiskFreeRate, newVolatility, process, 
				newPoissonIntensity, newJumpMean, newJumpVariance);    		
	}

	@Override
	public AssetModelMonteCarloSimulationInterface getCloneWithModifiedSeed(int seed) {
		// Create a corresponding MC process
		AbstractProcess process = new ProcessEulerScheme(new BrownianMotion(this.getTimeDiscretization(), 1 /* numberOfFactors */, this.getNumberOfPaths(), seed));
		return new MonteCarloMertonJumpDiffusionModel(initialValue, riskFreeRate, volatility, process, poissonIntensity, jumpMean, jumpVariance);
	}

	/**
	 * @return The number of paths.
	 * @see net.finmath.montecarlo.process.AbstractProcess#getNumberOfPaths()
	 */
	@Override
	public int getNumberOfPaths() {
		return getProcess().getNumberOfPaths();
	}
}
