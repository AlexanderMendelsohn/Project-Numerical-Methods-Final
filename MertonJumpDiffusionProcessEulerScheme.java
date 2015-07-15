package projectTest;



import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.process.AbstractProcess;

import net.finmath.stochastic.RandomVariableInterface;

	/**
	 * This class implements some numerical schemes for multi-dimensional multi-factor Ito process.
	 * 
	 * It features the standard Euler scheme and the standard predictor-corrector Euler scheme
	 * for <i>Y</i>, then applies the <i>state space transform</i> \( X = f(Y) \). For the standard Euler scheme
	 * the process Y is discretized as
	 * \[
	 * 	Y(t_{i+1}) = Y(t_{i}) + \mu(t_{i}) \Delta t_{i} + \sigma(t_{i}) \Delta W(t_{i}) \text{.} 
	 * \]

	 * 
	 * Hence, using the <i>state space transform</i>, it is possible to create a log-Eurler scheme, i.e.,
	 * \[
	 * 	X(t_{i+1}) = X(t_{i}) \cdot \exp\left( (\mu(t_{i}) - \frac{1}{2} sigma(t_{i})^2) \Delta t_{i} + \sigma(t_{i}) \Delta W(t_{i}) \right) \text{.} 
	 * \]
	 * 
	 * The dimension is called <code>numberOfComponents</code> here. The default for <code>numberOfFactors</code> is 1.
	 * 
	 * @author Christian Fries
	 * @see AbstractProcessInterface The interface definition contains more details.
	 * @version 1.4
	 */
	public class MertonJumpDiffusionProcessEulerScheme extends AbstractProcess {


		private BrownianMotionInterface brownianMotion;
		private PointProcessInterface compoundPoissonProcess;


		/*
		 * The storage of the simulated stochastic process.
		 */
		private transient RandomVariableInterface[][]	discreteProcess = null;
		private transient RandomVariableInterface[]		discreteProcessWeights;

		/**
		 * @param brownianMotion The Brownian driver of the process
		 */

		public MertonJumpDiffusionProcessEulerScheme(BrownianMotionInterface brownianMotion, PointProcessInterface compoundPoissonProcess) {
			// compare timediscr of brownian and compound
			super(brownianMotion.getTimeDiscretization());
			this.brownianMotion = brownianMotion;
			this.compoundPoissonProcess = compoundPoissonProcess;
		}
		/**
		 * This method returns the realization of the process at a certain time index.
		 * 
		 * @param timeIndex Time index at which the process should be observed
		 * @return A vector of process realizations (on path)
		 */
		@Override
		public RandomVariableInterface getProcessValue(int timeIndex, int componentIndex) {
			// Thread safe lazy initialization
			synchronized(this) {
				if (discreteProcess == null || discreteProcess.length == 0) {
					doPrecalculateProcess();
				}
			}

			if(discreteProcess[timeIndex][componentIndex] == null) {
				throw new NullPointerException("Generation of process component " + componentIndex + " at time index " + timeIndex + " failed. Likely due to out of memory");
			}
			
			// Return value of process
			return discreteProcess[timeIndex][componentIndex];
		}

		/**
		 * This method returns the weights of a weighted Monte Carlo method (the probability density).
		 * 
		 * @param timeIndex Time index at which the process should be observed
		 * @return A vector of positive weights
		 */
		@Override
		public RandomVariableInterface getMonteCarloWeights(int timeIndex) {
			// Thread safe lazy initialization
			synchronized(this) {
				if (discreteProcessWeights == null || discreteProcessWeights.length == 0) {
					doPrecalculateProcess();
				}
			}

			// Return value of process
			return discreteProcessWeights[timeIndex];
		}

		/**
		 * Calculates the whole (discrete) process.
		 */
		private void doPrecalculateProcess() {
			if (discreteProcess != null && discreteProcess.length != 0)	return;

			final int numberOfPaths			= this.getNumberOfPaths();
			final int numberOfFactors		= this.getNumberOfFactors();
			final int numberOfComponents	= this.getNumberOfComponents();

			// Allocate Memory
			discreteProcess			= new RandomVariableInterface[getTimeDiscretization().getNumberOfTimeSteps() + 1][getNumberOfComponents()];
			discreteProcessWeights	= new RandomVariableInterface[getTimeDiscretization().getNumberOfTimeSteps() + 1];

			// Set initial Monte-Carlo weights
			discreteProcessWeights[0] = brownianMotion.getRandomVariableForConstant(1.0 / numberOfPaths);

			// Set initial value
			RandomVariableInterface[] initialState = getInitialState();
			final RandomVariableInterface[] currentState = new RandomVariableInterface[numberOfComponents];
			for (int componentIndex = 0; componentIndex < numberOfComponents; componentIndex++) {
				currentState[componentIndex] = initialState[componentIndex];
				discreteProcess[0][componentIndex] = applyStateSpaceTransform(componentIndex, currentState[componentIndex]);
			}

			/*
			 * Evolve the process using an Euler scheme.
			 */


			// Evolve process
			for (int timeIndex2 = 1; timeIndex2 < getTimeDiscretization().getNumberOfTimeSteps()+1; timeIndex2++) {
				final int timeIndex = timeIndex2;
				// Generate process from timeIndex-1 to timeIndex
				final double deltaT = getTime(timeIndex) - getTime(timeIndex - 1);

				// Fetch drift vector
				RandomVariableInterface[] drift = getDrift(timeIndex - 1, discreteProcess[timeIndex - 1], null);
				RandomVariableInterface[] jump = new RandomVariableInterface[numberOfComponents];
				for (int comp = 0; comp < numberOfComponents;comp++){
					jump[comp] = discreteProcess[timeIndex-1][comp].mult(compoundPoissonProcess.getProcessIncrements(timeIndex-1));
					
				}

				// Calculate new realization
				for (int componentIndex2 = 0; componentIndex2 < numberOfComponents; componentIndex2++) {
					//System.out.println(componentIndex2);
					final int componentIndex = componentIndex2;

					final RandomVariableInterface	driftOfComponent	= drift[componentIndex];
					final RandomVariableInterface 	jumpOfComponent 	= jump[componentIndex];

					// Check if the component process has stopped to evolve
					if (driftOfComponent == null) continue;

					RandomVariableInterface[]	factorLoadings		= getFactorLoading(timeIndex - 1, componentIndex, discreteProcess[timeIndex - 1]);

							// Check if the component process has stopped to evolve
					
							// Temp storage for variance and diffusion
					RandomVariableInterface diffusionOfComponent	= brownianMotion.getRandomVariableForConstant(0.0);

							// Generate values for diffusionOfComponent and varianceOfComponent 
					for (int factor = 0; factor < numberOfFactors; factor++) {
						RandomVariableInterface factorLoading		= factorLoadings[factor];
						RandomVariableInterface brownianIncrement	= brownianMotion.getBrownianIncrement(timeIndex - 1, factor);

						diffusionOfComponent = diffusionOfComponent.addProduct(factorLoading, brownianIncrement);
					}

					RandomVariableInterface increment = diffusionOfComponent;
					if(driftOfComponent != null) increment = increment.addProduct(driftOfComponent, deltaT).add(jumpOfComponent);
					
					//if (componentIndex2 == 0) System.out.println(increment.get(0));

							// Add increment to state and applyStateSpaceTransform
					currentState[componentIndex] = currentState[componentIndex].add(increment);
							
							// Transform the state space to the value space and return it.
					//applyStateSpaceTransform(componentIndex, currentState[componentIndex]);
					discreteProcess[timeIndex][componentIndex] = applyStateSpaceTransform(componentIndex, currentState[componentIndex]);
					//System.out.println(currentState[componentIndex].get(0));
				}



				
				// Set Monte-Carlo weights
				discreteProcessWeights[timeIndex] = discreteProcessWeights[timeIndex - 1];



		// Set Monte-Carlo weights
		discreteProcessWeights[timeIndex] = discreteProcessWeights[timeIndex - 1];
	 // End for(timeIndex)
	}


}


		/**
		 * Reset all precalculated values
		 */
		private synchronized void reset() {
			this.discreteProcess = null;
			this.discreteProcessWeights = null;
		}

		/**
		 * @return Returns the numberOfPaths.
		 */
		@Override
		public int getNumberOfPaths() {
			return this.brownianMotion.getNumberOfPaths();
		}

		/**
		 * @return Returns the numberOfFactors.
		 */
		@Override
		public int getNumberOfFactors() {
			return this.brownianMotion.getNumberOfFactors();
		}

		/**
		 * @param seed The seed to set.
		 * @deprecated The class will soon be changed to be immutable
		 */
		@Deprecated
		public void setSeed(int seed) {
			// Create a new Brownian motion
			this.setBrownianMotion(new net.finmath.montecarlo.BrownianMotion(
					brownianMotion.getTimeDiscretization(), brownianMotion
					.getNumberOfFactors(), brownianMotion
					.getNumberOfPaths(), seed));
			// Force recalculation of the process
			this.reset();
		}

		/**
		 * @return Returns the Brownian motion used in the generation of the process
		 */
		@Override
		public BrownianMotionInterface getBrownianMotion() {
			return brownianMotion;
		}

		/**
		 * @param brownianMotion The brownianMotion to set.
		 * @deprecated Do not use anymore. Processes should be immutable.
		 */
		@Deprecated
		public void setBrownianMotion(
				net.finmath.montecarlo.BrownianMotion brownianMotion) {
			this.brownianMotion = brownianMotion;
			// Force recalculation of the process
			this.reset();
		}

		public PointProcessInterface getCompoundPoissonProcess(){
			return compoundPoissonProcess;
		}
	
		/* (non-Javadoc)
		 * @see net.finmath.montecarlo.process.AbstractProcess#clone()
		 */
		@Override
		public MertonJumpDiffusionProcessEulerScheme clone() {
			return new MertonJumpDiffusionProcessEulerScheme(getBrownianMotion(), getCompoundPoissonProcess());
		}

		/* (non-Javadoc)
		 * @see net.finmath.montecarlo.process.AbstractProcess#getCloneWithModifiedSeed(int)
		 */
		@Override
		public Object getCloneWithModifiedSeed(int seed) {
			return new MertonJumpDiffusionProcessEulerScheme(getBrownianMotion(), getCompoundPoissonProcess());
		}
		

	}