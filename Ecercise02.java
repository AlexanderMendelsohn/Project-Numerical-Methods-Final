package projectTest;

import java.text.DecimalFormat;


import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.assetderivativevaluation.AssetModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.assetderivativevaluation.MonteCarloBlackScholesModel;
import net.finmath.montecarlo.assetderivativevaluation.products.AbstractAssetMonteCarloProduct;
import net.finmath.montecarlo.assetderivativevaluation.products.EuropeanOption;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;


/**-----------------------------------------------------------
 * ------------------ Exercise 02 ----------------------------
 * -----------------------------------------------------------
 */

public class Ecercise02 {
	
	static final DecimalFormat formatterSci6 = new DecimalFormat (" 0.00000  ; -0.00000 ");
	
	public static void main(String[] args) throws CalculationException {
		
		/**====================================================
		 *================INITIALIZATION=======================
		  =====================================================*/
		
		//parameters of the black-scholes model
		double riskFreeRate = 0.05;
		double volatility = 0.3;
		
		/*parameters for the Merton-Jump-Diffusion-Model*/
		double[] intensityVector = {0.0, 0.4, 1.0, 1.4, 2.0};
		double nu = 0.15;
		double mu = -1.0/2.0 * nu * nu;
		double sigma = nu;
			
		double initialValue = 100;
		
		int numberOfPaths = 10000;
		
		/*initialize time discretization*/
		double maturity = 2.0;
		double t0 = 0.0;
		double deltaT = 0.001;
		TimeDiscretizationInterface times = new TimeDiscretization(t0, (int)(maturity/deltaT) , deltaT);
		
		/*strike for atm european option*/
		double strike = 100.0;
		
		/*european option*/
		AbstractAssetMonteCarloProduct product = new EuropeanOption(maturity, strike);
		
		/**====================================================
		 *================Computation & Printing===============
		  =====================================================*/
		
		System.out.println("Printing results for: riskFreeRate = " + riskFreeRate
				+ ", volatility = " + volatility
				+ ", nu = " + nu
				+ ", strike = " + strike
				+ ", maturity = " + maturity);
    	System.out.println("---------------------------------------------"
    			+ "---------------------");
    	System.out.println("Merton-Jump-Diffusion-Model:");
		System.out.println("intensity | Option price");
		for (int i=0; i<intensityVector.length;i++){
			AssetModelMonteCarloSimulationInterface jump = new MonteCarloMertonJumpDiffusionModel(times,
					numberOfPaths, initialValue, riskFreeRate, volatility, intensityVector[i], mu, sigma);
			System.out.println(intensityVector[i] + "        "   + formatterSci6.format(product.getValue(maturity,jump).getAverage()));

		}
    	System.out.println("---------------------------------------------"
    			+ "---------------------");
		AssetModelMonteCarloSimulationInterface black = new MonteCarloBlackScholesModel(times,
				numberOfPaths, initialValue, riskFreeRate, volatility);
		System.out.println("Option Price for the Black-Scholes-Model: " + formatterSci6.format(product.getValue(maturity,black).getAverage()));
		
	}

}
