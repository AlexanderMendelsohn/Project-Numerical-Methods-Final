/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christianfries.com.
 *
 * Created on 05.12.2013
 */

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
 * ------------------ Exercise 4 -----------------------------
 * -----------------------------------------------------------
 */


public class Exercise04 {
	static final DecimalFormat formatterSci6 = new DecimalFormat (" 0.0000  ; -0.0000 ");

	public static void main(String[] args) throws CalculationException {
		
		/**====================================================
		 *================INITIALIZATION=======================
		  =====================================================*/
		
		//parameters of the black-scholes model
		double riskFreeRate = 0.05;
		double volatility = 0.3;
		
		// parameters of the jump process
		double nu = 0.15;
		double intensity = 0.4	;
		double mu = -1.0/2.0 * nu * nu;
		double sigma = nu;

		
		double initialValue = 100;
		
		int numberOfPaths = 10000;
		
		//initialize time discretization
		double maturity = 2.0;
		double t0 = 0.0;
		double deltaT = 0.001;
		TimeDiscretizationInterface times = new TimeDiscretization(t0, (int)(maturity/deltaT) , deltaT);
		
		
		AssetModelMonteCarloSimulationInterface jump = new MonteCarloMertonJumpDiffusionModel(times,
				numberOfPaths, initialValue, riskFreeRate, volatility, intensity, mu, sigma);
		
		AssetModelMonteCarloSimulationInterface black = new MonteCarloBlackScholesModel(times,
				numberOfPaths, initialValue, riskFreeRate, volatility);
		
		//initialize strike vector
		double[] strikes = new double[12];
		for (int i = 0; i < 12; i++){
			strikes[i] = 30 + i*10;
		}
		
		/**====================================================
		 *===================Computation=======================
		  =====================================================*/
		
		//calculate prices of the european option at time t0 = 0
		double[] priceOfBlackScholes = new double[12];
		double[] priceOfJumpProcess = new double[12];
		for (int i = 0; i < 12; i++){
			AbstractAssetMonteCarloProduct product = new EuropeanOption(maturity, strikes[i]);
			priceOfBlackScholes[i] = product.getValue(black)
					/*net.finmath.functions.AnalyticFormulas.blackScholesGeneralizedOptionValue(initialValue * Math.exp(riskFreeRate * (maturity - t0)), volatility, maturity, strikes[i], Math.exp(-riskFreeRate * maturity))*/;
			priceOfJumpProcess[i] = product.getValue(jump);

		}

		//calculate the implied volatility of the jump process and black-scholes model by plugging price
		double[] impliedVolatilityBlackScholesByPrice = new double[12];
		double[] impliedVolatilityJumpProcessByPrice = new double[12];
		double forward = initialValue * Math.exp(riskFreeRate * (maturity - t0));
		for (int i = 0; i<12; i++){
			impliedVolatilityBlackScholesByPrice[i] = net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(
					forward,
					maturity,
					strikes[i],
					Math.exp(-riskFreeRate * (maturity - t0)),
					priceOfBlackScholes[i]);
			impliedVolatilityJumpProcessByPrice[i] = net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(
					forward,
					maturity,
					strikes[i],
					Math.exp(-riskFreeRate * (maturity - t0)),
					priceOfJumpProcess[i]);

		}
		
		/**====================================================
		 *===================Printing==========================
		  =====================================================*/
		
		//print results
		System.out.println("Printing results for: riskFreeRate = " + riskFreeRate
				+ ", volatility = " + volatility +", intensity = " + intensity
				+ ", nu = " + nu
				+ ", maturity = " + maturity);
    	System.out.println("---------------------------------------------"
    			+ "---------------------");
    	System.out.println("Merton-Jump-Diffusion-Model:");
		System.out.println("Strike |" + " " + "Option Price |" + " "  + "implied volatility");
		for (int i = 0; i < 12; i++){
			System.out.println(strikes[i]
					+ "  " + formatterSci6.format(priceOfJumpProcess[i])
					+ "  " + formatterSci6.format(impliedVolatilityJumpProcessByPrice[i]));
		}
    	System.out.println("---------------------------------------------"
    			+ "---------------------");
	}
}