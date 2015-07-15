package projectTest;
import net.finmath.stochastic.RandomVariableInterface;







import java.text.DecimalFormat;

 




import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.assetderivativevaluation.AssetModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.assetderivativevaluation.MonteCarloBlackScholesModel;
import net.finmath.montecarlo.assetderivativevaluation.products.AbstractAssetMonteCarloProduct;
import net.finmath.montecarlo.assetderivativevaluation.products.EuropeanOption;
import net.finmath.montecarlo.assetderivativevaluation.products.BlackScholesHedgedPortfolio.HedgeStrategy;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.montecarlo.assetderivativevaluation.products.BlackScholesHedgedPortfolio;
import net.finmath.montecarlo.assetderivativevaluation.products.BlackScholesHedgedPortfolio.HedgeStrategy;

public class Exercise06 {
	static final DecimalFormat formatterSci6 = new DecimalFormat (" 0.000000000  ; -0.000000000 ");
	public static void main(String[] args) throws CalculationException {
		
		/**====================================================
		 *================INITIALIZATION=======================
		  =====================================================*/
		
		/*parameters for the Black-Scholes-Model*/
		double riskFreeRate = 0.05;
        double volatility = 0.3;
        double initialValue = 100;
        
        /*parameters for the Merton-Jump-Diffusion-Model*/
        double nu = 0.15;
        double mu = -1.0/2 * nu * nu;
		double sigma = nu;
		double intensity = 0.4;
		
        /*parameters for the european option*/
        double maturity=2.0;
        double strike=100;
        
        int numberOfPaths = 3500;      
        
        /*time discretization*/
        double deltaT = 0.001;
        double t0 = 0.0;
        TimeDiscretizationInterface times = new TimeDiscretization(t0, (int)(maturity/deltaT), deltaT);
        
        /*Jump Diffusion Model*/
        AssetModelMonteCarloSimulationInterface jump = new MonteCarloMertonJumpDiffusionModel(times,
				numberOfPaths, initialValue, riskFreeRate, volatility, intensity, mu, sigma);
        /*Black Scholes Model*/
		AssetModelMonteCarloSimulationInterface black = new MonteCarloBlackScholesModel(times, numberOfPaths, initialValue, riskFreeRate, volatility);
		
		/*European Option*/
		AbstractAssetMonteCarloProduct product = new EuropeanOption(maturity,strike);
		
		/*Price of European Option on Jump-Diffusion-Model at time 0*/
  		double europeanPriceJump = product.getValue(jump);
  		/*Price of European Option on Black-Scholes-Model at time 0*/
  		double europeanPriceBlackScholes = product.getValue(black);
  		
  		/*Value of European Option in Jump-Process at maturity*/
  		RandomVariableInterface valueJumpEuropeanOptionAtMaturity = product.getValue(maturity, jump);
  		/*Value of European Option on BlackScholesModel at maturity*/
  		RandomVariableInterface valueBlackScholesEuropeanOptionAtMaturity = product.getValue(maturity, black);
  		
  		/*number of hedging times (number of timesteps is assumed to be 2000 here)*/
  		int[] numberOfHedgingTimesVector = {10,20,100,200,1000,2000};
  		int sizeOfHedgingTimesVector = numberOfHedgingTimesVector.length;
  		
  		/*vectors used to calculate the relative profit and loss*/
  		RandomVariableInterface[] relativePandLJump = new RandomVariableInterface[sizeOfHedgingTimesVector]; 
  		RandomVariableInterface[] relativePandLBlackScholes = new RandomVariableInterface[sizeOfHedgingTimesVector]; 
  		RandomVariableInterface[] differencePortfolioToOptionPriceJump = new RandomVariableInterface[sizeOfHedgingTimesVector]; 
  		RandomVariableInterface[] differencePortfolioToOptionPriceBlackScholes = new RandomVariableInterface[sizeOfHedgingTimesVector]; 
  		
		/**====================================================
		 *================Computation & Printing===============
		  =====================================================*/
  		
  		/*loop over numberOfHedgingTimes to get pathwise relative profit and loss, the printing is done in an inner loop over numberOfPaths*/
        for (int i=0;i<sizeOfHedgingTimesVector;i++){
        	
        	BlackScholesHedgedPortfolioWithModifiedTimeDiscretization hedgingPortfolioValue = 
        			new BlackScholesHedgedPortfolioWithModifiedTimeDiscretization(maturity,strike,riskFreeRate,volatility,numberOfHedgingTimesVector[i]);
        	
        	RandomVariableInterface portfolioValueJump = hedgingPortfolioValue.getValue(maturity, jump);
        	RandomVariableInterface portfolioValueBlackScholes = hedgingPortfolioValue.getValue(maturity, black);
        	
        	differencePortfolioToOptionPriceJump[i] = portfolioValueJump.sub(valueJumpEuropeanOptionAtMaturity);
        	differencePortfolioToOptionPriceBlackScholes[i] = portfolioValueBlackScholes.sub(valueBlackScholesEuropeanOptionAtMaturity);
        	
        	relativePandLJump[i] = differencePortfolioToOptionPriceJump[i].div(europeanPriceJump).mult(Math.exp(-riskFreeRate * (maturity-t0)));
        	relativePandLBlackScholes[i] = differencePortfolioToOptionPriceBlackScholes[i].div(europeanPriceBlackScholes).mult(Math.exp(-riskFreeRate * (maturity-t0)));
        	System.out.println("Printing for the hedgingtime: " + numberOfHedgingTimesVector[i]);
        	System.out.println("==================================================================");
        	/*System.out.println("Paths of the P&L for the Black-Scholes-Model: ");
        	for (int j = 0; j<numberOfPaths;j++){
        		System.out.println(relativePandLBlackScholes[i].get(j));
        	}     
           	System.out.println("------------------------------------------------------------------");
        	System.out.println("Paths of the P&L for the Merton-Jump-Diffusion-Model: ");
        	for (int j = 0; j<numberOfPaths;j++){
        		System.out.println(relativePandLJump[i].get(j));
        	}       	*/
        	System.out.println("------------------------------------------------------------------");
        	System.out.println("| Black Scholes Model | Merton-Jump-Diffusion Model:");
        	System.out.println("Mean: " + formatterSci6.format(relativePandLBlackScholes[i].getAverage()) + " " + formatterSci6.format(relativePandLJump[i].getAverage()));
        	System.out.println("Variance: " + formatterSci6.format(relativePandLBlackScholes[i].getVariance()) + " " +formatterSci6.format(relativePandLJump[i].getVariance()));
        	System.out.println("==================================================================");
        }
        
        
		/**====================================================
		 *======INITIALIZATION for the Delta-Gamma-Hedge=======
		  =====================================================*/
        //calculate implied volatility of the Merton-Jump-Diffusion-Process
        double forward = initialValue * Math.exp(riskFreeRate * (maturity - t0));
        double impliedVolatilityMertonJumpDiffusionProcess = net.finmath.functions.AnalyticFormulas.blackScholesOptionImpliedVolatility(
        		forward,
				maturity,
				strike,
				Math.exp(-riskFreeRate * (maturity - t0)),
				europeanPriceJump);
        
        //initialize the paramters for the european call option used as additional hedging instrument for the delta gamma hedge
        double additionalEuropeanCallMaturity = 5.0;
        double additionalEuropeanCallStrike = 110;
        
       	BlackScholesHedgedPortfolio deltaGamma = new BlackScholesHedgedPortfolio(maturity, strike, riskFreeRate, impliedVolatilityMertonJumpDiffusionProcess,
       			additionalEuropeanCallMaturity, additionalEuropeanCallStrike, HedgeStrategy.deltaGammaHedge);
       	
  		RandomVariableInterface portfolioValueJump = deltaGamma.getValue(maturity, jump);
  		RandomVariableInterface differencePortfolioToOptionPriceJumpDeltaGamma = portfolioValueJump.sub(valueJumpEuropeanOptionAtMaturity);
       	RandomVariableInterface relativePandLBlackJumpDeltaGamma = differencePortfolioToOptionPriceJumpDeltaGamma.div(europeanPriceJump).mult(Math.exp(-riskFreeRate * maturity));
    	System.out.println("==================================================================");
    	/*System.out.println("Paths of the P&L for the Merton-Jump-Diffusion-Model via Delta-Gamma-Hedging: ");
    	for (int j=0; j<numberOfPaths;j++){
    		System.out.println(relativePandLBlackJumpDeltaGamma.get(j));
    	}*/
    	System.out.println("------------------------------------------------------------------");
    	System.out.println("| Merton-Jump-Diffusion-Model via Delta-Gamma-Hedging ");
  		System.out.println("Mean" + formatterSci6.format(relativePandLBlackJumpDeltaGamma.getAverage()));
  		System.out.println("Variance" + formatterSci6.format(relativePandLBlackJumpDeltaGamma.getVariance()));
        
    }
	

}
