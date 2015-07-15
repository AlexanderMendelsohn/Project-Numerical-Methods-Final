package projectTest;

import net.finmath.stochastic.RandomVariableInterface;

public interface PointProcessInterface {
	
	public RandomVariableInterface getProcess(int timeIndex);
	
	public RandomVariableInterface getProcessIncrements(int timeIndex);
	
	

}
