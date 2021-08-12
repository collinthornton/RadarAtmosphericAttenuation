package propagationfactor;

public class PropagationFactorTest {
	public static void main(String [] args) {
		//passed

	
		PropagationFactorCalc initialTestOne = new PropagationFactorCalc(3, 0.001, 0 ,0, 8500000, 0.03, 10, 60, 40000);
		initialTestOne.computePropagation();
		
		//passed


		PropagationFactorCalc blockPathTestOne = new PropagationFactorCalc(3, 0.001, 0 ,0, 8500000, 0.03, 10, 60, 50000);
		blockPathTestOne.computePropagation();
		
		
//		//passed

		
		PropagationFactorCalc tryShortRangesTestOne = new PropagationFactorCalc(3, 0.001, 0 ,0, 8500000, 0.03, 10, 60, 60000);
		tryShortRangesTestOne.computePropagation();
		
//		//passed

	
		PropagationFactorCalc clearPathTestOne = new PropagationFactorCalc(3, 0.001, 0 ,0, 8500000, 0.03, 10, 60, 30000);
		clearPathTestOne.computePropagation();
//		//passed

		PropagationFactorCalc nearFirstNullTestOne = new PropagationFactorCalc(3, 0.001, 0 ,0, 8500000, 0.03, 10, 60, 25000);
		nearFirstNullTestOne.computePropagation();
		
//		//passed

		PropagationFactorCalc otherSideOfNullTestOne = new PropagationFactorCalc(3, 0.001, 0 ,0, 8500000, 0.03, 10, 60, 20000);
		otherSideOfNullTestOne.computePropagation();
		
	}
}
