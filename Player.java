import java.util.ArrayList;
import java.util.List;


class Player {
	private static final boolean SHOOTING = true;
	private static final boolean GUESSING = true;
	
	private static final int RUN = 45;
	private static final int FLY_PATH = 6;
	private static final double MIN_GUESS = 1000;
	
	private int tryGuess=0, wrongGuess=0, rightGuess=0;
	private int currentRound = 0, turn=0;
	private int lastGuess[] = new int[20];
	private double weightSpecies[] = new double[Constants.COUNT_SPECIES];
	private double[][][] species = new double[Constants.COUNT_SPECIES][Constants.COUNT_MOVE][Constants.COUNT_MOVE];
	private boolean[] knowSpecies = new boolean[Constants.COUNT_SPECIES];
	private boolean startBlackStork = false;
	
	public Player() {
		
	}
	
	public Action shoot(GameState pState, Deadline pDue) {
	    if(currentRound != pState.getRound()){
	        turn = 0;
	        currentRound = pState.getRound();
	    }
	    turn ++;
	    
	    if(turn <= RUN || !SHOOTING)     return cDontShoot;

	    int numBirds = pState.getNumBirds();
	    
	    double birdsProb[] = new double[numBirds];
	    int birdNextMove[] = new int[numBirds];
	    boolean birdsConf[] = new boolean[numBirds];
	    
	    int[] lGuess = guess(pState, pDue);
	    
	    for(int i=0; i<numBirds;i++){
	        Bird b = pState.getBird(i);
	        if(b.isAlive() && (!startBlackStork || lGuess[i] == 0 || lGuess[i] == 2 || lGuess[i] == 4 )){
	            int[] sequence = birdSequenceToArray(b);
	            
	            HMM hmm = new HMM(FLY_PATH, Constants.COUNT_MOVE, sequence.length, sequence);
	            hmm.modelEstimation();
	            birdsConf[i] = hmm.calculateNextMove(birdsProb, birdNextMove, i);
	            
	        }

	    }

	    int nextMove = Constants.COUNT_MOVE;
	    int indexBird = -1;
	    double m = 0.0;
	    for(int i=0; i<numBirds;i++){
	        if(birdsProb[i] > m && pState.getBird(i).isAlive() && birdsConf[i]){
	            m = birdsProb[i];
	            nextMove = birdNextMove[i];
	            indexBird = i;
	        }
	    }

	    if(nextMove == Constants.COUNT_MOVE)  return cDontShoot;

	    return new Action(indexBird, nextMove);

	}

	private int[] birdSequenceToArray(Bird b)
	{
		int l = b.getSeqLength();
		int n = 0;
	    for(int j=0; j<l && b.getObservation(j)!=-1; j++){
	    	n++;
	    }
	    
	    int[] sequence = new int[n];

	    for(int j=0; j<n; j++){
	        sequence[j]=b.getObservation(j);
	    }
	
	    return sequence;
	}
	
	public int[] guess(GameState pState, Deadline pDue) {
		System.err.println("-----GUESS----");
		
		int[] lGuess = new int[pState.getNumBirds()];
		for (int i = 0; i < pState.getNumBirds(); ++i){
			lGuess[i] = Constants.SPECIES_UNKNOWN;
			lastGuess[i] = lGuess[i];
		}
		if(!GUESSING) return lGuess;
		
		
		if(pState.getRound()==0){
			for (int i = 0; i < pState.getNumBirds(); ++i){
				lGuess[i] = Constants.SPECIES_PIGEON;
				lastGuess[i] = lGuess[i];
			}
			return lGuess;
		}
		
		int nBirds = pState.getNumBirds();
		
		for(int i=0; i<nBirds; i++){
			double[][] freq = countFrequence(pState.getBird(i));
	        
			double min = 1000;
			int iMin = -1;
			for(int j=0; j<Constants.COUNT_SPECIES; j++){
				double diff = calcDiff(species[j], freq);
				if(diff < min){
					min = diff;
					iMin = j;
				}
			}
			System.err.println(min);
			
			if(min < MIN_GUESS){
			    lGuess[i] = iMin;
			    tryGuess++;
			}
            else
                lGuess[i] = Constants.SPECIES_UNKNOWN;
		}
		
		return lGuess;
	}

	
	private double calcDiff(double[][] specie, double[][] freq) {
		double diff = 0.0;
		
		for(int i=0; i<Constants.COUNT_MOVE; i++){
			for(int j=0; j<Constants.COUNT_MOVE; j++){
				double d = (specie[i][j] - freq[i][j])*(specie[i][j] - freq[i][j]);
				diff +=  d;
			}
		}
		return diff;
	}

	public void hit(GameState pState, int pBird, Deadline pDue) {
		System.err.println("HIT BIRD!!!");
	}


	public void reveal(GameState pState, int[] pSpecies, Deadline pDue) {
		
		for(int i=0; i<pSpecies.length; i++){
			if(pSpecies[i]!=-1){
				knowSpecies[pSpecies[i]] = true;
				if(knowSpecies[0] && knowSpecies[1] && knowSpecies[2] && knowSpecies[3] && knowSpecies[4]) startBlackStork=true;
				
				double[][] freq = countFrequence(pState.getBird(i));
		        
				updateSpecies(weightSpecies[pSpecies[i]], species[pSpecies[i]], freq);
				weightSpecies[pSpecies[i]] += 1.0;
			}
		}
	}
	
	private void updateSpecies(double w, double[][] spec, double[][] freq) {
		for(int i=0; i<Constants.COUNT_MOVE; i++){
			for(int j=0; j<Constants.COUNT_MOVE; j++){
				spec[i][j] = ((spec[i][j]*w)+freq[i][j])/(w+1);
			}
		}
	}

	public static final Action cDontShoot = new Action(-1, -1);

	private double[][] countFrequence(Bird b){
		int l = b.getSeqLength();
		double[][] f = new double[Constants.COUNT_MOVE][Constants.COUNT_MOVE];
		
		for(int i=0; i<l-1 && b.getObservation(i+1)!=-1; i++){
			f[b.getObservation(i)][b.getObservation(i+1)] ++;
		}
		
		return f;
	}
}
