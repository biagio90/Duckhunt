

public class HMM {
	static private final int LIMIT = 50;
	static private final double CONFIDENZIAL = 0.79;
	static private final int PERC = 80;
	
	private double[][] a, b;
	private double[] pi;
	private int N, M, T;
	private int[] O;
	

	public HMM(double a[][], double b[][], double pi[][], int N, int M, int T, int[] seq){
		this.N = N;
		this.M = M;
		this.T = T;
		this.O = seq;
		
		this.a = a;
		this.b = b;
		this.pi = pi[0];
		
		for(double[] row : this.a){
			for(double ele : row){
				System.err.print(ele+" ");
			}
			System.err.println();
		}
		System.err.println();
		for(double[] row : this.b){
			for(double ele : row){
				System.err.print(ele+" ");
			}
			System.err.println();
		}
		System.err.println();
		for(double ele : this.pi){
			System.err.print(ele+" ");
		}
		System.err.println();
		for(int ele : this.O){
			System.err.print(ele+" ");
		}
		
		System.err.println();
		
	}
	
	
	HMM(int N, int M, int T, int[] seq){
		this.N = N;
		this.M = M;
		this.T = T;
		this.O = seq.clone();
		
		a = new double[N][N];
		b = new double[N][M];
		pi = new double[N];
		
		for(int i = 0; i<N; i++){
			double sum = 0.0;
			for(int j = 0; j<N; j++){
				double f = N;
	            f = 1/f;
	            double rn = (Math.random()*PERC)-(PERC/2);
	            rn = rn/100;
				a[i][j] = f+f*rn;
				sum += a[i][j];
			}
			for(int j=0; j<N; j++) {
	            a[i][j] = a[i][j] / sum;
	        }
		}
		

		for(int i = 0; i<N; i++){
			double sum = 0.0;
			for(int j = 0; j<M; j++){
				double f = M;
	            f = 1/f;
	            double rn = (Math.random()*PERC)-(PERC/2);
	            rn = rn/100;
				b[i][j] = f+f*rn;
				sum += b[i][j];
			}
			for(int j=0; j<N; j++) {
	            b[i][j] = b[i][j] / sum;
	        }
		}

        double sum = 0.0;
	    for(int i=0; i<N; i++){
	        double f = N;
            f = 1/f;
            double rn = (Math.random()*PERC)-(PERC/2);
            rn = rn/100;
			pi[i] = f+f*rn;
			sum += pi[i];
	    }
        for(int j=0; j<N; j++) {
            pi[j] = pi[j] / sum;
        }
	}
	
	HMM(){
		
	}

	public void modelEstimation()
	{
	    int i=0;
	    double logProb = -Double.MAX_VALUE, oldLogProb = -Double.MAX_VALUE;
	    double c[] = new double[T];
	    
	    do{
	        updateHMM(c);
	        i++;
	
	        oldLogProb = logProb;
	        logProb = 0.0;
	        for(int t=0; t<T; t++){
	            logProb += Math.log(c[t]);
	        }
	        logProb = -logProb;
	    }while(logProb > oldLogProb && i<LIMIT);
	}
	

	private void updateHMM(double c[])
	{
		double alpha[][] = new double[T][N];
		double beta[][] = new double[T][N];
		double Gamma[][] = new double[T][N];
		double diGamma[][][] = new double[T][N][N];
		
	    calculateAlphaBeta(alpha, beta, c);
	    calculateGamma(Gamma, diGamma, alpha, beta);
	
	
	    for(int i=0; i<N; i++){
	        pi[i] = Gamma[0][i];
	        for(int j=0; j<N; j++){
	            double num = 0.0, denom = 0.0;
	            for(int t=0; t<T-1; t++){
	                num += diGamma[t][i][j];
	                denom += Gamma[t][i];
	            }
	            a[i][j] = num / denom;
	        }
	        for(int j=0; j<M; j++){
	            double num = 0.0, denom = 0.0;
	            for(int t=0; t<T-1; t++){
	                if(O[t] == j) num += Gamma[t][i];
	                denom += Gamma[t][i];
	            }
	            b[i][j] = num / denom;
	        }
	    }
	
	}
	
		
	private void calculateGamma(double Gamma[][], double diGamma[][][], double alpha[][], double beta[][])
	{
	    for(int t=0; t<T-1; t++){
	        double denom = 0.0;
	        for(int i=0; i<N; i++){
	            for(int j=0; j<N; j++){
	                denom += alpha[t][i] * a[i][j] * b[j][O[t+1]] * beta[t+1][j];
	            }
	        }
	        for(int i=0; i<N; i++){
	            Gamma[t][i] = 0.0;
	            for(int j=0; j<N; j++){
	                diGamma[t][i][j] = (alpha[t][i] * a[i][j] * b[j][O[t+1]] * beta[t+1][j]) / denom;
	                Gamma[t][i] += diGamma[t][i][j];
	            }
	        }
	    }
	}
	
	private void calculateAlphaBeta(double[][] alpha, double[][] beta, double c[]){
		//a pass
	    c[0]=0.0;
	    for(int i=0; i<N; i++){
	        alpha[0][i] = pi[i] * b[i][O[0]];
	        c[0] += alpha[0][i];
	    }

	    c[0]=1/c[0];
	    for(int i=0; i<N; i++){
	        alpha[0][i] = c[0] * alpha[0][i];
	    }

	    for(int t=1; t<T; t++){
	        c[t]=0.0;
	        for(int i=0; i<N; i++){
	            alpha[t][i] = 0.0;
	            for(int j=0; j<N; j++){
	                alpha[t][i] += alpha[t-1][j] * a[j][i];
	            }
	            alpha[t][i] = alpha[t][i] * b[i][O[t]];
	            c[t] += alpha[t][i];
	        }

	        c[t]=1/c[t];
	        for(int i=0; i<N; i++){
	            alpha[t][i] = c[t] * alpha[t][i];
	        }
	    }

	//b pass
	    for(int i=0; i<N; i++){
	        beta[T-1][i] = c[T-1];
	    }

	    for(int t=T-2; t >= 0; t--){
	        for(int i=0; i<N; i++){
	            beta[t][i] = 0.0;
	            for(int j=0; j<N; j++){
	                beta[t][i] += a[i][j] * beta[t+1][j] * b[j][O[t+1]];
	            }
	            beta[t][i] = c[t] * beta[t][i];
	        }

	    }
	}


	public void printA() {
		System.err.print(N+" ");
		System.err.print(N+" ");
		for(double[] row : a){
			for(double ele : row){
				if(ele < 0.001) System.err.print("0.0 ");
				else			System.err.print(ele+" ");
			}
		}
	}


	public void printB() {
		System.err.print(N+" ");
		System.err.print(M+" ");
		for(double[] row : b){
			for(double ele : row){
				if(ele < 0.001) System.err.print("0.0 ");
				else			System.err.print(ele+" ");
			}
		}
	}
	
	
	private double probabilityEmission(int k, double[][] alpha, double[][] beta)
	{
	
	    double p = 0.0;
	    for(int i=0; i<N; i++){
	        double sum = 0.0, denom = 0.0;
	        for(int j=0; j<N; j++){
	            sum += a[j][i]*alpha[T-1][j];
	            denom += alpha[T-1][j];
	        }
	        sum = sum / denom;
	        p += b[i][k] * sum;
	    }
	
	    return p;
	}
	
	public boolean calculateNextMove(double[] pNextMove, int[] nextMove, int index)
	{
		double alpha[][] = new double[T][N];
		double beta[][]  = new double[T][N];
	    double c[] 		 = new double[T];
		
	    calculateAlphaBeta(alpha, beta, c);
	    
		double p[] = new double[M];
		
	    for(int k=0; k<M; k++) p[k] = probabilityEmission(k, alpha, beta);
	
	    double m = -1.0;
	    int im = -1;
	    for(int k=0; k<M; k++){
	        if(p[k] > m){
	            m = p[k];
	            im = k;
	        }
	        //std::cerr << p << "\n";
	    }
	
	    double m2 = -1.0;
	    for(int k=0; k<M; k++){
	        if(p[k] > m2 && k != im)   m2 = p[k];
	    }
	
	    pNextMove[index] = m;
	    nextMove[index] = im;
	
	    if( m > CONFIDENZIAL) 		return true;
	    else                        return false;
	
	}


	public double getProbabilitySequence(int[] seq) {
		O = seq.clone();
		T = seq.length;
		
		double alpha[][] = new double[T][N];
		double beta[][]  = new double[T][N];
	    double c[] 		 = new double[T];
		
	    for(int i=0; i<N; i++) pi[i]=1/((double)N);
	    
	    calculateAlphaBeta(alpha, beta, c);
	    
	    double ret = 0;
	    for(int i=0; i<N; i++){
	        ret += alpha[T-1][i];
	    }

	    return ret;
	}

}
