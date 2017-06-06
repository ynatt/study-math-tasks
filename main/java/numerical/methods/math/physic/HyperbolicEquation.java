package numerical.methods.math.physic;

import java.io.FileWriter;
import java.io.IOException;

public class HyperbolicEquation {
    public static void main(String[] args) {
        getSolution(50,50);
    }

    public static void getSolution(int stepNumber,int layNumber){
        double h = (1d/stepNumber);
        double tau = (1d/layNumber);
        double[][] solve = new double[layNumber+1][stepNumber+1];
        double gamma=(tau*tau)/(h*h);
        double A = gamma*gamma*0.5;
        double C = A;
        double B = -1-gamma*gamma;
        double[] alpha = new double[stepNumber+1];
        double[] beta = new double[layNumber+1];
        double[] F = new double[stepNumber];

        for(int i=0;i<=stepNumber;i++){
            solve[0][i]=firstStartCondition(i*h);
            solve[1][i]=secondStartCondition(i*h,tau);
        }

        for(int k =1;k<layNumber;k++){
            for(int i=1;i<stepNumber;i++){
                F[i] =
                        (2*solve[k][i]-solve[k-1][i])+0.5*tau*tau*lambda(solve,h,k-1,i)+tau*tau*f(i*h,k*tau);
            }
            solve[k+1][0]=firstEdgeCondition((k+1)*tau);
            solve[k+1][stepNumber]=secondEdgeCondition((k+1)*tau);
            alpha[1]=0;
            beta[1]=solve[k+1][0];
            for(int i=1;i<stepNumber;i++){
                alpha[i+1]= -C / (A*alpha[i] + B);
                beta[i+1]=(-F[i] - A*beta[i])/(A*alpha[i] + B);
            }
            for(int i=stepNumber-1;i>0;i--){
                solve[k+1][i]= alpha[i+1]*solve[k+1][i+1] + beta[i+1];
            }

        }
        showGraph(solve,layNumber,stepNumber,h,tau);
    }

    public static void showGraph(double[][] solve, int layNumber, int
            stepNumber, double h, double tau){
        try {
            FileWriter writer = new FileWriter("solve.txt", false);
            for(int k=0;k<layNumber+1;k++){
                for(int i=0;i<stepNumber+1;i++){
                    writer.write(i*h+" "+k*tau+" "+solve[k][i]+"\n");
                }
                writer.write("\n");
            }
            writer.flush();
        } catch(IOException ex){

            System.out.println(ex.getMessage());
        }
    }
    public static double lambda(double[][] solve,double h,int k,int i){
        return (solve[k][i-1]-2*solve[k][i]+solve[k][i+1])/(h*h);
    }

    public static double f(double x,double t){
        return x+t;
    }
    public static double firstStartCondition(double x){
        return 1.5*(x*x+0.8)*Math.sin(Math.PI*x);
    }

    public static double secondStartCondition(double x,double tau){
        return 1.5*(x*x+0.8)*Math.sin(Math.PI*x) + tau*0.1*x
                +
                0.5*tau*tau*(3-0.8*Math.PI*Math.PI*Math.sin(Math.PI*x)+f(x,0));
    }

    public static double firstEdgeCondition(double t){
        return 0;
    }
    public static double secondEdgeCondition(double t){
        return 0;
    }
}
