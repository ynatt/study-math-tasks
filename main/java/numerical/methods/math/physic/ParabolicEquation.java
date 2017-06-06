package numerical.methods.math.physic;


import java.io.FileWriter;
import java.io.IOException;

public class ParabolicEquation {

    public static void findSolve(int stepNumber,int layNumber){
        double h = (Math.PI/(2*stepNumber));
        double tau = (0.9/layNumber);
        double[][] solve = new double[layNumber+1][stepNumber+1];
        double gamma = tau/(h*h);
        double A = gamma/2;
        double C = A;
        double B = -1-gamma;
        double[] alpha = new double[stepNumber+1];
        double[] beta = new double[stepNumber+1];
        double[] F = new double[stepNumber];

        for(int i=0;i<=stepNumber;i++){
            solve[0][i]=startCondition(i*h);
        }
        for(int k =0;k<layNumber;k++){
            for(int i=1;i<stepNumber;i++){
                F[i] = A*solve[k][i-1] +(1-gamma)*solve[k][i] +
                        A*solve[k][i+1] + tau*f(i*h,(k+0.5)*tau);
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
        showGraf(solve,layNumber,stepNumber,h,tau);
    }

    public static void showGraf(double[][] solve,int layNumber,int
            stepNumber,double h,double tau){
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

    public static double f(double x,double t){
        return 2*Math.cos(x+t);
    }
    public static double startCondition(double x){
        return Math.cos(x) + Math.sin(x);
    }
    public static double firstEdgeCondition(double t){
        return Math.cos(t)+Math.sin(t);
    }
    public static double secondEdgeCondition(double t){
        return Math.cos(t)-Math.sin(t);
    }

    public static void main(String[] args) {
        findSolve(50,50);
    }
}
