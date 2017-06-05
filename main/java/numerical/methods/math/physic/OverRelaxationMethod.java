package numerical.methods.math.physic;

import java.io.FileWriter;
import java.io.IOException;

import static java.lang.Math.*;

public  class OverRelaxationMethod {
    private double lengthX;
    private double lengthY;
    private int stepNumberX;
    private int stepNumberY;
    private double stepX;
    private double stepY;
    private double delta;
    private double omega;
    private double A;
    private double B;
    private double C;
    private double epsilon;
    private double[][] solve;

    public OverRelaxationMethod(double lengthX, double lengthY, int stepNumberX, int stepNumberY, double factor) {
        this.lengthX = lengthX;
        this.lengthY = lengthY;
        this.stepNumberX = stepNumberX + 1;
        this.stepNumberY = stepNumberY + 1;
        solve = new double[this.stepNumberX][this.stepNumberY];
        stepX = lengthX / stepNumberX;
        stepY = lengthY / stepNumberY;
        delta = (2 * pow(stepY,2) * pow(sin((PI * pow(stepX,2)) / (2 * lengthX)),2)) / (pow(stepX,2) + pow(stepY,2)) +
                (2 * pow(stepX,2) * pow(sin((PI * pow(stepY,2)) / (2 * lengthY)),2)) / (pow(stepX,2) + pow(stepY,2));
        omega = 2 / (1 + sqrt(delta * (2 - delta)));
        A = -1/pow(stepX,2);
        B = -1/pow(stepY,2);
        C = 2/pow(stepX,2) + 2/pow(stepY,2);
        epsilon = factor * (pow(stepX,2) + pow(stepY,2));
        setEdgeValues();
        setInnerValues(0.5);
    }

    private void setEdgeValues(){
        for(int i = 0; i < stepNumberY ; i++ ){
            solve[stepNumberX-1][i] = 40;
            solve[0][i] = 40 * pow(i * stepY,2);
        }
        for(int i = 0; i < stepNumberX; i++ ){
            solve[i][0] = 40;
            solve[i][stepNumberY-1] = 40 * sin((PI * i * stepX)/2);
        }
    }

    private void setInnerValues(double innerValue){
        for (int i = 1; i < stepNumberX - 1; i++ ){
            for (int j = 1; j < stepNumberY - 1; j++ ){
                solve[i][j]=innerValue;
            }
        }
    }

    public void showGraf(){
        for (int i = 0; i < stepNumberX; i++ ){
            for (int j = 0; j < stepNumberY; j++ ){
                System.out.print(solve[i][j] + " ");
            }
            System.out.println();
        }
    }

    private double countDiscrepancy( int i, int j){
        return A * solve[i-1][j] + B * solve[i][j-1] + C * solve[i][j] + A * solve[i+1][j] + B * solve[j][i+1];
    }

    public double countNorm(){
        double r = 0;
        for(int i = 1; i < stepNumberX - 1; i++ ){
            for(int j = 1; j < stepNumberY - 1; j++ ){
                r= r + pow(countDiscrepancy(i,j),2) * stepX * stepY;
           }
        }
        return r;
    }

    public void getSolve(){
        double discrepancy;
        double norm = epsilon;
        int iterations = 0;
        while(norm>=epsilon) {
            norm=0;
            for (int i = 1; i < stepNumberX - 1; i++) {
                for (int j = 1; j < stepNumberY - 1; j++) {
                    discrepancy = A * solve[i-1][j] + B * solve[i][j-1] + C * solve[i][j] + A * solve[i+1][j] + B * solve[i][j+1];
                    solve[i][j] = solve[i][j] - omega/C * discrepancy;
                    norm += pow(discrepancy,2) * stepX * stepY;
                }
            }
            iterations++;
        }
        System.out.println("Итераций: " + iterations );
    }

    public void exportSolve(){
        try(FileWriter writer = new FileWriter("solve.txt", false)) {
            for(int i = 0; i < stepNumberX;i++){
                for(int j=0;j<stepNumberY;j++){
                    writer.write(i * stepX + " " + j * stepY + " " + solve[i][j]+"\n");
                }
                writer.write("\n");
            }
            writer.flush();
        } catch(IOException ex){
            System.out.println(ex.getMessage());
        }
    }

    @Override
    public String toString() {
        return "OverRelaxationMethod{" +
                "lengthX=" + lengthX +
                ", lengthY=" + lengthY +
                ", stepNumberX=" + stepNumberX +
                ", stepNumberY=" + stepNumberY +
                ", stepX=" + stepX +
                ", stepY=" + stepY +
                ", delta=" + delta +
                ", omega=" + omega +
                ", A=" + A +
                ", B=" + B +
                ", C=" + C +
                ", epsilon=" + epsilon +
                '}';
    }

    public static void main(String[] args) {
        OverRelaxationMethod orm = new OverRelaxationMethod(1,1,40,40,1);
        System.out.println(orm);
        orm.showGraf();
        orm.getSolve();
        orm.exportSolve();
    }


}
