import java.util.Scanner;

public class Main {
    public static Scanner sc = new Scanner(System.in);
    public static void main(String[] args) {
        int ch = 0;
        while(true) {
            System.out.println("1. Add matrices");
            System.out.println("2. Multiply matrix by a constant");
            System.out.println("3. Multiply matrices");
            System.out.println("4. Transpose matrix");
            System.out.println("5. Calculate a determinant");
            System.out.println("6. Inverse matrix");
            System.out.println("0. Exit");
            System.out.print("Your choice: ");
            ch = sc.nextInt();
            switch(ch) {
                case 1: {
                    add();
                    break;
                }
                case 2: {
                    matrixbyconstant();
                    break;
                }
                case 3: {
                    multiplymatrix();
                    break;
                }
                case 4: {
                    transpose();
                    break;
                }
                case 5: {
                    determinant();
                    break;
                }
                case 6: {
                    inverse();
                    break;
                }
                case 0: {
                    System.exit(0);
                }
            }
        }
    }

    public static void inverse() {
        int m = sc.nextInt();
        int n = sc.nextInt();
        double[][] array = new double[m][n];
        input(array, m, n);
        double [][] tarr = invert(array);
        double determinant = findDeterminant(array);
        if(determinant == 0) {
            System.out.println("This matrix doesn't have an inverse.");
        }
        else {
            System.out.println("The result is:");
            for(int i = 0; i < m; i++) {
                for(int j = 0; j < n; j++) {
                    System.out.print(tarr[i][j] + " ");
                }
                System.out.println();
            }
        }
    }

    public static void determinant() {
        int m = sc.nextInt();
        int n = sc.nextInt();
        double[][] array = new double[m][n];
        input(array, m, n);
        double ans = findDeterminant(array);
        System.out.println("The result is:\n" + ans);
    }

    public static double findDeterminant(double[][] matrix) {
        int m = matrix.length;
        if (m == 1) {
            return matrix[0][0];
        }
        if (m == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }
        double determinant = 0;
        for (int i = 0; i < m; i++) {
            int sign = i % 2 == 0 ? 1 : -1;
            determinant += matrix[0][i] * sign * findDeterminant(getMatrix(i, matrix));
        }
        return determinant;
    }

    public static double[][] invert(double a[][]) 
    {
        int n = a.length;
        double x[][] = new double[n][n];
        double b[][] = new double[n][n];
        int index[] = new int[n];
        for (int i=0; i<n; ++i) 
            b[i][i] = 1;

        // Transform the matrix into an upper triangle
        gaussian(a, index);

        // Update the matrix b[i][j] with the ratios stored
        for (int i=0; i<n-1; ++i)
            for (int j=i+1; j<n; ++j)
                for (int k=0; k<n; ++k)
                    b[index[j]][k]
                    -= a[index[j]][i]*b[index[i]][k];

        // Perform backward substitutions
        for (int i=0; i<n; ++i) 
        {
            x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
            for (int j=n-2; j>=0; --j) 
            {
                x[j][i] = b[index[j]][i];
                for (int k=j+1; k<n; ++k) 
                {
                    x[j][i] -= a[index[j]][k]*x[k][i];
                }
                x[j][i] /= a[index[j]][j];
            }
        }
        return x;
    }

    // Method to carry out the partial-pivoting Gaussian
    // elimination.  Here index[] stores pivoting order.

    public static void gaussian(double a[][], int index[]) 
    {
        int n = index.length;
        double c[] = new double[n];

        // Initialize the index
        for (int i=0; i<n; ++i) 
            index[i] = i;

        // Find the rescaling factors, one from each row
        for (int i=0; i<n; ++i) 
        {
            double c1 = 0;
            for (int j=0; j<n; ++j) 
            {
                double c0 = Math.abs(a[i][j]);
                if (c0 > c1) c1 = c0;
            }
            c[i] = c1;
        }

        // Search the pivoting element from each column
        int k = 0;
        for (int j=0; j<n-1; ++j) 
        {
            double pi1 = 0;
            for (int i=j; i<n; ++i) 
            {
                double pi0 = Math.abs(a[index[i]][j]);
                pi0 /= c[index[i]];
                if (pi0 > pi1) 
                {
                    pi1 = pi0;
                    k = i;
                }
            }

            // Interchange rows according to the pivoting order
            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i=j+1; i<n; ++i) 	
            {
                double pj = a[index[i]][j]/a[index[j]][j];

                // Record pivoting ratios below the diagonal
                a[index[i]][j] = pj;

                // Modify other elements accordingly
                for (int l=j+1; l<n; ++l)
                    a[index[i]][l] -= pj*a[index[j]][l];
            }
        }
    }

    public static double[][] getMatrix(int j, double[][] matrix) {
        double[][] newMatrix = new double[matrix.length - 1][matrix.length - 1];
        for (int i = 1; i < matrix.length; i++) {
            for (int k = 0; k < j ; k++) {
                newMatrix[i - 1][k] = matrix[i][k];
            }
        }
        for (int i = 1; i < matrix.length; i++) {
            for (int k = j + 1; k < matrix.length; k++) {
                newMatrix[i - 1][k - 1] = matrix[i][k];
            }
        }
        return newMatrix;
    }

    public static void add() {
        System.out.print("Enter size of first matrix: ");
        int m1 = sc.nextInt();
        int n1 = sc.nextInt();
        double[][] ar1 = new double [m1][n1];
        System.out.println("Enter first matrix:");
        input(ar1,m1,n1);
        System.out.print("Enter size of second matrix: ");
        int m2 = sc.nextInt();
        int n2 = sc.nextInt();
        System.out.println("Enter second matrix:");
        double[][] ar2 = new double [m2][n2];
        input(ar2,m2,n2);
        if((m1 == m2) && (n1 == n2)) {
            for(int i = 0; i < m1; i++) {
                for(int j = 0; j < n1; j++) {
                    System.out.print(ar1[i][j] + ar2[i][j] +" ");
                }
                System.out.println();
            }
        }
        else {
            System.out.println("The operation cannot be performed.");
        }
    }

    public static void matrixbyconstant() {
        System.out.print("Enter size of matrix: ");
        int m1 = sc.nextInt();
        int n1 = sc.nextInt();
        double[][] ar1 = new double [m1][n1];
        System.out.println("Enter matrix:");
        input(ar1,m1,n1);
        System.out.print("Enter constant: ");
        int c = sc.nextInt();
        System.out.println("The result is:");
        for (int i = 0; i < m1; i++) {
            for (int j = 0; j < n1; j++) {
                System.out.print((c * ar1[i][j]) + " ");
            }
            System.out.println();
        }    
    }

    public static void multiplymatrix() {
        System.out.print("Enter size of first matrix: ");
        int m1 = sc.nextInt();
        int n1 = sc.nextInt();
        double[][] ar1 = new double [m1][n1];
        input(ar1,m1,n1);
        System.out.print("Enter size of second matrix: ");
        int m2 = sc.nextInt();
        int n2 = sc.nextInt();
        double[][] ar2 = new double [m2][n2];
        input(ar2,m2,n2);
        if(n1 == m2){
            System.out.println("The result is:");
            double sum;
            for(int i = 0; i < m1; i++) {
                for(int j = 0; j < n2; j++) {
                    sum = 0;
                    for(int k = 0; k < n1; k++){
                        sum += ar1[i][k] * ar2[k][j];
                    }
                    System.out.print(sum +" ");
                }
                System.out.println();
            }
        }  
    }

    public static void input(double[][] array,int m,int n) {
        for(int i = 0; i < m; i++) {
            for(int j = 0; j < n; j++) {
                array[i][j] = sc.nextDouble();
            }
        }
    }

    public static void transpose() {
        System.out.println("1. Main diagonal");
        System.out.println("2. Side diagonal");
        System.out.println("3. Vertical line");
        System.out.println("4. Horizontal line");
        int choice = sc.nextInt();
        switch(choice) {
            case 1: {
                int m = sc.nextInt();
                int n = sc.nextInt();
                double[][] arr = new double [m][n];
                input(arr,m,n);
                double [][] tarr = new double[m][n];
                for(int i = 0; i < m; i++) {
                    for(int j = 0; j < n; j++) {
                        tarr[i][j] = arr[j][i];
                    }
                }
                System.out.println("The result is:");
                for(int i = 0; i < m; i++) {
                    for(int j = 0; j < n; j++) {
                        System.out.print(tarr[i][j]+" ");
                    }
                    System.out.println();
                }
                break;
            }
            case 2: {
                int m = sc.nextInt();
                int n = sc.nextInt();
                double[][] arr = new double [m][n];
                input(arr,m,n);
                double [][] tarr = new double [m][n];
                for(int i = 0; i < m; i++) {
                    for(int j = 0; j < n; j++) {
                        tarr[i][j] = arr[m-j-1][n-i-1];
                    }
                }
                System.out.println("The result is:");
                for(int i = 0; i < m; i++) {
                    for(int j = 0; j < n; j++) {
                        System.out.print(tarr[i][j]+" ");
                    }
                    System.out.println();
                }
                break;
            }
            case 3: {
                int m = sc.nextInt();
                int n = sc.nextInt();
                double[][] arr = new double [m][n];
                input(arr,m,n);
                double[][] tarr = new double [m][n];
                for(int i = 0; i < m; i++) {
                    for(int j = 0; j < n; j++) {
                        tarr[j][i] = arr[j][m-i-1];   
                    }
                }
                for(int i = 0; i < m; i++) {
                    for(int j = 0; j < n; j++) {
                        System.out.print(tarr[i][j]+" ");
                    }
                    System.out.println();
                }
                break;
            }
            case 4: {
                int m = sc.nextInt();
                int n = sc.nextInt();
                double[][] arr = new double [m][n];
                input(arr,m,n);
                double[][] tarr = new double [m][n];
                for(int i = 0; i < m; i++) {
                    for(int j = 0; j < n; j++) {
                        tarr[i][j] = arr[n-i-1][j]; 
                    }
                }
                System.out.println("The result is:");
                for(int i = 0; i < m; i++) {
                    for(int j = 0; j < n; j++) {
                        System.out.print(tarr[i][j]+" ");
                    }
                    System.out.println();
                }
                break;
            }
        }
    }
}