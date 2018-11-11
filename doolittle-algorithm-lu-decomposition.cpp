#include <bits/stdc++.h> 

using namespace std;

void getRandomMatrix(int n,vector<vector<double>> &matrix) {
    srand(0);
    for (int i=0 ; i<n ; i++) { 
        for (int j=0 ; j<n ; j++) {
            matrix[i][j] = rand()%150;
        }
    }
}

void LU_Decomposition(vector<vector<double>>& matrix) 
{ 

    ofstream opfile;
    opfile.open ("doolittle.txt");
    // Matrix Size
    int n = matrix.size();

    // The upper matrix
    vector<vector<double>> lower(n,vector<double>(n,0));
    // The Lower matrix
    vector<vector<double>> upper(n,vector<double>(n,0));
  

    double exec_time = clock();

    // Decomposing matrix into Upper and Lower 
    // triangular matrix 
    for (int i=0 ; i<n; i++) { 
  
        // Upper Triangular 
        for (int k=i; k<n ; k++) { 
  
            // Summation of L(i, j) * U(j, k) 
            double sum = 0; 
            for (int j=0 ; j<i ; j++) 
                sum += (lower[i][j] * upper[j][k]); 
  
            // Evaluating U(i, k) 
            upper[i][k] = matrix[i][k] - sum; 
        } 
  
        // Lower Triangular 
        for (int k=i ; k<n ; k++) { 
            if (i == k) 
                lower[i][i] = 1; // Diagonal as 1 
            else { 
  
                // Summation of L(k, j) * U(j, i) 
                double sum = 0; 
                for (int j=0 ; j<i ; j++) 
                    sum += (lower[k][j] * upper[j][i]); 
  
                // Evaluating L(k, i) 
                lower[k][i] = (matrix[k][i] - sum) / upper[i][i]; 
            } 
        } 
    } 
    
    exec_time = clock() - exec_time;

    // Output
    cout<<"Lower Triangular : "<<endl; 
  
    for (int i=0 ; i<n ; i++) { 
        for (int j=0 ; j<n ; j++) 
            cout<<lower[i][j]<<"\t";
        cout<<endl;  
    }

    cout<<endl<<endl;

    cout<<"Upper Triangular : "<<endl; 
  
    for (int i=0 ; i<n ; i++) { 
        for (int j=0 ; j<n ; j++) 
            cout<<upper[i][j]<<"\t"; 
        cout<<endl;  
    } 

    // cout<<"TIME TAKEN : "<<((double)exec_time/CLOCKS_PER_SEC)<<" milliseconds.";
} 

// Driver code 
int main() 
{ 
    int n=(2<<3);
    vector<vector<double>> matrix(n,vector<double>(n,0));
    getRandomMatrix(n,matrix);
  
    for (int i=0 ; i<n ; i++) { 
        for (int j=0 ; j<n ; j++) 
            cout<<matrix[i][j]<<"\t"; 
        cout<<endl;  
    } 

    LU_Decomposition(matrix); 
    return 0; 
} 